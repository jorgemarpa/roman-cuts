"""
roman_cutouts.cutout
~~~~~~~~~~~~~~~~~~~~
Pixel cutout extraction from Roman WFI L2 ASDF files using GWCS.
"""

import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D, NoOverlapError
from astropy.time import Time
from astropy.wcs.wcsapi import HighLevelWCSWrapper, SlicedLowLevelWCS
from tqdm.auto import tqdm

from .io import _detect_platform, _open_asdf

logger = logging.getLogger(__name__)

# ASDF metadata keys retained per cutout
_KEEP_CUTOUT = [
    "input_position_original",
    "origin_original",
    "input_position_cutout",
    "bbox_original",
]
_KEEP_EXP = [
    "effective_exposure_time",
    "exposure_time",
    "frame_time",
    "start_time",
    "end_time",
    "mid_time",
]
_KEEP_OBS = ["exposure", "observation", "observation_id", "visit"]
_KEEP_INSTR = ["detector", "optical_element"]



def _is_int(x):
    return isinstance(x, (int, np.integer))


def _normalize_targets(ra, dec, size):
    """
    Normalise target coordinates and sizes to equal-length per-target lists.

    Parameters
    ----------
    ra, dec : float or sequence of float
        Target coordinates (deg).  Scalars are treated as one target.
    size : int, tuple, or list
        ``int`` → square shared by all targets; ``(ny, nx)`` tuple → shared
        rectangle; list of ints/tuples → one size per target (a list of
        length one is broadcast).

    Returns
    -------
    ras, decs : ndarray, shape (n_targets,)
    sizes : list of (ny, nx)
    """
    ras = np.atleast_1d(np.asarray(ra, dtype=float))
    decs = np.atleast_1d(np.asarray(dec, dtype=float))
    if ras.ndim != 1 or ras.shape != decs.shape:
        raise ValueError(
            f"ra and dec must be scalars or 1-D sequences of equal length "
            f"(got shapes {ras.shape} and {decs.shape})."
        )
    n = len(ras)

    def _one(s):
        if _is_int(s):
            return (int(s), int(s))
        t = tuple(s)
        if len(t) == 2 and all(_is_int(x) for x in t):
            return (int(t[0]), int(t[1]))
        raise ValueError(f"Invalid cutout size {s!r}; expected int or (ny, nx).")

    if _is_int(size) or (
        isinstance(size, tuple) and len(size) == 2 and all(_is_int(x) for x in size)
    ):
        sizes = [_one(size)] * n
    else:
        seq = list(size)
        if len(seq) == 1:
            seq = seq * n
        if len(seq) != n:
            raise ValueError(
                f"size has {len(seq)} entries but there are {n} targets. Use an "
                "int, one (ny, nx) tuple, or a list with one entry per target."
            )
        sizes = [_one(s) for s in seq]
    return ras, decs, sizes


def _cut_plane(arr, flx, size):
    """Cut *arr* with the slices of an existing `Cutout2D`, NaN-padded to *size*."""
    out = np.full(size, np.nan)
    out[flx.slices_cutout] = arr[flx.slices_original]
    return out


def _placeholder(base_row, ra, dec):
    """All-NaN entry (out-of-bounds or failed): arrays are filled in later."""
    row = dict(base_row)
    row["ra_target"] = ra
    row["dec_target"] = dec
    row["in_bounds"] = False
    return {"flx": None, "err": None, "dq": None, "dumo": None, "gwcs": None, "row": row}


def _extract_file(idx, fname, ras, decs, sizes, lite, verbose, platform, mission, strict):
    """
    Extract the cutouts of *all* targets from one file (opened once).

    Module-level so it can be pickled for the process backend.

    Returns
    -------
    dict
        ``{"idx", "mid_jd", "targets"}`` where ``targets`` is a list (one per
        target) of dicts with ``flx``, ``err``, ``dq``, ``dumo``, ``gwcs``
        and ``row``.  Out-of-bounds targets or failed files give placeholders.
    """
    try:
        with _open_asdf(fname, mission=mission, platform=platform) as af:
            roman = af["roman"]
            meta = roman["meta"]
            gwcs_obj = meta["wcs"]

            # Exposure-level metadata is shared by all targets: read it once
            exp_node = meta["exposure"]
            mid_t = exp_node.get("mid_time")
            if mid_t is None:
                mid_t = Time(
                    (Time(exp_node["start_time"]).jd + Time(exp_node["end_time"]).jd)
                    / 2.0,
                    format="jd",
                )
            if not isinstance(mid_t, Time):
                mid_t = Time(mid_t)
            base_row = {k: exp_node[k] for k in _KEEP_EXP if k in exp_node}
            obs_node, instr_node = meta["observation"], meta["instrument"]
            base_row |= {k: obs_node[k] for k in _KEEP_OBS if k in obs_node}
            base_row |= {k: instr_node[k] for k in _KEEP_INSTR if k in instr_node}
            base_row["mid_time"] = mid_t

            # Sky → detector pixels for every target in one vectorised call (APE 14)
            coords = SkyCoord(ra=ras * u.deg, dec=decs * u.deg, frame="icrs")
            x_pix, y_pix = (np.atleast_1d(a) for a in gwcs_obj.world_to_pixel(coords))

            targets = []
            for t, (ra_t, dec_t, size) in enumerate(zip(ras, decs, sizes)):
                x, y = x_pix[t], y_pix[t]
                if not (np.isfinite(x) and np.isfinite(y)):
                    logger.debug("Skipping %s — target %d outside image bounds.", fname, t)
                    targets.append(_placeholder(base_row, ra_t, dec_t))
                    continue
                if verbose:
                    logger.info("  %s  target %d  →  pixel (%.2f, %.2f)", fname, t, x, y)
                try:
                    # Cutout2D reads the needed slice; .copy() outlives the file
                    flx = Cutout2D(
                        roman["data"],
                        position=(float(x), float(y)),
                        size=size,
                        mode="partial",
                        fill_value=np.nan,
                    )
                except NoOverlapError:
                    targets.append(_placeholder(base_row, ra_t, dec_t))
                    continue

                entry = {"flx": flx.data.copy(), "err": None, "dq": None, "dumo": None}
                if not lite:
                    for key in ("err", "dq", "dumo"):
                        entry[key] = _cut_plane(roman[key], flx, size)
                yslice, xslice = flx.slices_original
                entry["gwcs"] = HighLevelWCSWrapper(
                    SlicedLowLevelWCS(gwcs_obj, (yslice, xslice))
                )
                row = {k: flx.__dict__[k] for k in _KEEP_CUTOUT if k in flx.__dict__}
                row |= base_row
                row["ra_target"] = ra_t
                row["dec_target"] = dec_t
                row["in_bounds"] = True
                entry["row"] = row
                targets.append(entry)
            return {"idx": idx, "mid_jd": float(mid_t.jd), "targets": targets}
    except Exception:
        if strict:
            raise
        logger.exception("Failed to extract cutouts from %s; filling with NaN.", fname)
        return {
            "idx": idx,
            "mid_jd": np.nan,
            "targets": [_placeholder({}, r, d) for r, d in zip(ras, decs)],
        }


def make_cutout(
    file_names,
    ra,
    dec,
    size=(10, 10),
    mission=None,
    verbose=False,
    lite=False,
    n_workers=8,
    backend="thread",
    strict=False,
):
    """
    Create pixel cutout cubes for one or more targets from Roman WFI L2 ASDF files.

    The GWCS stored in each ASDF file maps (ra, dec) to detector pixel
    coordinates.  Each cutout is extracted with `~astropy.nddata.Cutout2D`
    (``mode='partial'``), which fills out-of-bounds pixels with NaN and
    handles coordinates near the detector edge gracefully.  Targets that fall
    outside a file's footprint are kept as NaN entries.

    Files are processed in parallel (one task per file; every target is cut
    from a file opened once) and the results are sorted by ``mid_time``.

    File access is dispatched by `~roman_cutouts.io._open_asdf`, which
    auto-detects the execution platform (RRN, SMDC, standalone) via
    `~roman_cutouts.io._detect_platform`.  Supply *mission* to force the
    RRN/MAST path from any platform.

    Parameters
    ----------
    file_names : str or list of str
        MAST product filenames (RRN), local/mounted-S3 paths (SMDC), plain
        local paths, or ``s3://`` URIs.
    ra : float or sequence of float
        Right Ascension(s) of the cutout centre(s) in degrees (ICRS).
    dec : float or sequence of float
        Declination(s) of the cutout centre(s) in degrees (ICRS).
    size : int, tuple of int, or list
        Cutout size in pixels.  An ``int`` gives a square cutout shared by
        all targets; a ``(ny, nx)`` tuple a shared rectangle; a list of ints
        or tuples gives one size per target (length one is broadcast).
    mission : `~astroquery.mast.MastMissions`, optional
        Authenticated mission object for RRN/MAST access.  With
        ``backend="process"`` it is not shared; workers create their own from
        ``MAST_API_TOKEN``.
    verbose : bool
        Print per-file pixel-coordinate messages when *True*.
    lite : bool
        Extract only the data plane (cube axis 0 has length 1).
    n_workers : int
        Number of parallel workers (capped at the number of files).  ``1``
        runs serially.  Default 8.
    backend : ``'thread'`` | ``'process'``
        ``'thread'`` (default) suits the I/O-bound workload and is safe in
        Jupyter and the CLI.  ``'process'`` uses a ``spawn`` process pool.
    strict : bool
        Re-raise per-file errors instead of logging them and returning NaN.

    Returns
    -------
    cutout_cubes : list of ndarray
        One per target, shape ``(n_types, n_files, ny, nx)`` where the first
        axis is ``[data, err, dq, dumo]`` (``n_types=4``, or 1 if *lite*).
    metadatas : list of `~pandas.DataFrame`
        One per target, ``n_files`` rows with exposure, observation and
        instrument columns, ``mid_time``, ``ra_target``, ``dec_target`` and
        Cutout2D geometry fields.  The boolean ``in_bounds`` column is
        *False* for NaN entries (target outside the file's footprint, or the
        file failed) and *True* otherwise, including partial overlaps, which
        are NaN-padded.
    gwcsss : list of list of `~astropy.wcs.wcsapi.HighLevelWCSWrapper`
        One per target, ``n_files`` sliced GWCS objects (*None* for NaN
        entries).  ``gwcsss[t][k].world_to_pixel(SkyCoord(...))`` returns
        pixel coordinates in the cutout frame.

    All outputs are sorted by ``mid_time``; entries with unknown time go last.
    A scalar *ra*/*dec* is one target, so each returned list has length 1.
    """
    if isinstance(file_names, str):
        file_names = [file_names]
    file_names = list(file_names)
    n_files = len(file_names)
    if backend not in ("thread", "process"):
        raise ValueError(f"backend must be 'thread' or 'process', got {backend!r}.")

    ras, decs, sizes = _normalize_targets(ra, dec, size)
    n_targets = len(ras)
    if verbose:
        for r, d in zip(ras, decs):
            logger.info("Target: RA=%.6f  Dec=%.6f", r, d)

    # Detect the platform once (it may probe the network); a mission forces RRN
    platform = "rrn" if mission is not None else _detect_platform()
    n_workers = max(1, min(int(n_workers), n_files))

    args = (ras, decs, sizes, lite, verbose, platform)
    results = [None] * n_files
    desc = "Extracting cutouts"
    if n_workers == 1:
        for i, fname in enumerate(tqdm(file_names, desc=desc, unit="file")):
            results[i] = _extract_file(i, fname, *args, mission, strict)
    else:
        if backend == "thread":
            executor = ThreadPoolExecutor(max_workers=n_workers)
            worker_mission = mission  # one shared session
        else:
            executor = ProcessPoolExecutor(
                max_workers=n_workers, mp_context=mp.get_context("spawn")
            )
            worker_mission = None  # not picklable; workers build their own
        with executor:
            futures = [
                executor.submit(_extract_file, i, fname, *args, worker_mission, strict)
                for i, fname in enumerate(file_names)
            ]
            for fut in tqdm(as_completed(futures), total=n_files, desc=desc, unit="file"):
                res = fut.result()
                results[res["idx"]] = res

    # Stable time sort: unknown times last, input order as tie-breaker
    mid_jd = np.array([r["mid_jd"] for r in results], dtype=float)
    order = sorted(
        range(n_files),
        key=lambda i: (np.isnan(mid_jd[i]), 0.0 if np.isnan(mid_jd[i]) else mid_jd[i], i),
    )

    n_types = 1 if lite else 4
    planes = ("flx",) if lite else ("flx", "err", "dq", "dumo")
    cutout_cubes, metadatas, gwcsss = [], [], []
    for t in range(n_targets):
        ny, nx = sizes[t]
        cube = np.full((n_types, n_files, ny, nx), np.nan)
        rows, gwcss = [], []
        for k, i in enumerate(order):
            entry = results[i]["targets"][t]
            for p, name in enumerate(planes):
                if entry[name] is not None:
                    cube[p, k] = entry[name]
            rows.append(entry["row"])
            gwcss.append(entry["gwcs"])
        cutout_cubes.append(cube)
        metadatas.append(pd.DataFrame(rows))
        gwcsss.append(gwcss)

    return cutout_cubes, metadatas, gwcsss
