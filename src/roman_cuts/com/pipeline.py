"""
roman_cutouts.pipeline
~~~~~~~~~~~~~~~~~~~~~~
High-level pipeline and ASDF serialisation helpers.
"""
import logging
import os
import re

import asdf
import numpy as np
from astropy.time import Time

from .cutout import make_cutout
from .query import query_mast_files

logger = logging.getLogger(__name__)


def prepare_dict_for_asdf(input_dict):
    """
    Convert a metadata dict to ASDF-serialisable NumPy arrays.

    - `~astropy.time.Time` objects → float JD arrays
    - Tuples / lists → 2-D ndarrays
    - Strings → str arrays
    - Everything else → ``np.asarray``
    """
    out = {}
    for key, value in input_dict.items():
        if isinstance(value[0], Time):
            out[key] = np.array([x.jd for x in value], dtype=float)
        elif isinstance(value[0], (tuple, list)):
            out[key] = np.array([list(x) for x in value])
        elif isinstance(value[0], str):
            out[key] = np.asarray(value).astype(str)
        else:
            out[key] = np.asarray(value)
    return out


def _default_filename(metadata):
    """
    Build ``cutouts_{pid}_{detector}_{filter}_{ra:.6f}_{dec:.6f}.asdf`` in the cwd.

    The program number ``pid`` is the leading ``PPPPP`` of the ``observation_id``
    metadata column (``PPPPPCCAAASSSOOOVVV`` — program, execution plan, pass,
    segment, observation, visit), tolerating an optional leading ``r``.
    """
    def first(col):
        if col not in metadata:
            return "unknown"
        valid = metadata[col].dropna()
        return valid.iloc[0] if len(valid) else "unknown"

    match = re.match(r"^[A-Za-z]?(\d{5})", str(first("observation_id")))
    pid = match.group(1) if match else "unknown"
    return (
        f"cutouts_{pid}_{first('detector')}_{first('optical_element')}"
        f"_{float(first('ra_target')):.6f}_{float(first('dec_target')):.6f}.asdf"
    )


def _save_one(cutout_cube, metadata, gwcss, filename):
    """Write a single target's cube, metadata and GWCS list to *filename*."""
    meta_dict = {col: metadata[col].tolist() for col in metadata.columns}
    meta_dict["df_index"] = metadata.index.to_numpy()
    # meta_serial = prepare_dict_for_asdf(meta_dict)

    if cutout_cube.ndim == 4:
        data = cutout_cube[0]
        has_all = cutout_cube.shape[0] >= 4
        err = cutout_cube[1] if has_all else None
        dq = cutout_cube[2] if has_all else None
        dumo = cutout_cube[3] if has_all else None
    elif cutout_cube.ndim == 3:
        data = cutout_cube
        err = None
        dq = None
        dumo = None
    else:
        raise ValueError(
            f"cutout_cube must be 3-D or 4-D, got {cutout_cube.ndim}-D."
        )

    tree = {
        "roman": {
            "meta": meta_dict,
            "wcs": gwcss,
            "data": data,
            "err": err,
            "dq": dq,
            "dumo": dumo,
            "time": [
                float(x.jd) if isinstance(x, Time) else np.nan
                for x in metadata["mid_time"].values
            ],
            "quality": np.zeros(len(metadata), dtype=np.int16),
        }
    }
    asdf.AsdfFile(tree).write_to(filename)
    logger.info("Wrote %s", filename)


def save_to_asdf(cutout_cube, metadata, gwcss, filename=None):
    """
    Serialise cutout cubes, metadata, and sliced GWCS lists to ASDF files.

    One file is written **per target**.  Pass the lists returned by
    `make_cutout` (``N`` targets → ``N`` files) or the arrays of a single
    target.

    Parameters
    ----------
    cutout_cube : ndarray or list of ndarray
        Per-target cube(s), shape ``(n_types, n_files, ny, nx)``
        (``[data, err, dq, dumo]``, or data only when ``n_types == 1``) or
        ``(n_files, ny, nx)``.  A bare array is treated as one target.
    metadata : `~pandas.DataFrame` or list of `~pandas.DataFrame`
        Per-exposure metadata as returned by `make_cutout`.
    gwcss : list or list of list
        Per-cutout sliced GWCS objects (one list per target when *cutout_cube*
        is a list).
    filename : str or list of str, optional
        Output file path(s) (``*.asdf``).  With ``N > 1`` targets it must be a
        list of length ``N``.  When *None*, each file is written to the
        current directory as
        ``cutouts_{pid}_{detector}_{filter}_{ra:.6f}_{dec:.6f}.asdf``, with
        ``pid`` parsed from the leading ``PPPPP`` of the ``observation_id``
        metadata column (``PPPPPCCAAASSSOOOVVV``), and detector
        (``detector``), filter (``optical_element``) and coordinates
        (``ra_target``, ``dec_target``) also taken from the metadata.
        Missing values are written as ``unknown``.

    Returns
    -------
    list of str
        The paths written, one per target.
    """
    if isinstance(cutout_cube, (list, tuple)):
        cubes, metas, gwcsss = list(cutout_cube), list(metadata), list(gwcss)
    else:
        cubes, metas, gwcsss = [cutout_cube], [metadata], [gwcss]
    n = len(cubes)
    if not (len(metas) == len(gwcsss) == n):
        raise ValueError(
            "cutout_cube, metadata and gwcss must have the same number of targets "
            f"(got {n}, {len(metas)}, {len(gwcsss)})."
        )

    if filename is None:
        filenames = [_default_filename(m) for m in metas]
    else:
        filenames = [filename] if isinstance(filename, (str, os.PathLike)) else list(filename)
        if len(filenames) != n:
            raise ValueError(
                f"filename must have one entry per target ({n}), got {len(filenames)}."
            )

    for cube, meta, gw, fname in zip(cubes, metas, gwcsss, filenames):
        _save_one(cube, meta, gw, fname)
    return [str(f) for f in filenames]


def run_cutout_pipeline(
    ra,
    dec,
    program_id,
    detector,
    filter_name,
    cutout_size,
    file_type="l2",
    execution_plan=None,
    pass_num=None,
    segment=None,
    observation=None,
    visit=None,
    max_files=None,
    output_path=None,
    mission=None,
    verbose=False,
    n_workers=8,
    backend="thread",
):
    """
    End-to-end pipeline: query MAST → extract cutouts → optionally save.

    Combines `query_mast_files`, `make_cutout`, and `save_to_asdf` into a
    single callable that works from an interactive session, an import in
    another script, or the CLI via ``roman-make-cutout``.

    Observation hierarchy parameters (``execution_plan`` … ``visit``) are
    forwarded to `query_mast_files`; see that function's docstring for the
    full naming-convention table.  Specifying ``observation`` is the minimum
    needed to pin a single field.

    Parameters
    ----------
    ra, dec : float or sequence of float
        ICRS coordinates of the cutout centre(s) (degrees).
    program_id : int
        Roman program ID.
    detector : str or int
        WFI detector (``'WFI05'``, ``5``, ``'05'`` …).
    filter_name : str
        Optical element (``'F146'``, ``'f106'`` …).
    cutout_size : int, tuple of int, or list
        Cutout size in pixels (square, ``(ny, nx)``, or one per target).
    file_type : ``'l1'`` | ``'l2'`` | ``'l4'``
        Data level.  Default ``'l2'``.
    execution_plan : int, optional
        Execution plan within the program.
    pass_num : int, optional
        Pass (epoch) within the execution plan.
    segment : int, optional
        Segment within the pass.
    observation : int, optional
        Observation number — primary sky-field identifier.
    visit : int, optional
        Visit (dither position) within the observation.
    max_files : int, optional
        Cap on the number of files sent to `make_cutout`.
    output_path : str, optional
        Directory in which to write the output ASDF file.  *None* skips
        saving.
    mission : `~astroquery.mast.MastMissions`, optional
        Authenticated mission object; auto-created from ``MAST_API_TOKEN``
        when *None*.
    verbose : bool
        Forward verbose flag to `make_cutout`.
    n_workers : int
        Parallel workers used by `make_cutout`.  Default 8.
    backend : ``'thread'`` | ``'process'``
        Parallel backend used by `make_cutout`.  Default ``'thread'``.

    Returns
    -------
    cutout_cubes : list of ndarray or None
        One cube per target, shape ``(4, n_files, ny, nx)``.
    metadatas : list of `~pandas.DataFrame` or None
    gwcsss : list of list or None
        All three are *None* when MAST returns no products.
    """
    products, mission = query_mast_files(
        program_id,
        detector,
        filter_name,
        file_type=file_type,
        execution_plan=execution_plan,
        pass_num=pass_num,
        segment=segment,
        observation=observation,
        visit=visit,
        mission=mission,
    )

    if len(products) == 0:
        logger.error("No products returned from MAST; aborting pipeline.")
        return None, None, None

    file_names = products["filename"].tolist()
    if max_files is not None:
        file_names = file_names[:max_files]
        logger.info("Capped to %d files (max_files=%d).", len(file_names), max_files)

    logger.info("Running make_cutout on %d files.", len(file_names))
    cutout_cubes, metadatas, gwcsss = make_cutout(
        file_names, ra, dec, cutout_size, mission=mission, verbose=verbose,
        n_workers=n_workers, backend=backend,
    )

    if output_path is not None:
        det_str = str(detector).upper()
        if not det_str.startswith("WFI"):
            det_str = f"WFI{int(det_str):02d}"

        hier = ""
        if execution_plan is not None:
            hier += f"_ep{execution_plan:02d}"
        if pass_num is not None:
            hier += f"_pass{pass_num:03d}"
        if segment is not None:
            hier += f"_seg{segment:03d}"
        if observation is not None:
            hier += f"_obs{observation:03d}"
        if visit is not None:
            hier += f"_vis{visit:03d}"

        # One output file per target; RA/Dec in the name keeps them distinct
        filenames = [
            os.path.join(
                output_path,
                f"roman_cutout_p{program_id:05d}{hier}"
                f"_{det_str}_{filter_name.upper()}"
                f"_ra{ra_t:.5f}_dec{dec_t:.5f}.asdf",
            )
            for ra_t, dec_t in zip(np.atleast_1d(ra), np.atleast_1d(dec))
        ]
        save_to_asdf(cutout_cubes, metadatas, gwcsss, filenames)

    return cutout_cubes, metadatas, gwcsss
