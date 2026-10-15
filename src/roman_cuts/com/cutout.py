"""
roman_cutouts.cutout
~~~~~~~~~~~~~~~~~~~~
Pixel cutout extraction from Roman WFI L2 ASDF files using GWCS.
"""
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.time import Time
from astropy.wcs.wcsapi import SlicedLowLevelWCS, HighLevelWCSWrapper

from .io import _open_asdf

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


def make_cutout(file_names, ra, dec, size, mission=None, verbose=False):
    """
    Create a pixel cutout (or cube of cutouts) from Roman WFI L2 ASDF files.

    The GWCS stored in each ASDF file maps (ra, dec) to detector pixel
    coordinates.  The cutout is extracted with `~astropy.nddata.Cutout2D`
    (``mode='partial'``), which fills out-of-bounds pixels with NaN and
    handles coordinates near the detector edge gracefully.

    File access is dispatched by `~roman_cutouts.io._open_asdf`, which
    auto-detects the execution platform (RRN, SMDC, standalone) via
    `~roman_cutouts.io._detect_platform`.  Supply *mission* to force the
    RRN/MAST path from any platform.

    Parameters
    ----------
    file_names : str or list of str
        MAST product filenames (RRN), local/mounted-S3 paths (SMDC), plain
        local paths, or ``s3://`` URIs.  When a list is supplied all files
        must cover the requested sky position.
    ra : float
        Right Ascension of the cutout centre in degrees (ICRS).
    dec : float
        Declination of the cutout centre in degrees (ICRS).
    size : int or tuple of int
        Cutout size in pixels.  An ``int`` gives a square cutout; a
        ``(ny, nx)`` tuple gives a rectangular one.
    mission : `~astroquery.mast.MastMissions`, optional
        Authenticated mission object for RRN/MAST access.
    verbose : bool
        Print per-file pixel-coordinate messages when *True*.

    Returns
    -------
    cutout_cube : ndarray
        Pixel-flux array, shape ``(N, ny, nx)`` for a list or ``(ny, nx)``
        for a single file.
    metadata : `~pandas.DataFrame`
        Per-exposure metadata with columns from exposure, observation, and
        instrument nodes, plus ``mid_time``, ``ra_target``, ``dec_target``,
        and Cutout2D geometry fields (``origin_original``, ``bbox_original``
        …).  One row per successfully extracted cutout.
    gwcss : list of `~astropy.wcs.wcsapi.HighLevelWCSWrapper`
        Per-cutout sliced GWCS objects (``SlicedLowLevelWCS`` wrapped in
        ``HighLevelWCSWrapper``).  ``gwcss[k].world_to_pixel(SkyCoord(...))``
        returns pixel coordinates in the cutout frame.
    """
    single = isinstance(file_names, str)
    if single:
        file_names = [file_names]

    target = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame="icrs")
    if verbose:
        logger.info("Target: RA=%.6f  Dec=%.6f", ra, dec)

    cutouts, meta_rows, gwcss = [], [], []

    for fname in tqdm(file_names, desc="Extracting cutouts", unit="file"):
        with _open_asdf(fname, mission=mission) as af:
            gwcs_obj = af["roman"]["meta"]["wcs"]

            # Convert sky coordinates to detector pixels via GWCS (APE 14)
            x_pix, y_pix = gwcs_obj.world_to_pixel(target)
            if np.isnan(x_pix) or np.isnan(y_pix):
                logger.debug("Skipping %s — target outside image bounds.", fname)
                continue
            if verbose:
                logger.info("  %s  →  pixel (%.2f, %.2f)", fname, x_pix, y_pix)

            # Cutout2D reads the needed slice during __init__;
            # .copy() keeps the data valid after the ASDF file closes
            cutout = Cutout2D(
                af["roman"]["data"],
                position=(float(x_pix), float(y_pix)),
                size=size,
                mode="partial",
                fill_value=np.nan,
            )
            cutouts.append(cutout.data.copy())

            # Build sliced GWCS for the cutout frame
            yslice, xslice = cutout.slices_original
            gwcss.append(
                HighLevelWCSWrapper(SlicedLowLevelWCS(gwcs_obj, (yslice, xslice)))
            )

            # Mid-time: use stored value; fall back to (start + end) / 2
            exp_node = af["roman"]["meta"]["exposure"]
            mid_t = exp_node.get("mid_time")
            if mid_t is None:
                mid_t = Time(
                    (Time(exp_node["start_time"]).jd + Time(exp_node["end_time"]).jd) / 2.0,
                    format="jd",
                )
            if not isinstance(mid_t, Time):
                mid_t = Time(mid_t)

            row = {k: cutout.__dict__[k] for k in _KEEP_CUTOUT if k in cutout.__dict__}
            row |= {k: af["roman"]["meta"]["exposure"][k] for k in _KEEP_EXP
                    if k in af["roman"]["meta"]["exposure"]}
            row |= {k: af["roman"]["meta"]["observation"][k] for k in _KEEP_OBS
                    if k in af["roman"]["meta"]["observation"]}
            row |= {k: af["roman"]["meta"]["instrument"][k] for k in _KEEP_INSTR
                    if k in af["roman"]["meta"]["instrument"]}
            row["mid_time"] = mid_t
            row["ra_target"] = ra
            row["dec_target"] = dec
            meta_rows.append(row)

    cutout_cube = np.array(cutouts)
    metadata = pd.DataFrame(meta_rows)

    if single:
        return cutout_cube[0], metadata.iloc[0], gwcss[0]

    return cutout_cube, metadata, gwcss
