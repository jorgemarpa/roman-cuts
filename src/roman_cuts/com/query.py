"""
roman_cutouts.query
~~~~~~~~~~~~~~~~~~~
MAST archive query helpers for Roman WFI data products.
"""
import os
import logging

logger = logging.getLogger(__name__)

# Map user-facing file-type labels to MAST productLevel and file suffix
FILE_TYPE_MAP = {
    "l1": {"product_level": 1, "file_suffix": "_uncal"},
    "l2": {"product_level": 2, "file_suffix": "_cal"},
    "l4": {"product_level": 4, "file_suffix": "_cat"},
}

QUERY_COLUMNS = [
    "program",
    "execution_plan",
    "pass",
    "segment",
    "visit",
    "observation",
    "optical_element",
    "exposure_type",
    "instrument_name",
    "detector",
    "productLevel",
    "product_type",
    "exposure_time",
    "exposure_start_time",
    "exposure_end_time",
    "fileSetName",
]


def query_mast_files(
    program_id,
    detector,
    filter_name,
    file_type="l2",
    execution_plan=None,
    pass_num=None,
    segment=None,
    observation=None,
    visit=None,
    mission=None,
):
    """
    Query the MAST archive for Roman WFI files at a specific sky field.

    Roman file names encode the full observation hierarchy as
    ``r{PPPPP}{CC}{AAA}{SSS}{OOO}{VVV}_{eeee}_{wfiss}_{filter}_{suffix}``:

    ======= ===  ================================================================
    Level   Len  Meaning
    ======= ===  ================================================================
    Program   5  Observing program (``program_id``)
    ExecPlan  2  Execution plan within the program (``execution_plan``)
    Pass      3  Pass within the execution plan (``pass_num``); repeat visits
                 of the same field across epochs share the same pass
    Segment   3  Segment within a pass (``segment``)
    Observ.   3  Observation within a segment (``observation``); **this is the
                 primary sky-field identifier** — all exposures with the same
                 observation number point to the same sky patch
    Visit     3  Visit (dither position) within an observation (``visit``)
    Exposure  4  Individual detector exposure
    ======= ===  ================================================================

    Parameters
    ----------
    program_id : int
        Roman program (proposal) ID, e.g. ``1054``.
    detector : str or int
        WFI detector.  Any of ``'WFI05'``, ``'wfi05'``, ``5``, ``'05'`` are
        accepted and normalised to ``'WFI##'``.
    filter_name : str
        Optical element, e.g. ``'F146'`` or ``'f106'``.  Case-insensitive.
    file_type : ``'l1'``, ``'l2'``, or ``'l4'``
        Data level:

        - ``'l1'`` – uncalibrated ramps (``_uncal.asdf``)
        - ``'l2'`` – calibrated rate images (``_cal.asdf``)  *(default)*
        - ``'l4'`` – source catalogs (``_cat.parquet``)
    execution_plan : int, optional
        Execution plan number within the program (``CC`` in the filename).
    pass_num : int, optional
        Pass number within the execution plan (``AAA``).  Different passes
        of the same observation cover the same sky at different epochs.
    segment : int, optional
        Segment number within the pass (``SSS``).
    observation : int, optional
        Observation number within the segment (``OOO``).  **This is the key
        parameter for selecting a specific field**: all exposures sharing the
        same observation number are co-pointed.
    visit : int, optional
        Visit (dither position) within the observation (``VVV``).  Omit to
        retrieve all dithers at the field; supply to pin a single dither.
    mission : `~astroquery.mast.MastMissions`, optional
        Authenticated mission object.  When *None*, one is created from the
        ``MAST_API_TOKEN`` environment variable.

    Returns
    -------
    products : `~astropy.table.Table`
        MAST product table.  Each row is one file; ``products['filename']``
        can be passed directly to `make_cutout`.
    """
    from astroquery.mast import MastMissions

    file_type = file_type.lower()
    if file_type not in FILE_TYPE_MAP:
        raise ValueError(
            f"file_type must be one of {list(FILE_TYPE_MAP)}, got {file_type!r}"
        )

    # Normalise detector to 'WFI##' (1 → 'WFI01', 'wfi5' → 'WFI05', etc.)
    det_str = str(detector).upper().lstrip("WFI").lstrip("0") or "0"
    detector_norm = f"WFI{int(det_str):02d}"

    filter_norm = filter_name.upper()

    if mission is None:
        mission = MastMissions(mission="roman")
        token = os.environ.get("MAST_API_TOKEN")
        if token:
            mission.login(token=token)

    level_info = FILE_TYPE_MAP[file_type]

    # Build search criteria; only include hierarchy levels when specified
    criteria = {
        "program": program_id,
        "detector": detector_norm,
        "optical_element": filter_norm,
        "productLevel": level_info["product_level"],
    }
    if execution_plan is not None:
        criteria["execution_plan"] = int(execution_plan)
    if pass_num is not None:
        criteria["pass"] = int(pass_num)
    if segment is not None:
        criteria["segment"] = int(segment)
    if observation is not None:
        criteria["observation"] = int(observation)
    if visit is not None:
        criteria["visit"] = int(visit)

    results = mission.query_criteria(**criteria, select_cols=QUERY_COLUMNS)

    if len(results) == 0:
        logger.warning(
            "No observations found for %s",
            ", ".join(f"{k}={v}" for k, v in criteria.items()),
        )
        return results

    logger.info(
        "Found %d observations — program=%s %s %s level=%s",
        len(results), program_id, detector_norm, filter_norm, file_type.upper(),
    )

    product_list = mission.get_unique_product_list(results)
    products = mission.filter_products(
        product_list, file_suffix=level_info["file_suffix"]
    )

    logger.info("Returning %d %s products.", len(products), file_type.upper())
    return products
