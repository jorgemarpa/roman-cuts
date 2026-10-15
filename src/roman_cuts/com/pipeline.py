"""
roman_cutouts.pipeline
~~~~~~~~~~~~~~~~~~~~~~
High-level pipeline and ASDF serialisation helpers.
"""
import os
import asdf
import logging
import numpy as np
from astropy.time import Time

from .query import query_mast_files
from .cutout import make_cutout

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


def save_to_asdf(cutout_cube, metadata, gwcss, filename):
    """
    Serialise a cutout cube, metadata, and sliced GWCS list to an ASDF file.

    Parameters
    ----------
    cutout_cube : ndarray
        Shape ``(N, ny, nx)``.
    metadata : `~pandas.DataFrame`
        Per-exposure metadata as returned by `make_cutout`.
    gwcss : list
        Per-cutout sliced GWCS objects.
    filename : str
        Output file path (``*.asdf``).
    """
    meta_dict = {col: metadata[col].tolist() for col in metadata.columns}
    meta_serial = prepare_dict_for_asdf(meta_dict)

    tree = {
        "roman": {
            "meta": meta_serial,
            "wcs": gwcss,
            "data": cutout_cube,
            "time": np.array([t.jd for t in metadata["mid_time"]], dtype=float),
            "quality": np.zeros(len(metadata), dtype=np.int16),
        }
    }
    asdf.AsdfFile(tree).write_to(filename)
    logger.info("Wrote %s", filename)


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
    ra, dec : float
        ICRS coordinates of the cutout centre (degrees).
    program_id : int
        Roman program ID.
    detector : str or int
        WFI detector (``'WFI05'``, ``5``, ``'05'`` …).
    filter_name : str
        Optical element (``'F146'``, ``'f106'`` …).
    cutout_size : int or tuple of int
        Cutout size in pixels (square or ``(ny, nx)``).
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

    Returns
    -------
    cutout_cube : ndarray or None
    metadata : `~pandas.DataFrame` or None
    gwcss : list or None
        All three are *None* when MAST returns no products.
    """
    products = query_mast_files(
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
    cutout_cube, metadata, gwcss = make_cutout(
        file_names, ra, dec, cutout_size, mission=mission, verbose=verbose,
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

        out_file = os.path.join(
            output_path,
            f"roman_cutout_p{program_id:05d}{hier}"
            f"_{det_str}_{filter_name.upper()}"
            f"_ra{ra:.5f}_dec{dec:.5f}.asdf",
        )
        save_to_asdf(cutout_cube, metadata, gwcss, out_file)

    return cutout_cube, metadata, gwcss
