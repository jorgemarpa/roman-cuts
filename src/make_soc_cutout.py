#!/usr/bin/env python3
import argparse
import logging
import multiprocessing
import os
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob
from typing import Any, Dict, Tuple

import asdf
import numpy as np
import pandas as pd
import s3fs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.time import Time
from astropy.wcs import WCS
from scipy import sparse
from tqdm import tqdm

asdf_dir_uri = "s3://stpubdata/roman/nexus/soc_simulations/r00340/l2"

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# Metadata retention keys provided in your specification
keep_keys_cutout = [
    "input_position_original",
    "origin_original",
    "input_position_cutout",
    "bbox_original",
]
keep_keys_exp = [
    "effective_exposure_time",
    "exposure_time",
    "frame_time",
    "ma_table_number",
    "start_time",
    "end_time",
    "mid_time",
]
keep_keys_obs = ["exposure", "observation", "observation_id", "visit"]


def get_exp_cutout(
    file_name: str,
    sky_coord: SkyCoord,
    cutout_size: int = 25,
    lite: bool = True,
    local: bool = False,
) -> Tuple[np.ndarray, WCS, Dict[str, Any]]:
    """
    Extracts an image cutout from a single file using the Roman Data Model structure.
    """
    # import roman_datamodels as rdm
    from asdf.exceptions import AsdfWarning

    warnings.filterwarnings("ignore", category=AsdfWarning)

    # Init file system inside the worker function so its multiprocessing safe
    if local:
        file_context = open(file_name, "rb")
    else:
        import s3fs

        fs = s3fs.S3FileSystem(anon=True)
        file_context = fs.open(file_name, "rb")

    # Using the specified filesystem stream handler
    with file_context as f:
        dm = asdf.open(f, lazy_load=True)

        # Find pixel coordinates from the input sky coordinates using the GWCS object
        gwcs = dm["roman"]["meta"]["wcs"]
        # pix, wcs_pix = get_center_pixel(gwcs, sky_coord.ra, sky_coord.dec)
        pix = (2553.652, 2249.9031)
        wcs_pix = WCS(gwcs.to_fits_sip())

        # Do cutout for data array
        cutout_data = Cutout2D(
            dm["roman"]["data"],
            position=pix,
            wcs=wcs_pix,
            size=cutout_size,
            mode="partial",
        )
        flux = cutout_data.data

        if isinstance(dm["roman"]["dq"], dict):
            dq = sparse.csr_matrix(
                (
                    dm["roman"]["dq"]["data"],
                    dm["roman"]["dq"]["indices"],
                    dm["roman"]["dq"]["indptr"],
                ),
                shape=tuple(dm["roman"]["dq"]["shape"]),
            ).toarray()
        else:
            dq = dm["roman"]["dq"]

        if not lite:
            cutout_err = Cutout2D(
                dm["roman"]["err"], position=pix, size=cutout_size, mode="partial"
            )
            cutout_dq = Cutout2D(dq, position=pix, size=cutout_size, mode="partial")
            # Stack components into a [3, rows, cols] matrix if lite=False
            flux = np.array([flux, cutout_err.data, cutout_dq.data])

        # Safely extract and merge specified dictionary metadata metadata
        metadata = {
            k: cutout_data.__dict__[k]
            for k in keep_keys_cutout
            if k in cutout_data.__dict__
        }
        metadata |= {
            k: dm["roman"]["meta"]["exposure"][k]
            for k in keep_keys_exp
            if k in dm["roman"]["meta"]["exposure"]
        }
        metadata |= {
            k: dm["roman"]["meta"]["observation"][k]
            for k in keep_keys_obs
            if k in dm["roman"]["meta"]["observation"]
        }

    return flux, cutout_data.wcs, metadata


def worker_wrapper(
    file_path: str,
    target_coord: SkyCoord,
    cutout_size: int,
    lite: bool,
    local: bool,
) -> Dict[str, Any]:
    """
    Wrapper function designed to isolate error boundaries for the parallel processing pool.
    """
    # try:
    if True:
        flux, wcs, metadata = get_exp_cutout(
            file_path, target_coord, cutout_size, lite, local
        )
        return {
            "status": "success",
            "flux": flux,
            "metadata": metadata,
            "wcs": wcs,
            "file": os.path.basename(file_path),
        }
    # except Exception as e:
    #     return {
    #         "status": "failed",
    #         "file": os.path.basename(file_path),
    #         "error": str(e)
    #     }

def prepare_dict_for_asdf(input_dict: Dict[str, Any]) -> Dict[str, np.ndarray]:
    """
    Converts dictionary values into ASDF-serializable NumPy arrays based on type.

    - Tuples are converted to ndarrays.
    - Astropy Time objects are converted to arrays of JD floats.
    - Scalar numbers/strings are converted to standard NumPy scalar arrays.
    """
    asdf_ready_dict = {}

    for key, value in input_dict.items():
        # 1. Handle Astropy Time objects
        if isinstance(value[0], Time):
            asdf_ready_dict[key] = np.array([x.jd for x in value], dtype=float)

        # 2. Handle Tuples (like bounding boxes or shapes)
        elif isinstance(value[0], (tuple, list)):
            asdf_ready_dict[key] = np.array([list(x) for x in value])

        # 4. Handle numerical scalars or strings
        elif isinstance(value[0], str):
            asdf_ready_dict[key] = np.asarray(value).astype(str)
        else:
            # np.asarray handles standard ints, floats, and strings smoothly
            asdf_ready_dict[key] = np.asarray(value)

    return asdf_ready_dict


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parallel extraction of Roman Space Telescope data model image cutouts."
    )
    parser.add_argument(
        "-ra",
        "--ra",
        type=float,
        required=True,
        help="Right Ascension of target center (deg).",
    )
    parser.add_argument(
        "-dec",
        "--dec",
        type=float,
        required=True,
        help="Declination of target center (deg).",
    )
    parser.add_argument(
        "-fld", "--field", type=int, required=True, help="Field number."
    )
    parser.add_argument(
        "-sca", "--sca", type=int, required=True, help="WFI SCA number."
    )
    parser.add_argument("-f", "--filter", type=str, default="F146", help="WFI filter .")
    parser.add_argument(
        "-s",
        "--cutout_size",
        type=int,
        default=25,
        help="Cutout dimension size in pixels. Default is 25.",
    )
    parser.add_argument(
        "-o", "--output_path", default="./", help="Prefix name for saving outputs."
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=4,
        help="Number of concurrent processes. Default is 4.",
    )
    parser.add_argument(
        "-tf",
        "--total_frames",
        type=int,
        default=100,
        help="Total number of frames to extract.",
    )
    parser.add_argument(
        "--lite",
        action="store_true",
        help="If flag is set, only returns flux array instead of [flux, err, dq].",
    )
    parser.add_argument(
        "--local", action="store_true", help="Load files from local directory."
    )

    args = parser.parse_args()

    if args.local:
        file_list = sorted(glob("/home/jorgemarpa/trexs/data/soc_sim/lite3/*.asdf"))
    else:
        fs = s3fs.S3FileSystem(anon=True)
        file_list = fs.glob(
            f"{asdf_dir_uri}/r*_{args.field:04}_wfi{args.sca:02}_{args.filter.lower()}_cal.asdf"
        )[: args.total_frames]
    total_files = len(file_list)
    if total_files == 0:
        raise FileNotFoundError("None files found for field/sca/filter.")
    logging.info(
        f"Total exposures ({'local' if args.local else 's3'}) for Field {args.field} SCA {args.sca} in filter {args.filter}: {total_files}"
    )

    target_coord = SkyCoord(ra=args.ra, dec=args.dec, unit=(u.deg, u.deg), frame="icrs")
    logging.info(
        f"Querying image cutout of size {args.cutout_size} with coordinates {target_coord}"
    )
    logging.info(
        f"Processing {total_files} files using {args.jobs} workers. Mode Lite = {args.lite}"
    )

    # Track results dynamically using lists to accommodate variable dictionary or matrix payloads
    results_list = [None] * total_files

    # Force the process pool executor context to use 'spawn'
    ctx = multiprocessing.get_context("spawn")

    # Execute workloads in parallel across process pools
    with ProcessPoolExecutor(max_workers=args.jobs, mp_context=ctx) as executor:
        future_to_index = {
            executor.submit(
                worker_wrapper,
                path,
                target_coord,
                args.cutout_size,
                args.lite,
                args.local,
            ): idx
            for idx, path in enumerate(file_list)
        }

        for future in tqdm(
            as_completed(future_to_index),
            total=total_files,
            desc="Extracting Cutouts",
            unit="frame",
        ):
            idx = future_to_index[future]
            result = future.result()

            if result["status"] == "success":
                results_list[idx] = {
                    "status": "success",
                    "data": result["flux"],
                    "metadata": result["metadata"],
                    "file": result["file"],
                    "wcs": result["wcs"],
                }
            else:
                logging.error(
                    f"Error handling entry {idx} [{result['file']}]: {result['error']}"
                )
                results_list[idx] = {
                    "status": "failed",
                    "file": result["file"],
                    "error": result["error"],
                }

    # Filter out failures and build separate structures
    clean_flux_list = []
    metadata_rows = []
    wcs_list = []

    for item in results_list:
        if item["status"] == "success":
            flat_meta = item["metadata"].copy()
            flat_meta["file"] = item["file"]
            clean_flux_list.append(item["data"])
            metadata_rows.append(flat_meta)
            wcs_list.append(item["wcs"])

    # Cast to final analytical data types
    if len(clean_flux_list) == 0:
        logging.error("No exposures were successfully processed.")

    # Converts Python lists of arrays directly into a single ndarray block
    flux_ndarray = np.array(clean_flux_list)
    metadata_df = pd.DataFrame(metadata_rows)

    # sorting data by visits
    idx_sort = np.argsort(metadata_df["visit"].values)
    flux_ndarray = flux_ndarray[idx_sort]
    metadata_df = metadata_df.sort_values("visit").reset_index(drop=True)

    # serialize df to save it in the asdf
    metadata_dict = {col: metadata_df[col].to_numpy() for col in metadata_df.columns}
    metadata_dict["df_index"] = metadata_df.index.to_numpy()

    metadata_serial = prepare_dict_for_asdf(metadata_dict)

    # Saving into a ASDF file
    tree = {
        "roman": {
            "meta": metadata_serial,
            "wcs": wcs_list,
            "data": flux_ndarray[:, 0],
            "err": flux_ndarray[:, 1],
            "dq": flux_ndarray[:, 2],
            "time": np.array([x.jd for x in metadata_df["mid_time"]]),
            "exposureno": metadata_df["visit"].values,
            "quality": np.zeros(len(metadata_df)),
        }
    }
    ff = asdf.AsdfFile(tree)
    out_file = (
        f"{args.output_path}/roman_cutout_{args.field:04}_wfi{args.sca:02}_{args.filter.lower()}"
        f"_ra{target_coord.ra.deg:.5f}_dec{target_coord.dec.deg:.5f}.asdf"
    )
    ff.write_to(out_file)

    logging.info("Task processing finished completely.")


if __name__ == "__main__":
    main()
