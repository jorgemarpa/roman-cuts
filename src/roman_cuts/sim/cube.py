import bz2
import functools
import json
import os
from typing import Any, Optional, Tuple

import asdf
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from tqdm import tqdm

from . import PACKAGEDIR, log
from .utils import extract_all_WCS, extract_average_WCS

RMIN = 0
RMAX = 4088
CMIN = 0
CMAX = 4088


class RomanCuts:
    """
    A class to create cutouts from Roman WFI simulated images created by TRExS group
    using the `RImTimSim` package.

    The class provides access to:
        - Per frame WCS
        - Season average WCS
        - Cutout cubes (ntime, npix, npix) from the simulated FFI stack
        - Save cubes to disk as ASDF
    """

    def __init__(
        self,
        field: int,
        sca: int,
        filter: str = "F146",
        file_list: list = [],
        file_format: str = "fits",
        file_version: str = "1.0",
    ):
        """
        Initializes the class with field, scs, filter, and file_list.

        Parameters
        ----------
        field : int
            The field number.
        sca : int
            The instrument Sensor Chip Assembly number.
        filter : str, optional
            The filter string (e.g., "F146"). Default is "F146".
        file_list : list, optional
            A list of file paths. Default is an empty list.

        """
        self.field = field
        self.sca = sca
        self.filter = filter
        self.file_format_in = file_format
        self.file_version = file_version

        if len(file_list) == 0:
            raise ValueError("Please provide a list of FFI files in `file_list`")
        if not isinstance(file_list, (list, np.ndarray)):
            file_list = np.sort([file_list])
        self.file_list = file_list

        self._check_input_files()

        log.info("Getting 1d arrays data...")
        self._get_arrays()
        log.info("Getting metadata")
        self._get_metadata()
        self.nt = self.metadata["NTIMES"]
        if self.file_format_in == "asdf":
            self.image_size = self.metadata["IMGSIZE"]
        else:
            self.image_size = [RMAX, CMAX]

    def __repr__(self):
        return f"Roman WFI Field {self.field} SCA {self.sca} Filter {self.filter} Frames {self.nt}"

    def _check_input_files(self):
        if self.file_format_in in ["asdf", "ASDF"]:
            if len(self.file_list) == 1:
                self._check_file_list()
            else:
                raise NotImplementedError
        elif self.file_format_in in ["fits", "FITS"]:
            if len(self.file_list) == 1:
                self._check_file_list()
            else:
                self._check_file_list()
        else:
            raise ValueError(
                f"File format {self.file_format_in} not supported. "
                "Please use 'fits' or 'asdf'."
            )

    def _check_file_list(self):
        """
        Hidden method to check that all files in `file_list` exist and are of
        Field/SCA/Filter.
        """
        if not any([os.path.isfile(x) for x in self.file_list]):
            raise ValueError("One of files in `file_list` do not exist in")

        field, sca, filter = [], [], []
        for f in self.file_list:
            if self.file_format_in == "fits":
                hdr = fits.getheader(f)
                sca.append(hdr["DETECTOR"])
                filter.append(hdr["FILTER"])
            if self.file_format_in == "asdf":
                self.asdf_model = asdf.open(f, lazy_tree=True, lazy_load=True)
                sca.append(self.asdf_model["roman"]["meta"]["DETECTOR"])
                field.append(self.asdf_model["roman"]["meta"]["FIELD"])
                filter.append(self.asdf_model["roman"]["meta"]["FILTER"])
                self.mmap_flux = self.asdf_model.tree["roman"]["data"]["flux"]
                self.mmap_flux_err = self.asdf_model.tree["roman"]["data"]["flux_err"]

        if len(set(field)) > 1:
            raise ValueError("File list contains more than one field")
        if len(set(sca)) > 1:
            raise ValueError("File list contains more than one detector")
        if len(set(filter)) > 1:
            raise ValueError("File list contains more than one filter")
        return

    def get_average_wcs(self):
        """
        Computes an average WCS from all available frames
        """
        dir = f"{PACKAGEDIR}/data/wcs/"
        filename = f"{dir}Roman_WFI_wcs_field{self.field:03}_sca{self.sca:02}_{self.filter}_v{self.file_version}.json.bz2"
        if not os.path.isfile(filename):
            self.wcs = extract_average_WCS(self.file_list)
            wcs_dict = {k: v for k, v in self.wcs.to_header().items()}
            os.makedirs(dir, exist_ok=True)
            with bz2.open(filename, "wt", encoding="utf-8") as f:
                f.write(json.dumps(wcs_dict))
        else:
            with bz2.open(filename, "rt", encoding="utf-8") as f:
                loaded_dict = json.load(f)
            self.wcs = WCS(loaded_dict, relax=True)
        return

    def get_all_wcs(self):
        """
        Extracts WCS information from all FFI files.
        """
        dir = f"{PACKAGEDIR}/data/wcs/"
        filename = f"{dir}Roman_WFI_wcss_field{self.field:03}_sca{self.sca:02}_{self.filter}_v{self.file_version}.json.bz2"
        if not os.path.isfile(filename):
            wcss_df = extract_all_WCS(self.file_list)
            wcss_df.to_json(filename, orient="index", compression="bz2")
        else:
            wcss_df = pd.read_json(filename, orient="index", compression="bz2")
        self.wcss = [
            WCS(wcs_dict, relax=True)
            for key, wcs_dict in wcss_df.to_dict(orient="index").items()
            if key in self.exposureno
        ]
        return

    def make_cutout(
        self,
        radec: Tuple = (None, None),
        rowcol: Tuple[int, int] = (0, 0),
        size: Tuple[int, int] = (15, 15),
        dithered: bool = False,
    ):
        """
        Creates a cutout from the data.

        Parameters
        ----------
        radec : tuple of floats or None, optional
            Right ascension and declination coordinates (ra, dec).
        rowcol : tuple of ints or None, optional
            Row and column pixel coordinates (row, col).
        size : tuple of ints, optional
            Size of the cutout in pixels (rows, columns). Default is (15, 15).
        dithered : bool, optional
            If True, recentres each frame using per-frame WCS to track the
            dithered pointing.  Default is False.
        """
        if not hasattr(self, "wcs"):
            self.get_average_wcs()
        if not hasattr(self, "wcss"):
            self.get_all_wcs()
        self.dithered = dithered
        if (
            radec != (None, None)
            and isinstance(radec[0], float)
            and isinstance(radec[1], float)
        ):
            self.ra = radec[0]
            self.dec = radec[1]
            if dithered:
                log.info(
                    "Using Ra, Dec coordinates and WCS per frame to center the cutout"
                )
                row, col = np.vstack(
                    [x.all_world2pix(self.ra, self.dec, 0) for x in self.wcss]
                ).T
                self.target_pixel = np.array([row, col]).T
                row = np.round(row).astype(int)
                col = np.round(col).astype(int)
            else:
                log.info(
                    "Using Ra, Dec coordinates and average WCS to center the cutout"
                )
                row, col = self.wcs.all_world2pix(self.ra, self.dec, 0)
                self.target_pixel = np.array([row, col])
                row = np.array([int(np.round(row))])
                col = np.array([int(np.round(col))])
        elif isinstance(rowcol[0], int) and isinstance(rowcol[1], int):
            row, col = rowcol[0], rowcol[1]
            self.ra, self.dec = self.wcs.all_pix2world(row, col, 0)
        else:
            raise ValueError("Please provide valid `radec` or `rowcol` values")

        log.info("Getting 3d data...")
        if dithered:
            center = tuple([(a, b) for a, b in np.vstack([row, col]).T])
            self._get_cutout_cube_dithered(center=center, size=size)
        else:
            origin = (int(row - size[0] / 2), int(col - size[1] / 2))
            self._get_cutout_cube_static(size=size, origin=origin)

        self._get_metadata()
        return

    @functools.lru_cache(maxsize=6)
    def _get_cutout_cube_static(
        self, size: Tuple[int, int] = (15, 15), origin: Tuple[int, int] = (0, 0)
    ):
        """
        Extracts a static cutout cube from the FFI files. Does not account for
        dithered observations; the cutout is fixed to the pixel grid.
        """
        rmin = RMIN + origin[0]
        cmin = CMIN + origin[1]
        rmax = rmin + size[0]
        cmax = cmin + size[1]

        if (
            (rmin < self.row_min_data)
            | (cmin < self.column_min_data)
            | (rmin > self.row_max_data)
            | (cmin > self.column_max_data)
        ):
            raise ValueError("`cutout_origin` must be within the image.")

        if (rmax > self.row_max_data) | (cmax > self.column_max_data):
            raise ValueError("Cutout exceeds image limits, reducing size.")

        if self.file_format_in == "fits":
            flux = []
            flux_err = []
            for f in tqdm(self.file_list, desc="Extracting cutout"):
                aux = fits.open(f)
                flux.append(aux[0].data[rmin:rmax, cmin:cmax])
                flux_err.append(aux[1].data[rmin:rmax, cmin:cmax])
                aux.close()
        elif self.file_format_in == "asdf":
            rmin -= self.row_min_data
            rmax -= self.row_min_data
            cmin -= self.column_min_data
            cmax -= self.column_min_data
            row_range = np.arange(rmin, rmax)
            col_range = np.arange(cmin, cmax)
            flux = self.mmap_flux[:][:, row_range[row_range >= 0]][
                :, :, col_range[col_range >= 0]
            ]
            flux_err = self.mmap_flux_err[:][:, row_range[row_range >= 0]][
                :, :, col_range[col_range >= 0]
            ]
        else:
            raise ValueError("File format not supported")

        self.flux = np.array(flux)
        self.flux_err = np.array(flux_err)
        self.row = np.arange(rmin, rmax)
        self.row = self.row[self.row >= 0]  # type: ignore
        self.row += self.row_min_data
        self.column = np.arange(cmin, cmax)
        self.column = self.column[self.column >= 0]  # type: ignore
        self.column += self.column_min_data
        self.target_pixel = np.array(
            [
                (self.row.min() + self.row.max()) / 2,
                (self.column.min() + self.column.max()) / 2,
            ]
        )
        return

    @functools.lru_cache(maxsize=6)
    def _get_cutout_cube_dithered(self, center: Any, size: Tuple[int, int] = (15, 15)):
        """
        Extracts a dithered cutout cube: the cutout centre shifts per-frame
        to track celestial coordinates through dithered pointings.
        """
        if not isinstance(center, (np.ndarray)):
            center = np.array(center)
        if len(center.shape) != 2:
            raise ValueError("`center` must be a 2D array with shape (nt, 2) or (1, 2)")
        if center.shape[0] != self.nt:
            raise ValueError(
                "The number of rows in `center` must match the number of files in `file_list`"
            )
        if center.shape[0] == 1:
            log.info("Using the same center for all frames, dithering not accounted.")
            center = np.tile(center, (len(self.file_list), 1))

        row0 = center[:, 0] - int(size[0] / 2)
        col0 = center[:, 1] - int(size[1] / 2)
        rmin = RMIN + row0
        cmin = CMIN + col0
        rmax = rmin + size[0]
        cmax = cmin + size[1]

        if (
            (center[:, 0] - int(size[0]) < self.row_min_data).any()
            | (center[:, 1] - int(size[1]) < self.column_min_data).any()
            | (center[:, 0] + int(size[0]) >= self.row_max_data).any()
            | (center[:, 1] + int(size[1]) >= self.column_max_data).any()
        ):
            log.warning(
                "Cutout out of CCD limits. This is due to the dithered observations"
                " and the size of the cutout. Some pixels will have nan values."
            )

        if self.file_format_in == "fits":
            flux = np.zeros((self.nt, size[0], size[1])) * np.nan
            flux_err = np.zeros((self.nt, size[0], size[1])) * np.nan

            for i, f in tqdm(
                enumerate(self.file_list),
                total=len(self.file_list),
                desc="Extracting cutout",
            ):
                aux = fits.open(f)
                img_rmax = self.row_max_data if rmax[i] > self.row_max_data else rmax[i]
                img_cmax = self.column_max_data if cmax[i] > self.column_max_data else cmax[i]
                img_rmin = self.row_min_data if rmin[i] < self.row_min_data else rmin[i]
                img_cmin = self.column_min_data if cmin[i] < self.column_min_data else cmin[i]

                cutout_rmin = int(img_rmin - row0[i])
                cutout_rmax = int(img_rmax - row0[i])
                cutout_cmin = int(img_cmin - col0[i])
                cutout_cmax = int(img_cmax - col0[i])

                if img_rmin >= img_rmax or img_cmin >= img_cmax:
                    log.warning(
                        f"Cutout for frame {i} is completely out of CCD limits. Skipping."
                    )
                    continue

                flux[i, cutout_rmin:cutout_rmax, cutout_cmin:cutout_cmax] = aux[0].data[
                    img_rmin:img_rmax, img_cmin:img_cmax
                ]
                flux_err[i, cutout_rmin:cutout_rmax, cutout_cmin:cutout_cmax] = aux[1].data[
                    img_rmin:img_rmax, img_cmin:img_cmax
                ]
                aux.close()

        elif self.file_format_in == "asdf":
            rmin -= self.row_min_data
            rmax -= self.row_min_data
            cmin -= self.column_min_data
            cmax -= self.column_min_data
            flux = []  # type: ignore
            flux_err = []  # type: ignore

            for i in range(self.nt):
                row_range = np.arange(rmin[i], rmax[i])
                col_range = np.arange(cmin[i], cmax[i])
                row_in_mask = (row_range >= 0) & (row_range < self.image_size[0])
                col_in_mask = (col_range >= 0) & (col_range < self.image_size[1])
                mask = row_in_mask[:, None] * col_in_mask[None, :]
                aux = np.zeros((size[0], size[1])) * np.nan
                aux[np.where(mask)] = self.mmap_flux[i][row_range[row_in_mask]][
                    :, col_range[col_in_mask]
                ].ravel()
                flux.append(aux)  # type: ignore
                aux = np.zeros((size[0], size[1])) * np.nan
                aux[np.where(mask)] = self.mmap_flux_err[i][row_range[row_in_mask]][
                    :, col_range[col_in_mask]
                ].ravel()
                flux_err.append(aux)  # type: ignore
        else:
            raise ValueError("File format not supported")

        self.flux = np.array(flux)
        self.flux_err = np.array(flux_err)
        self.row = np.vstack([np.arange(rn, rx) for rn, rx in zip(rmin, rmax)])  # type: ignore
        self.row += self.row_min_data
        self.column = np.vstack([np.arange(cn, cx) for cn, cx in zip(cmin, cmax)])  # type: ignore
        self.column += self.column_min_data
        if np.isnan(self.flux).all():
            log.warning("All pixels in the cutout are out of CCD limits.")
        return

    def _get_arrays(self):
        if self.file_format_in == "fits":
            self._fits_arrays()
        elif self.file_format_in == "asdf":
            self._asdf_arrays()
        return

    def _asdf_arrays(self):
        """
        Extracts time, exposureno, and quality arrays from the ASDF file.
        """
        self.time = self.asdf_model["roman"]["data"]["time"].copy()
        self.exposureno = self.asdf_model["roman"]["data"]["exposureno"].copy()
        self.quality = self.asdf_model["roman"]["data"]["quality"].copy()
        self.row_min_data = self.asdf_model["roman"]["data"]["row"]
        self.column_min_data = self.asdf_model["roman"]["data"]["column"]
        self.row_max_data = (
            self.asdf_model["roman"]["data"]["row"]
            + self.asdf_model["roman"]["meta"]["IMGSIZE"][0]
        )
        self.column_max_data = (
            self.asdf_model["roman"]["data"]["column"]
            + self.asdf_model["roman"]["meta"]["IMGSIZE"][1]
        )
        return

    def _fits_arrays(self):
        """
        Extracts time, exposureno, and quality arrays from the FFI files.
        """
        time, exposureno, quality = [], [], []
        for k, f in enumerate(self.file_list):
            hdu = fits.getheader(f)
            time.append((hdu["TSTART"] + hdu["TEND"]) / 2.0)
            exposureno.append(int(os.path.basename(f).split("_")[8]))
            quality.append(0)
        self.time = np.array(time)
        self.exposureno = np.array(exposureno)
        self.quality = np.array(quality)
        self.row_min_data = RMIN
        self.column_min_data = CMIN
        self.row_max_data = hdu["NROW"]
        self.column_max_data = hdu["NCOL"]
        return

    def _get_metadata(self):
        if self.file_format_in == "fits":
            self._fits_metadata()
        elif self.file_format_in == "asdf":
            self._asdf_metadata()
        return

    def _asdf_metadata(self):
        self.metadata = self.asdf_model["roman"]["meta"].copy()
        return

    def _fits_metadata(self):
        hdus = fits.getheader(self.file_list[0])
        hduf = fits.getheader(self.file_list[-1])
        self.metadata = {
            "MISSION": "Roman-Sim",
            "TELESCOP": "Roman",
            "CREATOR": "TRExS-roman-cuts",
            "SOFTWARE": hdus["SOFTWARE"],
            "FILEVER": self.file_version,
            "RADESYS": hdus["RADESYS"],
            "EQUINOX": hdus["EQUINOX"],
            "FILTER": hdus["FILTER"],
            "FIELD": int(os.path.basename(self.file_list[0]).split("_")[5][-2:]),
            "SCA": int(os.path.basename(self.file_list[0]).split("_")[4][-2:]),
            "DETECTOR": hdus["DETECTOR"],
            "EXPOSURE": hdus["EXPOSURE"],
            "READMODE": os.path.basename(self.file_list[0]).split("_")[6],
            "TSTART": hdus["TSTART"],
            "TEND": hduf["TEND"],
            "RA_CEN": float(self.ra) if hasattr(self, "ra") else None,
            "DEC_CEN": float(self.dec) if hasattr(self, "dec") else None,
            "DITHERED": self.dithered if hasattr(self, "dithered") else None,
            "NTIMES": len(self.file_list),
            "IMGSIZE": self.flux.shape[1:] if hasattr(self, "flux") else None,
        }
        return

    def save_cutout(self, output: Optional[str] = None, format: str = "asdf"):
        """
        Saves the cutout to a file.

        Parameters
        ----------
        output : str, optional
            Output file path.  A default filename is generated when *None*.
        format : str, optional
            File format (``"asdf"``).  Default is ``"asdf"``.
        """
        if output is None:
            cutout_str = f"{self.ra:.4f}_{self.dec:.4f}_s{self.flux.shape[1]}x{self.flux.shape[2]}"
            output = f"./roman_cutout_field{self.metadata['FIELD']:02}_{self.metadata['DETECTOR']:02}_{cutout_str}.{format}"
            log.info(f"Saving data to {output}")

        if isinstance(output, str) and not output.endswith(format):
            raise ValueError(
                "Use a valid and matching extension in `output` and `format`"
            )
        if self.dithered:
            save_row = self.row[:, 0]
            save_col = self.column[:, 0]
        else:
            save_row = self.row[0]
            save_col = self.column[0]

        if format in ["asdf", "ASDF"]:
            wcs = self.wcss if hasattr(self, "wcss") else self.wcs
            tree = {
                "roman": {
                    "meta": self.metadata,
                    "wcs": wcs,
                    "data": {
                        "flux": self.flux,
                        "flux_err": self.flux_err,
                        "time": self.time,
                        "exposureno": self.exposureno,
                        "quality": self.quality,
                        "row": save_row,
                        "column": save_col,
                    },
                }
            }
            asdf.AsdfFile(tree).write_to(output)
        elif format in ["fits", "FITS"]:
            raise NotImplementedError
        else:
            raise ValueError("Provide a valid format [FITS or asdf]")
