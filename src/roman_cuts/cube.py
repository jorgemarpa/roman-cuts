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
                # single ASDF file, this is a data cube
                self._check_file_list()
            else:
                # currently not supporting list of ASDF files
                raise NotImplementedError
                # list of ASDF files, each is a frame
        elif self.file_format_in in ["fits", "FITS"]:
            if len(self.file_list) == 1:
                # currently not supporting data cubes in FITS
                # this only works for single FFI
                self._check_file_list()
                # raise NotImplementedError
            else:
                self._check_file_list()
            # list of ASDF files, each is a frame
        else:
            raise ValueError(
                f"File format {self.file_format_in} not supported. "
                "Please use 'fits' or 'asdf'."
            )

    def _check_file_list(self):
        """
        HIdden method to check that all files in `file_list` exist and are of
        Field/SCA/Filter.
        """

        # check files exist
        if not any([os.path.isfile(x) for x in self.file_list]):
            raise ValueError("One of files in `file_list` do not exist in")

        field, sca, filter = [], [], []
        # check all files are same Field/SCA/Filter
        for f in self.file_list:
            if self.file_format_in == "fits":
                hdr = fits.getheader(f)
                # field.append(hdr["FIELD"])
                sca.append(hdr["DETECTOR"])
                filter.append(hdr["FILTER"])
            if self.file_format_in == "asdf":
                datamodel = asdf.open(f, lazy_tree=True, lazy_load=True)
                sca.append(datamodel["roman"]["meta"]["DETECTOR"])
                field.append(datamodel["roman"]["meta"]["FIELD"])
                filter.append(datamodel["roman"]["meta"]["FILTER"])
                datamodel.close()

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
        # check if wcs is in disk
        dir = f"{PACKAGEDIR}/data/wcs/"
        filename = f"{dir}Roman_WFI_wcs_field{self.field:03}_sca{self.sca:02}_{self.filter}.json.bz2"
        if not os.path.isfile(filename):
            # if not compute a new one and save it to disk
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
        # check if wcs is in disk
        dir = f"{PACKAGEDIR}/data/wcs/"
        filename = f"{dir}Roman_WFI_wcss_field{self.field:03}_sca{self.sca:02}_{self.filter}.json.bz2"
        if not os.path.isfile(filename):
            # if not compute a new one and save it to disk
            wcss_df = extract_all_WCS(self.file_list)
            wcss_df.to_json(filename, orient="index", compression="bz2")
        else:
            # if exist, load from disk
            wcss_df = pd.read_json(filename, orient="index", compression="bz2")
        # convert to list of WCS objects
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
            If None, rowcol is used. Default is (None, None).
        rowcol : tuple of ints or None, optional
            Row and column pixel coordinates (row, col). If None, radec is used.
            Default is (0, 0).
        size : tuple of ints, optional
            Size of the cutout in pixels (rows, columns). Default is (15, 15).
        """
        # check we have the wcs loaded
        if not hasattr(self, "wcs"):
            self.get_average_wcs()
        if not hasattr(self, "wcss"):
            self.get_all_wcs()
        self.dithered = dithered
        # use radec if given
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
                # use the each WCS to get the row col
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

        # if not use the rowcol
        elif isinstance(rowcol[0], int) and isinstance(rowcol[1], int):
            row, col = rowcol[0], rowcol[1]
            self.ra, self.dec = self.wcs.all_pix2world(row, col, 0)
        # raise error if values are not valid
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
        Extracts a static cutout cube from the FFI files. It does not account for
        dithered observations, therefore the cutout is fixed to the pixel grid.

        Parameters
        ----------
        size : tuple of ints, optional
            Size of the cutout in pixels (rows, columns). Default is (15, 15).
        origin : tuple of ints, optional
            Pixel coordinates of the origin (row, col). Default is (0, 0).
        """
        # set starting pixel
        rmin = RMIN + origin[0]
        cmin = CMIN + origin[1]
        # set ending pixels
        rmax = rmin + size[0]
        cmax = cmin + size[1]

        # make sure requested data is in cutout
        if (
            (rmin < self.row_min_data)
            | (cmin < self.column_min_data)
            | (rmin > self.row_max_data)
            | (cmin > self.column_max_data)
        ):
            raise ValueError("`cutout_origin` must be within the image.")

        if (rmax > self.row_max_data) | (cmax > self.column_max_data):
            raise ValueError("Cutout exceeds image limits, reducing size.")

        # get data from FITS assuming is the FFI
        if self.file_format_in == "fits":
            flux = []
            flux_err = []
            for f in tqdm(self.file_list, desc="Extracting cutout"):
                aux = fits.open(f)
                flux.append(aux[0].data[rmin:rmax, cmin:cmax])
                flux_err.append(aux[1].data[rmin:rmax, cmin:cmax])
                aux.close()
        # get data from ASDF, this could be an FFI or a cutout
        elif self.file_format_in == "asdf":
            rmin -= self.row_min_data
            rmax -= self.row_min_data
            cmin -= self.column_min_data
            cmax -= self.column_min_data
            row_range = np.arange(rmin, rmax)
            col_range = np.arange(cmin, cmax)

            cont = asdf.open(self.file_list[0], lazy_tree=True, lazy_load=True)
            flux = cont["roman"]["data"]["flux"][:][:, row_range[row_range >= 0]][
                :, :, col_range[col_range >= 0]
            ]
            flux_err = cont["roman"]["data"]["flux_err"][:][
                :, row_range[row_range >= 0]
            ][:, :, col_range[col_range >= 0]]
        else:
            raise ValueError("File format not supported")

        self.flux = np.array(flux)
        self.flux_err = np.array(flux_err)
        self.row = np.arange(rmin, rmax)
        self.row = self.row[self.row >= 0]
        self.row += self.row_min_data
        self.column = np.arange(cmin, cmax)
        self.column = self.column[self.column >= 0]
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
        Extracts a static cutout cube from the FFI files. The cutout is centered on
        the pixel coordinates equivalent to the celestial coordinates, therefore
        it accounts for dithered observations.

        Parameters
        ----------
        center : ndarray
            2D array of shape (nt, 2) with pixel coordinates of the center (row, col).
            If the shape is (1, 2), the same center is used for all frames.
        size : tuple of ints, optional
            Size of the cutout in pixels (rows, columns). Default is (15, 15).
        """
        if not isinstance(center, (np.ndarray)):
            center = np.array(center)
        # set starting pixel
        if len(center.shape) != 2:
            raise ValueError("`center` must be a 2D array with shape (nt, 2) or (1, 2)")
        if center.shape[0] != self.nt:
            raise ValueError(
                "The number of rows in `center` must match the number of files in `file_list`"
            )
        if center.shape[0] == 1:
            log.info("Using the same center for all frames, dithering not accounted.")
            center = np.tile(center, (len(self.file_list), 1))

        # get range of requested row/column pixels
        row0 = center[:, 0] - int(size[0] / 2)
        col0 = center[:, 1] - int(size[1] / 2)
        rmin = RMIN + row0
        cmin = CMIN + col0
        rmax = rmin + size[0]
        cmax = cmin + size[1]

        # make sure at least 50% of the requested data is on detector
        if (
            (center[:, 0] - int(size[0] / 4) < self.row_min_data).any()
            | (center[:, 1] - int(size[1] / 4) < self.column_min_data).any()
            | (center[:, 0] + int(size[0] / 4) >= self.row_max_data).any()
            | (center[:, 1] + int(size[1] / 4) >= self.column_max_data).any()
        ):
            raise ValueError(
                "Cutout out of CCD limits. This is due to the dithered observations"
                " and the size of the cutout. Please reduce the size or change the center."
            )
        # get data from FITS assuming is the FFI
        if self.file_format_in == "fits":
            flux = []
            flux_err = []
            for i, f in tqdm(
                enumerate(self.file_list),
                total=len(self.file_list),
                desc="Extracting cutout",
            ):
                aux = fits.open(f)
                flux.append(aux[0].data[rmin[i] : rmax[i], cmin[i] : cmax[i]])
                flux_err.append(aux[1].data[rmin[i] : rmax[i], cmin[i] : cmax[i]])
                aux.close()
        # get data from ASDF, this could be an FFI or a cutout
        elif self.file_format_in == "asdf":
            # reset range to array index in the data
            rmin -= self.row_min_data
            rmax -= self.row_min_data
            cmin -= self.column_min_data
            cmax -= self.column_min_data
            flux = []
            flux_err = []

            cont = asdf.open(self.file_list[0], lazy_tree=True, lazy_load=True)
            for i in range(self.nt):
                # find which requested row/column are in data range
                row_range = np.arange(rmin[i], rmax[i])
                col_range = np.arange(cmin[i], cmax[i])
                row_in_mask = (row_range >= 0) & (row_range < self.image_size[0])
                col_in_mask = (col_range >= 0) & (col_range < self.image_size[1])
                mask = row_in_mask[:, None] * col_in_mask[None, :]
                # requested pixels outside the image range are filled in with nans
                # to keep the cutout size consistent.
                # we accept up to 25% nan row/column in each edge
                aux = np.zeros((size[0], size[1])) * np.nan
                aux[np.where(mask)] = cont["roman"]["data"]["flux"][i][
                    row_range[row_in_mask]
                ][:, col_range[col_in_mask]].ravel()
                flux.append(aux)
                aux = np.zeros((size[0], size[1])) * np.nan
                aux[np.where(mask)] = cont["roman"]["data"]["flux_err"][i][
                    row_range[row_in_mask]
                ][:, col_range[col_in_mask]].ravel()
                flux_err.append(aux)
        else:
            raise ValueError("File format not supported")

        self.flux = np.array(flux)
        self.flux_err = np.array(flux_err)
        self.row = np.vstack([np.arange(rn, rx) for rn, rx in zip(rmin, rmax)])
        self.row += self.row_min_data
        self.column = np.vstack([np.arange(cn, cx) for cn, cx in zip(cmin, cmax)])
        self.column += self.column_min_data
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
        with asdf.open(self.file_list[0], lazy_tree=True, lazy_load=True) as cont:
            self.time = cont["roman"]["data"]["time"].copy()
            self.exposureno = cont["roman"]["data"]["exposureno"].copy()
            self.quality = cont["roman"]["data"]["quality"].copy()
            self.row_min_data = cont["roman"]["data"]["row"]
            self.column_min_data = cont["roman"]["data"]["column"]
            self.row_max_data = (
                cont["roman"]["data"]["row"] + cont["roman"]["meta"]["IMGSIZE"][0]
            )
            self.column_max_data = (
                cont["roman"]["data"]["column"] + cont["roman"]["meta"]["IMGSIZE"][1]
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
            # replace these two to corresponding keywords in future simulations
            exposureno.append(int(f.split("_")[-2]))
            quality.append(0)
        self.time = np.array(time)
        self.exposureno = np.array(exposureno)
        self.quality = np.array(quality)
        self.row_min_data = RMIN
        self.column_min_data = CMIN
        self.row_max_data = RMAX
        self.column_max_data = CMAX
        return

    def _get_metadata(self):
        if self.file_format_in == "fits":
            self._fits_metadata()
        elif self.file_format_in == "asdf":
            self._asdf_metadata()
        return

    def _asdf_metadata(self):
        """
        Extracts metadata from the ASDF file.
        """
        with asdf.open(self.file_list[0], lazy_tree=False, lazy_load=True) as cont:
            self.metadata = cont["roman"]["meta"].copy()
        return

    def _fits_metadata(self):
        """
        Extracts metadata from the first FFI file.
        """
        hdus = fits.getheader(self.file_list[0])
        hduf = fits.getheader(self.file_list[-1])
        self.metadata = {
            "MISSION": "Roman-Sim",
            "TELESCOP": "Roman",
            "CREATOR": "TRExS-roman-cuts",
            "SOFTWARE": hdus["SOFTWARE"],
            "RADESYS": hdus["RADESYS"],
            "EQUINOX": hdus["EQUINOX"],
            "FILTER": hdus["FILTER"],
            "FIELD": int(self.file_list[0].split("_")[-5][-2:]),
            "DETECTOR": hdus["DETECTOR"],
            "EXPOSURE": hdus["EXPOSURE"],
            "READMODE": self.file_list[0].split("_")[-4],
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
            The output file path. If None, a default filename is generated.
        format : str, optional
            The file format ("asdf" or "fits"). Default is "asdf".
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
            ff = asdf.AsdfFile(tree)
            ff.write_to(output)
        elif format in ["fits", "FITS"]:
            raise NotImplementedError
        else:
            raise ValueError("Provide a valid formate [FITS or asdf]")
