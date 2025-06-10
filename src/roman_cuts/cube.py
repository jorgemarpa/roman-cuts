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
        self, field: int, sca: int, filter: str = "F146", file_list: list = []
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

        if len(file_list) == 0:
            raise ValueError("Please provide a list of FFI files in `file_list`")
        if not isinstance(file_list, (list, np.ndarray)):
            file_list = np.sort([file_list])
        self.file_list = file_list
        self.nt = len(file_list)

        self._check_file_list()

        log.info("Getting 1d arrays data...")
        self._get_arrays()
        log.info("Getting metadata")
        self._get_metadata()

    def __repr__(self):
        return f"Roman WFI Field {self.field} SCA {self.sca} Filter {self.filter} Frames {self.nt}"

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
            hdr = fits.getheader(f)
            # field.append(hdr["FIELD"])
            sca.append(hdr["DETECTOR"])
            filter.append(hdr["FILTER"])

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
                row = np.round(row).astype(int)
                col = np.round(col).astype(int)
                self.target_pixel = np.array([row, col]).T
            else:
                log.info(
                    "Using Ra, Dec coordinates and average WCS to center the cutout"
                )
                row, col = self.wcs.all_world2pix(self.ra, self.dec, 0)
                row = np.array([int(np.round(row))])
                col = np.array([int(np.round(col))])
                self.target_pixel = np.array([row, col])
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

        if (rmin > RMAX) | (cmin > CMAX):
            raise ValueError("`cutout_origin` must be within the image.")

        # set ending pixels
        rmax = rmin + size[0]
        cmax = cmin + size[1]

        if (rmax > RMAX) | (cmax > CMAX):
            log.warning("Cutout exceeds image limits, reducing size.")
            rmax = np.min([rmax, RMAX])
            cmax = np.min([cmax, CMAX])

        flux = []
        flux_err = []
        for f in tqdm(self.file_list):
            aux = fits.open(f)
            flux.append(aux[0].data[rmin:rmax, cmin:cmax])
            flux_err.append(aux[1].data[rmin:rmax, cmin:cmax])
            aux.close()

        self.flux = np.array(flux)
        self.flux_err = np.array(flux_err)
        self.row = np.arange(rmin, rmax)
        self.column = np.arange(cmin, cmax)
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
        if center.shape[0] != len(self.file_list) and not center.shape[0] == 1:
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
            (rmin > RMAX).any()
            | (cmin > CMAX).any()
            | (rmax > RMAX).any()
            | (cmax > CMAX).any()
        ):
            raise ValueError(
                "Cutout out of CCD limits. This is due to the dithered observations"
                " and the size of the cutout. Please reduce the size or change the center."
            )

        flux = []
        flux_err = []
        for i, f in tqdm(enumerate(self.file_list), total=len(self.file_list)):
            aux = fits.open(f)
            flux.append(aux[0].data[rmin[i] : rmax[i], cmin[i] : cmax[i]])
            flux_err.append(aux[1].data[rmin[i] : rmax[i], cmin[i] : cmax[i]])
            aux.close()

        self.flux = np.array(flux)
        self.flux_err = np.array(flux_err)
        self.row = np.vstack([np.arange(rn, rx) for rn, rx in zip(rmin, rmax)])
        self.column = np.vstack([np.arange(cn, cx) for cn, cx in zip(cmin, cmax)])
        return

    def _get_arrays(self):
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

    def _get_metadata(self):
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
            "NTIMES": self.nt,
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
