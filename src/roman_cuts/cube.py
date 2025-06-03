import bz2
import json
import os
from typing import Optional, Tuple

import asdf
import numpy as np
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
        scs : int
            The scs number.
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
            wcss_df = pd.read_json(filename, oreint="index", compression="bz2")
        # convert to list of WCS objects
        self.wcss = [
            WCS(wcs_dict, relax=True)
            for wcs_dict in wcss_df.to_dict(orient="index").values()
        ]
        return

    def make_cutout(
        self,
        radec: Tuple = (None, None),
        rowcol: Tuple[int, int] = (0, 0),
        size: Tuple[int, int] = (15, 15),
        load_all_wcss: bool = False,
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
        # use radec if given
        if (
            radec != (None, None)
            and isinstance(radec[0], float)
            and isinstance(radec[1], float)
        ):
            log.info("Using Ra, Dec coordinates and average WCS to center the cutout")
            self.ra = radec[0]
            self.dec = radec[1]
            row, col = self.wcs.all_world2pix(self.ra, self.dec, 0)
            row = int(np.round(row))
            col = int(np.round(col))
        # if not use the rowcol
        elif isinstance(rowcol[0], int) and isinstance(rowcol[1], int):
            row, col = rowcol[0], rowcol[1]
        # raise error if values are not valid
        else:
            raise ValueError("Please provide valid `radec` or `rowcol` values")

        origin = (int(row - size[0] / 2), int(col - size[1] / 2))

        log.info("Getting 3d data...")
        self._get_cutout_cube(size=size, origin=origin)

        print()
        if not hasattr(self, "ra"):
            self.ra, self.dec = self.wcs.all_pix2world(row, col, 0)

        log.info("Getting 1d arrays data...")
        self._get_arrays()
        log.info("Getting metadata")
        self._get_metadata()
        if load_all_wcss:
            log.info("Getting WCSs per frame")
            self.get_all_wcs()
        return

    def _get_cutout_cube(self, size: Tuple = (15, 15), origin: Tuple = (0, 0)):
        """
        Extracts a cutout cube from the FFI files.

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

        rmax = rmin + size[0]
        cmax = cmin + size[1]

        if (rmax > RMAX) | (cmax > CMAX):
            log.warning("Cutout exceeds image limits, reducing size.")
            rmax = np.min([rmax, RMAX])
            cmax = np.min([cmax, CMAX])

        flux = []
        flux_err = []
        for f in tqdm(self.file_list):
            with fits.open(f, lazy_load_hdus=True) as aux:
                flux.append(aux[0].data[rmin:rmax, cmin:cmax])
                flux_err.append(aux[1].data[rmin:rmax, cmin:cmax])

        self.flux = np.array(flux)
        self.flux_err = np.array(flux_err)
        self.row = np.arange(rmin, rmax)
        self.column = np.arange(cmin, cmax)
        return

    def _get_arrays(self):
        """
        Extracts time, cadenceno, and quality arrays from the FFI files.
        """
        time, cadenceno, quality = [], [], []
        for k, f in enumerate(self.file_list):
            hdu = fits.getheader(f)
            time.append((hdu["TSTART"] + hdu["TEND"]) / 2.0)
            # replace these two to corresponding keywords in future simulations
            cadenceno.append(k)
            quality.append(0)
        self.time = np.array(time)
        self.cadenceno = np.array(cadenceno)
        self.quality = np.array(quality)

    def _get_metadata(self):
        """
        Extracts metadata from the first FFI file.
        """
        hdu = fits.getheader(self.file_list[0])
        self.metadata = {
            "MISSION": "Roman",
            "TELESCOP": "Roman",
            "CREATOR": "TRExS",
            "SOFTWARE": hdu["SOFTWARE"],
            "RADESYS": hdu["RADESYS"],
            "EQUINOX": hdu["EQUINOX"],
            "FILTER": hdu["FILTER"],
            "FIELD": int(self.file_list[0].split("_")[-5][-2:]),
            "DETECTOR": hdu["DETECTOR"],
            "EXPOSURE": hdu["EXPOSURE"],
            "READMODE": self.file_list[0].split("_")[-4],
            "RA_CEN": float(self.ra),
            "DEC_CEN": float(self.dec),
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
            cutout_str = (
                f"{self.ra:.4f}_{self.dec:.4f}_s{len(self.row)}x{len(self.column)}"
            )
            output = f"./roman_cutout_field{self.metadata['FIELD']:02}_{self.metadata['DETECTOR']:02}_{cutout_str}.{format}"
            log.info(f"Saving data to {output}")

        if isinstance(output, str) and not output.endswith(format):
            raise ValueError(
                "Use a valid and matching extension in `output` and `format`"
            )

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
                        "cadenceno": self.cadenceno,
                        "quality": self.quality,
                        "row": self.row,
                        "column": self.column,
                    },
                }
            }
            ff = asdf.AsdfFile(tree)
            ff.write_to(output)
        elif format in ["fits", "FITS"]:
            raise NotImplementedError
        else:
            raise ValueError("Provide a valid formate [FITS or asdf]")
