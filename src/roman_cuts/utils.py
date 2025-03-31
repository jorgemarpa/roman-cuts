"""Utilities to help work with cube data"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

WCS_ATTRS_STARTS = [
    "CTYPE",
    "CRVAL",
    "CRPIX",
    "CUNIT",
    "NAXIS",
    "CD1",
    "CD2",
    "CDELT",
    "WCS",
    "PC",
    "RADE",
    # "1P",
    # "2P",
    # "A_",
    # "AP_",
    # "B_",
    # "BP_",
]


def WCS_ATTRS(hdu, sip=True):
    wcs_attrs = np.hstack(
        [
            *[
                [key for key in hdu.keys() if key.startswith(keystart)]
                for keystart in [WCS_ATTRS_STARTS if sip else WCS_ATTRS_STARTS[:-4]][0]
            ],
        ]
    ).tolist()
    return wcs_attrs


def extract_average_WCS(file_list: list = []):
    """
    Reads WCS keys and values from `file_list` FITS and computes an average WCS
    """
    # read WCS keywords and values into a list of dictionaries for all times
    frame_wcs = []
    for f in file_list:
        hdu = fits.getheader(f)
        aux_wcs = {}
        for attr in WCS_ATTRS(hdu):
            aux_wcs[attr] = hdu[attr]
        frame_wcs.append(aux_wcs)

    # take the median value of every keyword with numeric values
    df = pd.DataFrame(frame_wcs)
    numeric = df.columns[[x in [float, int] for x in df.dtypes.values]].values
    strings = df.columns[[x not in [float, int] for x in df.dtypes.values]].values
    df_avg = df.loc[:, numeric].median()
    # add non-numeric items to the dataframe
    for k, v in df.loc[0, strings].items():  # noqa
        df_avg[k] = v

    df_avg.WCSAXES = int(df_avg.WCSAXES)
    df_avg.NAXIS = int(df_avg.NAXIS)
    # make HDU with WCS keys and values
    wcs_hdu = fits.PrimaryHDU()
    for attr in df_avg.index:
        wcs_hdu.header[attr] = df_avg[attr]
    # return a WCS object
    return WCS(wcs_hdu.header)


def extract_all_WCS(file_list: list = []):
    """
    Reads WCS keys and values from `file_list` FITS and return them as a list
    """
    # read WCS for each frame
    wcss = []
    for f in file_list:
        hdu = fits.getheader(f)
        wcss.append(WCS(hdu))

    # return a list of WCS object
    return wcss

