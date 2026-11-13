"""Utilities to help work with cube data"""

import os

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
    "1P",
    "2P",
    "A_",
    "AP_",
    "B_",
    "BP_",
]


def WCS_ATTRS(hdu, sip=True):
    """
    Extract WCS-related attribute keys from a FITS header.

    Parameters
    ----------
    hdu : astropy.io.fits.Header
    sip : bool, optional
        If True, include SIP-related WCS keywords. Default is True.

    Returns
    -------
    wcs_attrs : list of str
    """
    wcs_attrs = np.hstack(
        [
            *[
                [key for key in hdu.keys() if key.startswith(keystart)]
                for keystart in [WCS_ATTRS_STARTS if sip else WCS_ATTRS_STARTS[:-4]][0]
            ],
        ]
    ).tolist()
    return wcs_attrs


def extract_average_WCS(file_list: list):
    """
    Compute the average WCS from a list of FITS files.

    Parameters
    ----------
    file_list : list of str

    Returns
    -------
    wcs : astropy.wcs.WCS
    """
    df = extract_all_WCS(file_list).reset_index(drop=False)
    numeric = df.columns[[x in [float, int] for x in df.dtypes.values]].values
    strings = df.columns[[x not in [float, int] for x in df.dtypes.values]].values
    df_avg = df.loc[:, numeric].median()
    for k, v in df.loc[0, strings].items():  # noqa
        df_avg[k] = v
    df_avg.WCSAXES = int(df_avg.WCSAXES)
    df_avg.NAXIS = int(df_avg.NAXIS)
    wcs_hdu = fits.PrimaryHDU()
    for attr in df_avg.index:
        wcs_hdu.header[attr] = df_avg[attr]
    return WCS(wcs_hdu.header)


def extract_all_WCS(file_list: list):
    """
    Extract WCS keys and values from a list of FITS files.

    Parameters
    ----------
    file_list : list of str

    Returns
    -------
    df : pandas.DataFrame
        Indexed by exposure number.

    Raises
    ------
    ValueError
        If the file list is empty.
    """
    if len(file_list) == 0:
        raise ValueError("File list is empty. Cannot extract WCS.")
    frame_wcs = []
    frame_no = []
    for f in file_list:
        hdu = fits.getheader(f)
        aux_wcs = {}
        for attr in WCS_ATTRS(hdu):
            aux_wcs[attr] = hdu[attr]
        frame_wcs.append(aux_wcs)
        if "EXP_NO" in hdu.keys():
            exp_no = hdu["EXP_NO"]
        else:
            exp_no = int(os.path.basename(f).split("_")[8])
        frame_no.append(exp_no)
    df = pd.DataFrame(frame_wcs, index=frame_no)
    return df
