"""
roman_cutouts.io
~~~~~~~~~~~~~~~~
Platform detection and ASDF file-opening helpers.
"""
import logging
import os
import warnings
from contextlib import contextmanager
from functools import lru_cache

import asdf
import fsspec
import pandas as pd
import requests
import s3fs

logger = logging.getLogger(__name__)


def _detect_platform():
    """
    Return the current execution platform: ``'rrn'``, ``'smdc'``, or ``'standalone'``.

    Resolution order:

    1. ``ROMAN_ENV`` environment variable — set it explicitly to skip detection.
    2. JupyterHub env vars or a ``/teams`` mount → ``'rrn'``
       (STScI Roman Research Network).
    3. AWS EC2 instance-metadata endpoint reachable → ``'smdc'``
       (NASA Science Managed Data Center on AWS).
    4. Default → ``'standalone'``.
    """
    override = os.environ.get("ROMAN_ENV", "").lower()
    if override in ("rrn", "smdc", "standalone"):
        return override

    # RRN: STScI JupyterHub sets JUPYTERHUB_USER; /teams is the shared filesystem
    if os.environ.get("JUPYTERHUB_USER") or os.path.isdir("/teams"):
        return "rrn"

    # SMDC: running on an AWS EC2/ECS instance — link-local metadata is reachable
    try:
        import urllib.request
        urllib.request.urlopen("http://169.254.169.254/latest/meta-data/", timeout=0.3)
        return "smdc"
    except Exception:
        pass

    return "standalone"


@lru_cache(maxsize=None)
def _get_default_mission():
    """Build (once per process) a ``MastMissions`` object from ``MAST_API_TOKEN``."""
    from astroquery.mast import MastMissions
    mission = MastMissions(mission="roman")
    token = os.environ.get("MAST_API_TOKEN")
    if token:
        mission.login(token=token)
    return mission


@lru_cache(maxsize=None)
def _get_s3fs(anon):
    """Return a cached `s3fs.S3FileSystem` (one per process and credential mode)."""
    return s3fs.S3FileSystem(anon=anon)


@contextmanager
def _open_asdf(file_name, mission=None, platform=None):
    """
    Open an ASDF file from the appropriate backend for the current platform.

    **RRN** (STScI JupyterHub / Roman Research Network)
        Restricted L2 data must be streamed through the MAST API.  Pass an
        authenticated ``MastMissions`` instance as *mission*, or set the
        ``MAST_API_TOKEN`` environment variable to let this function create
        one automatically.

    **SMDC** (NASA Science Managed Data Center, AWS)
        Private S3 buckets are mounted as a POSIX filesystem; files are opened
        as ordinary local paths.  S3 URIs (``s3://``) are also supported and
        use IAM role credentials automatically.

    **Standalone**
        Public ``s3://`` URIs via anonymous `s3fs`, or plain local paths.

    Platform is resolved by `_detect_platform` (``ROMAN_ENV`` env var →
    JupyterHub indicators → AWS metadata endpoint → default).

    Parameters
    ----------
    file_name : str
        Local path, mounted S3 path, ``s3://`` URI, or MAST product filename.
    mission : `~astroquery.mast.MastMissions`, optional
        Logged-in mission object.  Supplying one forces the RRN/MAST path
        regardless of the detected platform.
    platform : {'rrn', 'smdc', 'standalone'}, optional
        Precomputed platform.  Pass it when opening many files so that
        detection (which may probe the network) runs only once.
    """
    from asdf.exceptions import AsdfWarning
    warnings.filterwarnings("ignore", category=AsdfWarning)

    if platform is None:
        platform = _detect_platform()

    if mission is not None or platform == "rrn":
        # RRN: restricted data via the MAST streaming API
        if mission is None:
            mission = _get_default_mission()
        product = mission.read_product(file_name)
        yield product

    elif str(file_name).startswith("s3://"):
        # SMDC: IAM role credentials; standalone: public S3
        fs = _get_s3fs(platform != "smdc")
        with fs.open(file_name, "rb") as f:
            with asdf.open(f, lazy_load=True) as af:
                yield af

    else:
        # Local path — covers both SMDC-mounted S3 and ordinary local files
        with asdf.open(str(file_name), lazy_load=True) as af:
            yield af


def open_parquet(filename, token):
    """
    Download one or more Parquet catalogs from the MAST S3 bucket using signed URLs.

    Fetches a pre-signed redirect URL from the MAST retrieve-product endpoint
    using *token* for authentication, then streams each Parquet file into a
    :class:`pandas.DataFrame` without writing to disk.

    Parameters
    ----------
    filename : str or list of str
        MAST product filename(s) (e.g. ``"roman/path/to/catalog.parquet"``).
        When a list is supplied a list of DataFrames is returned.
    token : str
        MAST API token used for the ``Authorization`` header.

    Returns
    -------
    pandas.DataFrame or list of pandas.DataFrame
        Catalog contents.  A single DataFrame when *filename* is a string;
        a list of DataFrames when *filename* is a list.
    """
    def _fetch_one(fname):
        resp = requests.get(
            f"https://mast.stsci.edu/search/roman/api/v0.1/retrieve_product?product_name={fname}",
            headers={"Authorization": f"token {token}"},
            allow_redirects=False,
        )
        resp.raise_for_status()
        fs = fsspec.filesystem("https")
        with fs.open(resp.headers["Location"], "rb") as f:
            return pd.read_parquet(f, engine="pyarrow")

    if isinstance(filename, str):
        return _fetch_one(filename)
    return [_fetch_one(fname) for fname in filename]
