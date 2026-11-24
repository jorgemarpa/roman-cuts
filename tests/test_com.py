"""
Unit tests for roman_cuts.com.

Integration tests that require MAST credentials or actual ASDF files are not
run in CI.  Each is marked with ``pytest.mark.integration`` and skipped by
default; run them locally with::

    pytest -m integration
"""

import pytest

from roman_cuts.com import (
    make_cutout,
    query_mast_files,
    run_cutout_pipeline,
    save_to_asdf,
)
from roman_cuts.com.io import _detect_platform
from roman_cuts.com.query import FILE_TYPE_MAP


def test_public_api_importable():
    """All public com symbols are importable from roman_cuts.com."""
    assert callable(query_mast_files)
    assert callable(make_cutout)
    assert callable(run_cutout_pipeline)
    assert callable(save_to_asdf)


def test_file_type_map_keys():
    """FILE_TYPE_MAP contains l1, l2, and l4 entries."""
    assert set(FILE_TYPE_MAP.keys()) == {"l1", "l2", "l4"}


def test_query_mast_files_bad_file_type():
    """query_mast_files raises ValueError for an unknown file_type."""
    with pytest.raises(ValueError, match="file_type must be one of"):
        query_mast_files(1055, 5, "F146", file_type="l3")


def test_detect_platform_returns_known_value():
    """_detect_platform always returns a known string."""
    platform = _detect_platform()
    assert platform in ("rrn", "smdc", "standalone")


def test_detect_platform_env_override(monkeypatch):
    """ROMAN_ENV env var overrides auto-detection."""
    monkeypatch.setenv("ROMAN_ENV", "smdc")
    assert _detect_platform() == "smdc"
    monkeypatch.setenv("ROMAN_ENV", "rrn")
    assert _detect_platform() == "rrn"
    monkeypatch.setenv("ROMAN_ENV", "standalone")
    assert _detect_platform() == "standalone"


@pytest.mark.integration
def test_query_mast_files_live():
    """Live MAST query — requires MAST_API_TOKEN in the environment."""
    import os
    if not os.environ.get("MAST_API_TOKEN"):
        pytest.skip("MAST_API_TOKEN not set")
    products = query_mast_files(1055, 5, "F146", observation=1)
    assert len(products) > 0
    assert "filename" in products.colnames
