"""
Actual test won't be possible as copying image files to GitHub repo will use too much 
memory.
"""

import pytest

from roman_cuts import RomanCuts


def test_romancuts():
    """
    Dummy function to test error when no files are provided
    """
    fl = []
    with pytest.raises(ValueError) as excinfo:
        _ = RomanCuts(field=3, sca=1, filter="F146", file_list=fl)
    assert str(excinfo.value) == "Please provide a list of FFI files in `file_list`"