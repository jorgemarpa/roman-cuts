"""
Actual test won't be possible as copying image files to GitHub repo will use too much
memory.
"""

import pytest

from roman_cuts import RomanCuts           # top-level re-export (backwards compat)
from roman_cuts.sim import RomanCuts as RomanCutsSim  # explicit submodule path


def test_romancuts_empty_file_list():
    """RomanCuts raises ValueError when no files are provided."""
    with pytest.raises(ValueError) as excinfo:
        _ = RomanCuts(field=3, sca=1, filter="F146", file_list=[])
    assert str(excinfo.value) == "Please provide a list of FFI files in `file_list`"


def test_romancuts_sim_alias():
    """roman_cuts.sim.RomanCuts is the same class as roman_cuts.RomanCuts."""
    assert RomanCuts is RomanCutsSim
