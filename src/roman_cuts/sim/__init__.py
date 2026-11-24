"""
roman_cuts.sim
~~~~~~~~~~~~~~
Cutout utilities for Roman WFI simulated images produced by RImTimSim.

Data lives on local disk as FITS or pre-stacked ASDF files.

Public API
----------
>>> from roman_cuts.sim import RomanCuts
"""
import logging
import os

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

log = logging.getLogger(__name__)

from .cube import RomanCuts  # noqa: E402

__all__ = ["RomanCuts"]
