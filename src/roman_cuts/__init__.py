#!/usr/bin/env python
"""
roman_cuts
==========
Roman WFI cutout utilities — simulation data and commissioning/science data.

Two submodules are provided:

``roman_cuts.sim``
    Cutouts from local FITS/ASDF simulation files produced by RImTimSim.
    OOP interface centred on the :class:`~roman_cuts.sim.RomanCuts` class.

``roman_cuts.com``
    Cutouts from MAST-archived commissioning and science ASDF files accessed
    via the MAST API or S3.  Functional interface via :func:`make_cutout`,
    :func:`query_mast_files`, and :func:`run_cutout_pipeline`.

Quick start
-----------
**Simulation data (local FITS)**::

    from roman_cuts import RomanCuts           # or: from roman_cuts.sim import RomanCuts
    rcube = RomanCuts(field=3, sca=1, filter="F146", file_list=[...])
    rcube.make_cutout(radec=(268.46, -29.23), size=(21, 21), dithered=True)
    rcube.save_cutout()

**Commissioning data (MAST cloud)**::

    from roman_cuts import query_mast_files, make_cutout
    products = query_mast_files(1055, 5, "F146", observation=13, mission=mission)
    cube, meta, gwcss = make_cutout(products["filename"].tolist(), ra, dec, 100,
                                    mission=mission)
"""
import logging
from importlib import metadata

try:
    __version__ = metadata.version(__package__)
except metadata.PackageNotFoundError:
    __version__ = "0.0.0.dev"

log = logging.getLogger(__name__)

# ── sim submodule (simulation data) ──────────────────────────────────────────
# ── com submodule (commissioning / science data) ─────────────────────────────
from .com import (  # noqa: E402
    make_cutout,
    query_mast_files,
    run_cutout_pipeline,
    save_to_asdf,
)
from .sim import RomanCuts  # noqa: E402

__all__ = [
    # sim
    "RomanCuts",
    # com
    "query_mast_files",
    "make_cutout",
    "run_cutout_pipeline",
    "save_to_asdf",
]
