#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import

import os
from importlib import metadata

# reads version number from package metadata
__version__ = metadata.version(__package__)

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

import logging  # noqa: E402


def get_logger():
    """
    Creates a logger
    """
    logger = logging.getLogger(__name__)
    # logger.addHandler(logging.StreamHandler())
    return logger


log = get_logger()

from .cube import RomanCuts  # noqa

__all__ = ["RomanCuts"]
