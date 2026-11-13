"""
roman_cuts.schema
~~~~~~~~~~~~~~~~~
Shared ASDF output-tree schema used by both ``sim`` and ``com`` modules.

Both modules write a ``roman`` top-level node with the same structure so
that outputs are interchangeable for downstream consumers (light-curve
extraction, TPF creation, …).

Schema
------
::

    roman/
      meta/      – dict of serialisable scalars / arrays (exposure, obs info …)
      wcs/       – list of per-frame WCS objects (astropy WCS or GWCS)
      data/      – pixel flux cube, shape (N, ny, nx)
      time/      – float array of mid-exposure times (JD)
      quality/   – int16 array of per-frame quality flags, shape (N,)

The ``sim`` module additionally stores ``flux_err``, ``exposureno``, ``row``,
and ``column`` under ``data/`` for backwards compatibility with the original
``roman-cuts`` format.
"""

# Top-level tree key used by both modules
ROMAN_KEY = "roman"

# Sub-keys
META_KEY = "meta"
WCS_KEY = "wcs"
DATA_KEY = "data"
TIME_KEY = "time"
QUALITY_KEY = "quality"
