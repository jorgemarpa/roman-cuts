[![PyPI](https://img.shields.io/pypi/v/roman-cuts.svg)](https://pypi.org/project/roman-cuts)
[![pytest](https://github.com/jorgemarpa/roman-cuts/actions/workflows/pytest.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/pytest.yaml/) [![mypy](https://github.com/jorgemarpa/roman-cuts/actions/workflows/mypy.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/mypy.yaml) [![ruff](https://github.com/jorgemarpa/roman-cuts/actions/workflows/ruff.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/ruff.yaml)[![Docs](https://github.com/jorgemarpa/roman-cuts/actions/workflows/deploy-mkdocs.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/deploy-mkdocs.yaml)

# Roman-cuts

Lightweight package to create image cutouts from simulations made with `RImTimSim`

## Install

Easy install with PyPI
```
pip install roman-cuts
```

## Usage

For more details check out the notebook tutorial [here](./tutorial.ipynb).

```python
from roman_cuts import RomanCuts

# make a list of your local FITS files
fl = paths to local FITS files

rcube = RomanCuts(field=3, sca=1, filter="F146", file_list=fl)

# using sky coord coordinates
radec = (268.461687, -29.232092)
rcube.make_cutout(radec=radec, size=(21, 21), dithered=True)

# or using rowcol pixel coordinates
rowcol = (256, 256)
rcube.make_cutout(rowcol=rowcol, size=(11, 11), dithered=False)

# we can save to disk, default is ASDF
rcube.save_cutout()
```

## Examples
The figure shows a sequence of 21x21 pixel cutouts taken from the FFI simulations 
centered on the target RA, Decl = (268.5112137932491, -29.24473947250156).
This account for dithered observations

![cutouts](./docs/figures/cutout_example.png)

The data is saved into a ASDF file as shown below:

![asdf](./docs/figures/asdf_example.png)