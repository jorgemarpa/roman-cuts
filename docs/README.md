[![PyPI](https://img.shields.io/pypi/v/roman-cuts.svg)](https://pypi.org/project/roman-cuts)
[![pytest](https://github.com/jorgemarpa/roman-cuts/actions/workflows/pytest.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/pytest.yaml/) [![mypy](https://github.com/jorgemarpa/roman-cuts/actions/workflows/mypy.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/mypy.yaml) [![ruff](https://github.com/jorgemarpa/roman-cuts/actions/workflows/ruff.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/ruff.yaml)[![Docs](https://github.com/jorgemarpa/roman-cuts/actions/workflows/deploy-mkdocs.yaml/badge.svg)](https://github.com/jorgemarpa/roman-cuts/actions/workflows/deploy-mkdocs.yaml)

# Roman-cuts

Package to create image cutouts from simulations made with `RImTimSim` or from Roman Commissioning activities

## Install

Easy install with pip from your local clone of the repository
```bash
pip install -e .           # sim only (no cloud/MAST deps)
```
To use the tools with commissioning data, checkout this branch (`romcom`) and then run
```bash
pip install -e ".[com]"    # sim + com (all features)
```

## Usage

### Simulated images

For more details check out the notebook tutorial [here](./tutorial-sim.ipynb).


### ROMCOM data

See the tutorial [here](./tutorial-com.ipynb).

## Examples
The figure shows a sequence of 21x21 pixel cutouts taken from the FFI simulations 
centered on the target RA, Decl = (268.5112137932491, -29.24473947250156).
This account for dithered observations

![cutouts](./docs/figures/cutout_example.png)

The data is saved into a ASDF file as shown below:

![asdf](./docs/figures/asdf_example.png)