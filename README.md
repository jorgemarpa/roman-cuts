# Roman-cuts

Lightweight package to create image cutouts from simulations made with `RImTimSim`

## Install

```
pip install git+https://github.com/jorgemarpa/roman-cuts
```

## Usage

```python
from roman_cuts import make_cutout

rcube = RomanCuts(field=3, sca=1, filter="F146")
rcube.get_average_wcs()
data_cube = rcube.make_cutout(ra, dec, outfile="my_cube.asdf")
```