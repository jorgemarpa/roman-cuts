# Roman-cuts

Lightweight package to create image cutouts from simulations made with `RImTimSim`

## Install

```
pip install roman-cuts
```

## Usage

```python
from roman_cuts import RomanCuts

# make a list of your local FITS files
fl = paths to local FITS files

rcube = RomanCuts(field=3, sca=1, filter="F146", file_list=fl)

# using sky coord coordinates
radec = (268.461687, -29.232092)
rcube.make_cutout(radec=radec, size=(21, 21))
# or using rowcol pixel coordinates

rowcol = (256, 256)
rcube.make_cutout(rowcol=rowcol, size=(11, 11))

# we can save to disk, default is ASDF
rcube.save_cutout()
```