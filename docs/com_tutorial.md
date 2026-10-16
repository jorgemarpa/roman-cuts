# Commissioning data walkthrough

> A full interactive notebook tutorial will be added here.  In the meantime
> see `notebooks/cutout_devs.ipynb` in the
> [romcom repository](https://github.com/jorgemarpa/roman-cuts).

## Quick start

```python
import os
from astroquery.mast import MastMissions
from roman_cuts.com import query_mast_files, make_cutout, run_cutout_pipeline

# 1. Authenticate (RRN / MAST API token required)
mission = MastMissions(mission="roman")
mission.login(token=os.environ["MAST_API_TOKEN"])

# 2. Find L2 files for a specific field
products = query_mast_files(
    program_id=1055,
    detector=5,
    filter_name="F146",
    observation=13,       # primary sky-field identifier
    mission=mission,
)
print(f"{len(products)} files found")

# 3. Extract cutouts
cube, meta, gwcss = make_cutout(
    products["filename"].tolist(),
    ra=250.585,
    dec=36.451,
    size=100,
    mission=mission,
)
print(cube.shape)   # (N, 100, 100)

# 4. Or run the full pipeline in one call (also saves to disk)
cube, meta, gwcss = run_cutout_pipeline(
    ra=250.585, dec=36.451,
    program_id=1055, detector=5, filter_name="F146",
    cutout_size=100, observation=13,
    output_path="./output",
    mission=mission,
)
```
