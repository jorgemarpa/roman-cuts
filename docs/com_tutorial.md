# Commissioning data walkthrough

> A full interactive notebook tutorial will be added here.  In the meantime
> see `notebooks/cutout_devs.ipynb` in the
> [romcom repository](https://github.com/jorgemarpa/roman-cuts).

## Quick start

```python
import os
from astroquery.mast import MastMissions
from roman_cuts.com import (
    query_mast_files, make_cutout, run_cutout_pipeline, save_to_asdf,
)

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

# 3. Extract cutouts (files are processed in parallel, 8 workers by default)
cubes, metas, gwcsss = make_cutout(
    products["filename"].tolist(),
    ra=250.585,
    dec=36.451,
    size=100,
    mission=mission,
)
# One entry per target, even for a single target
cube, meta, gwcss = cubes[0], metas[0], gwcsss[0]
print(cube.shape)   # (4, n_files, 100, 100): [data, err, dq, dumo] x time
print(len(meta), len(gwcss))   # n_files each, sorted by mid_time
```

## Multiple targets

`ra` and `dec` accept lists. Every file is opened once and all targets are cut
from it. `size` can be an `int`, a `(ny, nx)` tuple shared by all targets, or a
list with one size per target.

```python
cubes, metas, gwcsss = make_cutout(
    products["filename"].tolist(),
    ra=[250.585, 250.600],
    dec=[36.451, 36.440],
    size=[(100, 100), (50, 80)],
    mission=mission,
    n_workers=8,          # 1 = serial
    backend="thread",     # or "process" (spawn); both work in Jupyter and the CLI
)
cubes[1].shape   # (4, n_files, 50, 80)
```

Notes:

- Outputs are sorted by `mid_time`; `cubes[t]`, `metas[t]` and `gwcsss[t]` share
  the same time axis (`n_files`).
- A target that falls outside a file's footprint is kept as NaN (and `None` in
  `gwcsss`). Files that fail to open are also NaN and sorted last, unless
  `strict=True`.
- `lite=True` extracts the data plane only, shape `(1, n_files, ny, nx)`.

## Saving

`save_to_asdf` writes **one file per target**. With more than one target,
`filename` must be a list of the same length. If omitted, files are written to
the current directory as
`cutouts_{pid}_{detector}_{filter}_{ra:.6f}_{dec:.6f}.asdf`, where `pid` is the
program number (the leading `PPPPP`) of the `observation_id` column in the
metadata (`PPPPPCCAAASSSOOOVVV`, see the
[WFI data levels and products](https://roman-docs.stsci.edu/data-handbook/wfi-data-levels-and-products)
file name conventions).

```python
save_to_asdf(cubes, metas, gwcsss)
save_to_asdf(cubes, metas, gwcsss, ["t0.asdf", "t1.asdf"])
```

## Full pipeline

```python
# Query + extract + save in one call (one ASDF file per target in output_path)
cubes, metas, gwcsss = run_cutout_pipeline(
    ra=[250.585, 250.600], dec=[36.451, 36.440],
    program_id=1055, detector=5, filter_name="F146",
    cutout_size=100, observation=13,
    output_path="./output",
    mission=mission,
    n_workers=8,
)
```

## Command line

```bash
roman-com-cutout --ra 250.585 250.600 --dec 36.451 36.440 \
    --program-id 1055 --detector 5 --filter F146 --observation 13 \
    --cutout-size 100 50,80 --n-workers 8 --backend thread \
    --output-path ./output
```

`--cutout-size` takes `N` (square) or `NY,NX`; give one value to share it, or one
per target.
