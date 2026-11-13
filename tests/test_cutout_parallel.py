"""Tests for parallel, multi-target, time-sorted make_cutout."""

import asdf
import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.time import Time
from gwcs import coordinate_frames as cf
from gwcs import wcs as gwcs_wcs

from roman_cuts.com.cutout import _normalize_targets, make_cutout

pytest.importorskip("asdf")
pytest.importorskip("gwcs")

RA0, DEC0 = 150.0, 2.0
SHAPE = (100, 100)
# Files written in this (shuffled) order; mid_time MJD = 60000 + TIMES[i]
TIMES = [3.0, 1.0, 4.0, 2.0]


def _make_gwcs():
    scale = 0.11 / 3600.0  # deg / pix
    det2sky = (
        models.Shift(-50) & models.Shift(-50)
        | models.Scale(scale) & models.Scale(scale)
        | models.Pix2Sky_TAN()
        | models.RotateNative2Celestial(RA0, DEC0, 180)
    )
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    sky = cf.CelestialFrame(
        reference_frame=coord.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    return gwcs_wcs.WCS([(detector, det2sky), (sky, None)])


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    d = tmp_path_factory.mktemp("asdf")
    w = _make_gwcs()
    paths = []
    for i, t in enumerate(TIMES):
        tree = {
            "roman": {
                "meta": {
                    "wcs": w,
                    "exposure": {"mid_time": Time(60000 + t, format="mjd")},
                    "observation": {"visit": i},
                    "instrument": {"detector": "WFI01"},
                },
                "data": np.full(SHAPE, 10.0 * t, dtype="f4"),
                "err": np.full(SHAPE, t, dtype="f4"),
                "dq": np.zeros(SHAPE, dtype="u4"),
                "dumo": np.zeros(SHAPE, dtype="u4"),
            }
        }
        p = d / f"f{i}.asdf"
        asdf.AsdfFile(tree).write_to(p)
        paths.append(str(p))
    return paths


def _check_sorted(cube, meta, gwcss, n_files):
    assert cube.shape[1] == n_files
    assert len(meta) == n_files and len(gwcss) == n_files
    mjd = np.array([t.mjd for t in meta["mid_time"]])
    assert np.all(np.diff(mjd) > 0)
    # data plane value encodes the time of the file it came from
    np.testing.assert_allclose(cube[0, :, 0, 0], 10.0 * (mjd - 60000))


@pytest.mark.parametrize("backend", ["thread", "process"])
@pytest.mark.parametrize("n_workers", [1, 3])
def test_sorted_and_aligned(files, backend, n_workers):
    cubes, metas, gwcsss = make_cutout(
        files, RA0, DEC0, size=(6, 8), n_workers=n_workers, backend=backend
    )
    assert len(cubes) == len(metas) == len(gwcsss) == 1
    assert cubes[0].shape == (4, 4, 6, 8)
    _check_sorted(cubes[0], metas[0], gwcsss[0], 4)


def test_single_file_keeps_dims(files):
    cubes, metas, gwcsss = make_cutout(files[0], RA0, DEC0, size=5)
    assert cubes[0].shape == (4, 1, 5, 5)
    assert len(metas[0]) == 1 and len(gwcsss[0]) == 1


def test_multi_target_sizes_and_oob(files):
    far_ra = RA0 + 5.0  # far outside the 100x100 detector
    ras, decs = [RA0, far_ra, RA0], [DEC0, DEC0, DEC0 + 0.001]
    cubes, metas, gwcsss = make_cutout(
        files, ras, decs, size=[(4, 4), (6, 5), 3], n_workers=2, lite=True
    )
    assert [c.shape for c in cubes] == [(1, 4, 4, 4), (1, 4, 6, 5), (1, 4, 3, 3)]
    assert np.all(np.isnan(cubes[1]))
    assert all(g is None for g in gwcsss[1])
    assert len(metas[1]) == 4
    assert not metas[1]["in_bounds"].any()
    assert metas[0]["in_bounds"].all() and metas[2]["in_bounds"].all()
    # in-bounds targets match single-target calls
    for t in (0, 2):
        ref, _, _ = make_cutout(
            files, ras[t], decs[t], size=cubes[t].shape[2:], lite=True
        )
        np.testing.assert_array_equal(cubes[t], ref[0])
        _check_sorted(cubes[t], metas[t], gwcsss[t], 4)


def test_size_normalization_errors():
    with pytest.raises(ValueError):
        _normalize_targets([1, 2], [1, 2], [(3, 3), (4, 4), (5, 5)])
    with pytest.raises(ValueError):
        _normalize_targets([1, 2], [1], 5)
    _, _, s = _normalize_targets([1, 2], [1, 2], (3, 4))
    assert s == [(3, 4), (3, 4)]


def test_failed_file_is_nan_and_last(files, tmp_path):
    bad = str(tmp_path / "missing.asdf")
    cubes, metas, _ = make_cutout([bad] + files, RA0, DEC0, size=4, n_workers=2)
    assert cubes[0].shape == (4, 5, 4, 4)
    assert np.all(np.isnan(cubes[0][:, -1]))
    with pytest.raises(Exception):
        make_cutout([bad], RA0, DEC0, size=4, strict=True)
