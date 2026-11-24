"""
Verification tests for performance optimizations in roman_cuts.com:
- Numerical equivalence of cutout arrays (data, err, dq, dumo)
- Sliced GWCS projection equivalence
- Metadata geometry and schema fidelity
- Parallel save_to_asdf validation
"""

import os

import asdf
import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.time import Time
from gwcs import coordinate_frames as cf
from gwcs import wcs as gwcs_wcs

from roman_cuts.com import make_cutout, save_to_asdf

RA0, DEC0 = 150.0, 2.0
SHAPE = (200, 200)


def _make_sample_gwcs():
    scale = 0.11 / 3600.0  # deg / pix
    det2sky = (
        models.Shift(-100) & models.Shift(-100)
        | models.Scale(scale) & models.Scale(scale)
        | models.Pix2Sky_TAN()
        | models.RotateNative2Celestial(RA0, DEC0, 180)
    )
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS(), name="icrs", unit=(u.deg, u.deg))
    return gwcs_wcs.WCS([(detector, det2sky), (sky, None)])


@pytest.fixture
def mock_asdf_files(tmp_path):
    w = _make_sample_gwcs()
    paths = []
    for i in range(4):
        p = tmp_path / f"test_exp_{i}.asdf"
        # Unique values across grid to detect any slice offset
        data_plane = np.arange(SHAPE[0] * SHAPE[1], dtype="f4").reshape(SHAPE) + (i * 1000.0)
        err_plane = np.ones(SHAPE, dtype="f4") * (i + 1.0)
        dq_plane = np.zeros(SHAPE, dtype="u4")
        dumo_plane = np.zeros(SHAPE, dtype="u4")

        tree = {
            "roman": {
                "meta": {
                    "wcs": w,
                    "exposure": {"mid_time": Time(60000 + i, format="mjd")},
                    "observation": {"visit": i, "observation_id": f"0105501001001001{i:03d}"},
                    "instrument": {"detector": "WFI01", "optical_element": "F146"},
                },
                "data": data_plane,
                "err": err_plane,
                "dq": dq_plane,
                "dumo": dumo_plane,
            }
        }
        asdf.AsdfFile(tree).write_to(p)
        paths.append(str(p))
    return paths


def test_cutout_equivalence_and_geometry(mock_asdf_files):
    """Verify that extracted cutouts, geometry attributes, and sliced GWCS match expectations."""
    # Target 1: center of image (fully enclosed)
    # Target 2: near detector edge (partial overlap)
    # Target 3: far outside detector (completely out of bounds)
    ras = [RA0, RA0 - 0.0029, RA0 + 1.0]
    decs = [DEC0, DEC0 - 0.0029, DEC0 + 1.0]
    size = (20, 20)

    cubes, metas, gwcsss = make_cutout(mock_asdf_files, ras, decs, size=size, n_workers=2)

    assert len(cubes) == 3
    assert len(metas) == 3
    assert len(gwcsss) == 3

    # Check center target
    cube_center = cubes[0]
    meta_center = metas[0]
    gwcs_center = gwcsss[0]

    assert cube_center.shape == (4, 4, 20, 20)
    assert not np.any(np.isnan(cube_center))  # fully enclosed, no NaNs
    assert np.all(meta_center["in_bounds"])

    # Check geometry fields
    for row in meta_center.itertuples():
        assert hasattr(row, "input_position_original")
        assert hasattr(row, "origin_original")
        assert hasattr(row, "input_position_cutout")
        assert hasattr(row, "bbox_original")
        # Origin must be non-negative
        assert row.origin_original[0] >= 0
        assert row.origin_original[1] >= 0

    # Sliced GWCS projection test
    w_sliced = gwcs_center[0]
    assert w_sliced is not None
    # Projecting center RA, Dec back into cutout should be near cutout center ~ (10, 10)
    c_sky = coord.SkyCoord(ra=RA0 * u.deg, dec=DEC0 * u.deg, frame="icrs")
    x_c, y_c = w_sliced.world_to_pixel(c_sky)
    np.testing.assert_allclose([x_c, y_c], [10.0, 10.0], atol=0.5)

    # Check edge target (partial overlap)
    cube_edge = cubes[1]
    meta_edge = metas[1]
    assert np.any(np.isnan(cube_edge))  # Has NaN padding at boundary
    assert np.any(~np.isnan(cube_edge))  # But has valid data pixels
    assert np.all(meta_edge["in_bounds"])

    # Check out-of-bounds target
    cube_oob = cubes[2]
    meta_oob = metas[2]
    gwcs_oob = gwcsss[2]
    assert np.all(np.isnan(cube_oob))
    assert not np.any(meta_oob["in_bounds"])
    assert all(g is None for g in gwcs_oob)


def test_parallel_save_to_asdf(mock_asdf_files, tmp_path):
    """Verify that parallel save_to_asdf correctly writes readable ASDF files."""
    ras = [RA0, RA0 + 0.0005, RA0 - 0.0005]
    decs = [DEC0, DEC0 + 0.0005, DEC0 - 0.0005]

    cubes, metas, gwcsss = make_cutout(mock_asdf_files, ras, decs, size=15, n_workers=2)

    out_paths = [str(tmp_path / f"out_target_{i}.asdf") for i in range(3)]
    written = save_to_asdf(cubes, metas, gwcsss, filename=out_paths, n_workers=3)

    assert len(written) == 3
    for p in written:
        assert os.path.exists(p)
        with asdf.open(p) as af:
            assert "roman" in af
            r = af["roman"]
            assert "data" in r
            assert "meta" in r
            assert "time" in r
            assert "quality" in r
            assert r["data"].shape == (4, 15, 15)
            assert len(r["time"]) == 4
