"""Tests of r.in.s1slc on a synthetic Sentinel-1 IW SLC product."""

import json
import os

import numpy as np
import pytest

import grass.script as gs

from conftest import (
    BURST_STEP,
    DT,
    FIRST_VALID_LINE,
    FIRST_VALID_SAMPLE,
    LAST_VALID_LINE,
    LAST_VALID_SAMPLE,
    LPB,
    NSAMPLES,
    cal_value,
    latitude,
    noise_azimuth,
    noise_range,
    pixel_i,
    pixel_q,
    read_map,
    run_module,
)


def expected_deburst():
    """Burst, line-in-burst and azimuth time of every debursted output row.

    Between bursts a and a+1 the switch happens at the middle of the full
    burst overlap, i.e. at t_a + (LPB - 1 + BURST_STEP) / 2 lines.
    """
    rows = []
    t = FIRST_VALID_LINE * DT
    t_end = 2 * BURST_STEP * DT + LAST_VALID_LINE * DT
    while t <= t_end + 1e-9:
        k = 0
        for a in range(2):
            mid = (a * BURST_STEP + 0.5 * (LPB - 1 + BURST_STEP)) * DT
            if t > mid + 1e-9:
                k = a + 1
        line = int(round(t / DT)) - k * BURST_STEP
        rows.append((k, line, t))
        t += DT
    return rows


def valid_columns():
    cols = np.arange(NSAMPLES)
    return (cols >= FIRST_VALID_SAMPLE) & (cols <= LAST_VALID_SAMPLE)


def test_list(xy_session, safe_product):
    proc = run_module(xy_session, "l", input=safe_product)
    assert proc.returncode == 0, proc.stderr
    lines = proc.stdout.strip().splitlines()
    assert lines[0].startswith("swath|burst|burst_id")
    assert len(lines) == 4
    assert lines[1].startswith("IW1|1|500|")


def test_deburst_complex(xy_session, safe_product):
    proc = run_module(xy_session, input=safe_product, output="s1")
    assert proc.returncode == 0, proc.stderr
    i = read_map(xy_session, "s1_iw1_vv_i", null=-99999)
    q = read_map(xy_session, "s1_iw1_vv_q", null=-99999)
    rows = expected_deburst()
    assert i.shape == (len(rows), NSAMPLES)
    valid = valid_columns()
    for r, (k, line, _t) in enumerate(rows):
        assert np.all(i[r, valid] == pixel_i(k, line)), (r, k, line, i[r, valid][0])
        assert np.all(q[r, valid] == pixel_q(np.arange(NSAMPLES))[valid])
        assert np.all(i[r, ~valid] == -99999)
    # Seams: no source line is duplicated or skipped across the overlap.
    info = gs.raster_info("s1_iw1_vv_i", env=xy_session.env)
    assert info["datatype"] == "CELL"
    assert int(info["rows"]) == 104


def test_bursts_separate(xy_session, safe_product):
    proc = run_module(xy_session, "b", input=safe_product, output="s1", bursts="2", measure="intensity")
    assert proc.returncode == 0, proc.stderr
    power = read_map(xy_session, "s1_iw1_vv_b02_intensity")
    assert power.shape == (LPB, NSAMPLES)
    valid = valid_columns()
    q = pixel_q(np.arange(NSAMPLES)).astype(float)
    for line in range(LPB):
        if FIRST_VALID_LINE <= line <= LAST_VALID_LINE:
            np.testing.assert_allclose(power[line, valid], pixel_i(1, line) ** 2 + q[valid] ** 2, rtol=1e-6)
        else:
            assert np.all(np.isnan(power[line]))


def read_description(session, name):
    env = gs.gisenv(env=session.env)
    path = os.path.join(
        env["GISDBASE"], env["LOCATION_NAME"], env["MAPSET"], "cell_misc", name, "description.json"
    )
    with open(path) as fd:
        return json.load(fd)


def test_deburst_segments(xy_session, safe_product):
    proc = run_module(xy_session, input=safe_product, output="s1")
    assert proc.returncode == 0, proc.stderr
    segments = read_description(xy_session, "s1_iw1_vv_i")["raster_geometry"]["segments"]
    rows = expected_deburst()
    assert sum(s["rows"] for s in segments) == len(rows)
    assert [s["burst"] for s in segments] == [1, 2, 3]
    assert [s["burst_id"] for s in segments] == [500, 501, 502]
    for s in segments:
        for r in range(s["first_row"], s["first_row"] + s["rows"]):
            k, line, _t = rows[r]
            assert (k + 1, line) == (s["burst"], s["first_line"] + r - s["first_row"])


def test_bursts_complex_for_interferometry(xy_session, safe_product):
    """All bursts kept whole, overlaps included, as needed by ESD."""
    proc = run_module(xy_session, "b", input=safe_product, output="s1")
    assert proc.returncode == 0, proc.stderr
    valid = valid_columns()
    for k in range(3):
        name = f"s1_iw1_vv_b{k + 1:02d}"
        i = read_map(xy_session, name + "_i", null=-99999)
        assert i.shape == (LPB, NSAMPLES)
        for line in range(FIRST_VALID_LINE, LAST_VALID_LINE + 1):
            assert np.all(i[line, valid] == pixel_i(k, line))
        geometry = read_description(xy_session, name + "_i")["raster_geometry"]
        assert geometry["debursted"] is False
        assert geometry["segments"] == [
            {"burst": k + 1, "burst_id": 500 + k, "first_row": 0, "rows": LPB, "first_line": 0}
        ]


def test_calibration_and_noise(xy_session, safe_product):
    proc = run_module(
        xy_session,
        "n",
        input=safe_product,
        output="s1",
        calibration="sigma0",
        measure="intensity,db",
    )
    assert proc.returncode == 0, proc.stderr
    sigma0 = read_map(xy_session, "s1_iw1_vv_intensity")
    db = read_map(xy_session, "s1_iw1_vv_db")
    cols = np.arange(NSAMPLES)
    valid = valid_columns()
    for r, (k, line, t) in enumerate(expected_deburst()):
        dn2 = pixel_i(k, line) ** 2 + pixel_q(cols).astype(float) ** 2
        noise = noise_range(cols) * noise_azimuth(k * LPB + line)
        # The source line time is the burst line time, not the output time.
        t_src = k * BURST_STEP * DT + line * DT
        expected = np.maximum(dn2 - noise, 0) / cal_value(t_src, cols) ** 2
        np.testing.assert_allclose(sigma0[r, valid], expected[valid], rtol=1e-5)
        with np.errstate(divide="ignore"):
            exp_db = np.where(expected > 0, 10 * np.log10(expected), np.nan)
        np.testing.assert_allclose(db[r, valid], exp_db[valid], rtol=1e-5, atol=1e-4)


def test_calibrated_complex(xy_session, safe_product):
    proc = run_module(xy_session, input=safe_product, output="s1", calibration="sigma0", bursts="1")
    assert proc.returncode == 0, proc.stderr
    i = read_map(xy_session, "s1_iw1_vv_i")
    assert gs.raster_info("s1_iw1_vv_i", env=xy_session.env)["datatype"] == "FCELL"
    line = 10
    r = line - FIRST_VALID_LINE
    t = line * DT
    cols = np.arange(NSAMPLES)
    valid = valid_columns()
    np.testing.assert_allclose(i[r, valid], (pixel_i(0, line) / cal_value(t, cols))[valid], rtol=1e-5)


def test_geometry_and_gcps(xy_session, safe_product):
    proc = run_module(
        xy_session, input=safe_product, output="s1", measure="amplitude", geometry="latitude,longitude"
    )
    assert proc.returncode == 0, proc.stderr
    lat = read_map(xy_session, "s1_iw1_latitude")
    rows = expected_deburst()
    times = np.array([t for _k, _l, t in rows])
    np.testing.assert_allclose(lat[:, 5], latitude(times), rtol=0, atol=1e-6)
    env = gs.gisenv(env=xy_session.env)
    group_dir = os.path.join(env["GISDBASE"], env["LOCATION_NAME"], env["MAPSET"], "group", "s1_iw1_vv")
    refs = open(os.path.join(group_dir, "REF")).read()
    assert "s1_iw1_vv_amplitude" in refs and "s1_iw1_latitude" in refs
    points = np.loadtxt(os.path.join(group_dir, "POINTS"), comments="#")
    assert points.shape[1] == 5
    # First GCP: image row 0 centre, sample 0 centre.
    e1, n1, lon, lat0, status = points[0]
    assert (e1, n1) == (0.5, len(rows) - 0.5)
    np.testing.assert_allclose(lat0, latitude(times[0]), atol=1e-9)
    np.testing.assert_allclose(lon, 55.0, atol=1e-9)


@pytest.mark.parametrize("burst", [1, 2, 3])
def test_gcps_distinct(xy_session, safe_product, burst):
    # Grid lines 1 us after a burst start must not duplicate the first row.
    proc = run_module(xy_session, "b", input=safe_product, output="g", bursts=str(burst), measure="amplitude")
    assert proc.returncode == 0, proc.stderr
    env = gs.gisenv(env=xy_session.env)
    points = np.loadtxt(
        os.path.join(
            env["GISDBASE"], env["LOCATION_NAME"], env["MAPSET"], "group", f"g_iw1_vv_b{burst:02d}", "POINTS"
        ),
        comments="#",
    )
    image_rows = np.unique(points[:, 1])
    assert np.diff(image_rows).min() > 0.5
    assert np.unique(points[:, :2], axis=0).shape[0] == points.shape[0]
    assert image_rows[0] == 0.5 and image_rows[-1] == LPB - 0.5


def test_metadata(xy_session, safe_product):
    proc = run_module(xy_session, input=safe_product, output="s1", measure="phase")
    assert proc.returncode == 0, proc.stderr
    info = gs.raster_info("s1_iw1_vv_phase", env=xy_session.env)
    assert info["semantic_label"] == "S1_VV_PHASE"
    assert info["units"] == "radians"
    stamp = gs.read_command("r.timestamp", map="s1_iw1_vv_phase", env=xy_session.env)
    assert stamp.startswith("12 Jan 2023 14:25:09")
    env = gs.gisenv(env=xy_session.env)
    path = os.path.join(
        env["GISDBASE"], env["LOCATION_NAME"], env["MAPSET"], "cell_misc", "s1_iw1_vv_phase", "description.json"
    )
    meta = json.load(open(path))
    assert meta["product"]["relative_orbit_start"] == "130"
    assert meta["product"]["ipf_version"] == "003.52"
    assert meta["swath"]["lines_per_burst"] == LPB
    assert len(meta["swath"]["orbit_state_vectors"]) == 2
    assert meta["raster_geometry"]["rows"] == 104
    assert meta["raster_geometry"]["first_line_time"] == "2023-01-12T14:25:09.020000"


def test_zip_and_bbox(xy_session, safe_zip):
    # The bbox touches the latitude span of burst 3 only.
    proc = run_module(
        xy_session, "b", input=safe_zip, output="z", measure="amplitude", bbox="54.9,25.0081,55.1,25.0085"
    )
    assert proc.returncode == 0, proc.stderr
    maps = gs.list_strings("raster", pattern="z_*", env=xy_session.env)
    assert [m.split("@")[0] for m in maps] == ["z_iw1_vv_b03_amplitude"]


@pytest.mark.parametrize(
    ("flags", "kwargs", "message"),
    [
        ("n", {"measure": "complex"}, "noise"),
        ("", {"bursts": "1,3"}, "contiguous"),
        ("b", {"bursts": "2-4"}, "out of range"),
    ],
)
def test_failures(xy_session, safe_product, flags, kwargs, message):
    args = (flags,) if flags else ()
    proc = run_module(xy_session, *args, input=safe_product, output="f", **kwargs)
    assert proc.returncode != 0
    assert message in proc.stderr


def test_projected_project_fails(tmp_path, safe_product):
    project = tmp_path / "ll"
    gs.create_project(project, epsg=4326)
    with gs.setup.init(project, env=os.environ.copy()) as session:
        proc = run_module(session, input=safe_product, output="p")
    assert proc.returncode != 0
    assert "XY" in proc.stderr
