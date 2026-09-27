"""Tests of r.in.s1slc against a real Sentinel-1 IW SLC product.

The outputs are checked against an independent reading of the product
(measurement TIFF through GDAL, annotation XML through ElementTree) and an
independent implementation of the ESA SNAP (microwave toolbox) algorithms:
TOPSARDeburstOp line selection, Sentinel1Calibrator and
Sentinel1RemoveThermalNoiseOp for TOPS SLC. No code of the module is reused.

The tests are skipped unless R_IN_S1SLC_PRODUCT gives the path of an IW SLC
product (.zip archive, .SAFE directory or manifest.safe). Three consecutive
bursts around the middle of the sub-swath are imported; this needs about
5 GB of disk under the pytest base temporary directory and 4 GB of memory:

    R_IN_S1SLC_PRODUCT=/data/S1A_IW_SLC__1SDV_....zip \\
        grass --tmp-project XY --exec python3 -m pytest tests/real_product_test.py \\
        --basetemp=/data/tmp/pytest

R_IN_S1SLC_SWATH (default IW2) and R_IN_S1SLC_POLARIZATION (default VV)
select the image.
"""

import json
import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import grass.script as gs

from conftest import run_module

PRODUCT = os.environ.get("R_IN_S1SLC_PRODUCT")
SWATH = os.environ.get("R_IN_S1SLC_SWATH", "IW2").lower()
POL = os.environ.get("R_IN_S1SLC_POLARIZATION", "VV").lower()

pytestmark = pytest.mark.skipif(
    not PRODUCT, reason="R_IN_S1SLC_PRODUCT is not set (path of a real IW SLC product)"
)

# Float32 outputs against float64 references.
RTOL = 1e-5


def tsec(text):
    """Seconds since midnight of the acquisition day, to the nanosecond."""
    t = np.datetime64(text.strip().rstrip("Z"), "ns")
    return (t - t.astype("datetime64[D]")).astype(np.int64) * 1e-9


def floats(text):
    return np.array(text.split(), dtype=np.float64)


class Reference:
    """Direct reading of one sub-swath/polarization of the product."""

    def __init__(self, path):
        from osgeo import gdal

        gdal.UseExceptions()
        self.gdal = gdal
        path = os.path.abspath(path)
        if path.lower().endswith(".zip"):
            safe = os.path.basename(path)[: -len(".zip")] + ".SAFE"
            self.root = "/vsizip/" + path + "/" + safe
        else:
            self.root = path[: -len("/manifest.safe")] if path.endswith("manifest.safe") else path
        name = self._find("annotation", "")
        self.ann = ET.fromstring(self.read("annotation/" + name))
        self.cal = ET.fromstring(self.read("annotation/calibration/calibration-" + name))
        self.noise = ET.fromstring(self.read("annotation/calibration/noise-" + name))
        info = self.ann.find("imageAnnotation/imageInformation")
        self.lpb = int(self.ann.findtext("swathTiming/linesPerBurst"))
        self.dt = float(info.findtext("azimuthTimeInterval"))
        self.nlines = int(info.findtext("numberOfLines"))
        self.nsamples = int(info.findtext("numberOfSamples"))
        self.first_time = tsec(info.findtext("productFirstLineUtcTime"))
        self.last_time = tsec(info.findtext("productLastLineUtcTime"))
        self.start_time = tsec(self.ann.findtext("adsHeader/startTime"))
        self.bursts = []
        for b in self.ann.findall("swathTiming/burstList/burst"):
            fvs = np.array(b.findtext("firstValidSample").split(), dtype=int)
            lvs = np.array(b.findtext("lastValidSample").split(), dtype=int)
            valid = np.flatnonzero(fvs != -1)
            self.bursts.append(
                {"t": tsec(b.findtext("azimuthTime")), "fvs": fvs, "lvs": lvs, "fvl": valid[0], "lvl": valid[-1]}
            )
        self.fvs = np.concatenate([b["fvs"] for b in self.bursts])
        self.lvs = np.concatenate([b["lvs"] for b in self.bursts])
        self.ds = gdal.Open(self.root + "/measurement/" + name[: -len(".xml")] + ".tiff")
        assert (self.ds.RasterXSize, self.ds.RasterYSize) == (self.nsamples, self.nlines)

    def _find(self, directory, prefix):
        for name in self.gdal.ReadDir(self.root + "/" + directory):
            if name.startswith(prefix + "s1") and f"-{SWATH}-slc-{POL}-" in name and name.endswith(".xml"):
                return name
        pytest.skip(f"No {SWATH.upper()} {POL.upper()} image in {PRODUCT}")

    def read(self, rel):
        gdal = self.gdal
        handle = gdal.VSIFOpenL(self.root + "/" + rel, "rb")
        gdal.VSIFSeekL(handle, 0, 2)
        size = gdal.VSIFTellL(handle)
        gdal.VSIFSeekL(handle, 0, 0)
        data = gdal.VSIFReadL(1, size, handle)
        gdal.VSIFCloseL(handle)
        return data

    def raw(self, lines):
        """Complex digital numbers of stacked source lines."""
        lines = np.asarray(lines)
        out = np.zeros((lines.size, self.nsamples), dtype=np.complex64)
        start = 0
        while start < lines.size:
            stop = start + 1
            while stop < lines.size and lines[stop] == lines[stop - 1] + 1:
                stop += 1
            if lines[start] >= 0:
                out[start:stop] = self.ds.GetRasterBand(1).ReadAsArray(
                    0, int(lines[start]), self.nsamples, stop - start
                )
            start = stop
        return out

    def valid(self, lines):
        """Annotated valid samples (firstValidSample..lastValidSample) of lines."""
        lines = np.asarray(lines)
        cols = np.arange(self.nsamples)
        fvs, lvs = self.fvs[lines], self.lvs[lines]
        return (cols >= fvs[:, None]) & (cols <= lvs[:, None]) & ((fvs >= 0) & (lines >= 0))[:, None]

    def snap_deburst_lines(self, selected):
        """Source line of every output row, as TOPSARDeburstOp selects it.

        After TOPSAR-Split the first line time is that of the first selected
        burst; a burst covers [first line, last line) and in an overlap the
        later burst is used after the middle of the two burst line spans.
        """
        dt = self.dt
        t_first = np.array([self.bursts[k]["t"] for k in selected])
        t_last = t_first + (self.lpb - 1) * dt
        t0 = t_first[0]
        # SNAP truncates; the offset only absorbs floating-point noise.
        ymin = int((self.bursts[selected[0]]["fvl"] * dt) / dt + 1e-6)
        ymax = int((t_first[-1] + self.bursts[selected[-1]]["lvl"] * dt - t0) / dt + 1e-6)
        src = []
        for y in range(ymin, ymax + 1):
            t = t0 + y * dt
            hits = [
                (j, selected[j] * self.lpb + int((t - t_first[j]) / dt + 0.5))
                for j in range(len(selected))
                if t_first[j] <= t < t_last[j]
            ]
            if len(hits) >= 2 and t > 0.5 * (t_last[hits[0][0]] + t_first[hits[1][0]]):
                src.append(hits[1][1])
            else:
                src.append(hits[0][1] if hits else -1)
        return np.array(src), t0 + ymin * dt

    def calibration_lut(self, tag, burst):
        """Calibration vector tag interpolated on the lines of one burst.

        Bilinear in pixel and azimuth time as Sentinel1Calibrator, but with
        the true zero-Doppler time of each burst line, which is what the
        module documents (SNAP uses a pseudo-continuous time over the burst
        stack; the difference is below 1e-4 dB).
        """
        vectors = self.cal.findall("calibrationVectorList/calibrationVector")
        times = np.array([tsec(v.findtext("azimuthTime")) for v in vectors])
        cols = np.arange(self.nsamples)
        values = np.stack([np.interp(cols, floats(v.findtext("pixel")), floats(v.findtext(tag))) for v in vectors])
        t = self.bursts[burst]["t"] + np.arange(self.lpb) * self.dt
        idx = np.clip(np.searchsorted(times, t, side="right") - 1, 0, len(vectors) - 2)
        mu = (t - times[idx]) / (times[idx + 1] - times[idx])
        return (1 - mu)[:, None] * values[idx] + mu[:, None] * values[idx + 1]

    def noise_power(self, burst):
        """Thermal noise power of one burst (populateNoiseMatrixForTOPSSLC)."""
        rvecs = self.noise.findall("noiseRangeVectorList/noiseRangeVector")
        rtimes = np.array([tsec(v.findtext("azimuthTime")) for v in rvecs])
        skip = int(np.argmin(np.abs(rtimes - self.start_time)))
        rv = rvecs[burst + skip]
        cols = np.arange(self.nsamples)
        rng = np.interp(cols, floats(rv.findtext("pixel")), floats(rv.findtext("noiseRangeLut")))
        az = self.noise.find("noiseAzimuthVectorList/noiseAzimuthVector")
        lines = burst * self.lpb + np.arange(self.lpb)
        azi = np.interp(lines, floats(az.findtext("line")), floats(az.findtext("noiseAzimuthLut")))
        return azi[:, None] * rng[None, :]


@pytest.fixture(scope="module")
def ref():
    return Reference(PRODUCT)


@pytest.fixture(scope="module")
def imported(ref, tmp_path_factory):
    """Import the reference bursts once in an XY project; yield (session, bursts)."""
    middle = len(ref.bursts) // 2
    selected = [middle - 1, middle, middle + 1]
    project = tmp_path_factory.mktemp("real") / "xy"
    gs.create_project(project)
    common = {"input": PRODUCT, "swath": SWATH.upper(), "polarization": POL.upper()}
    single = str(middle + 1)
    runs = [
        (
            (),
            {
                "output": "deb",
                "bursts": f"{selected[0] + 1}-{selected[-1] + 1}",
                "measure": "complex,amplitude,intensity,db,phase",
                "geometry": "latitude,longitude,height,incidence_angle,elevation_angle",
            },
        ),
        (("b",), {"output": "raw", "bursts": single, "measure": "complex"}),
        (("b",), {"output": "s0", "bursts": single, "measure": "complex,intensity", "calibration": "sigma0"}),
        (("b",), {"output": "b0", "bursts": single, "measure": "intensity", "calibration": "beta0"}),
        (("b",), {"output": "g0", "bursts": single, "measure": "intensity", "calibration": "gamma0"}),
        (("b", "n"), {"output": "s0n", "bursts": single, "measure": "intensity,db", "calibration": "sigma0"}),
        (("b", "n"), {"output": "n", "bursts": single, "measure": "intensity"}),
    ]
    with gs.setup.init(project, env=os.environ.copy()) as session:
        for flags, kwargs in runs:
            proc = run_module(session, *flags, **common, **kwargs)
            assert proc.returncode == 0, proc.stderr
        yield session, selected


def read(session, name):
    """Map as float64 with nulls as NaN (CELL maps included)."""
    import grass.script.array as garray

    gs.run_command("g.region", raster=name, env=session.env)
    return np.asarray(garray.array(name, null=np.nan, dtype=np.float64, env=session.env))


def mapset_path(session, *parts):
    env = gs.gisenv(env=session.env)
    return os.path.join(env["GISDBASE"], env["LOCATION_NAME"], env["MAPSET"], *parts)


def tag(name, burst):
    return f"{name}_{SWATH}_{POL}_b{burst + 1:02d}"


def assert_close(actual, expected, mask):
    """Same nulls as mask, values within RTOL, exact zeros where expected is 0."""
    assert np.array_equal(~np.isnan(actual), mask)
    nonzero = mask & (expected != 0)
    np.testing.assert_allclose(actual[nonzero], expected[nonzero], rtol=RTOL)
    assert (actual[mask & (expected == 0)] == 0).all()


def test_burst_digital_numbers(ref, imported):
    session, selected = imported
    burst = selected[1]
    lines = burst * ref.lpb + np.arange(ref.lpb)
    z = ref.raw(lines)
    mask = ref.valid(lines)
    for comp, expected in (("i", z.real), ("q", z.imag)):
        actual = read(session, f"{tag('raw', burst)}_{comp}")
        assert gs.raster_info(f"{tag('raw', burst)}_{comp}", env=session.env)["datatype"] == "CELL"
        assert np.array_equal(~np.isnan(actual), mask)
        assert np.array_equal(actual[mask], expected[mask])


def test_deburst_matches_snap(ref, imported):
    session, selected = imported
    src, _t0 = ref.snap_deburst_lines(selected)
    with open(mapset_path(session, "cell_misc", f"deb_{SWATH}_{POL}_i", "description.json")) as fd:
        geometry = json.load(fd)["raster_geometry"]
    rows = np.concatenate(
        [(s["burst"] - 1) * ref.lpb + np.arange(s["first_line"], s["first_line"] + s["rows"]) for s in geometry["segments"]]
    )
    assert geometry["rows"] == src.size
    assert np.array_equal(rows, src)
    z = ref.raw(src)
    mask = ref.valid(src)
    for comp, expected in (("i", z.real), ("q", z.imag)):
        actual = read(session, f"deb_{SWATH}_{POL}_{comp}")
        assert np.array_equal(~np.isnan(actual), mask)
        assert np.array_equal(actual[mask], expected[mask])


def test_derived_measures(imported):
    session, _selected = imported
    i = read(session, f"deb_{SWATH}_{POL}_i")
    q = read(session, f"deb_{SWATH}_{POL}_q")
    mask = ~np.isnan(i)
    power = i**2 + q**2
    with np.errstate(divide="ignore"):
        db = 10 * np.log10(power)
    assert_close(read(session, f"deb_{SWATH}_{POL}_intensity"), power, mask)
    assert_close(read(session, f"deb_{SWATH}_{POL}_amplitude"), np.sqrt(power), mask)
    assert_close(read(session, f"deb_{SWATH}_{POL}_db"), db, mask & (power > 0))
    phase = read(session, f"deb_{SWATH}_{POL}_phase")
    assert np.array_equal(~np.isnan(phase), mask)
    np.testing.assert_allclose(phase[mask], np.arctan2(q, i)[mask], rtol=0, atol=1e-6)


@pytest.mark.parametrize(
    ("output", "lut_tag"), [("s0", "sigmaNought"), ("b0", "betaNought"), ("g0", "gamma")]
)
def test_calibrated_intensity(ref, imported, output, lut_tag):
    session, selected = imported
    burst = selected[1]
    lines = burst * ref.lpb + np.arange(ref.lpb)
    power = np.abs(ref.raw(lines).astype(np.complex128)) ** 2
    lut = ref.calibration_lut(lut_tag, burst)
    actual = read(session, f"{tag(output, burst)}_intensity")
    assert_close(actual, power / lut**2, ref.valid(lines))


def test_calibrated_complex(ref, imported):
    session, selected = imported
    burst = selected[1]
    lines = burst * ref.lpb + np.arange(ref.lpb)
    z = ref.raw(lines).astype(np.complex128)
    lut = ref.calibration_lut("sigmaNought", burst)
    mask = ref.valid(lines)
    assert_close(read(session, f"{tag('s0', burst)}_i"), z.real / lut, mask)
    assert_close(read(session, f"{tag('s0', burst)}_q"), z.imag / lut, mask)


def test_thermal_noise_removal(ref, imported):
    session, selected = imported
    burst = selected[1]
    lines = burst * ref.lpb + np.arange(ref.lpb)
    power = np.abs(ref.raw(lines).astype(np.complex128)) ** 2
    noise = ref.noise_power(burst)
    # Zero noise marks missing noise data and leaves the power unchanged.
    denoised = np.where(noise == 0, power, np.maximum(power - noise, 0.0))
    mask = ref.valid(lines)
    actual = read(session, f"{tag('n', burst)}_intensity")
    assert_close(actual, denoised, mask)
    assert (actual[mask] >= 0).all()
    sigma0 = denoised / ref.calibration_lut("sigmaNought", burst) ** 2
    actual = read(session, f"{tag('s0n', burst)}_intensity")
    assert_close(actual, sigma0, mask)
    db = read(session, f"{tag('s0n', burst)}_db")
    positive = mask & (sigma0 > 0)
    assert np.array_equal(~np.isnan(db), ~np.isnan(actual) & (actual > 0))
    np.testing.assert_allclose(db[positive], 10 * np.log10(sigma0[positive]), rtol=0, atol=1e-4)


@pytest.mark.parametrize(
    ("layer", "grid_tag"),
    [
        ("latitude", "latitude"),
        ("longitude", "longitude"),
        ("height", "height"),
        ("incidence_angle", "incidenceAngle"),
        ("elevation_angle", "elevationAngle"),
    ],
)
def test_geometry_at_grid_nodes(ref, imported, layer, grid_tag):
    session, selected = imported
    _src, t0 = ref.snap_deburst_lines(selected)
    actual = read(session, f"deb_{SWATH}_{layer}")
    checked = 0
    for point in ref.ann.findall("geolocationGrid/geolocationGridPointList/geolocationGridPoint"):
        row = (tsec(point.findtext("azimuthTime")) - t0) / ref.dt
        col = int(point.findtext("pixel"))
        # Nodes falling on an output row are exact interpolation knots.
        if abs(row - round(row)) < 0.05 and 0 <= round(row) < actual.shape[0] and col < ref.nsamples:
            assert abs(actual[round(row), col] - float(point.findtext(grid_tag))) < 1e-4
            checked += 1
    assert checked > 0


@pytest.mark.parametrize("group", ["deb", "raw"])
def test_gcps_distinct(imported, group):
    session, selected = imported
    name = f"deb_{SWATH}_{POL}" if group == "deb" else tag("raw", selected[1])
    points = np.loadtxt(mapset_path(session, "group", name, "POINTS"), comments="#")
    assert np.unique(points[:, :2], axis=0).shape[0] == points.shape[0]
    assert np.diff(np.unique(points[:, 1])).min() > 0.5
