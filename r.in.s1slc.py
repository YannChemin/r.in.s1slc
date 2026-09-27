#!/usr/bin/env python3

# MODULE:    r.in.s1slc
# AUTHOR(S): Yann Chemin
# PURPOSE:   Imports Sentinel-1 TOPS (IW/EW) Single Look Complex (SLC)
#            SAFE products in radar geometry, with optional debursting,
#            radiometric calibration, thermal noise removal, geolocation
#            layers, ground control points and full product metadata.
# COPYRIGHT: (C) 2026 by Yann Chemin and the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later

# %module
# % description: Imports Sentinel-1 SLC (Single Look Complex) SAFE products in radar geometry.
# % keyword: raster
# % keyword: import
# % keyword: imagery
# % keyword: SAR
# % keyword: radar
# % keyword: Sentinel-1
# % keyword: SLC
# % keyword: TOPS
# % keyword: metadata
# %end

# %option G_OPT_F_BIN_INPUT
# % key: input
# % label: Sentinel-1 SLC product
# % description: SAFE product as .zip archive, .SAFE directory or its manifest.safe file
# %end

# %option G_OPT_R_BASENAME_OUTPUT
# % required: yes
# % description: Name prefix for output raster maps and imagery groups
# %end

# %option
# % key: swath
# % type: string
# % required: no
# % multiple: yes
# % options: IW1,IW2,IW3,EW1,EW2,EW3,EW4,EW5,S1,S2,S3,S4,S5,S6
# % label: Sub-swaths to import
# % description: Default: all sub-swaths present in the product
# % guisection: Subset
# %end

# %option
# % key: polarization
# % type: string
# % required: no
# % multiple: yes
# % options: VV,VH,HH,HV
# % label: Polarizations to import
# % description: Default: all polarizations present in the product
# % guisection: Subset
# %end

# %option
# % key: bursts
# % type: string
# % required: no
# % multiple: yes
# % key_desc: range
# % label: Bursts to import (1-based indices or ranges, e.g. 3-5)
# % description: Applied to each selected sub-swath
# % guisection: Subset
# %end

# %option
# % key: bbox
# % type: double
# % required: no
# % multiple: yes
# % key_desc: west,south,east,north
# % label: Import only bursts intersecting this WGS84 longitude/latitude box
# % description: Burst footprints come from the annotation geolocation grid
# % guisection: Subset
# %end

# %option
# % key: measure
# % type: string
# % required: no
# % multiple: yes
# % options: complex,amplitude,intensity,db,phase
# % answer: complex
# % label: Quantities to write for each sub-swath and polarization
# % description: complex writes two maps (_i, _q)
# % descriptions: complex;In-phase and quadrature components (maps _i and _q);amplitude;Amplitude |z|;intensity;Intensity |z|^2 (linear power);db;Intensity in decibels 10*log10(|z|^2);phase;Interferometric phase atan2(q, i) in radians
# % guisection: Output
# %end

# %option
# % key: calibration
# % type: string
# % required: no
# % multiple: no
# % options: none,sigma0,beta0,gamma0
# % answer: none
# % label: Radiometric calibration applied to the measures
# % description: Uses the calibration LUT of the product (ESA S1-TN-MDA-52-7448)
# % descriptions: none;Digital numbers;sigma0;Sigma nought;beta0;Beta nought;gamma0;Gamma (ellipsoid)
# % guisection: Output
# %end

# %option
# % key: geometry
# % type: string
# % required: no
# % multiple: yes
# % options: latitude,longitude,height,incidence_angle,elevation_angle
# % label: Geolocation layers to write in radar geometry
# % description: Bilinear interpolation of the annotation geolocation grid
# % guisection: Output
# %end

# %option
# % key: target
# % type: string
# % required: no
# % multiple: no
# % key_desc: name
# % label: Name of the target project for the ground control points
# % description: GCPs are reprojected to its CRS and it becomes the imagery group target (default: GCPs in WGS84 longitude/latitude, no target)
# % guisection: Output
# %end

# %flag
# % key: b
# % label: Keep bursts separate (do not deburst)
# % description: Writes one map per burst with the full burst line count
# %end

# %flag
# % key: n
# % label: Remove thermal noise
# % description: Subtracts the annotated noise power; not defined for the complex and phase measures
# %end

# %flag
# % key: l
# % description: List sub-swaths and bursts with their footprints and exit
# % suppress_required: yes
# %end

# %rules
# % exclusive: bursts,bbox
# %end

import atexit
import json
import os
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

import grass.script as gs

# Output rows processed at once; bounds memory to a few hundred MB per chunk.
CHUNK_ROWS = 512
INT_NULL = -2147483648

GEOMETRY_TAGS = {
    "latitude": ("latitude", "degrees"),
    "longitude": ("longitude", "degrees"),
    "height": ("height", "meters"),
    "incidence_angle": ("incidenceAngle", "degrees"),
    "elevation_angle": ("elevationAngle", "degrees"),
}

TMP_FILES = []


def cleanup():
    for path in TMP_FILES:
        try:
            os.remove(path)
        except OSError:
            pass


def local(tag):
    """Return the tag name without its XML namespace."""
    return tag.rsplit("}", 1)[-1]


def find_local(root, name):
    """Return the first element whose local tag name is name, or None."""
    for elem in root.iter():
        if local(elem.tag) == name:
            return elem
    return None


def text_local(root, name, default=None):
    elem = find_local(root, name)
    if elem is None or elem.text is None:
        return default
    return elem.text.strip()


def parse_time(text):
    """Parse an annotation UTC time with any number of fractional digits."""
    text = text.strip().rstrip("Z")
    if "." in text:
        head, frac = text.split(".")
        frac = (frac + "000000")[:6]
    else:
        head, frac = text, "000000"
    return datetime.strptime(head, "%Y-%m-%dT%H:%M:%S").replace(
        microsecond=int(frac)
    )


class TimeAxis:
    """Seconds since a reference epoch, exact to the microsecond."""

    def __init__(self, ref):
        self.ref = ref

    def sec(self, text_or_dt):
        dt = text_or_dt if isinstance(text_or_dt, datetime) else parse_time(text_or_dt)
        return (dt - self.ref) // timedelta(microseconds=1) * 1e-6

    def iso(self, sec):
        return (self.ref + timedelta(microseconds=round(sec * 1e6))).isoformat()

    def datetime(self, sec):
        return self.ref + timedelta(microseconds=round(sec * 1e6))


def floats(text):
    return [float(v) for v in text.split()]


def ints(text):
    return [int(v) for v in text.split()]


class SafeProduct:
    """Access to a SAFE product stored as a directory or a zip archive."""

    def __init__(self, path):
        from osgeo import gdal

        self.gdal = gdal
        gdal.UseExceptions()
        self.root = self._find_root(path)
        self.name = os.path.basename(self.root)
        self.manifest = ET.fromstring(self.read("manifest.safe"))

    def _find_root(self, path):
        gdal = self.gdal
        path = os.path.abspath(path)
        if os.path.isfile(path) and path.lower().endswith(".zip"):
            vsi = "/vsizip/" + path
            for entry in gdal.ReadDirRecursive(vsi) or []:
                parts = entry.strip("/").split("/")
                if len(parts) == 2 and parts[1] == "manifest.safe":
                    return vsi + "/" + parts[0]
            gs.fatal(_("No SAFE manifest found in <{}>").format(path))
        if os.path.basename(path) == "manifest.safe":
            path = os.path.dirname(path)
        if os.path.isdir(path) and os.path.isfile(os.path.join(path, "manifest.safe")):
            return path
        gs.fatal(
            _(
                "<{}> is not a SAFE product (zip archive, .SAFE directory or "
                "manifest.safe expected)"
            ).format(path)
        )

    def path(self, rel):
        return self.root + "/" + rel

    def read(self, rel):
        gdal = self.gdal
        handle = gdal.VSIFOpenL(self.path(rel), "rb")
        if handle is None:
            gs.fatal(_("Unable to read <{}>").format(self.path(rel)))
        try:
            gdal.VSIFSeekL(handle, 0, 2)
            size = gdal.VSIFTellL(handle)
            gdal.VSIFSeekL(handle, 0, 0)
            return gdal.VSIFReadL(1, size, handle)
        finally:
            gdal.VSIFCloseL(handle)

    def listdir(self, rel):
        return sorted(self.gdal.ReadDir(self.path(rel)) or [])

    def exists(self, rel):
        return self.gdal.VSIStatL(self.path(rel)) is not None

    def info(self):
        """Product-level metadata from manifest.safe."""
        m = self.manifest
        info = {
            "product_name": self.name,
            "mission": "{}{}".format(
                text_local(m, "familyName", ""), text_local(m, "number", "")
            ),
            "product_type": text_local(m, "productType"),
            "product_class": text_local(m, "productClass"),
            "mode": text_local(m, "mode"),
            "polarizations": [
                e.text.strip()
                for e in m.iter()
                if local(e.tag) == "transmitterReceiverPolarisation"
            ],
            "start_time": text_local(m, "startTime"),
            "stop_time": text_local(m, "stopTime"),
            "pass": text_local(m, "pass"),
            "timeliness": text_local(m, "productTimelinessCategory"),
            "mission_data_take_id": text_local(m, "missionDataTakeID"),
            "slice_number": text_local(m, "sliceNumber"),
            "total_slices": text_local(m, "totalSlices"),
            "cycle_number": text_local(m, "cycleNumber"),
            "footprint": text_local(m, "coordinates"),
        }
        for elem in m.iter():
            tag = local(elem.tag)
            if tag == "orbitNumber":
                info["absolute_orbit_" + elem.get("type", "start")] = elem.text
            elif tag == "relativeOrbitNumber":
                info["relative_orbit_" + elem.get("type", "start")] = elem.text
            elif tag == "software" and "IPF" in (elem.get("name") or ""):
                info["ipf_version"] = elem.get("version")
        return info


class RowLut:
    """Values given on knot rows (azimuth time) and knot columns (pixel).

    Each knot row is first interpolated linearly along range to every image
    column; values are then interpolated linearly in azimuth time, with
    linear extrapolation beyond the first and last knots (as SNAP does for
    the calibration vectors).
    """

    def __init__(self, knot_times, knot_pixels, knot_values, ncols):
        import numpy as np

        self.np = np
        order = np.argsort(knot_times)
        self.times = np.asarray(knot_times, dtype=np.float64)[order]
        cols = np.arange(ncols, dtype=np.float64)
        self.values = np.stack(
            [
                np.interp(cols, knot_pixels[k], knot_values[k]).astype(np.float64)
                for k in order
            ]
        )

    def weights(self, times):
        np = self.np
        times = np.asarray(times, dtype=np.float64)
        if len(self.times) == 1:
            idx = np.zeros(times.shape, dtype=int)
            return idx, idx, np.zeros(times.shape)
        idx = np.clip(np.searchsorted(self.times, times, side="right") - 1, 0, len(self.times) - 2)
        mu = (times - self.times[idx]) / (self.times[idx + 1] - self.times[idx])
        return idx, idx + 1, mu

    def rows(self, times, cols=None):
        i0, i1, mu = self.weights(times)
        v0 = self.values[i0] if cols is None else self.values[i0][:, cols]
        v1 = self.values[i1] if cols is None else self.values[i1][:, cols]
        return (1.0 - mu)[:, None] * v0 + mu[:, None] * v1


class Swath:
    """One sub-swath/polarization image: annotation, calibration and noise."""

    def __init__(self, product, ann_name, axis):
        import numpy as np

        self.np = np
        self.product = product
        self.axis = axis
        self.ann_name = ann_name
        self.stem = ann_name[: -len(".xml")]
        self.tiff = "measurement/" + self.stem + ".tiff"
        ann = ET.fromstring(product.read("annotation/" + ann_name))
        self.ann = ann

        header = ann.find("adsHeader")
        self.mission = header.findtext("missionId")
        self.product_type = header.findtext("productType")
        self.pol = header.findtext("polarisation")
        self.mode = header.findtext("mode")
        self.swath = header.findtext("swath")
        self.start = axis.sec(header.findtext("startTime"))
        self.stop = axis.sec(header.findtext("stopTime"))

        pinfo = ann.find("generalAnnotation/productInformation")
        iinfo = ann.find("imageAnnotation/imageInformation")
        self.radar_frequency = float(pinfo.findtext("radarFrequency"))
        self.range_sampling_rate = float(pinfo.findtext("rangeSamplingRate"))
        self.azimuth_steering_rate = float(pinfo.findtext("azimuthSteeringRate"))
        self.heading = float(pinfo.findtext("platformHeading"))
        self.pass_direction = pinfo.findtext("pass")
        self.nlines = int(iinfo.findtext("numberOfLines"))
        self.nsamples = int(iinfo.findtext("numberOfSamples"))
        self.dt = float(iinfo.findtext("azimuthTimeInterval"))
        self.first_line_time = axis.sec(iinfo.findtext("productFirstLineUtcTime"))
        self.last_line_time = axis.sec(iinfo.findtext("productLastLineUtcTime"))
        self.slant_range_time = float(iinfo.findtext("slantRangeTime"))
        self.range_pixel_spacing = float(iinfo.findtext("rangePixelSpacing"))
        self.azimuth_pixel_spacing = float(iinfo.findtext("azimuthPixelSpacing"))
        self.incidence_mid = float(iinfo.findtext("incidenceAngleMidSwath"))
        self.pixel_value = iinfo.findtext("pixelValue")
        prf = ann.findtext("generalAnnotation/downlinkInformationList/downlinkInformation/prf")
        self.prf = float(prf) if prf else None

        self._read_bursts()
        self._read_geolocation()
        self._cal = None
        self._noise = None

    def _read_bursts(self):
        np = self.np
        timing = self.ann.find("swathTiming")
        blist = timing.findall("burstList/burst") if timing is not None else []
        self.bursts = []
        if not blist:
            # Stripmap-like product: one pseudo-burst spanning the image.
            self.lpb = self.nlines
            full = np.full(self.nlines, 0, dtype=np.int32)
            self.bursts.append(
                {
                    "t": self.first_line_time,
                    "azimuth_time": self.axis.iso(self.first_line_time),
                    "burst_id": None,
                    "first_valid_sample": full,
                    "last_valid_sample": full + self.nsamples - 1,
                    "first_valid_line": 0,
                    "last_valid_line": self.nlines - 1,
                }
            )
        else:
            self.lpb = int(timing.findtext("linesPerBurst"))
            for elem in blist:
                fvs = np.array(ints(elem.findtext("firstValidSample")), dtype=np.int32)
                lvs = np.array(ints(elem.findtext("lastValidSample")), dtype=np.int32)
                valid = np.flatnonzero(fvs >= 0)
                bid = elem.find("burstId")
                self.bursts.append(
                    {
                        "t": self.axis.sec(elem.findtext("azimuthTime")),
                        "azimuth_time": elem.findtext("azimuthTime"),
                        "azimuth_anx_time": float(elem.findtext("azimuthAnxTime")),
                        "sensing_time": elem.findtext("sensingTime"),
                        "byte_offset": int(elem.findtext("byteOffset")),
                        "burst_id": int(bid.text) if bid is not None else None,
                        "burst_id_absolute": bid.get("absolute") if bid is not None else None,
                        "first_valid_sample": fvs,
                        "last_valid_sample": lvs,
                        "first_valid_line": int(valid[0]) if valid.size else -1,
                        "last_valid_line": int(valid[-1]) if valid.size else -1,
                    }
                )
        fvs = np.concatenate([b["first_valid_sample"] for b in self.bursts])
        lvs = np.concatenate([b["last_valid_sample"] for b in self.bursts])
        if fvs.size != self.nlines:
            gs.fatal(
                _("Burst valid-sample arrays of <{}> cover {} lines, image has {}").format(
                    self.ann_name, fvs.size, self.nlines
                )
            )
        self.first_valid_sample = fvs
        self.last_valid_sample = lvs

    def _read_geolocation(self):
        np = self.np
        points = self.ann.findall("geolocationGrid/geolocationGridPointList/geolocationGridPoint")
        rows = {}
        for p in points:
            line = int(p.findtext("line"))
            rows.setdefault(line, []).append(p)
        self.grid_lines = sorted(rows)
        self.grid_times = []
        self.grid_pixels = []
        self.grid = {tag: [] for tag, _unit in GEOMETRY_TAGS.values()}
        for line in self.grid_lines:
            pts = sorted(rows[line], key=lambda p: int(p.findtext("pixel")))
            self.grid_times.append(
                float(np.mean([self.axis.sec(p.findtext("azimuthTime")) for p in pts]))
            )
            self.grid_pixels.append([int(p.findtext("pixel")) for p in pts])
            for tag in self.grid:
                self.grid[tag].append([float(p.findtext(tag)) for p in pts])
        self._geo_luts = {}

    def geo_lut(self, tag):
        if tag not in self._geo_luts:
            self._geo_luts[tag] = RowLut(
                self.grid_times, self.grid_pixels, self.grid[tag], self.nsamples
            )
        return self._geo_luts[tag]

    def burst_line_time(self, k, line):
        return self.bursts[k]["t"] + line * self.dt

    def burst_footprint(self, k):
        """Longitude/latitude extent of the valid area of burst k."""
        np = self.np
        b = self.bursts[k]
        t = self.burst_line_time(k, np.linspace(b["first_valid_line"], b["last_valid_line"], 3))
        fvs = b["first_valid_sample"][b["first_valid_sample"] >= 0]
        lvs = b["last_valid_sample"][b["last_valid_sample"] >= 0]
        cols = np.linspace(fvs.min(), lvs.max(), 5).round().astype(int)
        lat = self.geo_lut("latitude").rows(t, cols)
        lon = self.geo_lut("longitude").rows(t, cols)
        return lon.min(), lat.min(), lon.max(), lat.max()

    # Calibration.
    def calibration_lut(self, kind):
        """RowLut of the calibration vector kind (sigmaNought, betaNought, gamma, dn)."""
        if self._cal is None:
            self._cal = ET.fromstring(
                self.product.read("annotation/calibration/calibration-" + self.ann_name)
            )
        vectors = self._cal.findall("calibrationVectorList/calibrationVector")
        if len(vectors) < 1:
            gs.fatal(_("No calibration vectors for <{}>").format(self.ann_name))
        times = [self.axis.sec(v.findtext("azimuthTime")) for v in vectors]
        pixels = [ints(v.findtext("pixel")) for v in vectors]
        values = [floats(v.findtext(kind)) for v in vectors]
        return RowLut(times, pixels, values, self.nsamples)

    def calibration_constant(self):
        if self._cal is None:
            self.calibration_lut("dn")
        return float(self._cal.findtext("calibrationInformation/absoluteCalibrationConstant"))

    # Thermal noise.
    def noise_model(self):
        """Return (range vectors per burst, azimuth blocks) of the noise annotation.

        Range vectors are interpolated to every column; each burst uses the
        vector closest in time to its first line, as SNAP does for TOPS SLC.
        Azimuth blocks (IPF >= 2.9) are (first_line, last_line, first_sample,
        last_sample, lines, lut) in the burst-stacked line numbering.
        """
        np = self.np
        if self._noise is not None:
            return self._noise
        name = "annotation/calibration/noise-" + self.ann_name
        if not self.product.exists(name):
            gs.fatal(_("Noise annotation <{}> not found").format(name))
        root = ET.fromstring(self.product.read(name))
        rvecs = root.findall("noiseRangeVectorList/noiseRangeVector")
        lut_tag = "noiseRangeLut"
        if not rvecs:
            # IPF < 2.9 products only have range noise vectors.
            rvecs = root.findall("noiseVectorList/noiseVector")
            lut_tag = "noiseLut"
        if not rvecs:
            gs.fatal(_("No noise vectors in <{}>").format(name))
        rtimes = np.array([self.axis.sec(v.findtext("azimuthTime")) for v in rvecs])
        cols = np.arange(self.nsamples, dtype=np.float64)
        per_burst = []
        for b in self.bursts:
            v = rvecs[int(np.argmin(np.abs(rtimes - b["t"])))]
            per_burst.append(
                np.interp(cols, ints(v.findtext("pixel")), floats(v.findtext(lut_tag)))
            )
        blocks = []
        for v in root.findall("noiseAzimuthVectorList/noiseAzimuthVector"):
            blocks.append(
                (
                    int(v.findtext("firstAzimuthLine")),
                    int(v.findtext("lastAzimuthLine")),
                    int(v.findtext("firstRangeSample")),
                    int(v.findtext("lastRangeSample")),
                    np.array(ints(v.findtext("line")), dtype=np.float64),
                    np.array(floats(v.findtext("noiseAzimuthLut")), dtype=np.float64),
                )
            )
        self._noise = (per_burst, blocks)
        return self._noise

    def noise_power(self, src_lines):
        """Noise power (DN^2) for stacked source lines, shape (rows, nsamples)."""
        np = self.np
        per_burst, blocks = self.noise_model()
        burst = np.clip(src_lines // self.lpb, 0, len(self.bursts) - 1)
        rng = np.stack([per_burst[k] for k in burst])
        if not blocks:
            return rng
        azi = np.ones_like(rng)
        for first_line, last_line, first_sample, last_sample, lines, lut in blocks:
            rows = (src_lines >= first_line) & (src_lines <= last_line)
            if not rows.any():
                continue
            values = np.interp(src_lines[rows], lines, lut) if lines.size > 1 else np.full(rows.sum(), lut[0])
            azi[rows, first_sample : last_sample + 1] = values[:, None]
        return rng * azi

    def metadata(self):
        """Annotation metadata needed downstream (e.g. for interferometry)."""
        ann = self.ann
        np = self.np

        def elem_dict(elem):
            return {local(c.tag): (c.text or "").strip() for c in elem if len(c) == 0}

        orbits = []
        for o in ann.findall("generalAnnotation/orbitList/orbit"):
            orbits.append(
                {
                    "time": o.findtext("time"),
                    "frame": o.findtext("frame"),
                    "position": [float(o.findtext("position/" + a)) for a in "xyz"],
                    "velocity": [float(o.findtext("velocity/" + a)) for a in "xyz"],
                }
            )
        dc = []
        for e in ann.findall("dopplerCentroid/dcEstimateList/dcEstimate"):
            dc.append(
                {
                    "azimuth_time": e.findtext("azimuthTime"),
                    "t0": float(e.findtext("t0")),
                    "geometry_dc_polynomial": floats(e.findtext("geometryDcPolynomial")),
                    "data_dc_polynomial": floats(e.findtext("dataDcPolynomial")),
                }
            )
        fm = []
        for e in ann.findall("generalAnnotation/azimuthFmRateList/azimuthFmRate"):
            poly = e.findtext("azimuthFmRatePolynomial")
            if poly is None:  # IPF < 2.8 layout
                poly = " ".join(e.findtext(c) for c in ("c0", "c1", "c2"))
            fm.append(
                {
                    "azimuth_time": e.findtext("azimuthTime"),
                    "t0": float(e.findtext("t0")),
                    "polynomial": floats(poly),
                }
            )
        proc = ann.find("imageAnnotation/processingInformation")
        swath_proc = proc.find("swathProcParamsList/swathProcParams") if proc is not None else None
        processing = elem_dict(proc) if proc is not None else {}
        if swath_proc is not None:
            for sub in ("rangeProcessing", "azimuthProcessing"):
                elem = swath_proc.find(sub)
                if elem is not None:
                    processing[sub] = elem_dict(elem)
        bursts = []
        for k, b in enumerate(self.bursts):
            bursts.append(
                {
                    "index": k + 1,
                    "burst_id": b.get("burst_id"),
                    "burst_id_absolute": b.get("burst_id_absolute"),
                    "azimuth_time": b["azimuth_time"],
                    "azimuth_anx_time": b.get("azimuth_anx_time"),
                    "sensing_time": b.get("sensing_time"),
                    "byte_offset": b.get("byte_offset"),
                    "first_valid_line": b["first_valid_line"],
                    "last_valid_line": b["last_valid_line"],
                    "first_valid_sample": b["first_valid_sample"].tolist(),
                    "last_valid_sample": b["last_valid_sample"].tolist(),
                }
            )
        grid = []
        for i, line in enumerate(self.grid_lines):
            for j, pixel in enumerate(self.grid_pixels[i]):
                point = {"line": line, "pixel": pixel, "azimuth_time": self.axis.iso(self.grid_times[i])}
                for tag in self.grid:
                    point[tag] = self.grid[tag][i][j]
                grid.append(point)
        return {
            "annotation": self.ann_name,
            "measurement": self.tiff,
            "mission": self.mission,
            "mode": self.mode,
            "swath": self.swath,
            "polarization": self.pol,
            "product_type": self.product_type,
            "pixel_value": self.pixel_value,
            "start_time": self.axis.iso(self.start),
            "stop_time": self.axis.iso(self.stop),
            "pass": self.pass_direction,
            "platform_heading": self.heading,
            "radar_frequency": self.radar_frequency,
            "wavelength": 299792458.0 / self.radar_frequency,
            "range_sampling_rate": self.range_sampling_rate,
            "azimuth_steering_rate": self.azimuth_steering_rate,
            "prf": self.prf,
            "number_of_lines": self.nlines,
            "number_of_samples": self.nsamples,
            "lines_per_burst": self.lpb,
            "azimuth_time_interval": self.dt,
            "product_first_line_utc_time": self.axis.iso(self.first_line_time),
            "product_last_line_utc_time": self.axis.iso(self.last_line_time),
            "slant_range_time": self.slant_range_time,
            "range_pixel_spacing": self.range_pixel_spacing,
            "azimuth_pixel_spacing": self.azimuth_pixel_spacing,
            "incidence_angle_mid_swath": self.incidence_mid,
            "terrain_height": float(
                np.mean([float(e.text) for e in ann.findall("generalAnnotation/terrainHeightList/terrainHeight/value")] or [0.0])
            ),
            "processing": processing,
            "bursts": bursts,
            "orbit_state_vectors": orbits,
            "doppler_centroid": dc,
            "azimuth_fm_rate": fm,
            "geolocation_grid": grid,
        }


class Layout:
    """Mapping of output rows to source lines and azimuth times."""

    def __init__(self, swath, bursts, deburst):
        np = swath.np
        self.swath = swath
        self.bursts = bursts
        self.deburst = deburst
        dt = swath.dt
        if not deburst:
            (k,) = bursts
            line = np.arange(swath.lpb)
            self.src = k * swath.lpb + line
            self.times = swath.burst_line_time(k, line)
        else:
            self._deburst(swath, bursts, dt)
        self.nrows = self.src.size
        self.ncols = swath.nsamples
        self.t0 = float(self.times[0])

    def _deburst(self, sw, bursts, dt):
        """Assign each output line to one burst like SNAP TOPSAR-Deburst.

        Output lines are regularly sampled in azimuth time between the first
        valid line of the first burst and the last valid line of the last
        burst. In the overlap of two bursts the line is taken from the earlier
        burst up to the middle of the overlap, from the later one after it.
        """
        np = sw.np
        first, last = sw.bursts[bursts[0]], sw.bursts[bursts[-1]]
        t0 = sw.burst_line_time(bursts[0], first["first_valid_line"])
        t1 = sw.burst_line_time(bursts[-1], last["last_valid_line"])
        n = int(round((t1 - t0) / dt)) + 1
        times = t0 + np.arange(n) * dt
        owner = np.full(n, bursts[0])
        for a, b in zip(bursts[:-1], bursts[1:]):
            mid = 0.5 * (sw.burst_line_time(a, sw.lpb - 1) + sw.bursts[b]["t"])
            owner[times > mid] = b
        src = np.full(n, -1, dtype=np.int64)
        for k in bursts:
            rows = np.flatnonzero(owner == k)
            line = np.floor((times[rows] - sw.bursts[k]["t"]) / dt + 0.5).astype(np.int64)
            ok = (line >= sw.bursts[k]["first_valid_line"]) & (line <= sw.bursts[k]["last_valid_line"])
            src[rows[ok]] = k * sw.lpb + line[ok]
        # A line rejected by its owner (invalid edge line) falls back to a
        # neighbouring selected burst when that one has it valid.
        for r in np.flatnonzero(src < 0):
            for k in bursts:
                line = int(np.floor((times[r] - sw.bursts[k]["t"]) / dt + 0.5))
                if sw.bursts[k]["first_valid_line"] <= line <= sw.bursts[k]["last_valid_line"]:
                    src[r] = k * sw.lpb + line
                    break
        self.src = src
        self.times = times

    def chunks(self):
        for r0 in range(0, self.nrows, CHUNK_ROWS):
            yield r0, min(r0 + CHUNK_ROWS, self.nrows)


class MeasurementReader:
    def __init__(self, product, swath):
        gdal = product.gdal
        self.np = swath.np
        self.gdal = gdal
        self.ds = gdal.Open(product.path(swath.tiff))
        if (self.ds.RasterXSize, self.ds.RasterYSize) != (swath.nsamples, swath.nlines):
            gs.fatal(
                _("Measurement <{}> is {}x{}, annotation says {}x{}").format(
                    swath.tiff,
                    self.ds.RasterXSize,
                    self.ds.RasterYSize,
                    swath.nsamples,
                    swath.nlines,
                )
            )

    def read(self, line0, nlines):
        """Return (i, q) float32 arrays for consecutive source lines."""
        np = self.np
        gdal = self.gdal
        ds = self.ds
        ncols = ds.RasterXSize
        if ds.RasterCount == 2:
            i = ds.GetRasterBand(1).ReadAsArray(0, line0, ncols, nlines).astype(np.float32)
            q = ds.GetRasterBand(2).ReadAsArray(0, line0, ncols, nlines).astype(np.float32)
            return i, q
        raw = ds.GetRasterBand(1).ReadRaster(
            0, line0, ncols, nlines, buf_type=gdal.GDT_CInt16
        )
        iq = np.frombuffer(raw, dtype=np.int16).reshape(nlines, ncols, 2)
        return iq[..., 0].astype(np.float32), iq[..., 1].astype(np.float32)


class BinaryRaster:
    """Raw row-major buffer written chunk by chunk, then imported by r.in.bin."""

    def __init__(self, name, is_float):
        self.name = name
        self.is_float = is_float
        self.path = gs.tempfile(create=False)
        TMP_FILES.append(self.path)
        self.handle = open(self.path, "wb")

    def write(self, array):
        np = sys.modules["numpy"]
        array.astype(np.float32 if self.is_float else np.int32).tofile(self.handle)

    def finish(self, nrows, ncols, title):
        self.handle.close()
        gs.run_command(
            "r.in.bin",
            flags="f" if self.is_float else "s",
            input=self.path,
            output=self.name,
            title=title,
            bytes=4,
            order="native",
            north=nrows,
            south=0,
            east=ncols,
            west=0,
            rows=nrows,
            cols=ncols,
            anull="nan" if self.is_float else INT_NULL,
            overwrite=gs.overwrite(),
            quiet=True,
        )
        os.remove(self.path)
        TMP_FILES.remove(self.path)


def parse_bursts(answers, nbursts):
    selected = set()
    for item in answers:
        item = item.strip()
        try:
            if "-" in item:
                a, b = (int(v) for v in item.split("-", 1))
            else:
                a = b = int(item)
        except ValueError:
            gs.fatal(_("Invalid burst selection <{}>").format(item))
        if a < 1 or b < a:
            gs.fatal(_("Invalid burst range <{}>").format(item))
        selected.update(range(a, b + 1))
    return sorted(k - 1 for k in selected if k <= nbursts)


def contiguous_runs(indices):
    runs = []
    for k in indices:
        if runs and k == runs[-1][-1] + 1:
            runs[-1].append(k)
        else:
            runs.append([k])
    return runs


def ensure_xy_project():
    if int(gs.region()["projection"]) != 0:
        gs.fatal(
            _(
                "Sentinel-1 SLC data are imported in radar geometry (lines, "
                "samples); the current project must be unprojected (XY). "
                "Create one with: grass -c XY <path_to_project> -e"
            )
        )


def gcp_transformer(target):
    """Return a function mapping lon/lat arrays to the target project CRS."""
    if not target:
        return lambda lon, lat: (lon, lat)
    from osgeo import osr

    env = gs.gisenv()
    tdir = os.path.join(env["GISDBASE"], target, "PERMANENT")
    if not os.path.isdir(tdir):
        gs.fatal(_("Target project <{}> not found in <{}>").format(target, env["GISDBASE"]))
    tenv = gs.create_environment(env["GISDBASE"], target, "PERMANENT")[1]
    wkt = gs.read_command("g.proj", flags="p", format="wkt", env=tenv)
    dst = osr.SpatialReference()
    if dst.ImportFromWkt(wkt) != 0:
        gs.fatal(_("Unable to read the CRS of target project <{}>").format(target))
    src = osr.SpatialReference()
    src.ImportFromEPSG(4326)
    for srs in (src, dst):
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(src, dst)

    def transform(lon, lat):
        pts = ct.TransformPoints([(float(x), float(y)) for x, y in zip(lon, lat)])
        return [p[0] for p in pts], [p[1] for p in pts]

    return transform


def write_gcps(group, layout, transform):
    """Write the imagery group POINTS file from the geolocation grid.

    Image coordinates follow r.in.gdal: x = sample + 0.5 (pixel centre),
    y = rows - (line + 0.5), in the unprojected project where the map spans
    0..cols and 0..rows.
    """
    sw = layout.swath
    np = sw.np
    rows = (np.asarray(sw.grid_times) - layout.t0) / sw.dt
    rows = rows[(rows >= 0) & (rows <= layout.nrows - 1)]
    rows = np.unique(np.concatenate([[0.0, layout.nrows - 1.0], rows]))
    cols = np.array(sw.grid_pixels[0])
    times = layout.t0 + rows * sw.dt
    lat = sw.geo_lut("latitude").rows(times, cols)
    lon = sw.geo_lut("longitude").rows(times, cols)
    ee, nn = transform(lon.ravel(), lat.ravel())
    img_e = np.tile(cols + 0.5, rows.size)
    img_n = np.repeat(layout.nrows - (rows + 0.5), cols.size)
    env = gs.gisenv()
    gdir = os.path.join(env["GISDBASE"], env["LOCATION_NAME"], env["MAPSET"], "group", group)
    with open(os.path.join(gdir, "POINTS"), "w") as fd:
        fd.write("# %7s %15s %15s %15s %9s status\n" % ("", "image", "", "target", ""))
        fd.write("# %15s %15s %15s %15s   (1=ok)\n" % ("east", "north", "east", "north"))
        fd.write("#\n")
        for e1, n1, e2, n2 in zip(img_e, img_n, ee, nn):
            fd.write("  %15f %15f %15f %15f %4d\n" % (e1, n1, e2, n2, 1))
    return img_e.size


def list_product(product, swaths):
    print(
        "swath|burst|burst_id|azimuth_time|lines|samples|west|south|east|north"
    )
    for sw in swaths:
        for k, b in enumerate(sw.bursts):
            w, s, e, n = sw.burst_footprint(k)
            print(
                "{}|{}|{}|{}|{}|{}|{:.5f}|{:.5f}|{:.5f}|{:.5f}".format(
                    sw.swath,
                    k + 1,
                    b.get("burst_id") or "",
                    b["azimuth_time"],
                    sw.lpb,
                    sw.nsamples,
                    w,
                    s,
                    e,
                    n,
                )
            )


def main():
    import numpy as np

    options, flags = gs.parser()
    atexit.register(cleanup)

    product = SafeProduct(options["input"])
    info = product.info()

    ann_names = [n for n in product.listdir("annotation") if n.lower().endswith(".xml")]
    if not ann_names:
        gs.fatal(_("No annotation files in <{}>").format(product.name))
    first_header = ET.fromstring(product.read("annotation/" + ann_names[0])).find("adsHeader")
    start = parse_time(first_header.findtext("startTime"))
    axis = TimeAxis(datetime(start.year, start.month, start.day))

    want_swaths = options["swath"].split(",") if options["swath"] else None
    want_pols = options["polarization"].split(",") if options["polarization"] else None
    swaths = []
    for name in ann_names:
        parts = name.split("-")
        if len(parts) < 4:
            continue
        if want_swaths and parts[1].upper() not in want_swaths:
            continue
        if want_pols and parts[3].upper() not in want_pols:
            continue
        swaths.append(Swath(product, name, axis))
    if not swaths:
        gs.fatal(_("No sub-swath/polarization of <{}> matches the selection").format(product.name))
    for sw in swaths:
        if sw.product_type != "SLC":
            gs.fatal(
                _("<{}> is a {} product; this module imports SLC products only").format(
                    product.name, sw.product_type
                )
            )
    swaths.sort(key=lambda s: (s.swath, s.pol))

    if flags["l"]:
        seen = set()
        list_product(product, [s for s in swaths if not (s.swath in seen or seen.add(s.swath))])
        return 0

    ensure_xy_project()

    measures = options["measure"].split(",") if options["measure"] else []
    geometry = options["geometry"].split(",") if options["geometry"] else []
    calibration = options["calibration"]
    denoise = flags["n"]
    deburst = not flags["b"]
    if not measures and not geometry:
        gs.fatal(_("Nothing to import: select at least one measure or geometry layer"))
    if denoise and ({"complex", "phase"} & set(measures)):
        gs.fatal(
            _(
                "Thermal noise removal is defined on power only; it cannot be "
                "applied to the complex or phase measures"
            )
        )
    cal_tag = {"sigma0": "sigmaNought", "beta0": "betaNought", "gamma0": "gamma"}.get(calibration)

    bbox = None
    if options["bbox"]:
        bbox = [float(v) for v in options["bbox"].split(",")]
        if len(bbox) != 4 or bbox[0] >= bbox[2] or bbox[1] >= bbox[3]:
            gs.fatal(_("bbox must be west,south,east,north with west<east and south<north"))

    # Burst selection per sub-swath (identical for its polarizations).
    selection = {}
    for sw in swaths:
        if sw.swath in selection:
            continue
        nb = len(sw.bursts)
        if options["bursts"]:
            sel = parse_bursts(options["bursts"].split(","), nb)
        elif bbox:
            sel = []
            for k in range(nb):
                w, s, e, n = sw.burst_footprint(k)
                if w <= bbox[2] and e >= bbox[0] and s <= bbox[3] and n >= bbox[1]:
                    sel.append(k)
        else:
            sel = list(range(nb))
        if not sel:
            gs.message(_("Sub-swath {}: no burst selected, skipped").format(sw.swath))
        selection[sw.swath] = sel
    swaths = [s for s in swaths if selection[s.swath]]
    if not swaths:
        gs.fatal(_("No burst selected in any sub-swath"))

    # Output plan: one layout per (swath, burst run) or per burst.
    jobs = []
    for sw in swaths:
        sel = selection[sw.swath]
        if deburst:
            runs = contiguous_runs(sel)
            if len(runs) > 1:
                gs.fatal(
                    _(
                        "Sub-swath {}: selected bursts {} are not contiguous; "
                        "debursting needs a contiguous range (or use -b)"
                    ).format(sw.swath, ",".join(str(k + 1) for k in sel))
                )
            jobs.append((sw, runs[0], ""))
        else:
            for k in sel:
                jobs.append((sw, [k], "_b{:02d}".format(k + 1)))

    prefix = options["output"]
    suffixes = []
    for m in measures:
        suffixes += ["i", "q"] if m == "complex" else [m]

    def names_for(sw, tag):
        base = "{}_{}_{}{}".format(prefix, sw.swath.lower(), sw.pol.lower(), tag)
        geo_base = "{}_{}{}".format(prefix, sw.swath.lower(), tag)
        return base, geo_base

    planned = []
    geo_done = set()
    for sw, _bursts, tag in jobs:
        base, geo_base = names_for(sw, tag)
        planned += ["{}_{}".format(base, s) for s in suffixes]
        if geo_base not in geo_done:
            planned += ["{}_{}".format(geo_base, g) for g in geometry]
            geo_done.add(geo_base)
    if not gs.overwrite():
        for name in planned:
            if gs.find_file(name, element="cell", mapset=".")["file"]:
                gs.fatal(
                    _("Raster map <{}> already exists, use --overwrite to replace it").format(name)
                )

    transform = gcp_transformer(options["target"])
    env = gs.gisenv()
    mapset_dir = os.path.join(env["GISDBASE"], env["LOCATION_NAME"], env["MAPSET"])
    geo_written = {}

    for sw, bursts, tag in jobs:
        base, geo_base = names_for(sw, tag)
        layout = Layout(sw, bursts, deburst)
        burst_label = "{}-{}".format(bursts[0] + 1, bursts[-1] + 1)
        gs.message(
            _("Importing {} {} bursts {} ({} lines x {} samples)...").format(
                sw.swath, sw.pol, burst_label, layout.nrows, layout.ncols
            )
        )
        reader = MeasurementReader(product, sw)
        cal_lut = sw.calibration_lut(cal_tag) if cal_tag else None

        outputs = {}
        for s in suffixes:
            is_float = not (s in ("i", "q") and cal_lut is None)
            outputs[s] = BinaryRaster("{}_{}".format(base, s), is_float)
        geo_out = {}
        if geo_base not in geo_written:
            for g in geometry:
                geo_out[g] = BinaryRaster("{}_{}".format(geo_base, g), True)

        cols = np.arange(layout.ncols)
        for r0, r1 in layout.chunks():
            gs.percent(r0, layout.nrows, 5)
            src = layout.src[r0:r1]
            times = layout.times[r0:r1]
            n = r1 - r0
            i = np.zeros((n, layout.ncols), dtype=np.float32)
            q = np.zeros_like(i)
            valid = np.zeros((n, layout.ncols), dtype=bool)
            good = np.flatnonzero(src >= 0)
            # Read runs of consecutive source lines in one call each.
            start = 0
            while start < good.size:
                stop = start + 1
                while stop < good.size and src[good[stop]] == src[good[stop - 1]] + 1 and good[stop] == good[stop - 1] + 1:
                    stop += 1
                rows = good[start:stop]
                bi, bq = reader.read(int(src[rows[0]]), rows.size)
                i[rows], q[rows] = bi, bq
                fvs = sw.first_valid_sample[src[rows]]
                lvs = sw.last_valid_sample[src[rows]]
                valid[rows] = (cols >= fvs[:, None]) & (cols <= lvs[:, None]) & (fvs >= 0)[:, None]
                start = stop

            power = None
            if {"amplitude", "intensity", "db"} & set(measures):
                power = i.astype(np.float64) ** 2 + q.astype(np.float64) ** 2
                if denoise:
                    noise = sw.noise_power(np.where(src >= 0, src, 0))
                    power = np.where(noise > 0, np.maximum(power - noise, 0.0), power)
            if cal_lut is not None:
                a = cal_lut.rows(times)
                if power is not None:
                    power = power / (a * a)
                i = i / a
                q = q / a

            for s, out in outputs.items():
                if s == "i":
                    data = i
                elif s == "q":
                    data = q
                elif s == "amplitude":
                    data = np.sqrt(power)
                elif s == "intensity":
                    data = power
                elif s == "db":
                    with np.errstate(divide="ignore"):
                        data = np.where(power > 0, 10.0 * np.log10(power), np.nan)
                elif s == "phase":
                    data = np.arctan2(q, i)
                if out.is_float:
                    data = np.where(valid, data, np.nan)
                else:
                    data = np.where(valid, data, INT_NULL)
                out.write(data)
            for g, out in geo_out.items():
                out.write(sw.geo_lut(GEOMETRY_TAGS[g][0]).rows(times))
        gs.percent(1, 1, 1)

        map_names = []
        product_meta = dict(info)
        swath_meta = sw.metadata()
        cal_const = sw.calibration_constant() if cal_tag else None
        t_first = layout.t0
        t_last = float(layout.times[-1])
        layout_meta = {
            "debursted": deburst,
            "bursts": [k + 1 for k in bursts],
            "rows": layout.nrows,
            "cols": layout.ncols,
            "first_line_time": sw.axis.iso(t_first),
            "last_line_time": sw.axis.iso(t_last),
            "azimuth_time_interval": sw.dt,
            "slant_range_time_first_sample": sw.slant_range_time,
            "range_pixel_spacing": sw.range_pixel_spacing,
            "source_first_line": int(layout.src[layout.src >= 0].min()),
            "image_coordinates": "x = sample + 0.5, y = rows - (line + 0.5)",
        }
        mission = info.get("mission") or sw.mission
        date_span = "{}/{}".format(
            sw.axis.datetime(t_first).strftime("%d %b %Y %H:%M:%S"),
            (sw.axis.datetime(t_last) + timedelta(seconds=1)).strftime("%d %b %Y %H:%M:%S"),
        )

        to_finish = [(s, out, False) for s, out in outputs.items()]
        to_finish += [(g, out, True) for g, out in geo_out.items()]
        for key, out, is_geo in to_finish:
            if is_geo:
                title = "{} {} {} {}".format(mission, sw.mode, sw.swath, key.replace("_", " "))
                units = GEOMETRY_TAGS[key][1]
                label = "S1_{}".format(key)
                measure_meta = {"geometry": key, "source": "annotation geolocation grid (bilinear)"}
            else:
                cal_text = calibration if cal_tag else "uncalibrated"
                title = "{} {} SLC {} {} {} ({}{})".format(
                    mission, sw.mode, sw.swath, sw.pol, key, cal_text,
                    ", noise removed" if denoise else "",
                )
                units = {
                    "i": "DN" if not cal_tag else "calibrated amplitude",
                    "q": "DN" if not cal_tag else "calibrated amplitude",
                    "amplitude": "DN" if not cal_tag else "calibrated amplitude",
                    "intensity": "DN^2" if not cal_tag else "linear power ratio",
                    "db": "dB",
                    "phase": "radians",
                }[key]
                label = "S1_{}_{}".format(sw.pol, key.upper())
                measure_meta = {
                    "measure": key,
                    "calibration": calibration,
                    "absolute_calibration_constant": cal_const,
                    "thermal_noise_removed": denoise,
                }
            out.finish(layout.nrows, layout.ncols, title)
            name = out.name
            map_names.append(name)
            gs.run_command(
                "r.support",
                map=name,
                title=title,
                units=units,
                source1=product.name,
                source2=sw.tiff,
                description="{} {} orbit {} (relative {}), {} pass, {}; bursts {}; "
                "metadata in cell_misc/{}/description.json".format(
                    mission,
                    info.get("product_type"),
                    info.get("absolute_orbit_start"),
                    info.get("relative_orbit_start"),
                    sw.pass_direction,
                    "debursted" if deburst else "single burst",
                    burst_label,
                    name,
                ),
                semantic_label=label,
                quiet=True,
            )
            gs.run_command("r.timestamp", map=name, date=date_span, quiet=True)
            gs.raster_history(name, overwrite=True)
            meta = {
                "product": product_meta,
                "swath": swath_meta,
                "raster_geometry": layout_meta,
            }
            meta.update(measure_meta)
            meta_dir = os.path.join(mapset_dir, "cell_misc", name)
            os.makedirs(meta_dir, exist_ok=True)
            with open(os.path.join(meta_dir, "description.json"), "w") as fd:
                json.dump(meta, fd, indent=1)

        if geo_out:
            geo_written[geo_base] = list(map_names[-len(geo_out):])
        group_maps = [m for m in map_names if m not in geo_written.get(geo_base, [])]
        group_maps += geo_written.get(geo_base, [])
        gs.run_command("i.group", group=base, input=group_maps, quiet=True)
        npts = write_gcps(base, layout, transform)
        if options["target"]:
            gs.run_command(
                "i.target", group=base, location=options["target"], mapset="PERMANENT", quiet=True
            )
        gs.message(
            _("Imagery group <{}>: {} maps, {} ground control points").format(
                base, len(group_maps), npts
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
