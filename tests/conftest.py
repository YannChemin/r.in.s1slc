"""Fixtures for r.in.s1slc tests: a tiny synthetic Sentinel-1 IW SLC product.

The synthetic product has one sub-swath (IW1, VV) of three bursts whose
annotations follow the SAFE layout read by the module. Its values are chosen
so that every processing step has an exact expected result:

- pixel I = 100 * burst + line_in_burst, Q = sample - 30,
- calibration LUT A(t, p) = 200 + 0.5 * p + 10 * (t - T0), linear in both
  azimuth time and range so bilinear interpolation is exact,
- noise = (5 + p) * (1 + 0.01 * stacked_line), also exactly interpolable,
- latitude = 25 + 0.01 * (t - T0), longitude = 55 + 0.001 * p.
"""

import os
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest

import grass.script as gs

T0 = datetime(2023, 1, 12, 14, 25, 9)
DT = 0.01
NBURSTS = 3
LPB = 40
NSAMPLES = 60
BURST_STEP = 34  # lines between burst starts, i.e. a 6-line overlap
FIRST_VALID_LINE = 2
LAST_VALID_LINE = 37
FIRST_VALID_SAMPLE = 3
LAST_VALID_SAMPLE = 56
NAME = "S1A_IW_SLC__1SDV_20230112T142509_20230112T142510_046752_059AD8_TEST.SAFE"
STEM = "s1a-iw1-slc-vv-20230112t142509-20230112t142510-046752-059ad8-004"
SCRIPT = str(Path(__file__).resolve().parent.parent / "r.in.s1slc.py")


def iso(seconds):
    return (T0 + timedelta(microseconds=round(seconds * 1e6))).isoformat(
        timespec="microseconds"
    )


def burst_time(k):
    return k * BURST_STEP * DT


def cal_value(t, p):
    return 200.0 + 0.5 * p + 10.0 * t


def noise_range(p):
    return 5.0 + p


def noise_azimuth(line):
    return 1.0 + 0.01 * line


def latitude(t):
    return 25.0 + 0.01 * t


def longitude(p):
    return 55.0 + 0.001 * p


def pixel_i(k, line):
    return 100 * k + line


def pixel_q(sample):
    return sample - 30


def annotation_xml():
    nlines = NBURSTS * LPB
    last = burst_time(NBURSTS - 1) + (LPB - 1) * DT
    fvs = [-1] * FIRST_VALID_LINE + [FIRST_VALID_SAMPLE] * (
        LAST_VALID_LINE - FIRST_VALID_LINE + 1
    )
    fvs += [-1] * (LPB - len(fvs))
    lvs = [-1 if v < 0 else LAST_VALID_SAMPLE for v in fvs]
    bursts = ""
    for k in range(NBURSTS):
        bursts += f"""
      <burst>
        <azimuthTime>{iso(burst_time(k))}</azimuthTime>
        <azimuthAnxTime>{400.0 + burst_time(k)}</azimuthAnxTime>
        <sensingTime>{iso(burst_time(k) + 1.0)}</sensingTime>
        <byteOffset>{1000 + k * LPB * NSAMPLES * 4}</byteOffset>
        <firstValidSample count="{LPB}">{" ".join(map(str, fvs))}</firstValidSample>
        <lastValidSample count="{LPB}">{" ".join(map(str, lvs))}</lastValidSample>
        <burstId absolute="{1000 + k}">{500 + k}</burstId>
      </burst>"""
    points = ""
    for k in range(NBURSTS + 1):
        line = min(k * LPB, nlines - 1)
        t = burst_time(k) if k < NBURSTS else last
        for p in (0, 20, 40, NSAMPLES - 1):
            points += f"""
      <geolocationGridPoint>
        <azimuthTime>{iso(t)}</azimuthTime>
        <slantRangeTime>5.33e-03</slantRangeTime>
        <line>{line}</line>
        <pixel>{p}</pixel>
        <latitude>{latitude(t)}</latitude>
        <longitude>{longitude(p)}</longitude>
        <height>{10.0 + p}</height>
        <incidenceAngle>{30.0 + 0.1 * p}</incidenceAngle>
        <elevationAngle>{27.0 + 0.1 * p}</elevationAngle>
      </geolocationGridPoint>"""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<product>
  <adsHeader>
    <missionId>S1A</missionId>
    <productType>SLC</productType>
    <polarisation>VV</polarisation>
    <mode>IW</mode>
    <swath>IW1</swath>
    <startTime>{iso(0.0)}</startTime>
    <stopTime>{iso(last)}</stopTime>
    <absoluteOrbitNumber>46752</absoluteOrbitNumber>
    <missionDataTakeId>367320</missionDataTakeId>
    <imageNumber>004</imageNumber>
  </adsHeader>
  <generalAnnotation>
    <productInformation>
      <pass>Ascending</pass>
      <timelinessCategory>Fast-24h</timelinessCategory>
      <platformHeading>-12.4</platformHeading>
      <projection>Slant Range</projection>
      <rangeSamplingRate>6.434523812571428e+07</rangeSamplingRate>
      <radarFrequency>5.405000454334350e+09</radarFrequency>
      <azimuthSteeringRate>1.590368784000000e+00</azimuthSteeringRate>
    </productInformation>
    <downlinkInformationList count="1">
      <downlinkInformation><prf>1717.128973878037</prf></downlinkInformation>
    </downlinkInformationList>
    <orbitList count="2">
      <orbit><time>{iso(-10.0)}</time><frame>Earth Fixed</frame>
        <position><x>1.0</x><y>2.0</y><z>3.0</z></position>
        <velocity><x>4.0</x><y>5.0</y><z>6.0</z></velocity></orbit>
      <orbit><time>{iso(0.0)}</time><frame>Earth Fixed</frame>
        <position><x>7.0</x><y>8.0</y><z>9.0</z></position>
        <velocity><x>1.5</x><y>2.5</y><z>3.5</z></velocity></orbit>
    </orbitList>
    <terrainHeightList count="1">
      <terrainHeight><azimuthTime>{iso(0.0)}</azimuthTime><value>50.0</value></terrainHeight>
    </terrainHeightList>
    <azimuthFmRateList count="1">
      <azimuthFmRate><azimuthTime>{iso(0.0)}</azimuthTime><t0>5.33e-03</t0>
        <azimuthFmRatePolynomial count="3">-2335.8 449443.7 -78230906.4</azimuthFmRatePolynomial>
      </azimuthFmRate>
    </azimuthFmRateList>
  </generalAnnotation>
  <imageAnnotation>
    <imageInformation>
      <productFirstLineUtcTime>{iso(0.0)}</productFirstLineUtcTime>
      <productLastLineUtcTime>{iso(last)}</productLastLineUtcTime>
      <slantRangeTime>5.330568075155056e-03</slantRangeTime>
      <pixelValue>Complex</pixelValue>
      <outputPixels>16 bit Signed Integer</outputPixels>
      <rangePixelSpacing>2.329562e+00</rangePixelSpacing>
      <azimuthPixelSpacing>1.399187e+01</azimuthPixelSpacing>
      <azimuthTimeInterval>{DT}</azimuthTimeInterval>
      <numberOfSamples>{NSAMPLES}</numberOfSamples>
      <numberOfLines>{nlines}</numberOfLines>
      <incidenceAngleMidSwath>33.97</incidenceAngleMidSwath>
    </imageInformation>
    <processingInformation>
      <thermalNoiseCorrectionPerformed>false</thermalNoiseCorrectionPerformed>
      <swathProcParamsList count="1"><swathProcParams>
        <swath>IW1</swath>
        <rangeProcessing><windowType>Hamming</windowType><processingBandwidth>5.6e+07</processingBandwidth></rangeProcessing>
        <azimuthProcessing><windowType>Hamming</windowType><processingBandwidth>327.0</processingBandwidth></azimuthProcessing>
      </swathProcParams></swathProcParamsList>
    </processingInformation>
  </imageAnnotation>
  <dopplerCentroid>
    <dcEstimateList count="1">
      <dcEstimate><azimuthTime>{iso(0.0)}</azimuthTime><t0>5.33e-03</t0>
        <geometryDcPolynomial count="3">2.09 -280.5 37703.9</geometryDcPolynomial>
        <dataDcPolynomial count="3">65.9 -121450.0 86145270.0</dataDcPolynomial>
      </dcEstimate>
    </dcEstimateList>
  </dopplerCentroid>
  <swathTiming>
    <linesPerBurst>{LPB}</linesPerBurst>
    <samplesPerBurst>{NSAMPLES}</samplesPerBurst>
    <burstList count="{NBURSTS}">{bursts}
    </burstList>
  </swathTiming>
  <geolocationGrid>
    <geolocationGridPointList count="{4 * (NBURSTS + 1)}">{points}
    </geolocationGridPointList>
  </geolocationGrid>
</product>
"""


def calibration_xml():
    pixels = list(range(0, NSAMPLES, 10)) + [NSAMPLES - 1]
    vectors = ""
    times = [-0.2 + 0.2 * j for j in range(8)]
    for t in times:
        sigma = " ".join(f"{cal_value(t, p):.9f}" for p in pixels)
        other = " ".join("237.0" for _p in pixels)
        vectors += f"""
    <calibrationVector>
      <azimuthTime>{iso(t)}</azimuthTime>
      <line>{round(t / DT)}</line>
      <pixel count="{len(pixels)}">{" ".join(map(str, pixels))}</pixel>
      <sigmaNought count="{len(pixels)}">{sigma}</sigmaNought>
      <betaNought count="{len(pixels)}">{other}</betaNought>
      <gamma count="{len(pixels)}">{sigma}</gamma>
      <dn count="{len(pixels)}">{other}</dn>
    </calibrationVector>"""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<calibration>
  <calibrationInformation><absoluteCalibrationConstant>1.0</absoluteCalibrationConstant></calibrationInformation>
  <calibrationVectorList count="{len(times)}">{vectors}
  </calibrationVectorList>
</calibration>
"""


def noise_xml():
    pixels = list(range(0, NSAMPLES, 10)) + [NSAMPLES - 1]
    rvec = ""
    for k in range(-1, NBURSTS + 1):
        rvec += f"""
    <noiseRangeVector>
      <azimuthTime>{iso(burst_time(k))}</azimuthTime>
      <line>{k * LPB}</line>
      <pixel count="{len(pixels)}">{" ".join(map(str, pixels))}</pixel>
      <noiseRangeLut count="{len(pixels)}">{" ".join(str(noise_range(p)) for p in pixels)}</noiseRangeLut>
    </noiseRangeVector>"""
    lines = list(range(0, NBURSTS * LPB, 10)) + [NBURSTS * LPB - 1]
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<noise>
  <noiseRangeVectorList count="{NBURSTS + 2}">{rvec}
  </noiseRangeVectorList>
  <noiseAzimuthVectorList count="1">
    <noiseAzimuthVector>
      <swath>IW1</swath>
      <firstAzimuthLine>0</firstAzimuthLine>
      <firstRangeSample>0</firstRangeSample>
      <lastAzimuthLine>{NBURSTS * LPB - 1}</lastAzimuthLine>
      <lastRangeSample>{NSAMPLES - 1}</lastRangeSample>
      <line count="{len(lines)}">{" ".join(map(str, lines))}</line>
      <noiseAzimuthLut count="{len(lines)}">{" ".join(str(noise_azimuth(v)) for v in lines)}</noiseAzimuthLut>
    </noiseAzimuthVector>
  </noiseAzimuthVectorList>
</noise>
"""


MANIFEST = """<?xml version="1.0" encoding="UTF-8"?>
<xfdu:XFDU xmlns:xfdu="urn:ccsds:schema:xfdu:1"
  xmlns:safe="http://www.esa.int/safe/sentinel-1.0"
  xmlns:s1sarl1="http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1">
  <metadataSection>
    <safe:processing name="SLC Processing">
      <safe:facility name="Test">
        <safe:software name="Sentinel-1 IPF" version="003.52"/>
      </safe:facility>
    </safe:processing>
    <safe:platform>
      <safe:familyName>SENTINEL-1</safe:familyName>
      <safe:number>A</safe:number>
      <safe:instrument><safe:extension><s1sarl1:instrumentMode>
        <s1sarl1:mode>IW</s1sarl1:mode></s1sarl1:instrumentMode></safe:extension></safe:instrument>
    </safe:platform>
    <safe:orbitReference>
      <safe:orbitNumber type="start">46752</safe:orbitNumber>
      <safe:relativeOrbitNumber type="start">130</safe:relativeOrbitNumber>
      <safe:cycleNumber>279</safe:cycleNumber>
      <safe:extension><s1:orbitProperties xmlns:s1="http://www.esa.int/safe/sentinel-1.0/sentinel-1">
        <s1:pass>ASCENDING</s1:pass></s1:orbitProperties></safe:extension>
    </safe:orbitReference>
    <safe:acquisitionPeriod>
      <safe:startTime>2023-01-12T14:25:09.000000</safe:startTime>
      <safe:stopTime>2023-01-12T14:25:10.070000</safe:stopTime>
    </safe:acquisitionPeriod>
    <s1sarl1:standAloneProductInformation>
      <s1sarl1:transmitterReceiverPolarisation>VV</s1sarl1:transmitterReceiverPolarisation>
      <s1sarl1:productType>SLC</s1sarl1:productType>
    </s1sarl1:standAloneProductInformation>
  </metadataSection>
</xfdu:XFDU>
"""


def write_safe(parent):
    from osgeo import gdal

    gdal.UseExceptions()
    safe = Path(parent) / NAME
    for sub in ("annotation/calibration", "measurement"):
        (safe / sub).mkdir(parents=True, exist_ok=True)
    (safe / "manifest.safe").write_text(MANIFEST)
    (safe / "annotation" / f"{STEM}.xml").write_text(annotation_xml())
    (safe / "annotation" / "calibration" / f"calibration-{STEM}.xml").write_text(
        calibration_xml()
    )
    (safe / "annotation" / "calibration" / f"noise-{STEM}.xml").write_text(noise_xml())

    nlines = NBURSTS * LPB
    data = np.zeros((nlines, NSAMPLES), dtype=np.complex64)
    for k in range(NBURSTS):
        for line in range(LPB):
            data[k * LPB + line] = pixel_i(k, line) + 1j * pixel_q(np.arange(NSAMPLES))
    ds = gdal.GetDriverByName("GTiff").Create(
        str(safe / "measurement" / f"{STEM}.tiff"), NSAMPLES, nlines, 1, gdal.GDT_CInt16
    )
    ds.GetRasterBand(1).WriteArray(data)
    ds = None
    return safe


@pytest.fixture(scope="session")
def safe_product(tmp_path_factory):
    return write_safe(tmp_path_factory.mktemp("safe"))


@pytest.fixture(scope="session")
def safe_zip(safe_product, tmp_path_factory):
    import shutil

    base = tmp_path_factory.mktemp("zip") / NAME.replace(".SAFE", "")
    return Path(shutil.make_archive(str(base), "zip", safe_product.parent, NAME))


@pytest.fixture
def xy_session(tmp_path):
    project = tmp_path / "xy"
    gs.create_project(project)
    with gs.setup.init(project, env=os.environ.copy()) as session:
        yield session


def run_module(session, *flags, **kwargs):
    """Run the script from the source tree; return the completed process."""
    import subprocess
    import sys

    args = [sys.executable, SCRIPT]
    args += ["-" + f for f in flags]
    args += [f"{k}={v}" for k, v in kwargs.items()]
    return subprocess.run(args, env=session.env, capture_output=True, text=True, check=False)


def read_map(session, name, null=np.nan):
    import grass.script.array as garray

    gs.run_command("g.region", raster=name, env=session.env)
    return np.asarray(garray.array(name, null="nan" if null is np.nan else null, env=session.env))
