## DESCRIPTION

*r.in.s1slc* imports Sentinel-1 TOPS (IW, EW) **Single Look Complex
(SLC)** products in their native radar geometry (azimuth lines ×
slant-range samples). The product can be given as the downloaded `.zip`
archive (read in place, without extraction), as the `.SAFE` directory or
as its `manifest.safe` file.

For every selected sub-swath and polarization the module can:

- **deburst** the TOPS bursts into one continuous image (default), or
  keep every burst as its own map (**-b** flag);
- write the complex signal (**measure=complex**, maps `_i` and `_q`),
  the amplitude, the intensity, the intensity in decibels or the phase;
- apply the **radiometric calibration** of the product
  (**calibration=sigma0|beta0|gamma0**);
- remove the **thermal noise** (**-n** flag) from the power measures;
- write **geolocation layers** (latitude, longitude, height, incidence
  and elevation angles) interpolated from the annotation geolocation grid;
- store **ground control points** in an imagery group for *i.rectify*;
- import the **product metadata** (see below).

Radar geometry has no map projection: the current project must be
unprojected (XY). Create one with `grass -c XY <path_to_project> -e`.
The module fails otherwise.

### Processing

The algorithms follow ESA SNAP (microwave toolbox), so that results can
be compared with a SNAP *TOPSAR-Split → Calibration → TOPSAR-Deburst*
chain:

- **Debursting** (*TOPSAR-Deburst*): output lines are sampled every
  azimuth time interval from the first valid line of the first burst to
  the last valid line of the last burst. In the overlap of two bursts a
  line is taken from the earlier burst up to the middle of the overlap and
  from the later burst after it; lines are nearest neighbours (no
  resampling), so the complex values are the original ones. Samples
  outside the annotated valid range of each line (`firstValidSample`,
  `lastValidSample`) are null.
- **Calibration**: `value / A` for the complex components and amplitude,
  `|DN|² / A²` for the intensity, where `A` is the sigmaNought,
  betaNought or gamma calibration vector, interpolated bilinearly in
  range (pixel) and azimuth time. The azimuth time used is the true
  zero-Doppler time of the source line in its burst.
- **Thermal noise removal**: the noise power `N` is the product of the
  range noise vector of the burst (the one closest in time to the burst
  start, as SNAP does) and the azimuth noise vector (IPF ≥ 2.9), both
  linearly interpolated. `|DN|² − N` is clipped at 0 before calibration.
  It is not defined for the complex and phase measures, which are refused.

### Output naming

Maps are named `{output}_{swath}_{pol}_{measure}`, e.g.
`s1_iw2_vv_i`, `s1_iw2_vv_q`, `s1_iw2_vv_intensity`. With the **-b**
flag a burst tag is added: `s1_iw2_vv_b04_i`. Geolocation layers do not
depend on the polarization: `s1_iw2_latitude`, `s1_iw2_b04_latitude`.

Uncalibrated complex components are integer maps (CELL) holding the
original 16-bit digital numbers; all other outputs are FCELL.

Each map covers `0..cols` east and `0..rows` north, one cell per sample
and line (like *r.in.gdal* for non-georeferenced data): use
`g.region raster=<map>` before working with it.

### Imagery groups and geocoding

One imagery group `{output}_{swath}_{pol}` (plus the burst tag with
**-b**) is created for each image, containing its maps and the
geolocation layers of the sub-swath. Its `POINTS` file holds ground
control points taken from the annotation geolocation grid, mapped to the
output lines through their azimuth time: image coordinates are cell
centres, target coordinates are WGS84 longitude/latitude, or the CRS of
the **target** project, which then also becomes the group target (see
*i.target*). *i.rectify* can then geocode the maps; the thin plate
spline method (**-t**) follows the grid best.

The geolocation grid is given on the ellipsoid plus the annotated
terrain height; in relief, terrain correction with a DEM is needed for
accurate geocoding.

### Metadata

Product metadata are imported with every map:

- title, units, source (product and measurement file), description and
  semantic label (e.g. `S1_VV_I`, `S1_VV_INTENSITY`, `S1_latitude`) with
  *r.support*;
- the azimuth time span of the image as timestamp (*r.timestamp*),
  in whole seconds;
- the command history;
- a JSON file `$MAPSET/cell_misc/<map>/description.json` with the
  manifest information (mission, absolute and relative orbit, pass, IPF
  version, footprint...), the sub-swath annotation (radar frequency,
  wavelength, range sampling rate, PRF, slant range time, pixel
  spacings, azimuth time interval, burst list with burst IDs, times and
  valid samples, orbit state vectors, Doppler centroid and azimuth FM
  rate polynomials, processing information, geolocation grid), the
  output raster geometry (first line azimuth time, source lines, bursts,
  debursting) and the measure, calibration and noise settings.

These are the parameters needed by downstream processing such as
coregistration and interferometry.

### Subsetting

A full IW SLC product holds three sub-swaths of nine bursts each,
about 1.2 GB per sub-swath and polarization. Use **swath**,
**polarization** and either **bursts** (1-based indices or ranges) or
**bbox** (bursts intersecting a longitude/latitude box) to import only
the area of interest. The **-l** flag lists the bursts of each
sub-swath with their burst ID, azimuth time and footprint, and exits.

## NOTES

Merging the sub-swaths (*TOPSAR-Merge*) is not done: each sub-swath has
its own slant range origin and is written separately.

The module requires the GDAL Python bindings and NumPy.

## EXAMPLES

List the bursts of a product:

```sh
r.in.s1slc -l input=S1A_IW_SLC__1SDV_20230112T142507_20230112T142534_046752_059AD8_A55D.zip
```

Import the debursted complex signal of the bursts over Sharjah (UAE)
from the IW2 VV image, with latitude and longitude layers:

```sh
grass -c XY $HOME/grassdata/s1slc_sharjah -e
grass $HOME/grassdata/s1slc_sharjah/PERMANENT --exec \
    r.in.s1slc input=S1A_IW_SLC__1SDV_20230112T142507_20230112T142534_046752_059AD8_A55D.zip \
        output=s1_20230112 swath=IW2 polarization=VV bbox=55.35,25.25,55.50,25.40 \
        measure=complex geometry=latitude,longitude
```

Calibrated, noise-free sigma0 in dB of burst 5 of each sub-swath:

```sh
r.in.s1slc -n -b input=S1A_IW_SLC__1SDV_20230112T142507_20230112T142534_046752_059AD8_A55D.zip \
    output=s1 bursts=5 measure=db calibration=sigma0
```

Geocode the result into a UTM project with the GCPs:

```sh
r.in.s1slc input=... output=s1 swath=IW2 bursts=4-6 measure=db \
    calibration=sigma0 target=sharjah_utm40n
i.rectify -t group=s1_iw2_vv extension=_utm resolution=10
```

## REFERENCES

- ESA, *Sentinel-1 Product Specification*, S1-RS-MDA-52-7441.
- ESA, *Radiometric Calibration of S-1 Level-1 Products Generated by
  the S-1 IPF*, ESA-EOPG-CSCOP-TN-0002.
- ESA, *Thermal Denoising of Products Generated by the S-1 IPF*,
  MPC-0392.
- ESA SNAP microwave toolbox:
  <https://github.com/senbox-org/microwave-toolbox>

## SEE ALSO

*[i.group](https://grass.osgeo.org/grass-stable/manuals/i.group.html),
[i.rectify](https://grass.osgeo.org/grass-stable/manuals/i.rectify.html),
[i.target](https://grass.osgeo.org/grass-stable/manuals/i.target.html),
[r.in.gdal](https://grass.osgeo.org/grass-stable/manuals/r.in.gdal.html),
[r.in.sentinel](r.in.sentinel.html),
[r.support](https://grass.osgeo.org/grass-stable/manuals/r.support.html),
[r.timestamp](https://grass.osgeo.org/grass-stable/manuals/r.timestamp.html)*

## AUTHORS

Yann Chemin
