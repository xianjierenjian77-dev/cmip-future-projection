# Variable and input recipes

Contents: 1. Scope; 2. Temperature; 3. Precipitation; 4. Circulation/flux variables;
5. Series; 6. Daily events and unsupported inputs.

## 1. Start from input semantics

The numerical kernels can operate on many continuous variables, but that does
not make every transformation scientifically appropriate. Establish units,
frequency, aggregation and physical bounds before choosing correction.

| Input | Bundled route | Additional work |
|---|---|---|
| Rectilinear NetCDF, one variable/level | Gridded route | Confirm calendar and units; remap if necessary. |
| Rectilinear GRIB | Gridded route with cfgrib/ecCodes | Resolve message groups and metadata before analysis. |
| Annual CSV/NPZ | Series route | Document how the annual index was produced. |
| WorldClim GeoTIFF climatologies | No direct reader | Preserve period/month/CRS/nodata semantics in an adapter; a 20-year climatology is not 20 annual observations. |
| Curvilinear ocean/ice grid | No native general remapper | Remap with cell-area/mask-aware tools or compute a verified annual series upstream. |
| Daily station CSV | No direct daily reader | Add a date-aware adapter or use an appropriate external daily analysis. |

Do not rerun bias adjustment on an already adjusted product unless a new,
explicitly justified adjustment is intended. Interpolation to a finer grid is
not dynamical or statistical downscaling.

## 2. Temperature

For temperature supplied in Kelvin and reported in Celsius:

~~~yaml
variable: {name: tas, units: degC, long_name: Near-surface air temperature}
data:
  reference:
    variable: tas
    units_in: K
    units_transform: {scale: 1.0, offset: -273.15}
    aggregate: {months: JJA, stat: mean, frequency: monthly, weight_by_days: true}
  models:
    variable: tas
    units_in: K
    units_transform: {scale: 1.0, offset: -273.15}
    aggregate: {months: JJA, stat: mean, frequency: monthly, weight_by_days: true}
correction: {method: mean, multiplicative: false}
extremes: {enabled: false}
~~~

This snippet is not a complete configuration: supply data paths, periods and a
region. Inspect `cell_methods` when distinguishing daily mean tas, daily maximum
tasmax and daily minimum tasmin. Monthly tasmax means are not the hottest daily
values. Use additive rather than ratio correction for Celsius temperature.

For seasonal maxima of a regional daily-maximum-temperature mean, use
`templates/config_regional_extremes.yaml`. A different event, such as the spatial
average of local hottest days, needs a different operator order and label.
A lower-tail cold event needs a defined sign/variable transformation, not merely
an upper-tail GEV applied to annual minima.

## 3. Precipitation

| Verified source quantity | Conversion before seasonal reduction |
|---|---|
| Mean flux in kg m−2 s−1 | ×86400 gives mm day−1; take day-weighted seasonal mean for a mean rate. |
| Mean flux in m s−1 | ×86400000 gives mm day−1. |
| Accumulated depth in m over each input interval | ×1000 gives mm per interval; sum nonoverlapping intervals for a total. |
| Monthly mean daily accumulation in m/day | ×1000 gives mm/day; weight monthly means by month length. |

ERA5 precipitation products differ in accumulation/averaging semantics. Inspect
the actual product documentation, units and time bounds; do not reuse a scale
factor merely because a field is named tp. To calculate a monthly total from a
mean rate in mm/day, multiply by the calendar month's days before summing.
The CLI does not automatically perform this rate integration.

`templates/config_precip_qdm.yaml` assumes both reference and model data are
verified monthly mean rates. It converts kg m−2 s−1 to mm/day and computes a
seasonal mean rate. QDM's ratio form is optional and should be checked around
zeros. A zero lower bound prevents negative corrected values; it does not model
wet-day occurrence or enforce water/energy-budget consistency.

Keep GEV disabled for seasonal totals or mean rates by default. For heavy-rain
events use daily data and a precisely defined daily regional or pointwise index.
A high quantile of monthly data is not a daily precipitation extreme.

## 4. Other continuous variables

| Variable | Method considerations |
|---|---|
| Wind speed | Nonnegative magnitude; bound at zero if justified. Vector components and direction require separate handling. |
| Geopotential/height | Distinguish geopotential m²/s² from height m: height = geopotential / 9.80665. Select the intended pressure level and verify its coordinate. |
| Mean sea-level pressure | Pa to hPa uses ×0.01. Avoid confusing surface pressure with sea-level pressure. |
| Relative humidity | Distinguish fraction from percent; 0–1 or 0–100 bounds can matter, but clipping distorts distributions. |
| Specific humidity | kg/kg to g/kg uses ×1000. Near-zero values complicate ratio correction. |
| Radiative/turbulent fluxes | Confirm direction/sign convention and whether values are W/m² or accumulated J/m². Accumulation-to-flux conversion needs actual interval duration. |
| SST | Apply temperature units and ocean masks consistently; native ocean grids often need prior remapping. |
| Sea-ice concentration | Distinguish fraction and percent, apply the ocean mask, and use true cell areas for ice area. A regional average concentration is not sea-ice area or volume. |

A circulation index can use `region_diff` or a fixed `projection_index`.
For the latter, keep longitude coverage adequate for the intended zonal mean,
set a fixed historical template, and report dimensionless units. Never infer
causal mechanisms from an index trend alone.

## 5. Existing annual series

Use a CSV with one row per year, for example:

~~~text
year,value
1985,1.20
1986,1.35
~~~

Keep the reference in one file and models in
`data/models/<MODEL>/historical.csv` and `data/models/<MODEL>/ssp245.csv`.
NPZ files use `years` and `values` arrays. All series must describe the same
quantity, units, annual/seasonal label and reference baseline.

Use `index.type: series` and both source kinds `series`. Apply any desired
anomaly, detrending or standardization upstream with historical parameters;
the configuration rejects unsupported transformations in this route. A compound
or custom index is supported as a precomputed series, not automatically derived
from multiple variables by this package.

## 6. Daily extremes and event definitions

Daily gridded inputs require a spatial bounding box. Use one consistent member
for each model's historical and scenario data. Combining members before finding
peaks changes the event distribution and should not be used as a routine
extremes shortcut.

- `kind: block`: maximum of the daily regional series per complete year/season.
- `kind: pot`: date-declustered exceedances of a stated threshold, fitted with GPD.
- Heatwave duration, consecutive wet/dry days, storm tracking and compound-event
  dependence: additional event-detection algorithms are required.

A fixed historical or physical threshold supports comparisons of the same event.
A per-series quantile moves with each climate and addresses a different question.
Record which definition was used, especially for projected exceedance frequency.

A file of twelve long-term monthly climatologies cannot support interannual
trend fitting, annual GEV, or daily POT. Report changes between climatological
periods through a dedicated adapter without manufacturing time samples.
