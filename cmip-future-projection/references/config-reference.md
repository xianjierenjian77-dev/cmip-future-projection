# Configuration reference

Every key is optional unless marked **required**. Unset keys fall back to the
default shown. Unknown keys are reported as warnings and ignored, so a config
written for a newer version still runs.

Placeholders: strings are expanded with `{name}` (project name), `{root}`
(`data.models.root`), `{outdir}` and `{variable}`. `{model}`, `{experiment}`,
`{member}`, `{year}` are **not** expanded at load time — they are resolved per
file by the loader.

---

## `project`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `name` | str | `cmip_projection` | run name, used in `outdir` and the report title |
| `outdir` | str | `runs/{name}` | output directory (absolute or relative to the CWD) |
| `description` | str | `""` | free text stored in the resolved config |

## `variable`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `name` | str | `variable` | short name used in file names and titles |
| `units` | str | `""` | units **after** `units_transform`; used on axis labels |
| `long_name` | str | `""` | descriptive label for figure axes |

`data.*.units_in` is informational (the units in the source files); the actual
conversion is always `units_transform`.

## `data.reference` — observations / reanalysis

| Key | Type | Default | Meaning |
|---|---|---|---|
| `kind` | `grid` \| `series` | `grid` | gridded files or a pre-computed annual series |
| `paths` | list[str] | `[]` | file paths or globs; `{year}` expands over `years` |
| `years` | `[start, end]` | `null` | period to load (also used to expand `{year}`) |
| `variable` | str | auto | variable name inside the file; omit to auto-detect |
| `units_in` | str | `null` | informational only |
| `units_transform` | `{scale, offset}` | `{1.0, 0.0}` | `x*scale + offset` applied on load |
| `plev` | number | `null` | pressure level in hPa (or Pa); nearest match is used |
| `bbox` | `{lat: [lo,hi], lon: [lo,hi]}` | `null` | spatial subset before regridding (saves memory) |
| `aggregate.months` | list/int/str | `ANN` | `[9,10]`, `"JJA"`, `"DJF"`, `"ANN"`, `"9-10"` |
| `aggregate.stat` | str | `mean` | `mean,sum,max,min,median,std,p05..p99,count_gt,count_ge,count_lt,count_le,max_abs,mean_abs` |
| `aggregate.threshold` | number | `null` | required by `count_*` (units of the converted data) |
| `series_file` | str | `null` | CSV/NPZ, used when `kind: series` |
| `series_columns` | `{year,value}` | `{year,value}` | column names for CSV input |

`kind: series` expects a CSV with `year,value` columns or an NPZ with arrays
`years` and `values`.

## `data.models` — the model ensemble

| Key | Type | Default | Meaning |
|---|---|---|---|
| `kind` | `grid` \| `series` | `grid` | |
| `root` | str | **required** (grid) | archive root |
| `patterns` | list[str] | `[]` | candidate patterns, tried in order; empty = built-in CMIP6 layouts |
| `pattern` | str | `null` | single pattern (legacy); prepended to the candidates |
| `variable` | str | `null` | variable name inside the files |
| `units_in` | str | `null` | informational |
| `units_transform` | `{scale, offset}` | `{1.0, 0.0}` | applied on load |
| `plev` | number | `null` | |
| `bbox` | `{lat, lon}` | `null` | |
| `historical_experiment` | str | `historical` | experiment name of the baseline run |
| `scenarios` | list[str] | `[ssp245]` | future experiments to process |
| `scenario_labels` | dict | `{}` | overrides the default `SSP2-4.5` style labels |
| `members` | list[str] | `[r1i1p1f1]` | tried in order; `"*"` discovers every member |
| `member_combine` | `first` \| `mean` \| `median` | `first` | `first`: one member only; `mean`/`median`: load all configured/discovered members and combine them per year (multi-member ensemble) |
| `models` | `all` \| list[str] | `all` | explicit list is recommended for reproducibility |
| `exclude` | list[str] | `[]` | models to drop |
| `aggregate.*` | as above | | season/stat for the model data |
| `series_dir` | str | `null` | for `kind: series`: `<model>/<scenario>.csv` or `<model>_<scenario>.npz` |

Built-in patterns (used when `patterns` is empty):

```
{root}/{experiment}/{model}/{variable}_Amon_{model}_{experiment}_{member}_*.nc
{root}/{experiment}/{model}/{variable}/{variable}_Amon_{model}_{experiment}_{member}_*.nc
```

Explicit patterns are tried first, the built-ins are kept as fallbacks. Other
placeholders available to `patterns`: `{root}`, `{model}`, `{experiment}`,
`{member}`, `{variable}`. Globs are allowed anywhere.

**Discovery**: when `models: all`, model names are captured from the `{model}`
position of the first matching pattern — this works whether `{model}` is a
directory level or part of the file name.

## `grid`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `lat` | `{start,stop,step}` or list | `null` (copied from reference) | target latitudes |
| `lon` | `{start,stop,step}` or list | `null` (copied from reference) | target longitudes |
| `method` | `linear` \| `nearest` | `linear` | interpolation to the target grid |
| `lon_360` | bool | `true` | target longitude convention |

Use `{start: 0, stop: 70, step: 1}`; the number of points is derived, the `stop`
value is inclusive.

## `periods`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `historical` | `[start,end]` | `[1950,2014]` | model historical run to process |
| `future` | `[start,end]` | `[2026,2100]` | projection period |
| `reference_climate` | `[start,end]` | `[1950,2014]` | climatology used for anomalies and map changes |
| `future_windows` | list of pairs | 4 × 20-year windows | averaging windows for the map figure |

`historical` must end before `future` starts. If the model data do not cover a
period, the available intersection is used and logged.

## `index`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `type` | str | `region_mean` | `region_mean`, `region_mean_anomaly`, `region_diff`, `projection_index`, `series` |
| `region` | `{lat,lon,name?}` | `null` | region for `region_mean*` and `projection_index` |
| `regions` | list of 2 regions | `null` | `[reference, target]` for `region_diff` |
| `remove_zonal_mean` | bool | `false` | subtract each latitude's zonal mean before reducing |
| `anomaly` | bool | `false` | subtract the baseline climatology |
| `detrend` | bool | `false` | remove the series' linear trend |
| `standardize` | bool | `false` | divide by the baseline standard deviation |
| `units_label` | str | `null` | overrides `variable.units` on index axes |
| `projection.template_source` | `reference` \| `model` \| `year` | `reference` | where the projection template comes from |
| `projection.template_period` | `[start,end]` | `reference_climate` | climatology period of the template |
| `projection.template_year` | int | `null` | single year template (`template_source: year`) |
| `projection.template_model` | str | `null` | model to use (`template_source: model`) |

Region bounds: `{lat: [25, 45], lon: [90, 120]}`. Longitudes may be 0–360 or
−180–180 and may wrap (`lon: [350, 20]`).

## `correction`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `method` | str | `mean_std` | `none`, `mean`, `mean_std`, `delta_change`, `quantile_mapping`, `quantile_delta_mapping` |
| `params_period` | `[start,end]` | `reference_climate` | baseline used to fit the transfer function |
| `apply_to` | list | `[historical, future]` | which fields get corrected |
| `multiplicative` | `auto` \| bool | `auto` | ratio-scale correction; `auto` = precipitation/wind keywords and non-negative data |
| `lower_bound` | number | `null` | clip after correction (use `0` for precipitation) |
| `upper_bound` | number | `null` | clip after correction |
| `n_quantiles` | int | `100` | quantile table size for the two mapping methods |
| `trend_correction.enabled` | bool | `false` | rescale future trends to the observed trend |
| `trend_correction.mode` | `difference` \| `ratio` | `difference` | |

The transfer function is fitted **per grid point** on the years common to the
reference and the model historical run inside `params_period`.

## `ensemble`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `statistic` | `mean` \| `median` | `mean` | ensemble centre |
| `bands` | list[int] | `[10,90,20,80]` | percentile bands (two pairs drawn) |
| `weights` | `equal` \| `skill` | `equal` | inverse-RMSE weights from the historical scores |
| `exclude_models` | list[str] | `[]` | dropped after discovery |
| `min_models` | int | `3` | below this the run still proceeds but logs degenerate bands (single-model runs are allowed) |

## `extremes`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `true` | run the GEV stage |
| `source` | `future` \| `historical` | `future` | `historical` = only the observed fit |
| `fit_window` | `{historical, future}` | `null` → `periods.*` | fit windows |
| `distribution` | `gev` \| `gumbel` | `gev` | |
| `fit_method` | `mle` \| `mm` | `mle` | falls back automatically |
| `return_periods` | list[float] | `[1.1,2,5,10,20,30,50,100,200,500,1000]` | written to `gev_return_levels.csv` and marked on the curves |
| `reference_value` | number | `null` | value whose return period is reported (defaults to the observed value at the reference year) |
| `bootstrap` | int | `0` | percentile-bootstrap resamples of the ensemble-median curve |
| `show_empirical` | bool | `false` | add empirical (Gringorten) points to the figure |
| `min_years` | int | `20` | minimum sample for a fit |
| `shape_comparison` | bool | `true` | Kruskal–Wallis across scenarios |

## `daily` — daily-scale extremes

The annual pipeline uses one value per year. This block adds a second, independent
analysis on **daily** data: it loads the daily field for a season, reduces it to a
daily regional index, bias-corrects that series, and fits either a peaks-over-
threshold GPD or annual maxima.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `false` | run the daily stage |
| `months` | list/int/str | `[6,7,8]` | season to keep from the daily files |
| `index` | str | `region_mean` | `region_mean`, `region_mean_anomaly`, `region_diff` (daily) |
| `baseline` | `[start,end]` | `reference_climate` | period used to fit the daily correction and the historical extreme fit |
| `correction` | str | `quantile_delta_mapping` | correction applied to the daily series (`none`/`mean`/`mean_std`/`quantile_mapping`/`quantile_delta_mapping`) |
| `extremes.kind` | `pot` \| `block` | `pot` | GPD peaks-over-threshold or annual-maxima GEV |
| `extremes.threshold_mode` | `observed` \| `quantile` \| `absolute` | `observed` | `observed`: one threshold from the reference baseline (comparable across models); `quantile`: each series' own quantile (keeps the exceedance rate constant under warming); `absolute`: a physical value such as 35 °C |
| `extremes.threshold_quantile` | float | `0.95` | quantile for the observed/quantile modes |
| `extremes.threshold_abs` | number | `null` | required for `threshold_mode: absolute` |
| `extremes.decluster_days` | int | `3` | runs of exceedances separated by fewer days form one event |
| `extremes.min_peaks` | int | `30` | minimum number of declustered peaks for a GPD fit |
| `extremes.min_years` | int | `20` | minimum years for a block-maxima fit |
| `extremes.distribution` | `gev` \| `gumbel` | `gev` | for `kind: block` |
| `extremes.return_periods` | list[float] | `[2,5,10,20,50,100]` | written to `daily_extremes.csv` |
| `plots` | bool | `true` | draw `fig_daily_return_levels.png` |

Daily input must be gridded (`data.*.kind: grid`). Set `data.*.bbox` to the study
region: the loader refuses fields larger than ~30 M values, which is exactly what
a global 0.25° daily field would be.

The threshold choice matters scientifically. With a warming variable, a fixed
`observed` threshold yields rapidly increasing exceedance rates (a legitimate
impact metric), whereas `quantile` keeps the rate constant and compares how the
tail level itself shifts. State which one you used.

## `plots`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `dpi` | int | `300` | output resolution |
| `formats` | list[str] | `[png]` | any matplotlib format (`pdf`, `svg`, …) |
| `figsize_scale` | float | `1.0` | global font/axes scaling |
| `font` | str | `null` | font family; `null` = DejaVu, or a CJK font when labels contain Chinese |
| `language` | `en` \| `zh` | `en` | only affects automatic labels |
| `colors.historical` | colour | `#1f4e79` | observed line |
| `colors.observation` | colour | `#111111` | secondary observed elements |
| `colors.models` | colour | `#9e9e9e` | individual model points |
| `colors.scenario.<name>` | colour | palette | per-scenario colour |
| `figures` | list[str] | `[timeseries, gev, combined, maps, kde, shape_box]` | which figures to draw |
| `timeseries.show_members` | bool | `true` | individual model scatter |
| `timeseries.show_bands` | bool | `true` | percentile bands |
| `timeseries.show_trend` | bool | `true` | trend lines + `/10a` annotations |
| `timeseries.reference_year` | int | `2025` | year marked between history and future |
| `timeseries.reference_value` | number | `null` | value at that year; auto-read from the reference series when available |
| `timeseries.per_model` | bool | `false` | also write `fig_timeseries_<model>.png`: that model's own projection, all scenarios in one panel, ensemble mean as dashed reference |
| `timeseries.ylim` / `xlim` | pair | `null` | axis limits |
| `gev.xlim` | pair | `[1.1, 1100]` | return-period axis |
| `gev.ylim` | pair | `null` | |
| `maps.variable` | `change` \| `absolute` | `change` | map quantity |
| `maps.cmap` | str | `RdBu_r` | colormap |
| `maps.levels` | `[vmin,vmax]` or list | `null` | explicit colour limits |
| `maps.n_levels` | int | `21` | number of levels when `levels` is a range |
| `maps.stipple` | bool | `true` | draw the model-agreement stipple |
| `maps.stipple_alpha` | float | `0.05` | dot transparency |
| `maps.agreement` | float | `0.8` | fraction of models that must agree on the sign |
| `maps.regions` | bool | `true` | draw `index.region(s)` boxes |
| `maps.extent` | `[lon0,lon1,lat0,lat1]` | `null` | map extent (default: data extent) |
| `maps.projection` | `PlateCarree` | `PlateCarree` | cartopy projection |
| `kde.grid` | `[rows,cols]` | `[4,6]` | maximum panel grid; adapts down to the number of models |
| `kde.max_models` | int | `24` | cap on panels |
| `kde.mode` | `auto` \| `field` \| `series` | `auto` | `field`: pool all grid points of the study region over the baseline (smooth PDF, the 5.19 style); `series`: use the annual index values |
| `kde.per_model` | bool | `true` | also write one `fig_pdf_<model>.png` per model (observed vs raw vs corrected + KS) |

## `runtime`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `verbose` | bool | `true` | DEBUG logging to stdout |
| `cache` | bool | `true` | reuse yearly fields in `<outdir>/cache` |
| `strict` | bool | `false` | abort the whole run when one model fails |

## Complete example

```yaml
project: {name: ea-zg-projection, outdir: /disk1/wsr/proj/ea-zg/run}

variable: {name: zg, units: gpm, long_name: 300 hPa geopotential height}

data:
  reference:
    paths: ["/disk3/era5-0.25/monthly/pressure/hgt.{year}.grib"]
    years: [1950, 2014]
    variable: z
    plev: 300
    units_transform: {scale: 0.1019716, offset: 0.0}   # m -> gpm
    aggregate: {months: [9, 10], stat: mean}
  models:
    root: /disk5/cmip6/monthly
    variable: zg
    plev: 300
    units_transform: {scale: 1.0, offset: 0.0}
    scenarios: [ssp245, ssp370, ssp585]
    models: [ACCESS-CM2, CESM2, MIROC6, MPI-ESM1-2-LR]
    aggregate: {months: [9, 10], stat: mean}

grid: {lat: {start: 20, stop: 80, step: 1}, lon: {start: 0, stop: 359, step: 1}}

periods:
  historical: [1950, 2014]
  future: [2026, 2100]
  reference_climate: [1950, 2014]
  future_windows: [[2026, 2045], [2046, 2065], [2066, 2085], [2086, 2100]]

index:
  type: projection_index
  remove_zonal_mean: true
  region: {lat: [25, 80], lon: [30, 150]}
  projection: {template_source: reference, template_period: [1995, 2014]}

correction: {method: mean_std, lower_bound: null}

ensemble: {statistic: mean, bands: [10, 90, 20, 80]}

extremes: {enabled: true, distribution: gev, return_periods: [2, 10, 50, 100]}

plots: {dpi: 300, figures: [combined, maps, kde, shape_box]}
```
