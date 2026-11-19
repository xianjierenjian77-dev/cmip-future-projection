# Per-variable recipes

The pipeline is variable-agnostic; what changes is the units, the aggregation,
whether the correction is additive or multiplicative, and whether a GEV fit is
meaningful. Copy the block that matches and adapt the paths.

General rules:

* Convert to the units you will report **before** correcting (`units_transform`),
  so the correction parameters are in reporting units.
* Seasonal precipitation: use `sum`, not `mean`, over monthly data.
* Anything non-negative: `correction.lower_bound: 0`.
* Extremes need one value per year that is (or is close to) a block maximum:
  `max`, `min`, `p95`, `count_gt`, or a composite index.

---

## 2 m temperature (`tas`, `t2m`)

```yaml
variable: {name: tas, units: degC, long_name: 2 m temperature}
data:
  reference:
    paths: ["/disk3/era5-1.0/monthly/surface/t2m.{year}.grib"]
    variable: t2m
    units_transform: {scale: 1.0, offset: -273.15}
    aggregate: {months: "JJA", stat: mean}
  models:
    root: /disk5/cmip6/monthly
    variable: tas
    units_transform: {scale: 1.0, offset: -273.15}
    aggregate: {months: "JJA", stat: mean}
correction: {method: mean_std, multiplicative: false}
```

For **TXx / TNn** (hottest day / coldest night) use daily data and
`stat: max` / `min`; then a GEV fit is a genuine block-maxima analysis. For
"days above 35 °C" use `stat: count_gt`, `threshold: 35`, `correction.method:
quantile_delta_mapping` (counts are badly served by mean corrections).

## Precipitation (`pr`, `tp`, `prec`)

```yaml
variable: {name: pr, units: mm/day, long_name: precipitation}
data:
  reference:
    paths: ["/disk3/era5-1.0/monthly/surface/precip/prec.{year}.grib"]
    variable: tp
    units_transform: {scale: 86400.0, offset: 0.0}   # m/s -> mm/day
    aggregate: {months: [9, 10], stat: sum}
  models:
    root: /disk5/cmip6/monthly
    variable: pr
    units_transform: {scale: 86400.0, offset: 0.0}   # kg m-2 s-1 -> mm/day
    aggregate: {months: [9, 10], stat: sum}
correction:
  method: quantile_delta_mapping
  multiplicative: true
  lower_bound: 0
extremes: {distribution: gev}      # precipitation maxima are usually Frechet (c<0)
```

Notes: ERA5 `prec` is accumulated, CMIP6 `pr` is a flux — check the magnitude
(`order 1e-5` means flux, `order 1e-4` m means accumulated). Many zeros break
multiplicative scaling; QDM handles them. For a GEV of seasonal totals, use
`stat: sum` per year — that is a block total, not a maximum; state it.

## Wind speed (`sfcWind`, `wsped`)

```yaml
variable: {name: sfcWind, units: m/s}
data:
  models: {variable: sfcWind, aggregate: {months: [12, 1, 2], stat: p95}}
correction: {method: quantile_delta_mapping, multiplicative: true, lower_bound: 0}
extremes: {distribution: gev}
```

`p95` of daily wind is a defensible extreme index; `max` is noisier. Wind speed is
non-negative — never use an additive correction without `lower_bound: 0`.

## Geopotential height (`zg`, `hgt`)

```yaml
variable: {name: zg, units: gpm, long_name: 300 hPa geopotential height}
data:
  reference:
    paths: ["/disk3/era5-0.25/monthly/pressure/hgt.{year}.grib"]
    variable: z
    plev: 300                       # hPa for ERA5
    units_transform: {scale: 0.1019716, offset: 0}   # m2/s2 -> gpm
    bbox: {lat: [20, 80], lon: [40, 160]}
  models:
    variable: zg
    plev: 300                       # the loader converts Pa -> hPa automatically
    units_transform: {scale: 1.0, offset: 0}
correction: {method: mean_std}
```

For a wave-train study use `index.type: projection_index` with a template from the
reference climatology, as in the 6.25 Fig. 5 workflow. For a zonal-anomaly
difference index use `index.type: region_diff` with two regions; the pipeline
removes the zonal mean automatically.

## Mean sea level pressure (`psl`, `msl`, `slp`)

`method: mean_std`, `multiplicative: false`. For storm-track indices use
`stat: std` (sub-seasonal variability) or a `count_lt` threshold (deep lows).

## Humidity (`hurs`, `hus`, `huss`)

Additive correction; `units_transform` for `%` (hurs is already %) or `g/kg`
(hus is kg/kg → ×1000). Avoid multiplicative correction on relative humidity near
100 % (it can exceed 100 — set `upper_bound: 100`).

## Radiation and fluxes (`rsds`, `rlds`, `nsr`, `ssrd`)

Additive correction, `stat: mean` or `sum` for accumulated fluxes. Fluxes can be
negative (net) or strictly positive (downward shortwave): set `lower_bound: 0`
only for the latter. Switching sign conventions between ERA5 and CMIP6 is a common
error — compare seasonal means of both before correcting.

## Sea surface temperature / sea ice

SST: `method: mean_std`, mask land (the loader keeps NaN). Sea ice concentration:
`method: quantile_mapping` on `%` with `lower_bound: 0`, `upper_bound: 100`; for
"ice-free days" use `count_lt` with `threshold: 15`.

## Pre-computed series (no grids)

Use `kind: series` when the index already exists (station records, published
indices, a previous run's `index_series.npz` exported to CSV):

```yaml
data:
  reference:
    kind: series
    series_file: /path/obs_index.csv        # columns year,value
    years: [1950, 2014]
  models:
    kind: series
    series_dir: /path/models/               # <model>/historical.csv, <model>/ssp245.csv
    scenarios: [ssp245, ssp585]
correction: {method: mean_std}
index: {type: series}
grid: {}                                    # not needed
plots: {figures: [timeseries, gev, combined, shape_box]}
```

Series mode exercises the same correction, ensemble, GEV and plotting code, so a
figure produced from series and from grids is directly comparable.

## Daily extremes (heatwaves, heavy-rain days, wind gusts)

The annual recipe above compresses the season to one number per year. When the
question is about *days* — a 35 °C threshold, a 3-day heatwave, a 50 mm/day rain
event — use the `daily` block on daily input:

```yaml
data:
  reference:
    paths: ["/disk3/era5-0.25/daily/t2m_max.{year}.nc"]
    years: [1985, 2014]
    bbox: {lat: [30, 48], lon: [100, 130]}
  models:
    root: /disk5/cmip6/daily
    patterns: ["{root}/{experiment}/{model}/tasmax_day_{model}_{experiment}_{member}_*.nc"]
    variable: tasmax
    members: [r1i1p1f1, r2i1p1f1, r3i1p1f1]
    member_combine: mean
    bbox: {lat: [30, 48], lon: [100, 130]}

daily:
  enabled: true
  months: [6, 7, 8]
  index: region_mean
  correction: quantile_delta_mapping
  extremes: {kind: pot, threshold_mode: absolute, threshold_abs: 35.0, decluster_days: 3}
```

* **Heatwaves** — `kind: pot`, `threshold_mode: absolute`, `threshold_abs` = the
  physical definition (35 °C), `decluster_days: 3`. The Poisson rate then reads
  as "events per summer", and return levels as "the temperature of a 1-in-20-year
  event".
* **TXx / TNn** — `kind: block`, or the annual pipeline with `stat: max` / `min`;
  both give annual maxima, the block route avoids a threshold choice.
* **Heavy precipitation** — daily `pr`, `threshold_mode: observed` (a fixed
  percentile of the historical distribution), `multiplicative` handling is
  automatic for precipitation; expect a heavy (negative shape) GPD tail.
* **Wind gusts** — same POT recipe with `threshold_mode: quantile` to keep the
  event rate stable across models.

Memory: a global 0.25° daily field is far too large; always set `bbox`. The stage
refuses fields above ~30 M values.

## Variables where GEV should be switched off

* Annual **means** of a smooth field — a GEV on an annual mean is a tail model for
  a non-extreme quantity; either switch `extremes.enabled: false` or describe it
  explicitly as "GEV fitted to the annual index".
* Series with fewer than `extremes.min_years` values (fits are skipped anyway).
* Indices bounded in a narrow range (e.g. a correlation-based index in [−1, 1]).
