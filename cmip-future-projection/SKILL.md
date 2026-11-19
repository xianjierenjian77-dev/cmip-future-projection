---
name: cmip-future-projection
description: "Bias-correct CMIP6/CMIP5 (or any multi-model) projections against an observational reference and turn them into a complete future-projection analysis: seasonal/annual index series, multi-model ensemble statistics, GEV return-period/extreme analysis, and publication figures (time series with spread, maps with agreement stippling, KDE skill diagnostics, shape-parameter comparison). Use when a task involves 模式偏差校正 / bias correction / 未来预估 / future projection / 情景预估 / SSP scenarios / CMIP6 数据处理 / 重现期 / return period / GEV / 极端事件 / multi-model ensemble / 集合平均 / 预估图 / projection figures, for any variable (temperature, precipitation, wind, geopotential height, humidity, radiation, custom index), any region, any season, gridded or pre-computed series."
---

# CMIP future-projection pipeline

One YAML file describes a run; one command produces bias-corrected fields, index
series, ensemble statistics, a GEV analysis, figures, CSV tables and a written
`report.md`. The code is deliberately generic: **any variable, any region, any
season, any index, any model ensemble** — you supply the data paths and the
scientific choices, the pipeline does the rest.

Skill root: `/disk1/wsr/.dsh/skills/cmip-future-projection`
(mirrored at `/disk1/wsr/.agents/skills/cmip-future-projection` for other agent
tools; `~/.dsh/skills` is the same directory on this server because `$HOME` is
`/disk1/wsr`)
Python (this server): `/disk1/wsr/miniconda3/envs/climate/bin/python`

If this skill does not appear in the skill catalog of an older session, read
`SKILL.md` directly from the path above and follow it — the scripts do not depend
on the catalog.

```bash
PY=/disk1/wsr/miniconda3/envs/climate/bin/python
SK=/disk1/wsr/.dsh/skills/cmip-future-projection
$PY $SK/scripts/run.py doctor            # environment check, run once per session
$PY $SK/scripts/run.py init  myrun       # writes myrun/config.yaml + README
$PY $SK/scripts/run.py run   myrun/config.yaml
```

Never re-implement the pipeline by hand in a notebook: edit the config, run the
CLI, read `runs/<name>/report.md`. If the config cannot express what is needed,
extend the package (see *Extending* at the end) instead of writing a one-off
script that the next agent cannot reuse.

## What you must establish before writing the config

Ask the user — or infer from the data — these eight things. Do not guess silently;
record the answers in the config comments.

| # | Question | Config key |
|---|---|---|
| 1 | Variable, target units, and the units in the files | `variable.name/units`, `data.*.units_in`, `data.*.units_transform` |
| 2 | Which season / months, and what aggregation (mean, sum, max, percentile, threshold count) | `data.*.aggregate` |
| 3 | Region(s) and what "one number per year" means (index type) | `index.*` |
| 4 | Observational reference: files, variable name, years, grid | `data.reference.*` |
| 5 | Model archive: root, layout, member(s) and whether to average them, experiments, models, scenarios | `data.models.*`, `data.models.member_combine` |
| 6 | Baseline period for fitting and for anomalies; projection period | `periods.*` |
| 7 | Bias-correction method | `correction.method` |
| 8 | Which analyses/figures are wanted (extremes? maps? KDE?) | `extremes.enabled`, `plots.figures` |

If a required path does not exist, stop and ask the user; never fabricate data or
silently substitute a different model/variable.

## The pipeline in one picture

```
reference files ─┐
                 ├─ load → unit convert → month/season aggregate → yearly field → common grid
model files ─────┘
                        │
        correction params fitted on a common baseline (obs vs model historical)
                        │
        corrected historical + corrected future fields  (netCDF per model/scenario)
                        │
        index series (region mean / zonal-anomaly region difference / projection index)
                        │
        ensemble statistics + trend tests (+ optional trend scaling)
                        │
        GEV fits (per model & scenario) → return levels, shape comparison, Kruskal–Wallis
                        │
        [daily] daily series → daily bias correction → POT/GPD or annual-maxima GEV
                        │
        figures + CSV tables + report.md
```

Stages: `prepare → correct → ensemble → extremes → daily → figures → report`.
Run a subset with `--stages correct,figures` (e.g. after editing the plotting
section only). A subset run still needs its upstream stages in the *same*
process: `--stages correct,figures` runs both.

## Outputs

```
<outdir>/
  config.resolved.yaml     fully merged config (reproducibility)
  manifest.json            versions, scenarios, models
  report.md                human-readable summary: skill scores, trends, GEV, figure list
  corrected/<model>_historical_corrected.nc
  corrected/<model>_<scenario>_corrected.nc
  params/<model>_params.nc          per-gridpoint transfer functions
  series/index_series.npz           obs + every model series
  ensemble/ensemble_<scenario>.npz  matrix, mean/median, percentile bands
  extremes/gev_<scenario>.npz       curves, shapes, model names
  tables/model_scores.csv, trends.csv, gev_fits.csv, gev_shape_summary.csv
  figures/fig_*.png                 ensemble figures
  figures/fig_pdf_<model>.png       per-model PDF (observed vs raw vs corrected)
  figures/fig_timeseries_<model>.png  per-model projection (optional)
  figures/fig_daily_return_levels.png daily POT/GPD or block-maxima return levels
  extremes/daily_<scenario>.npz     daily extreme curves
  tables/daily_extremes.csv         per-model GPD/GEV parameters and return levels
  logs/run.log
```

## Choose the index (`index.type`)

| Type | One value per year is | Use for |
|---|---|---|
| `region_mean` | area-weighted regional mean | temperature, precipitation, wind, any regional average |
| `region_mean_anomaly` | same, minus the baseline climatology | warming/anomaly series (pipeline adds `anomaly: true` anyway) |
| `region_diff` | mean(target region) − mean(reference region), both after removing the zonal mean | jet/ridge/wave-train indices, "north minus south" contrasts |
| `projection_index` | projection of the zonal-anomaly field onto a template pattern, `Σ(A′·T·w)/Σ(T²·w)` | wave-train / teleconnection indices (the 6.25 Fig. 5 workflow) |
| `series` | the input already *is* a series (CSV/NPZ, no grids) | station data, pre-computed indices, non-gridded products |

Add `index.anomaly: true` to remove the baseline climatology, `index.detrend`,
`index.standardize` for a z-score, and `index.region`/`index.regions` in
`{lat: [lo, hi], lon: [lo, hi]}` form (longitudes may be 0–360 or −180–180).

## Choose the correction method (`correction.method`)

| Method | What it fixes | Use when |
|---|---|---|
| `none` | nothing | control run, or the ensemble is already bias-corrected |
| `mean` | mean bias | simplest anchoring; keeps model variability and trend |
| `mean_std` | mean + variance/CV | **default**; the method used in the 6.25/5.19 workflow |
| `quantile_mapping` | full distribution (per grid point) | validating the historical period; **not** for future values |
| `quantile_delta_mapping` | full distribution + preserves the model's quantile change | **recommended for extremes/precipitation projections** |
| `delta_change` | same arithmetic as `mean` | when the user/papers call it "delta change" |

Set `correction.multiplicative: auto` (default) and the pipeline uses a
multiplicative correction for precipitation/wind-like variables that are
non-negative. For strictly positive variables also set `correction.lower_bound: 0`
so corrections cannot create negative rain.

`correction.trend_correction.enabled: true` rescales each model's future trend to
the observed historical trend (`mode: difference`, the Fig. 5 approach). Use it
only when the scientific question is explicitly "what if the observed trend
continues"; it overrides the model's own trend and must be stated in the methods.

## Choose the extreme analysis (`extremes.*`)

GEV is fitted to one value per year (block maxima or an annual index) per model
and per scenario. `extremes.enabled: true` produces:
return-level curves, a multi-model median + 10–90/20–80 bands, per-model shape
parameters, a Kruskal–Wallis test across scenarios, and the return period of a
reference value (e.g. the 2025 observed record).

* `distribution: gev` (default) or `gumbel`.
* `min_years: 20` — fits with fewer years are skipped (they are unstable).
* `reference_value` — the value whose return period is reported (defaults to
  `plots.timeseries.reference_value`, i.e. the observed value at the reference year).
* `fit_window: {historical: [1950, 2024], future: [2026, 2100]}` — fit windows.
* `source: future|historical` — `historical` reports only the observed fit.
* `show_empirical: true` adds Gringorten plotting positions to the figure.

## Daily extremes (optional)

The annual pipeline has one value per year, so a heatwave or a storm cannot be
represented. `daily.enabled: true` adds an independent daily analysis: read the
daily field for a season, reduce it to a daily regional index, bias-correct that
series, then

* `daily.extremes.kind: pot` — peaks over a threshold, declustered by runs,
  GPD fit with a Poisson rate; return levels in years.
* `daily.extremes.kind: block` — annual maxima of the daily series, GEV fit.

`threshold_mode` is the key scientific choice: `observed` (one fixed threshold
from the historical baseline, comparable across models), `quantile` (each series'
own quantile, constant exceedance rate under warming) or `absolute` (a physical
threshold such as 35 °C). Always set `data.*.bbox` for daily data — the loader
refuses global daily fields.

Full recipe and caveats: `references/variables.md`, `references/methods.md` §8.

## Multi-member ensembles

`data.models.members: [r1i1p1f1, r2i1p1f1, r3i1p1f1]` together with
`member_combine: mean` (or `median`) loads every member of a model, aligns the
years and averages them per grid point before any correction. With
`member_combine: first` (default) only one member is used. The report records how
many members went into each model. Averaging members reduces internal
variability but also narrows the spread between models — state which one you did.

## Figures

`plots.figures` accepts any of: `timeseries`, `gev`, `combined` (3×2 panel like
6.25 Fig. 4), `maps`, `kde`, `shape_box`. Set `plots.dpi: 300` for publication.
Chinese labels are switched on automatically when the labels contain CJK text
(WenQuanYi fonts are installed); override with `plots.font`.
Details and per-journal sizing: `references/plotting.md`.

## Verify before you report back

1. `report.md` — raw vs corrected bias should drop by an order of magnitude; KS
   should drop; `sd_ratio_corr` should be near 1.
2. `figures/fig_timeseries.png` — history joins the future at the reference year,
   ensemble mean inside its bands, trend annotations have units `/10a`.
3. `tables/gev_fits.csv` — every `ks_p > 0.05` and `shape` within a plausible
   range (|c| < 1); `n` equals the number of years.
4. `figures/fig_maps_*.png` — stippling only where models agree (default ≥80%).
5. State in the answer: variable, units, season, region, index, method, period,
   scenarios, number of models, and where the outputs are.

## Hard rules

* **Units and calendars first.** A factor-of-1000 precipitation error or a
  360-day calendar silently ruins everything. Print min/max/mean of the loaded
  field before correcting; if the values are implausible (e.g. tas in °C but
  configured as K), fix the config, do not "fix" the science.
* **Never apply `quantile_mapping` to future values** unless the user explicitly
  wants a stationary distribution; it removes the climate-change signal in the
  tails. Use `quantile_delta_mapping`.
* **Never mix members** (e.g. r1i1p1f1 and r101i1p1f1) in one model; the loader
  warns and keeps the member with the most files.
* **Do not delete the run directory** to "start clean" without telling the user;
  it contains the only record of the configuration used.
* Every number that goes into the final answer must come from `report.md` or the
  CSV tables — not from memory of a previous run.

## Troubleshooting (first aid)

| Symptom | Cause / fix |
|---|---|
| `no models discovered` | `data.models.root` or `pattern(s)` wrong; run with `runtime.verbose: true` and check the glob |
| `no input files matched` | member/scenario string mismatch; set `members: ["*"]` to auto-pick |
| empty/NaN series | the index baseline does not overlap the field; check `periods.reference_climate` |
| all values 1000× too large | `units_transform: {scale: 0.001}` for precipitation in mm vs m |
| `Reindexing only valid with uniquely valued Index` | duplicate longitudes (0 and 360); the loader now de-duplicates — update the skill if seen |
| GRIB fails to open | `cfgrib`/eccodes missing; convert to NetCDF first |
| cartopy downloads hang | offline server; coastlines are cached, do not clear `~/.local/share/cartopy` |
| figures empty | wrong `plots.figures` name; valid names are listed above |

Full list: `references/troubleshooting.md`.

## Where to read more

| Need | File |
|---|---|
| Step-by-step operating procedure with checkpoints | `references/workflow.md` |
| Equations, assumptions, when each method is valid | `references/methods.md` |
| Every configuration key, type and default | `references/config-reference.md` |
| Figure types, styling, panel composition | `references/plotting.md` |
| Per-variable recipes (tas, pr, wind, zg, humidity, radiation, ice) | `references/variables.md` |
| Errors and fixes | `references/troubleshooting.md` |
| Ready-to-edit configs (grid projection, regional extremes, precipitation QDM, series-only, daily POT) | `templates/*.yaml` |

## Extending the skill

* New data layout → add a pattern to `data.models.patterns`; no code change.
* New aggregation (e.g. heatwave duration) → add a `stat` branch in
  `cmipkit/io.py::_reduce` and one row in `references/variables.md`.
* New index → add a branch in `cmipkit/indices.py::build_index` and document it
  in `references/methods.md` and the table above.
* New extreme model (GPD variants, r-largest, non-stationary GEV) → extend
  `cmipkit/extremes.py` / `cmipkit/daily.py` and add a `daily.extremes.kind` branch.
* New figure → add a function in `cmipkit/plotting.py`, register the name in
  `cmipkit/config.py::validate_config` and call it from `Pipeline.figures`.
* After any change run **all three**: `run.py selftest`, `tests/test_methods.py`
  and `tests/test_daily.py`; a change is not done until all pass.
