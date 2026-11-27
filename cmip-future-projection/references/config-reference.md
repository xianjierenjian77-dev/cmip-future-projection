# Configuration reference

Contents: 1. Resolution and CLI; 2. Project and variables; 3. Inputs;
4. Grid and periods; 5. Indices; 6. Correction; 7. Ensemble; 8. Extremes;
9. Daily; 10. Plots and runtime; 11. Manifest and saved outputs; 12. Raw input quality.

The code's `DEFAULTS` in `scripts/cmipkit/config.py` is authoritative. Validate a
configuration before a production run. Defaults are a safe starting point, not
an inferred scientific design. The templates explicitly choose example regions,
periods and methods; replace them for the task.

## 1. Resolution and CLI

Resolve relative paths against the YAML file's directory. Expand `~` and
`$ENVIRONMENT_VARIABLE` in path values. Project placeholders `{name}`,
`{root}`, `{outdir}` and `{variable}` resolve during configuration loading;
file placeholders `{model}`, `{experiment}`, `{member}` and `{year}` resolve
when discovering inputs. Quote paths containing YAML-special characters.

~~~bash
python /path/to/skill/scripts/run.py init myrun --template grid
python /path/to/skill/scripts/run.py validate myrun/config.yaml --check-inputs
python /path/to/skill/scripts/run.py validate myrun/config.yaml --check-data --report-json myrun/input-quality.json
python /path/to/skill/scripts/run.py run myrun/config.yaml --dry-run
python /path/to/skill/scripts/run.py run myrun/config.yaml
~~~

Template choices: `grid`, `series`, `precip`, `regional-extremes`, `daily`.
`validate` checks the configuration; `--check-inputs` additionally checks discovery
and every selected file's header/time axis. `--check-data` includes those checks
and scans selected raw values in chunks. Both `validate` and `run` accept
`--report-json PATH`; this explicitly writes the inspection report, including
failed checks. A report's relative path resolves from the shell's working directory.
`run --dry-run` performs preflight without projection calculations; configured
full scans or `--check-data` can still read all selected input values.
`run --workers N` overrides `runtime.workers`; `--worker-threads N` overrides
`runtime.worker_threads` for that invocation. Both require positive integers.
Record the resolved configuration when comparing runs.
`--stages` includes required upstream stages in the same run; it does not reload
all saved analysis state. Use `--overwrite` only when intentionally replacing
outputs with a changed scientific configuration. The previous run directory is
moved to a sibling named with `.previous-<timestamp>` before the new run starts;
it is not silently deleted.

## 2. Project and variables

| Key | Default | Meaning |
|---|---|---|
| `project.name` | `cmip_projection` | Run label; avoid path separators. |
| `project.outdir` | `runs/{name}` | Output directory relative to this YAML, or an absolute path. |
| `project.description` | empty | Scientific purpose and decisions. |
| `variable.name` | `variable` | Output variable label. |
| `variable.units` | empty | Actual units after conversion and aggregation. |
| `variable.long_name` | empty | Descriptive label. |

## 3. Inputs

Use matching `kind` values for reference/models. Grid inputs need rectilinear
one-dimensional latitude/longitude coordinates. Native GeoTIFF, Zarr stores,
station daily CSV, rotated/curvilinear grids and unstructured meshes need adapters.

### Shared source options

| Key below `data.reference` / `data.models` | Default | Meaning |
|---|---|---|
| `kind` | `grid` | `grid` or precomputed annual `series`. |
| `variable` | null | Variable inside the input file; explicitly set for multi-variable datasets. |
| `units_in` | null | Expected source units; checked against available file metadata with a limited set of equivalent spellings. Conversion still requires `units_transform`. |
| `units_transform` | `{scale: 1, offset: 0}` | Affine conversion before aggregation. |
| `plev` | null | Exact pressure-level selection within numerical tolerance, with hPa/Pa conversion; missing levels fail. |
| `bbox` | null | `{lat: [south, north], lon: [west, east]}`; subset before loading/regridding. |
| `aggregate.months` | `ANN` | Ordered month list or `DJF/MAM/JJA/SON/ANN`. |
| `aggregate.stat` | `mean` | `mean,sum,max,min,median,std,p05,p10,p90,p95,p99,count_gt,count_ge,count_lt,count_le,max_abs,mean_abs`. |
| `aggregate.threshold` | null | Threshold in transformed input units for count statistics. |
| `aggregate.frequency` | `auto` | `auto,monthly,daily,annual`; declare if detection is ambiguous. |
| `aggregate.weight_by_days` | true | Weight monthly mean/mean_abs by month length. Does not integrate rates for sum. |
| `aggregate.require_complete` | true | Exclude seasons lacking required months/daily steps. |
| `aggregate.min_valid_fraction` | 1.0 | Required finite time-sample fraction per grid cell, between 0 and 1. |

Counts count timesteps, not automatically days. The maximum of monthly means
is not a daily maximum. Missing-value policy changes the sample and must be recorded.

### Reference

| Key below `data.reference` | Default | Meaning |
|---|---|---|
| `paths` | `[]` | Local NetCDF/GRIB paths or globs; `{year}` expands over years. |
| `years` | null | Requested reference range `[start,end]`. |
| `series_file` | null | Required CSV/NPZ file for series input. |
| `series_columns` | `{year: year, value: value}` | CSV column names. |

CSV input contains one row per integer year. NPZ input contains `years` and
`values` arrays. A saved multi-model `index_series.npz` is not automatically the
same schema; export the intended single series first. Units and event definitions
must be compatible across reference and model series.

### Models

| Key below `data.models` | Default | Meaning |
|---|---|---|
| `root` | null | Grid archive root. |
| `patterns` | `[]` | File-layout patterns; specify them for nonstandard layouts or daily tables. |
| `pattern` | null | Legacy single-pattern option. |
| `historical_experiment` | `historical` | Baseline experiment name. |
| `scenarios` | `[ssp245]` | Future experiments, with no assumed probability weights. |
| `scenario_labels` | `{}` | Display-name overrides. |
| `members` | `[r1i1p1f1]` | Member preference/list; wildcard discovery needs explicit review. |
| `member_combine` | `first` | `first,mean,median`. Averaging alters variability and extremes. |
| `models` | `all` | Discovery or an explicit reproducible list. |
| `exclude` | `[]` | Declared exclusions, supported by reasons. |
| `series_dir` | null | Model series archive. |
| `series_columns` | `{year: year, value: value}` | Model CSV column names; model `units_transform` is also applied. |

Typical monthly pattern:

~~~text
{root}/{experiment}/{model}/{variable}_Amon_{model}_{experiment}_{member}_*.nc
~~~

A common alternative adds a `{variable}` directory below `{model}`. Daily data
need a `day` pattern. Explicit patterns do not silently fall back to other
layouts, and an explicit member list does not select an unlisted substitute.
Historical/scenario members are kept consistent. Constrain globs to one data version and grid label; do not
concatenate replicas. Standard series layout is `<series_dir>/<model>/<experiment>.csv`
or the corresponding NPZ; flat `<model>_<experiment>` files are also recognized.

## 4. Grid and periods

| Key | Default | Meaning |
|---|---|---|
| `grid.lat`, `grid.lon` | null | Reference axes, explicit arrays, or `{start,stop,step}` with inclusive stop. |
| `grid.method` | `linear` | `linear` or `nearest`; neither is conservative remapping. |
| `grid.lon_360` | true | Longitude convention. |
| `periods.historical` | `[1950,2014]` | Model historical loading range. |
| `periods.future` | `[2026,2100]` | Scenario loading range; choose for the scientific task. |
| `periods.reference_climate` | `[1950,2014]` | Baseline for anomalies/map changes. |
| `periods.validation` | null | Optional `[start,end]` window for annual historical model-score diagnostics on common years. |
| `periods.future_windows` | four default windows | Averaging windows for maps; set explicitly. |

Do not fabricate observations or silently bridge 2015–2025 if the requested
projection begins in 2026. Historical and scenario experiments have distinct
coverage. Read the actual years used, especially after completeness filtering.

`periods.validation` selects the annual diagnostic sample; it does not change the
calibration period, fetch more observations, or automatically make diagnostics
out-of-sample. An overlap with calibration generates a warning. Scores require
enough common valid annual observations. Field-pooled KDE diagnostics still have
their own sampling definition; do not describe every figure as a held-out score
merely because this window was supplied.

## 5. Indices

| Key below `index` | Default | Meaning |
|---|---|---|
| `type` | `region_mean` | `region_mean,region_mean_anomaly,region_diff,projection_index,series`. |
| `region` | null | `{lat,lon,name?}` for the target region. |
| `regions` | null | `[reference_region,target_region]` for region_diff. |
| `remove_zonal_mean` | false | Remove the longitude mean before reduction. |
| `anomaly` | false | Subtract historical climatology, preserving future changes. |
| `detrend` | false | Deliberately remove the input segment's linear trend. |
| `standardize` | false | Reuse historical mean/std; do not standardize future periods independently. |
| `units_label` | null | Axis-unit override, e.g. dimensionless projection coefficient. |
| `projection.template_source` | `reference` | `reference,year,model` where supported by supplied data. |
| `projection.template_period` | baseline | Historical pattern period. |
| `projection.template_year` | null | Observed template year for `year`. |
| `projection.template_model` | null | Requested historical model pattern for `model`. |

Series mode requires its annual index to be preprocessed upstream; do not combine
it with `index.anomaly`, `standardize` or `detrend` and assume those are applied.
Use `methods.md` to choose the operator order and projection normalization.

## 6. Correction

| Key below `correction` | Default | Meaning |
|---|---|---|
| `method` | `none` | `none,mean,mean_std,delta_change,quantile_mapping,quantile_delta_mapping`. |
| `params_period` | baseline | Historical calibration interval. |
| `apply_to` | `[historical,future]` | Adjust selected segments; document inconsistent-scale choices. |
| `multiplicative` | `auto` | true/false or variable/data heuristic; set explicitly for final work. |
| `lower_bound`, `upper_bound` | null | Optional physical clipping; not a wet-day model. |
| `n_quantiles` | 100 | Quantile-table resolution; does not increase independent sample size. |
| `trend_correction.enabled` | false | Imposed index-trend sensitivity experiment. |
| `trend_correction.mode` | `difference` | `difference` or `ratio`. |

QDM estimates its future CDF from the configured future window. Mean_std scales
trends as well as variability. Neither is a universal default for every variable.

## 7. Ensemble

| Key below `ensemble` | Default | Meaning |
|---|---|---|
| `statistic` | `mean` | Mean or median over models. |
| `bands` | `[10,90,20,80]` | Model-spread percentiles; plotting expects two ordered lower/upper pairs. |
| `weights` | `equal` | `equal` or heuristic historical `skill` weights. |
| `exclude_models` | `[]` | Explicit exclusions after discovery. |
| `min_models` | 3 | Required contributing models per year; if no year meets it, the ensemble stage fails. Set explicitly for a pilot; a single-model run does not estimate multi-model spread. |

## 8. Annual extremes

| Key below `extremes` | Default | Meaning |
|---|---|---|
| `enabled` | false | Explicit opt-in for annual GEV/Gumbel analysis. |
| `sample_kind` | `annual_index` | `block_maxima` or explicitly justified empirical `annual_index`. |
| `source` | `future` | Future plus observed comparator, historical-only, or both as supported. |
| `fit_window` | null | `{historical: [a,b], future: [c,d]}`. |
| `distribution` | `gev` | `gev` or `gumbel`. |
| `fit_method` | `mle` | `mle` or `mm` (method of moments); inspect the recorded fallback. PWM is not implemented. |
| `return_periods` | configured list | Years for one complete block per year; all must exceed 1. |
| `reference_value` | null | Value whose return period is requested. |
| `bootstrap` | 0 | Nonzero pipeline bootstrap is rejected; no hierarchical fit interval is implemented. |
| `show_empirical` | false | Empirical plotting positions. |
| `min_years` | 20 | Minimum valid blocks; not a precision guarantee. |
| `shape_comparison` | true | Exploratory scenario comparison with dependence caveats. |

`ks_p_naive` is an uncalibrated fitted-sample diagnostic. `ks_p` is missing unless
a valid refitted-distribution test was actually performed. Do not filter to p > 0.05.

## 9. Daily extremes

| Key below `daily` | Default | Meaning |
|---|---|---|
| `enabled` | false | Read daily gridded inputs; provide a bounded study domain. |
| `months` | `[6,7,8]` | Selected seasonal days. |
| `index` | `region_mean` | `region_mean` or `region_diff`; daily anomaly climatology is not implemented. |
| `baseline` | baseline | Reference/model correction and historical extreme-fit window. |
| `correction` | `quantile_delta_mapping` | Daily series adjustment, separate from the annual route. |
| `extremes.kind` | `pot` | `pot` or `block`. |
| `extremes.threshold_mode` | `observed` | `observed,absolute,quantile`. |
| `extremes.threshold_quantile` | 0.95 | Probability for observed/quantile thresholds. |
| `extremes.threshold_abs` | null | Required for absolute threshold, in converted units. |
| `extremes.decluster_days` | 3 | Exceedance gaps at most this many days join a cluster only over continuous finite observations. |
| `extremes.min_peaks` | 30 | Minimum declustered excesses. |
| `extremes.min_years` | 20 | Minimum valid annual/seasonal blocks. |
| `extremes.distribution` | `gev` | GEV/Gumbel for block maxima. |
| `extremes.return_periods` | `[2,5,10,20,50,100]` | Return levels to report. |
| `plots` | true | Produce the daily return-level figure. |

Daily extrema always concern the configured regional series. They are not
heatwave durations or pointwise extremes. The daily route requires
`data.models.member_combine: first` and locks the selected historical member for
its scenario files. It rejects daily anomaly requests instead of silently
subtracting an unrelated climatology. Preserve actual timestamps for clustering
and use a justified completeness/exposure rule.

## 10. Plots and runtime

See `plotting.md` for figure choices and interpretation. Principal controls:
`plots.figures` (default `[timeseries,maps,kde]`), `dpi` (300), `formats`
(`[png]`), `figsize_scale` (1), `font` (null), `language` (`en`),
`colors`, `timeseries`, `gev`, `maps` and `kde`.

| Key | Default | Meaning |
|---|---|---|
| `plots.timeseries.reference_year` | null | Optional marker; do not invent a join or missing observed value. |
| `runtime.verbose` | true | Detailed logs. |
| `runtime.cache` | true | Reuse eligible prepared fields. |
| `runtime.workers` | 1 | Maximum simultaneous model workers for correction and the daily branch, excluding the coordinator and runtime helper processes. One worker or one selected model keeps serial execution; active workers cannot exceed the number of models. `run --workers N` overrides this value. |
| `runtime.worker_threads` | 1 | Limit supported BLAS/OpenMP numerical thread pools within each spawned worker. Positive integer; `run --worker-threads N` overrides it. This is separate from the number of model processes and is not a system-wide thread limit. |
| `runtime.max_input_values` | 30000000 | Guard selected loads before reading: one season-year or a cumulative daily field. Also bounds combined annual fields/cache reads. This is an array-size guard, not a total memory budget. |
| `runtime.strict` | true | Fail on processing errors instead of quietly losing models. |
| `runtime.strict_config` | true | Reject unknown keys; false logs warnings for deliberate compatibility review. |
| `runtime.overwrite` | false | Protect an existing run from a changed scientific configuration. |

Do not relax strictness just to obtain a successful run. Fix the inputs or state
and justify the intended exclusions. Preserve failed-run logs for diagnosis.
Read [parallel.md](parallel.md) for process scope, memory/I/O planning and safe
Python entry points. `max_input_values` applies within each process; multiplying
workers does not reserve memory or raise that limit. Changing worker/thread counts
is an operational change and does not require `--overwrite` by itself.


## 11. Manifest and saved outputs

`manifest.json` distinguishes `preflight_inputs` (candidate files discovered
before computation) from `actual_inputs` (files consumed by the run). The legacy
`inputs` entry mirrors the actual-input list. A preflight match is not proof that
a file participated in the result.

Use `artifacts_current_invocation` for outputs generated or changed by this
invocation and `artifacts_not_produced_in_this_invocation` for unchanged existing
files. The manifest also records status, completed stages, failures, loaded
models/members and any `previous_run_directory` made by an explicit overwrite.
Do not report an old figure as newly generated solely because it is present in
the output directory.

The `parallel` entry records requested concurrency and per-stage effective
workers, process IDs, task status/timing and worker log paths. Inspect each
stage's `start_method` (`serial` or `spawn`) for its actual execution mode;
the top-level `parallel.start_method: spawn` names the parallel strategy only.
Task records describe received results and can be incomplete after interruption.
See [parallel.md](parallel.md) for field meanings and timing scope.

`quality_report.json` records input inspection separately from statistical
results. The manifest's `quality` entry records the report path, scan level,
whether values were checked, pass status and summary counts. The Markdown report
links to the full quality report. Preserve failed/interrupted-run status and logs;
completed files do not imply that later stages succeeded.

## 12. Raw input quality

| Key below `quality` | Default | Meaning |
|---|---|---|
| `scan` | `headers` | Production preflight level: `headers` or `full`. `full` scans all values within the configured selection. |
| `coverage` | `warn` | Report missing requested months/days/years as `warn` or `error`; use the source's calendar. |
| `require_units` | false | Missing units metadata warns by default; true makes it an error for gridded input. |
| `max_missing_fraction` | null | Optional maximum nonfinite fraction across selected time/space values, including static masks and overlapping file records. |
| `valid_min`, `valid_max` | null | Optional lower/upper value limits in transformed input units, before temporal aggregation. Violations are errors. |
| `chunk_size` | 32 | Maximum time samples per value-scan chunk. |
| `max_chunk_values` | 2000000 | Maximum values per value-scan chunk; additional spatial splitting may be used. |

Setting a missingness or value limit forces a full scan in production, even with
`scan: headers`. `validate --check-inputs` remains a header check; use
`validate --check-data` to evaluate value limits without running projections.
No universal physical range is inferred. The checker does not overwrite inputs,
fill missing observations, or clip values. See [quality-control.md](quality-control.md)
for coverage rules, static-mask interpretation, examples and report scope.
