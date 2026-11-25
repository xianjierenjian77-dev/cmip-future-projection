---
name: cmip-future-projection
description: "Prepare, inspect, run, or repair reproducible CMIP6/CMIP5 future-projection analyses using local rectilinear NetCDF/GRIB fields or annual CSV/NPZ series. Use for raw CMIP quality checks, future scenarios, seasonal or regional indices, explicitly chosen bias adjustment, multi-model summaries, projection figures, and optional daily POT/GPD or block-maxima analysis. Also use when a user provides an existing CMIP configuration or asks to generalize a projection workflow. 支持原始 CMIP 数据质量检查、CMIP6 未来预估、情景比较、偏差订正、区域指数和极值分析；不将年度均值自动视为块极值。"
---

# CMIP future projections

Turn the user's scientific definition and existing data into a reproducible, inspectable analysis.
Use the bundled CLI where its supported methods fit. Keep the main task moving; request only missing choices that would change the scientific meaning.

## 1. Choose the smallest useful workflow

| User has / wants | Action |
| --- | --- |
| Existing configuration or failed run | Inspect the configuration, log and manifest; reproduce the failure on a small subset before editing. |
| Raw-file quality review | Use all-file input checks, then a chunked value scan where needed; return the report without starting projections. |
| A new analysis with local fields | Start with the grid template; establish variable, units, calendar, aggregation and index. |
| Precomputed yearly indices | Use the series template; do not invent gridded observations or apply index preprocessing again. |
| Daily regional extremes | Use a daily template; reduce each day's field to the requested regional series before fitting POT or block maxima. |
| Only a plan or method review | Review the inputs and assumptions; do not start a full analysis or download an archive. |
| Data not yet available | Prepare a concrete configuration and list the exact missing data. Discover available paths/sources first; never substitute synthetic data into a scientific result. |
| A different method or unsupported format | State the boundary, use a documented adapter or an isolated extension, and validate its scientific invariants before relying on it. |

Do not automatically enable bias correction, trend adjustment or extreme-value fitting.
Preserve explicit user choices and the conventions of an existing study; make changes to scientific definitions visible.

## 2. Establish the data contract

Recover information from the user's request, existing configuration, file headers and prior results before asking questions.

Record in the configuration:

- Variable name in each source, original units, requested output units and explicit conversion.
- Input frequency, calendar, selected months and aggregation. Distinguish a mean rate from an accumulated total.
- Region, spatial statistic and the order of temporal/spatial operations.
- Observation/reanalysis product; model, experiment, member and grid identifiers.
- Historical period, calibration period, anomaly baseline, future window and any held-out validation period.
- Adjustment method, physical bounds, member handling and requested outputs.
- For extremes: the actual sample definition, threshold, declustering, exposure and stationarity assumption.

Treat unknown units, ambiguous variables, mismatched members, incomplete seasons and absent data as actionable input problems.
Do not silently pick a different variable/member, drop failed models, or claim that an unknown format is supported.
Routine choices such as filenames and plot size can use sensible documented defaults.

For detailed choices, read [workflow.md](references/workflow.md) and [config-reference.md](references/config-reference.md).
For units and variable-specific recipes, read [variables.md](references/variables.md).

## 3. Locate and prepare the runtime

Resolve this skill's directory from the location of this SKILL.md. Use the project's active Python interpreter or a suitable environment found on the machine.
Do not assume a server name, home directory, conda environment, font or map-data cache.

Run the entry point relative to the resolved directory:

~~~bash
python <skill-dir>/scripts/run.py doctor
python <skill-dir>/scripts/run.py init <run-dir> --template grid
~~~

Choose one starter:

| Template | Intended data and question |
| --- | --- |
| grid | Gridded monthly temperature and seasonal/regional change |
| series | Already computed annual indices in CSV or NPZ |
| precip | Monthly precipitation-rate QDM, with explicit rate units |
| regional-extremes | Maxima of a regional daily series |
| daily | Threshold events in a regional daily series using POT/GPD |

Edit the generated config.yaml before running. Relative paths resolve from the configuration file's directory, regardless of the shell's working directory.
Paths support ~ and environment variables; undefined environment variables are errors.
Use the same interpreter for dependency checks, analysis and tests.

The core uses NumPy, SciPy, pandas, xarray, Matplotlib and PyYAML. NetCDF-4 reading needs netCDF4 or h5netcdf; CF calendars need cftime; GRIB needs cfgrib/ecCodes.
Cartopy is optional. Offline plain-axis maps are valid; do not trigger large downloads solely to draw a coastline.

## 4. Preflight, then execute

~~~bash
python <skill-dir>/scripts/run.py validate <run-dir>/config.yaml
python <skill-dir>/scripts/run.py validate <run-dir>/config.yaml --check-inputs
python <skill-dir>/scripts/run.py validate <run-dir>/config.yaml --check-data --report-json <run-dir>/input-quality.json
python <skill-dir>/scripts/run.py run <run-dir>/config.yaml --dry-run
python <skill-dir>/scripts/run.py run <run-dir>/config.yaml
~~~

- Configuration validation checks settings, types and incompatible combinations.
- Input preflight checks every selected file's header and time coordinates, including cross-file consistency and requested seasonal coverage. It does not read every data value.
- Use --check-data to scan the selected time/level/region in bounded chunks. Review missingness, infinities and explicit value limits before accepting the input.
- Use quality.coverage: error when complete requested coverage is required; the default reports missing coverage as warnings. Define any value limits in transformed, pre-aggregation units.
- A dry run prints the input inventory and execution plan without creating analysis outputs.
- An explicit --report-json writes the inspection report even when input checks fail. Production runs record quality_report.json and link its scope and findings from the manifest/report.
- Runtime loading checks actual coordinates, calendars, time coverage and data values.
- Unknown configuration keys are errors by default. runtime.strict_config: false retains warnings for migration.
- Runtime data failures stop the run by default. If deliberately using runtime.strict: false, report every skipped dataset and the resulting ensemble size.
- A changed scientific configuration requires a new project.outdir or explicit --overwrite. Never delete the user's run directory as a routine troubleshooting step.

Stages are prepare, correct, ensemble, extremes, daily, figures, report.
A request such as --stages figures automatically includes the upstream computations needed by the requested figures.
This is dependency planning within one process, not a promise to reload all previous products.
The loader cache is keyed by source-file metadata, preprocessing settings, grid and implementation version. It is not a cryptographic archive checksum.
Read [quality-control.md](references/quality-control.md) for scan levels, thresholds, report interpretation and memory limits.
Preserve raw inputs. Do not automatically fill missing dates, interpolate missing values or clip outliers to make checks pass.

## 5. Apply methods with their actual meaning

Read [methods.md](references/methods.md) before selecting adjustment or extreme-value methods.

### Temporal and spatial processing

Use aggregate.frequency to state auto, monthly, daily or annual.
Monthly means use calendar-aware month-length weights by default.
Cross-year seasons such as DJF use the ending year and include the previous December.
Drop incomplete seasons by default; retain missing data as missing, including for counts and sums.

The bundled interpolation supports one-dimensional latitude/longitude grids.
Cosine-latitude regional weighting requires regularly spaced coordinates; use a validated cell-area workflow for other grids.
Curvilinear/rotated grids, GeoTIFF collections, conservative remapping and downscaling require external preprocessing.
Record that preprocessing and do not call interpolation bias correction or downscaling.

For precipitation, multiplying kg m-2 s-1 by 86400 yields mm/day.
Summing monthly mean rates does not yield a seasonal total; integrate over the actual represented durations first.

### Adjustment and indices

| Method | Interpretation |
| --- | --- |
| none | Preserve the supplied values; appropriate for controls or already adjusted input. |
| mean / delta_change | Shift additive means or use an explicit multiplicative mean ratio. |
| mean_std | Adjust mean and spread; scaling generally also changes modeled trends. |
| quantile_mapping | Historical distribution mapping; may distort projected quantile changes. |
| quantile_delta_mapping | Preserve the selected additive/multiplicative quantile changes under the implemented finite-sample assumptions. |

Fit historical calibration parameters only on the declared baseline.
QDM additionally estimates a CDF from the supplied application window: record the future window, and do not claim it is a sliding-window implementation.
Dry-day frequency adjustment is not implemented. A multiplicative ratio with a zero modeled denominator and a nonzero observed quantile requires a different, justified method.
Physical clipping changes moments and tails; examine the result.

Use region_mean, region_mean_anomaly, region_diff, projection_index or precomputed series as supported.
Apply future standardization with fixed historical parameters. Detrending removes the projected trend and belongs to explicit variability diagnostics.
For series inputs, preprocess anomalies or standardized values upstream and document that definition.

Keep historical and future members matched. Averaging members before correction/extremes changes variability and can suppress extremes.
Treat observed-trend rescaling as a separately labeled sensitivity experiment.

### Ensembles and extremes

Count available models by year; mask summaries below ensemble.min_models.
Report model percentiles as model spread. Do not label them confidence intervals or calibrated probabilities without an additional justified analysis.
Skill weighting here uses historical RMSE and does not establish independent skill or account for model relatedness.

Annual/seasonal means are not block maxima.
A GEV fit to an annual index is an explicit descriptive distribution assumption; declare extremes.sample_kind: annual_index.
Regional daily maxima require the daily route: the mean of grid-cell maxima is not the maximum of the regional mean.

POT declustering uses actual dates, and exposure uses valid days relative to the selected season's calendar.
Block-maxima fits default to complete seasonal blocks.
POT waiting times and GEV block-exceedance return periods have different definitions.
Inspect stationarity, threshold stability, tail fit and sample size before interpreting extrapolation.

SciPy GEV c equals -xi in the common extreme-value convention.
ks_p_naive is not a calibrated goodness-of-fit p-value after fitting parameters to those same data.
Do not reject or retain models by a universal p > 0.05 or |c| < 1 rule.
Formal refitted goodness-of-fit tests, nonstationary extremes and ensemble uncertainty bootstrap are not supplied by this pipeline.

## 6. Inspect and deliver

Read the current run's quality_report.json, report.md, manifest.json, tables and generated figures.
Check:

1. Units, calendars, baseline overlap, member matching and the list of completed/skipped inputs.
2. Calibration diagnostics separately from held-out evaluation. Free-running climate simulations need not reproduce observed year-to-year weather.
3. Historical/future continuity without shifting a projection merely to make a plot join neatly.
4. Per-year model counts, missing-year alignment, actual uncertainty labels and sign-agreement masks.
5. Return-level monotonicity, finite support, exposure and extrapolation. Treat a successful fit as a numerical result, not scientific acceptance.
6. Plot labels, scales, units, coverage, legends and readable panel layout.

Use [plotting.md](references/plotting.md) for figure choices and [troubleshooting.md](references/troubleshooting.md) for failures.
Report the variable, region, season, index, methods, periods, scenarios, contributing model/member counts and artifact locations.
Use numbers from this run's recorded outputs; distinguish completed analysis, synthetic verification and proposed work.

## 7. Maintain and extend

Keep changes local to the needed modules. Preserve public interfaces when practical and document unsupported behavior.
Add targeted regression cases for scientific invariants affected by a change.
Run the bundled checks with the same interpreter:

~~~bash
python <skill-dir>/tests/test_methods.py
python <skill-dir>/tests/test_io_contracts.py
python <skill-dir>/tests/test_raw_io.py
python <skill-dir>/tests/test_quality.py
python <skill-dir>/tests/test_quality_workflow.py
python <skill-dir>/tests/test_statistical_contracts.py
python <skill-dir>/tests/test_workflow.py
python <skill-dir>/tests/test_daily.py
python <skill-dir>/scripts/run.py selftest
~~~

Synthetic tests establish implementation behavior; they do not validate a real CMIP projection or replace validation on the user's archive.
For method provenance and official documentation, read [sources.md](references/sources.md).
