# Operating procedure

Contents: 1. Choose the route; 2. Inspect inputs; 3. Fix the analysis contract;
4. Configure and preflight; 5. Pilot and run; 6. Assess and deliver.

## 1. Choose the smallest useful route

Use the user's existing data, decisions and analysis code when they are available.
Do not turn a request to inspect a configuration or redraw a figure into a new
CMIP study. Infer routine choices from verified metadata; ask only about missing
choices that change the scientific question.

| Request | Route |
|---|---|
| Inspect or repair a configuration | Validate, inspect metadata, explain necessary edits; stop before production computation. |
| Project an existing annual index | Use the series template; do not invent a gridded reconstruction. |
| Project a seasonal mean or spatial pattern | Use the grid template and an explicit index definition. |
| Estimate daily regional extremes | Use daily block maxima or POT with daily gridded input. |
| Use WorldClim GeoTIFF, station daily CSV, curvilinear ocean grids, or custom event durations | Identify the required adapter; the bundled CLI does not ingest these directly. |
| Draw from existing results | Reuse verified outputs with a small plotting script if suitable; do not claim the CLI can resume all computed stages from saved NPZ files. |

The CLI operates on local files. Dataset discovery, authenticated downloads,
conservative remapping, nonstationary extreme-value models and downscaling are
separate tasks. Do not claim they happened merely because a configuration exists.

## 2. Inspect inputs

Use the active Python environment and the actual location of this skill. No
server-specific interpreter or data directory is required:

~~~bash
python /absolute/path/to/cmip-future-projection/scripts/run.py doctor
~~~

Run the environment check once when setup is uncertain. Use representative files
to understand each layout, then let `validate --check-inputs` inspect all selected
files. Use `xarray.open_dataset` with the appropriate backend for a targeted manual
inspection. Confirm:

| Evidence | Decision |
|---|---|
| Variable name, units, `cell_methods`, time bounds | Affine conversion, rate versus accumulated quantity, sampling frequency. |
| Time coordinates and calendar | Seasonal labels, missing/duplicate dates, valid historical and future coverage. |
| Latitude/longitude dimensions and extents | Rectilinear grid compatibility, longitude convention, subset and interpolation. |
| Model, experiment, variant, grid label, version | Files belonging to one consistent simulation; historical/scenario pairing. |
| Missing values and regional coverage | Validity threshold and a stable spatial support. |
| Minimum/maximum and a few concrete values | Obvious unit or scale problems before correction. |

CF metadata defines the distinction between point values and cell means/sums;
consult the [CF conventions](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html)
when the source's time semantics are ambiguous. Do not infer a precipitation
accumulation interval from the variable name alone.

The grid reader expects one-dimensional latitude and longitude axes after level
selection. Bilinear interpolation does not create resolved fine-scale climate
information and does not conserve integral quantities. Pre-remap irregular grids
with an appropriate method before using the skill.

## 3. Fix the analysis contract

Record the variable, output units, input frequency, season, region, index,
reference, periods, experiments, members, correction and requested products.
Keep scientific choices in the YAML and its comments, not only in the chat.

Separate these periods:

- `periods.historical`: model historical data to load.
- `correction.params_period`: calibration data for the adjustment.
- `periods.reference_climate`: historical baseline for anomalies and map changes.
- `periods.future`: scenario data and the future CDF window for QDM.
- `extremes.fit_window`: stationary fit windows for the annual extremes route.

Do not force historical simulations to follow the observed year-to-year weather.
They are not forecasts initialized to reproduce individual observed years.
Evaluate distributions and relevant processes. Calibration improvement alone is
not evidence of future predictive skill.

Decide whether correction is needed. Use `none` for an unadjusted benchmark or
already adjusted input. Use `mean` when an offset is the intended operation.
Use `mean_std` only if scaling variability and the trend amplitude is acceptable.
Use QDM when preservation of quantile changes is part of the design, with explicit
future windows and checks around zeros. Read `methods.md` before extreme analysis.

Keep annual GEV off for seasonal means, seasonal totals and arbitrary circulation
indices unless the user explicitly requests and justifies an empirical GEV model.
A fitted curve alone does not supply extreme-value asymptotics.

## 4. Configure and preflight

~~~bash
python /absolute/path/to/skill/scripts/run.py init myrun --template grid
python /absolute/path/to/skill/scripts/run.py validate myrun/config.yaml
python /absolute/path/to/skill/scripts/run.py validate myrun/config.yaml --check-inputs
python /absolute/path/to/skill/scripts/run.py validate myrun/config.yaml --check-data --report-json myrun/input-quality.json
python /absolute/path/to/skill/scripts/run.py run myrun/config.yaml --dry-run
~~~

Edit paths, units, periods and region before the input check. Relative data and
output paths resolve against the configuration file's directory, not the shell's
current directory. `~` and environment variables are supported. Brace placeholders
such as `{model}`, `{experiment}` and `{member}` are resolved by file discovery.

| Template option | File | Starting point |
|---|---|---|
| `grid` | `templates/config_grid_projection.yaml` | Monthly regional temperature projection and maps. |
| `series` | `templates/config_series.yaml` | Existing annual index CSV/NPZ, no grids. |
| `precip` | `templates/config_precip_qdm.yaml` | Monthly precipitation rate, seasonal average, QDM. |
| `regional-extremes` | `templates/config_regional_extremes.yaml` | Annual/seasonal maxima of a daily regional temperature series. |
| `daily` | `templates/config_daily_extremes.yaml` | Daily regional POT/GPD with a fixed historical threshold. |

Input preflight checks every selected header and time axis. Review the reported
gaps, overlaps, units and source/member/grid identities. Missing requested coverage
is a warning by default; set `quality.coverage: error` when it must stop the run.
Use `--check-data` for the selected raw values, including later files and chunks.
Keep the JSON report: header checks and full value scans have distinct scopes.
Neither is a test of scientific model skill. If a path is unavailable, complete
the configuration and state exactly which input is missing.

Read [quality-control.md](quality-control.md) before setting missingness or value
limits. Bounds apply after the configured unit conversion and before aggregation.
Do not copy a temperature bound into a precipitation or anomaly analysis.

## 5. Pilot and run

Create a separate pilot configuration/output directory. Use one model, one
scenario, a small region and enough years for the selected method. A two-year
pilot can test I/O, but cannot validate a 30-year GEV or QDM analysis.

~~~bash
python /absolute/path/to/skill/scripts/run.py run pilot/config.yaml
python /absolute/path/to/skill/scripts/run.py run myrun/config.yaml
~~~

`--stages` selects requested outputs and includes necessary upstream computations
in the same process. For example, `--stages figures,report` still needs calculated
indices and ensemble statistics. This is dependency resolution, not a promise of
resuming every stage from previously written products. The cache can reuse eligible
prepared fields. Preserve separate run directories for changed scientific designs.

Inspect logs before scaling up. Strict processing is on by default. In permissive mode, failed models can be skipped;
a successful exit is not proof that every requested model/scenario survived.
Prefer `runtime.strict: true` for final analyses when all requested members are
required. Freeze the discovered model list for a reproducible production run.

Production runs use `quality.scan: headers` by default. Set `quality.scan: full`
or pass `run --check-data` for a value scan; explicit missingness/value thresholds
also require full scanning. The scan reads source data, so allow for additional
I/O even when prepared-field caches already exist. Preserve the generated
`quality_report.json` beside the resolved configuration and manifest.

## 6. Assess and deliver

Check unit ranges, calendars, period completeness, missing-data fractions,
selected members, per-scenario model counts and map coverage. Compare a few
hand-computable values to the output. Review every requested figure for labels,
baselines, sample definitions and layout.

Read the quality report's scan level and findings before describing input quality.
Report unresolved warnings and any deliberately relaxed coverage/valid-fraction
settings. A successful headers-only run does not establish exhaustive value checks.

For correction, compare raw/adjusted diagnostics on calibration data and an
independent period when available. Do not require an arbitrary tenfold bias
reduction or remove models until p-values pass. For GEV/GPD, assess the support,
scale, sample count, QQ/return-level behaviour, threshold/window sensitivity and
stationarity; a naive KS p-value is not an acceptance gate.

Report what was actually computed: scenario-conditioned changes, model counts,
spread, correction choices, fit assumptions and exact output locations. Distinguish
completed computation from configuration-only work, missing-data exclusions and
unimplemented extensions. Keep numerical claims traceable to tables/arrays.
