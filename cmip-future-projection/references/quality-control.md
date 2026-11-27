# Raw CMIP input quality

Contents: 1. Choose a scan; 2. Interpret checks; 3. Set limits;
4. Read the report; 5. Keep processing stable; 6. Act on findings.

Use this reference when accepting a new CMIP archive, diagnosing inconsistent
files, or preparing a production run. Preserve raw inputs and record any external
repair or exclusion. These checks establish input properties for the configured
analysis; they do not certify a model's scientific skill or full CF compliance.

## 1. Choose a scan

Run from any working directory with an explicit skill/configuration path:

~~~bash
python /path/to/skill/scripts/run.py validate myrun/config.yaml --check-inputs
python /path/to/skill/scripts/run.py validate myrun/config.yaml --check-data --report-json myrun/input-quality.json
~~~

| Mode | Reads and checks | Output |
|---|---|---|
| `validate CONFIG` | Configuration types, options and compatible combinations. | Console/JSON summary; no analysis. |
| `validate CONFIG --check-inputs` | Discovery and every selected file's header, coordinates, time axis and compatible identifiers. | Findings and requested coverage; raw values are not exhaustively scanned. |
| `validate CONFIG --check-data` | All input checks plus all selected raw values in bounded chunks. | Missingness, finite range, infinities and configured value-limit findings. |
| `run CONFIG` | Input checks, then configured analysis. | `quality_report.json`, resolved configuration, manifest and requested products. |

`--check-data` implies input discovery. "Selected" means the configured variable,
level, source bbox, season and requested period; it does not mean every variable
or all years in the archive. A configured index region may be smaller than the
source bbox; choose `data.*.bbox` deliberately.

Both `validate` and `run` accept `--report-json PATH`, including for failed input
checks. This explicit option writes a report even during a dry run. The report
path resolves from the shell's working directory; data/output paths inside YAML
resolve from the YAML directory. Do not point a report at a source or configuration
file. Without this option, validation does not create analysis outputs.
An existing destination is replaced only if it is a previously generated
`cmipkit.input_quality` JSON report; other existing files are protected.

Production runs use the configured `quality.scan` level. `run --check-data`
forces a full scan. Explicit missingness/value limits also force a full scan in
production, so a limit cannot silently pass without reading values. The explicit
`validate --check-inputs` command still performs only header/time checks.

## 2. Interpret checks

### Files, variable and grid

Inspect all discovered file paths, not only the first file in a series. Check
that files are readable and nonempty, the requested variable and level exist,
and the grid is supported. Incompatible coordinates or conflicting source,
experiment, member or grid identifiers must be resolved before concatenation.
Metadata absent from a file cannot establish its simulation identity: retain
the filename/path evidence and inspect absent fields in the reported identity.

Each model may use its own calendar. Files being combined into one time series
must have compatible calendars and units. Keep historical/scenario members
matched; do not independently choose the member with the most files in each era.
Restrict input globs to one version and grid label. File modification times and
sizes identify inputs for this run; they are not cryptographic integrity checks.

### Time axes and coverage

Read actual time coordinates. Inspect nonmonotonic ordering, duplicate timestamps,
overlapping files and missing requested samples. Report overlaps separately from
conflicting data: header inspection alone cannot prove two copies have equal
values. Value comparisons and the runtime loader reject conflicting overlaps;
identical duplicate timestamps can be collapsed with an explicit warning.

Evaluate required samples in the source calendar. A `360_day` year has 360 days;
do not insert Gregorian leap days or coerce its timestamps to Gregorian dates.
For DJF, the season-year includes December of the preceding calendar year. Check
missing months/days within the selected season, not the intentionally excluded
parts of the year. Declare `aggregate.frequency` if sparse sampling makes automatic
frequency inference ambiguous.

`quality.coverage: warn` reports gaps without treating every shorter archive as
fatal. Use `error` when all requested samples are necessary. This setting is
separate from `aggregate.require_complete`: aggregation can discard incomplete
seasons even when the quality inspection reported only a warning. Inspect the
years actually retained by the analysis.

### Units and decoded values

`units_in` is an assertion about input metadata. The checker recognizes a limited
set of equivalent spellings; it does not perform dimensional analysis or infer
an arbitrary conversion. Keep `units_transform` explicit. A Kelvin/Celsius
mismatch is not resolved merely by changing a label.

CF packing and fill-value decoding are handled by the file reader. Apply the
configured affine conversion once to the decoded field, then evaluate explicit
quality bounds before temporal aggregation. This differs from `correction` bounds
applied to adjusted products.

The full scan reports nonfinite values, finite minima/maxima, missingness and
all-empty selections/time slices. Infinities are errors. Ordinary NaNs, including
static land/ocean masks, are not automatically a failed dataset; a wholly empty
selected dataset cannot support an analysis. No universal physical temperature,
precipitation or anomaly range is guessed.

## 3. Set limits for the actual analysis

When the daily branch selects different months or a different reference period,
inspect its separately labelled audit as well as the annual-analysis audit. Only
identical daily scopes are reused. Sampling-frequency metadata and time-coordinate
inference must agree with the frequency used by the loader; declare the frequency
explicitly when sparse metadata-free data are ambiguous.

Check the discovered model count against `ensemble.min_models`. The preflight
warns when the requested ensemble minimum cannot be reached with the available
models. Passing raw-input checks does not guarantee that a later statistical
stage has enough samples or models for its own requirements.

Start with an inspection and derive limits from the variable and selected domain:

~~~yaml
quality:
  scan: full
  coverage: warn
  require_units: false
  max_missing_fraction: null
  valid_min: null
  valid_max: null
  chunk_size: 32
  max_chunk_values: 2000000

runtime:
  max_input_values: 30000000
~~~

Use `coverage: error` or `require_units: true` when those requirements apply.
Set `valid_min`/`valid_max` only after verifying the configured conversion; the
limits use transformed input units. Setting only one bound is supported. Bounds
detect problems and stop the run; they do not clip values or alter the source.

The missing fraction uses selected time/space values as its denominator and
includes static masks and overlapping file records before runtime deduplication.
A land-only variable on a bbox containing ocean may have
a large legitimate missing fraction. Inspect masks before choosing a threshold.
Missing timestamps/files are diagnosed through coverage, not counted as observed
NaN values. Infinities contribute to the nonfinite fraction and are also reported
as errors. A missingness limit therefore does not replace a coverage check.

## 4. Read the report before scaling up

Inspect the report's top-level `ok`, `errors`, `warnings` and discovered datasets.
Its `quality` section records:

- `level`: configuration-only, headers or full value inspection.
- `data_values_checked`: whether the requested value scan was performed.
- `summary`: inspected dataset/file and finding counts, plus scanned values.
- `datasets`: per-dataset headers, identity, coverage, findings and available
  value statistics.

Each finding has a severity, machine-readable code and message; relevant file
paths, counts and limited examples help locate the problem. Keep the full JSON
when console output is shortened. Distinguish requested scan level from how much
could be inspected when files were unreadable; read the errors and per-file data.

Production records `quality_report.json`. The manifest links to it and records
its scope and summary; `report.md` includes a brief quality section. Header-only
success means raw values were not exhaustively scanned. Full-scan success means
the implemented checks and explicit limits passed within the configured selection.
Neither statement establishes validity outside that selection or guarantees a
scientific conclusion.

## 5. Keep processing stable

Crop the requested level/region/time before loading numeric arrays. Annual
aggregation reads selected season-years; daily processing has a selected-data
size guard. `runtime.max_input_values` is checked before those loads, on combined
annual output size and when accepting a prepared-field cache. Estimate
memory for the dtype, temporary arrays, interpolation, fitting and final products
before increasing it. The guard is not a total-process RAM guarantee.
With model parallelism, each worker has its own loaded arrays and temporary
allocations. Budget for the coordinating process plus all active workers; the
per-load guard is not a shared memory reservation. See [parallel.md](parallel.md).

Value scanning has separate `chunk_size` and `max_chunk_values` controls. Smaller
chunks reduce individual array loads at the cost of more I/O. Coordinates and
reporting structures also consume memory. These controls do not automatically
turn the pipeline into a Dask/distributed or arbitrary-size workflow.

Read required values before closing each dataset. Write caches and run metadata
through temporary files before replacing the complete target. Rebuild unreadable
caches from source with a diagnostic rather than treating them as valid inputs.
Record implementation-version changes in cache identity; preserve source files.

After failure or interruption, read the manifest status and log. Completed stages
and existing plots are evidence only for those products. Stage selection still
resolves upstream work in one run; it is not a general resume-from-artifacts
system.

## 6. Act on findings

Fix the first identified cause: correct a path/variable, replace a damaged file,
select consistent versions/members, or supply the missing boundary period. If the
scientific design deliberately accepts incomplete coverage, change the relevant
policy explicitly and report the retained sample. Preserve missing data and
outliers until their cause is understood. Do not automatically interpolate gaps,
drop failed models, overwrite source NetCDFs or weaken all checks to obtain a pass.

Use [troubleshooting.md](troubleshooting.md) for targeted fixes and
[sources.md](sources.md) for CF and xarray references.
