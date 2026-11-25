# Troubleshooting

Contents: 1. Preflight; 2. Input integrity; 3. Scientific semantics;
4. Correction; 5. Extremes; 6. Runtime and figures.

Diagnose the first concrete failure, preserve the run log and change the smallest
relevant setting. Do not weaken all validation, delete output directories or
substitute other data just to obtain a successful exit.

## 1. Configuration and preflight

| Symptom | Action |
|---|---|
| Unknown configuration key | Fix spelling or remove a genuinely obsolete key. Strict config is on by default. |
| Relative path works only from one shell directory | Paths resolve from the YAML directory; inspect the resolved path and use a correct relative/absolute path. |
| No input files matched | Inspect the precise model/experiment/member pattern, table name, extension, root and environment variables. |
| No models discovered | Confirm the layout contains the configured model placeholder; try one verified explicit model for diagnosis. |
| Environment check fails | Use the intended environment and install only required dependencies through the environment's normal process. |
| Existing run refuses changed configuration | Give the new design a new output directory, or use explicit overwrite only when replacing that run is intended. |
| Header preflight succeeds but numerical processing fails | Headers/time coordinates were checked, but raw values may not have been scanned. Run `validate --check-data` and inspect the implicated file and `logs/run.log`. |

Use `validate CONFIG --check-inputs`, followed by `validate CONFIG --check-data
--report-json input-quality.json` when a value scan is needed, before an expensive
run. Missing-path errors do not authorize invented data; provide the ready config
and the exact missing inputs if data cannot be obtained within the request.

## 2. Input integrity

**A variable has unexpected extra dimensions.** Select the intended level and
member upstream or through documented source options. Do not average ensemble,
forecast lead-time or vertical dimensions merely to force a three-dimensional
array. Check variable names in multi-variable files.

**GRIB messages do not form one hypercube.** Inspect the product/message groups
with cfgrib/ecCodes. Select one compatible field or convert it to a documented
NetCDF intermediate. Merely changing a file extension is not conversion.

**Calendar cannot be decoded.** Install cftime for supported non-Gregorian dates;
retain the native calendar. Do not rewrite a 360-day calendar as Gregorian.
Inspect malformed time units against the source specification.

**Seasons disappear.** Completeness checking may remove boundary or missing-date
seasons. Supply the previous December for DJF and other required boundary data.
Inspect missingness before considering a lower valid-fraction threshold. Record
any deliberate relaxation and its effect on results.

**Duplicate times or overlapping files.** Constrain patterns to one version,
grid and variant, and distinguish exact replicas from conflicting overlapping
values. Header inspection reports the overlap; value comparison is required to
establish whether copies are identical. Do not silently average them into a new
simulation.

**A later file fails while the first file looks correct.** All-file preflight
reports each incompatible or unreadable file. Compare units, calendar, variable,
grid and source/member/experiment identifiers across the reported pair. Replace
a damaged download or narrow an overbroad glob; do not skip it without documenting
the resulting coverage.

**Missing coverage is only a warning.** Set `quality.coverage: error` if complete
requested periods are required. Check the previous December for DJF and the native
calendar before changing dates or the coverage policy. See
[quality-control.md](quality-control.md) for inspection scope.

**Missingness exceeds a configured limit.** The denominator includes the selected
time/space domain, including static land/ocean masks. Examine whether the threshold
matches that domain. An absent file/time coordinate is a coverage issue rather
than a NaN value. Do not replace missing data with zeros to pass the limit.

**Infinity or a value-limit error.** Trace the affected variable and file, inspect
decoded values and the affine conversion, and check for corrupt input or double
unit conversion. Value limits use transformed, pre-aggregation units. Fix the
cause; raising the limit or clipping values requires a scientific reason.

**Empty regional subset.** Inspect latitude order, longitude convention and
actual overlap. A region crossing 0°/180° needs wrap-aware bounds. Daily input
still needs a bounded domain; removing the bbox entirely is not a general fix.

**Curvilinear latitude/longitude or GeoTIFF input.** The native reader does not
cover these formats/grids. Build an explicit tested adapter or precompute a
verified series. Preserve CRS, nodata, masks, time semantics and spatial weights.

**Different member across historical/scenario.** Select a member present in both
experiments and document any exceptional lineage; do not use the member with the
most files independently for each segment.

## 3. Quantity and units

**Rainfall is off by factors of 30, 1000 or 86400.** Check whether each value is a
rate, a monthly mean daily accumulation or an interval total. The affine converter
cannot infer interval durations. A sum of monthly mean rates is not a seasonal
total. Use the mean-rate template or integrate correctly upstream.

**A temperature conversion looks unreasonable.** Compare metadata and a handful
of actual raw values. Kelvin-to-Celsius is an offset, not a multiplicative ratio.
Do not apply it again to a file already stored in Celsius.

**`units_in` disagrees with the header.** Correct the declared source units or
select consistent source files, then verify `units_transform`. Recognized aliases
do not turn the checker into a general unit converter. Missing metadata warns by
default; `quality.require_units: true` makes it an error.

**A future anomaly is near zero despite warming.** Verify a fixed historical
baseline and fixed historical standardization parameters. Never subtract each
future window's own mean if the goal is to measure change.

**An index differs from a published one.** Check spatial/temporal operator order,
longitude coverage for zonal means, mask, pressure level, template normalization,
season-year labeling and units. A daily regional maximum differs from an average
of cellwise maxima.

## 4. Bias correction and ensembles

**Correction fails on constant/all-missing data.** Inspect the baseline sample.
Do not turn absent data into zeros or repair zero variance with random noise.
Use a justified simpler correction, omit unsupported cells/models with a reason,
or obtain adequate baseline data.

**QDM gives large ratios at dry locations.** Inspect zero/near-zero historical
quantiles, dry-day frequency and future-window definition. A wet-day model may
be required; clipping to zero does not solve unstable positive ratios.

**Mean_std changes the projected slope.** This follows from its variance-scaling
factor. Use additive mean adjustment or another justified method if preserving
the raw additive change is the requirement.

**Distributional diagnostics improve but a model is still physically poor.**
That is possible: adjustment can fit calibration statistics while leaving
circulation/feedback biases. Use independent relevant process evidence.

**A band is narrow or changes suddenly.** Inspect valid model counts and whether
the model subset changes over time/scenarios. A single-model spread is not a
confidence interval. Member averaging also suppresses variability.

## 5. GEV and daily extremes

**No GEV figure.** Annual GEV is opt-in and distinct from daily extremes. Confirm
that the requested sample is a block maximum or an explicitly justified empirical
index fit, then check valid sample count, variance and support.

**KS p-value is missing.** A fitted-sample naive KS result is retained only as
`ks_p_naive`. A calibrated test has not been performed unless explicitly reported.
Do not replace missing `ks_p` with a naive p-value or require p > 0.05 as a gate.

**Hundred-year levels are unstable.** Inspect fit-window stationarity, sample
length, outliers, support and parameter/threshold sensitivity. Short-sample
extrapolation cannot be repaired by increasing plot resolution or trying fits
until one passes a p-value test.

**POT yields too few peaks.** Confirm units, threshold, selected days, coverage
and declustering. Choose a threshold using stability diagnostics; do not lower it
solely to clear a sample-count check. Explain when the data do not support a fit.

**Daily gaps affect clustering.** Preserve timestamps. Missing days are unknown,
and successive selected-season samples can be months apart. Inspect the exposure
and completeness rule before interpreting event rates in years.

**A requested heatwave duration is absent.** The daily stage fits intensity peaks
or block maxima. Consecutive-day duration is an additional event definition and
algorithm, not a plot option.

## 6. Runtime and figures

**Memory pressure.** Subset the region, level, model list and required periods;
use a pilot to estimate resources. `runtime.max_input_values` guards selected
loads before reading, per season-year for annual aggregation and cumulatively
for daily fields. `quality.max_chunk_values` limits value-scan chunks separately.
These controls do not bound every downstream temporary array or total process RAM.
Increase a limit only after estimating available memory. Changing calendars or
dropping data silently is not a memory-management method.

**Stale prepared fields suspected.** Disable `runtime.cache` for a small diagnosis
or use a new run directory. Preserve the old products and logs rather than
deleting the only record of a previous result. Corrupt cache files are rebuilt
from the source with a warning; temporary writes are not accepted as complete
caches. Input metadata and implementation-version keys reduce accidental reuse,
but do not replace archive checksums.

**A run was interrupted.** Inspect the manifest's status and completed stages,
the quality report and the log. Retain completed products; do not describe the
run as completed or assume every stage can resume from saved artifacts. Atomic
metadata/cache writes protect complete files from partial replacements.

**Only some scenarios appear.** Read the model/scenario failures. Strict processing
is on by default; if permissive mode was requested, report every omission and
actual count. A zero exit code alone does not guarantee complete coverage.

**Fonts or coastlines unavailable.** Select an installed font; use available
cartographic data or a clearly labelled plain-coordinate map. Check the final
rendering for missing glyphs. Do not assume the original server's font/cache.

**Figures requested from an earlier run.** The CLI's stage selection resolves
upstream computations in-process; it is not a full artifact-resume engine. Use a
small plotting adapter on verified saved arrays if that best satisfies the task.
