# Troubleshooting

Errors are grouped by stage. Each entry: message → cause → fix. When a fix
requires a code change, re-run `run.py selftest` and `tests/test_methods.py`.

---

## Configuration / startup

**`invalid configuration: ...`**
The validator lists every problem at once. Fix them all, do not iterate one by one.

**`missing required configuration key: data.models.root`**
`root` is required for `kind: grid`. For `kind: series` set `series_dir` instead.

**`unknown configuration key (ignored): plots.maps.period`**
A typo (or a key from an older version). The key is ignored, so the setting has no
effect. Compare with `references/config-reference.md`. Warnings are not fatal.

**`data.models.scenarios must not contain 'historical'`**
`historical` is the baseline experiment (`historical_experiment`), not a scenario.

---

## Data loading

**`no models discovered; check data.models.pattern / models / series_dir`**
The discovery glob matched nothing. Check:
1. `data.models.root` exists and is the directory that *contains* the experiments;
2. for `models: all`, discovery uses the `historical` experiment — if only future
   files exist, list `models:` explicitly;
3. the archive layout: some archives keep files flat
   (`<model>/<var>_Amon_<model>_<exp>_<member>_*.nc`), others put them in a
   per-variable sub-directory (`<model>/<var>/...`). Both built-ins are tried;
   anything else needs `data.models.patterns`.

**`no input files matched the configured pattern`**
The pattern resolved but the member/experiment string did not match. Typical cause:
the archive uses `r101i1p1f1` while the config says `r1i1p1f1`. Set
`members: ["*"]` — the loader then keeps the member with the most files.

**`variable 'tas' not in <file>; available: [...]`**
Wrong `variable`. Use the exact name printed in the message. For GRIB, ERA5 uses
`t2m`, `tp`, `z`, `msl`, `wsped`; CMIP6 uses `tas`, `pr`, `zg`, `psl`, `sfcWind`.

**`skipping <file> (no time dimension)`**
A static file (orography, area fractions) matched the glob. Narrow the pattern or
set an explicit `variable`.

**`cannot open GRIB file ...`**
`cfgrib`/eccodes missing, or a GRIB edition the local eccodes cannot read. Convert
to NetCDF first (`cdo -f nc copy in.grib out.nc`) and point the config at the
NetCDF.

**Values are ~×1000 too large / too small**
Unit conversion. ERA5 precipitation is m (×1000 → mm); CMIP6 `pr` is kg m⁻² s⁻¹
(×86400 → mm/day); `zg` in m² s⁻² needs ÷9.80665; temperature in K needs −273.15.
Print `mean` of the loaded field and compare with climatology before running.

**Values are plausible but the region is empty / all NaN**
Longitude convention or `bbox`. ERA5 is 0–360, some models are −180–180. The
loader normalises both, but a `bbox` given in the wrong convention can select
nothing. Remove `bbox` and rely on `index.region` for the spatial selection.

**`Reindexing only valid with uniquely valued Index objects`**
Duplicate longitudes (a 0–360 *inclusive* grid has both 0 and 360). The loader
de-duplicates; if this appears, the file was opened outside the loader or the
skill is outdated.

**Polar/high-latitude grid points are NaN**
Some models have missing corner cells or a different pole treatment. The index is
NaN-safe (it weights available points), but map panels will show holes. Consider a
`bbox` excluding the pole.

---

## Correction

**`fewer than 5 common years between reference and model for bias correction`**
`data.reference.years` and the model historical run do not overlap. Check the
model's historical time axis (CMIP6 ends 2014; some models start in 1850).

**`baseline climatology is all NaN`**
`periods.reference_climate` does not overlap the field being converted. This
happens when a future field is asked for an anomaly without the historical
climatology — the pipeline passes it automatically, so seeing this means a custom
script used `build_index` directly without `climatology=`.

**Corrected bias is 0 but raw bias is also 0**
The diagnostic index is already an anomaly (mean ≈ 0). The pipeline scores on the
*absolute* index, so this should not happen; if it does, check that
`index.type` is not `series` with pre-centred values.

**Corrected values are negative for precipitation**
Set `correction.lower_bound: 0`. Multiply-corrected precipitation can also be
zero-inflated; use `quantile_delta_mapping`.

**KS does not improve**
Possible causes: (a) the correction is fitted on a different quantity than the one
diagnosed (e.g. absolute vs anomaly — the pipeline handles this, a custom script
may not); (b) the model has a strong trend bias that a stationary correction
cannot fix; (c) too few years. Inspect `figures/fig_diagnostics_kde.png`.

**Trend after correction points the wrong way**
The correction is stationary: it cannot change the sign of a model trend. If the
projected trend is unphysical, the model or the index is wrong, not the
correction. Do not "fix" it with `trend_correction` unless the study design is
explicitly trend-scaling.

---

## Ensemble

**`only N model(s), skipping ensemble`**
Fewer than `ensemble.min_models` (default 3) models produced a series. Either add
models or lower `min_models` (and say so in the report).

**Bands are very wide / narrow**
Model spread is real; with <5 models percentile bands are not meaningful. Report
the model count and prefer min–max or individual lines.

**`RuntimeError: cannot reindex or align along dimension`**
Mismatched year axes between models (e.g. different future start years). The
alignment pads with NaN; if it errors, check that all series have integer years.

---

## Extremes / GEV

**`no valid GEV fits`**
Every fit failed or had < `min_years` values. Check the future period length
(2026–2050 = 25 years ≥ 20, OK; 2026–2040 = 15, too short), or whether a series is
constant.

**`shape` is huge (|c| > 2) or NaN**
Unstable fit from too few years or a bounded index. Compare across models; drop
models with extreme shapes from the *shape* comparison only, not from the ensemble.

**Return level curve stops before T = 1000**
Correct behaviour for c > 0 (upper-bounded distribution). Do not "extend" it.

**Kruskal–Wallis p is NaN**
Fewer than two scenarios with fits. Run at least two scenarios or disable
`extremes.shape_comparison`.

**`return_period` is `inf`**
The value is at or above the fitted upper bound. Report "≥ the fitted upper bound"
rather than a number.

---

## Daily extremes

**`daily: no valid POT fit` / `only N peaks above ...`**
Too few declustered peaks. Lower `threshold_quantile` (0.95 → 0.90), relax
`min_peaks`, extend the season or the period, or switch `threshold_mode` to
`quantile`. With `threshold_mode: observed` a warming variable can leave almost
no days above a fixed historical threshold — that is a design signal, not a bug.

**`daily field would hold ...M values`**
The memory guard rail. Set `data.reference.bbox` and `data.models.bbox` to the
study region, or shorten `periods.future`.

**daily stage skipped**
`daily.enabled` is false, or the input is `kind: series` (daily extremes need
gridded sub-daily values).

**daily correction changes the threshold meaning**
For a strongly warming variable, a fixed `observed` threshold produces rapidly
increasing exceedance counts. Use `threshold_mode: absolute` for an impact metric
("days above 35 °C") or `quantile` to keep the rate constant; say which in the
report.

**`daily.correction` unknown**
Must be one of `none, mean, mean_std, quantile_mapping, quantile_delta_mapping`.

## Figures

**No figures produced, no error**
`plots.figures` is empty or contains names not in the catalogue. Validation would
normally catch this; check the log for `figures` warnings.

**Figure is blank / axes empty**
The stage data was empty (e.g. `maps` requested but the model source is series).
Check the warning in the log.

**Chinese text renders as boxes**
Set `plots.font: WenQuanYi Micro Hei` (installed). Do not install new fonts into
the shared environment without asking.

**cartopy hangs at "Downloading"**
The machine is offline and the Natural Earth cache is missing. The cache at
`~/.local/share/cartopy` is complete for 110 m/50 m — do not delete it. If it is
genuinely missing, maps still render without coastlines.

**`OverflowError` / `ValueError` in a GEV plot**
Non-finite return levels from an unstable fit. Filter with
`extremes.min_years` or switch to `gumbel`.

---

## Runtime

**A run takes hours**
Cost is dominated by opening many small files (yearly chunks). Options: use the
cache (`runtime.cache: true`, default) so a re-run is fast; reduce `grid`
resolution; restrict `data.models.bbox`; reduce `periods.future` for the pilot.

**Memory blow-up**
The pipeline streams corrected fields to disk, so memory is dominated by one
model/scenario field at a time. If a single field is too large (0.25° global,
daily), set `data.*.bbox` to the study domain.

**A model fails but the run continues**
By default failures are logged and the model is skipped. With
`runtime.strict: true` the run aborts on the first failure — use it when
completeness matters.

**Two runs in the same `outdir`**
They overwrite each other's files. Use a different `project.name`/`outdir`, or
wait.

---

## Reporting mistakes to avoid

* Quoting a trend without units or per-decade conversion.
* Calling the 10–90 % band a confidence interval.
* Presenting a return period without saying which distribution and fit window.
* Comparing scenarios produced with different `params_period` or correction methods.
* Reporting `n` models without saying that some scenarios had fewer (check
  `ensemble.n_models` in the NPZ).
