# Methods and interpretation

Contents: 1. Time and units; 2. Spatial operators and indices; 3. Bias adjustment;
4. Ensembles and trends; 5. GEV; 6. Daily POT and block maxima; 7. Scope.

Use this document to choose methods and interpret outputs. For implementation
parameters, see `config-reference.md`. Sources are collected in `sources.md`.

## 1. Time and units

Apply `x * scale + offset` before aggregation. `units_in` records the source
units; changing that string does not convert the numbers. Do not apply a Kelvin
conversion twice. The output label must match the actual transformed quantity.

For complete monthly means, compute a seasonal mean with calendar-aware month
lengths:

~~~text
seasonal_mean = sum(days_in_month * monthly_mean) / sum(days_in_month)
~~~

`aggregate.weight_by_days: true` performs this weighting for monthly mean and
mean_abs. Daily samples are equally weighted. Set `frequency` explicitly when
metadata or sparse sampling makes automatic detection ambiguous. Keep
`require_complete: true` to reject incomplete seasons, and specify
`min_valid_fraction` for finite values at each grid cell. All-missing cells must
remain missing, including sums and threshold counts.

`sum` sums the input values; it is not a flux integrator. For monthly mean
precipitation flux p in kg m−2 s−1, the monthly water-equivalent total is
`p × 86400 × days_in_month` in mm. Monthly accumulated totals are simply summed.
Fluxes in m s−1 need `1000 × 86400` to become mm day−1. Confirm the source
accumulation interval separately for ERA5 products. The precipitation template
avoids this ambiguity by calculating a seasonal mean rate in mm day−1.

Cross-year seasons must use consistent year labels: winter 2000 includes
December 1999 and January–February 2000. Include the boundary files needed for
complete seasons. Retain each model's calendar for aggregation instead of
pretending that a 360-day calendar is Gregorian. Compare annual/seasonal values
on their labelled years; this does not require an artificial daily date match.

## 2. Spatial operators and indices

For a regular latitude–longitude grid, a regional mean uses finite grid cells
and latitude weights proportional to cos(latitude). Use actual cell areas or an
appropriate remapping step for irregular grids. Missing cells alter the sampled
area; assess the coverage instead of assuming a temporal completeness check
also guarantees constant spatial support.

The order of spatial averaging and nonlinear temporal reduction matters:

~~~text
max_t(area_mean[x(t)]) != area_mean[max_t(x(t))]
~~~

The annual gridded route aggregates each cell in time before constructing the
regional index. The daily route constructs a regional series first, then takes
its block maxima or POT peaks. Choose according to the event being studied.

| Index | Meaning and caution |
|---|---|
| `region_mean` | Area-weighted regional value after the requested temporal aggregation. |
| `region_mean_anomaly` | Regional anomaly relative to the specified historical baseline. |
| `region_diff` | Region B minus region A after subtracting the zonal mean. This is not an ordinary unadjusted north–south difference. |
| `projection_index` | Projection on a fixed historical zonal-anomaly pattern; dimensionless for compatible field/template units. |
| `series` | Use the supplied annual series; its upstream spatial and temporal definition remains the user's responsibility. |

For a projection index, with zonal anomaly A′, fixed pattern T and spatial
weights w:

~~~text
I(t) = sum[A′(t) T w] / sum[T² w]
~~~

The pattern and denominator must be valid on the selected region. A projection
coefficient of one does not imply a perfect pattern correlation: orthogonal
structure can remain. A regional longitude crop changes the meaning of the
zonal mean; retain a full longitude circle when a true zonal anomaly is intended.

Subtract each model's historical climatology for its future anomalies; do not
subtract the future mean. Standardization must likewise reuse historical mean
and standard deviation. Detrending deliberately removes part of the signal and
should be disabled for an ordinary future-change calculation.

## 3. Bias adjustment

Let μo, σo and μm, σm be historical reference/model mean and standard deviation
on the calibration data. Fit the adjustment separately for each grid cell or
annual series. These are statistical adjustments, not repairs of incorrect
circulation dynamics or spatial dependence.

### Mean and variance adjustment

~~~text
additive mean:        x* = x + (μo − μm)
multiplicative mean:  x* = x μo / μm
additive mean_std:    x* = μo + (x − μm) σo / σm
~~~

Additive mean adjustment preserves the raw additive changes. Multiplicative mean
adjustment preserves ratios where the denominator is meaningful. Mean_std
rescales variability and any linear trend by σo/σm; it does not generally
preserve climate-change amplitude.

The implementation's multiplicative mean_std combines mean scaling and a
coefficient-of-variation adjustment. In the nondegenerate positive-mean case,
its affine result is algebraically the same as the mean_std expression above.
Do not describe it as simultaneous preservation of relative change and variance.
Zeros, near-zero denominators and clipping need explicit diagnostics.

`delta_change` is an alias of mean adjustment in this package. It is not a
separate weather-generator or observed-series perturbation method.

### Empirical quantile mapping

With historical CDF Fm and reference quantile function Qo:

~~~text
QM: x* = Qo(Fm(x))
~~~

The finite empirical tables interpolate between estimated quantiles. This
implementation clips outside the trained support. Applied to a future climate,
QM can distort changes and clip novel tail values. It does not mathematically
force every possible future distribution to equal the historical distribution.
Prefer a change-preserving method when that property is required.

### Quantile delta mapping

For a future sample at probability p = Ffuture(x):

~~~text
additive QDM:        x* = x + Qo(p) − Qm(p)
multiplicative QDM:  x* = x Qo(p) / Qm(p)
~~~

Historical Qo and Qm remain fitted to the baseline; the future empirical CDF is
estimated from the future sample being adjusted. In the CLI this is the entire
configured future period for that scenario, not a rolling window. Keep the
window fixed in comparisons or run explicitly documented windows separately.

The design preserves quantile-wise additive or relative changes, subject to
sampling, interpolation, zero handling and bounds. It does not promise exact
preservation of a fitted temporal trend. The bundled implementation does not
provide a full precipitation occurrence/wet-day-frequency model. All-zero or
nearly dry baselines need an alternative documented treatment; a bound of zero
alone does not solve this.

### Calibration, validation and optional trend manipulation

Assess calibration and out-of-sample behaviour separately. For free-running
climate simulations, same-year RMSE/correlation mixes forced-response skill with
unmatched internal variability. Prefer physical process checks and distributional
metrics aligned with the scientific question. Do not select models to produce a
preferred projected signal or pass a desired significance threshold.

Leave `correction.trend_correction.enabled` false for normal projections. The
optional operation changes an index using an observed/model historical slope
comparison. It is an imposed sensitivity experiment; it is not a general
calibration of future trends and is not guaranteed to make the future slope
exactly equal to the observed historical slope.

## 4. Ensembles and trends

Construct one model value per year, then compute the configured mean/median and
percentile bands over finite model values. Report counts by scenario and year.
A changing set of models can create artificial changes in the ensemble; compare
a common model subset when making scenario differences.

The percentile bands are model spread, not confidence intervals or calibrated
probabilities. They omit or incompletely sample observational uncertainty,
parameter-fit uncertainty, internal variability and model dependence. Models
sharing components are not independent simply because one member was selected.

`member_combine: mean|median` combines members before downstream correction and
index analysis. That can help estimate a forced mean, but suppresses variability
and changes extreme statistics. For extreme-event distributions, analyse members
separately and combine summaries with an explicit hierarchy; the bundled flat
multi-model workflow does not implement that full hierarchy.

Equal model weights are the transparent starting point. Inverse historical-RMSE
weights, where requested, are heuristic and require independent validation;
they do not implement model-independence weighting. Percentile bands must not be
labelled weighted probabilities when only the centre is weighted.

OLS slopes are descriptive trends, with default p-values relying on usual error
assumptions. Theil–Sen and the supplied Mann–Kendall function are alternatives,
but neither automatically removes serial-correlation effects. Use a separately
justified autocorrelation-aware uncertainty method when required. Map stippling
based on agreement indicates the fraction sharing a sign, not a significance
level or detection/attribution claim.

## 5. Annual GEV analysis

Enable this route explicitly and describe the sample. A defensible block-maxima
sample contains one maximum per approximately independent complete year/season.
A seasonal mean, total or arbitrary index is not a block maximum. If an empirical
GEV is deliberately fitted to such an index, record `sample_kind: annual_index`,
justify it against alternatives, and avoid claims from extreme-value asymptotics.

For z = (x−μ)/σ and SciPy's shape c = −ξ:

~~~text
GEV: F(x) = exp(−[1 − c z]^(1/c)),  1 − c z > 0, c != 0
Gumbel: F(x) = exp(−exp(−z)),       c = 0
~~~

Positive SciPy c implies an upper bound; negative c implies an unbounded heavy
upper tail. SciPy's GPD shape has the opposite convention: its c equals ξ.
Do not copy a GEV sign interpretation onto a GPD table.

For one maximum per year, the upper-tail return level is F−1(1−1/T), T > 1.
The return period of x is 1/SF(x), using the survival function for numerical
stability. A finite T-quantile approaches a bounded endpoint smoothly; it does
not become invalid simply because the fitted distribution is bounded.

All included fits are stationary. A fit over a strongly warming century pools
different climates. It is a window-conditional summary, not a return period for
every year in that century. Prefer a scientifically meaningful window; shorter
windows trade reduced nonstationarity for greater parameter uncertainty.

Use sample count, support, scale, fit-method/fallback records, diagnostic plots
and sensitivity analysis. The naive one-sample KS p-value from the same data used
to estimate the parameters is uncalibrated. `ks_p_naive` is diagnostic only;
`ks_p` is left missing when no fitted-distribution Monte Carlo test was performed.
A valid Monte Carlo test must re-fit each simulated sample. A large p-value does
not establish that extrapolated tails are accurate.

Fit each model separately and summarize its return levels across models.
Do not fit an ensemble-mean time series and call it a typical model's extremes.
Model spread and a single fit's sampling interval answer different questions.
A bootstrap of an ensemble-median annual series would concern that derived
series, not the median of model return-level curves. The CLI therefore rejects
nonzero extremes.bootstrap until an appropriate multi-model method is implemented.
The separate single-series bootstrap helper does not change this limitation.

Scenario shape comparisons are exploratory: the same models recur across
scenarios, violating an independent-groups reading of Kruskal–Wallis. Use paired
or hierarchical inference for a formal scenario test. Do not treat shape signs
from short samples as certain changes in physical tail behaviour.

## 6. Daily POT and block maxima

The daily route reads daily fields within a bounding box, constructs a regional
series, adjusts it on the chosen scale and fits extremes. It does not calculate
heatwave duration, storm tracks, spell length or compound-event dependence.

### Block maxima

`daily.extremes.kind: block` takes annual/seasonal maxima of the daily regional
series, then fits GEV/Gumbel. Preserve the calendar and complete seasonal blocks.
A minimum temperature index needs an explicit lower-tail/sign transformation;
setting `kind: block` alone always targets upper maxima.

### Peaks over threshold

Choose threshold u. Join exceedances whose date separation is at most `decluster_days` over
continuous finite daily observations; retain one peak per cluster. Missing or
noncontinuous steps split a cluster rather than being assumed below threshold. Missing days and season gaps
are not known below-threshold days. Estimate λ as cluster count divided by the
valid calendar-year/season exposure used in the analysis, and document the
coverage rule. In the date-aware implementation, each season contributes its
finite observed days divided by its calendar-specific expected days. Partial
coverage therefore gives fractional exposure. This assumes missingness does not
preferentially remove extremes and cannot reconstruct lost peaks. A season-only
rate means events during that season per year.

For excess y = x−u, the GPD survival function is

~~~text
S(y) = (1 + ξ y/σ)^(−1/ξ),  y >= 0, 1 + ξ y/σ > 0
S(y) = exp(−y/σ),          ξ = 0
~~~

The mean waiting time between events above x is

~~~text
T_wait(x) = 1 / [λ S(x−u)]
z_T = u + σ/ξ [(λT)^ξ − 1]     (ξ != 0)
z_T = u + σ log(λT)            (ξ = 0)
~~~

These are event-rate return periods. Under a Poisson assumption, the probability
of at least one exceedance in a year is 1−exp[−λS(x−u)], so the inverse annual
exceedance probability differs from T_wait except in the rare-event limit.
Return levels below u lie outside the fitted excess model and should be missing.

| Threshold mode | Interpretation |
|---|---|
| `observed` | One threshold estimated from the reference baseline and reused for comparability. |
| `absolute` | One physical threshold in the converted variable units. |
| `quantile` | Each series' own quantile; thresholds move, so exceedance counts do not compare the same event. Clustering can still change the peak rate. |

Check threshold stability and declustering sensitivity. Thirty peaks or twenty
years are computational cutoffs, not guarantees of precise hundred-year levels.

## 7. What this skill does not infer

A successful pipeline run does not establish causality, future forecast skill,
stationarity, calibrated risk probabilities, physical consistency among separately
adjusted variables, or a publication-ready scientific conclusion. Missing
capabilities should become a named, tested extension, not an undocumented
change to the meaning of an existing option.
