# Methods, equations and assumptions

Everything below is what `cmipkit` actually computes. Notation: `x(t, i, j)` is a
model field for year `t` at grid point `(i, j)`; `o(t, i, j)` the observational
reference; `μ`, `σ` the temporal mean and sample standard deviation over a common
baseline; `Q(p)` a quantile function; `F` a CDF.

---

## 1. From files to one value per year

For each file the loader selects the configured months, applies the unit
transform (`x·scale + offset`), and reduces within each year:

| `aggregate.stat` | reduction over the selected months/days |
|---|---|
| `mean` | `mean_t x` |
| `sum` | `Σ_t x` |
| `max` / `min` | `max_t x` / `min_t x` |
| `median`, `std` | `median_t x`, `std_t x` |
| `p95`, `p99`, … | `quantile_t x` at that percentile |
| `count_gt` / `count_ge` / `count_lt` / `count_le` | number of timesteps above/below `aggregate.threshold` |
| `max_abs`, `mean_abs` | `max_t |x|`, `mean_t |x|` |

Years follow the calendar of the file. **DJF** is special: December is assigned to
the *following* year (`Dec 1999` → winter 2000), so a DJF series never mixes two
winters in one value.

All fields are then interpolated bilinearly (`grid.method: linear`) or
nearest-neighbour to `grid.lat`/`grid.lon`; longitudes are converted to the
configured convention and the seam at 0/360 is closed before interpolation.

---

## 2. Indices

### 2.1 Area-weighted regional mean

```
I(t) = Σ_ij w_i x(t,i,j) m_ij / Σ_ij w_i m_ij ,      w_i = cos φ_i
```

with `m_ij` the region mask (rectangles may wrap across 0°). Missing values are
excluded from both sums, so a partly-missing field still yields a value — check
`n` in the tables if that is a concern.

### 2.2 Zonal-anomaly region difference

```
x'(t,i,j) = x(t,i,j) − mean_j x(t,i,·)
I(t) = mean_{region B} x' − mean_{region A} x'
```

This removes the zonally symmetric component (e.g. the global-mean warming) and
isolates the longitudinal contrast — the "ridge/trough" or "region2 − region1"
index used in the reference workflow.

### 2.3 Projection (wave-train) index

Template `T` is the zonal anomaly of the reference climatology over the template
period (or a chosen year/model):

```
T(i,j) = mean_template[ x_ref(t,i,j) ] − zonal mean
I(t) = Σ_ij x'(t,i,j) T(i,j) w_i / Σ_ij T(i,j)² w_i
```

`I = 1` means "the anomaly field looks exactly like the template". `I` is
dimensionless and linear in the field amplitude, which makes it suitable for
trend and GEV analysis of circulation patterns. The template is fixed from the
reference (historical) climate, so future values measure similarity to that
pattern, not to the future mean state.

### 2.4 Anomaly / standardisation

`anomaly: true` subtracts the **reference climatology** (not the field's own
mean): for historical fields the baseline period of the same field; for future
fields the historical corrected climatology of the same model. This keeps all
models and scenarios on one common baseline and makes `I = 0` mean "baseline
climate". `standardize: true` additionally divides by the baseline standard
deviation (z-score); `detrend` removes the series' linear trend.

---

## 3. Bias correction

Fitted per grid point on the years common to the reference and the model
historical run inside `correction.params_period`, then applied unchanged to the
historical and future fields.

### 3.1 `mean` (first moment, additive)

```
x_corr = x + (μ_o − μ_m)
```

### 3.2 `mean` (first moment, multiplicative)

```
x_corr = x · μ_o / μ_m
```

Used when `correction.multiplicative: auto` resolves to true (precipitation/wind
like, non-negative). Guarded against zero/negative denominators.

### 3.3 `mean_std` (additive, the 6.25 / 5.19 method)

```
x_corr = (x − μ_m) · σ_o/σ_m + μ_o
```

Corrects mean and variance; the historical distribution then matches the
observed mean and standard deviation exactly. It assumes the model's *bias* is
stationary and that its future trend is already realistic — the correction does
not change the trend because it is a linear rescaling.

### 3.4 `mean_std` (multiplicative, mean + CV anchored)

For positive ratio-scale variables a single multiplicative factor cannot fix mean
and variance simultaneously, so two steps are used:

```
x'  = x · μ_o/μ_m
x'' = μ_o + (x' − μ_o) · (CV_o / CV_m),      CV = σ/μ
```

`x''` has mean `μ_o` and standard deviation `σ_o` on the baseline. Outside the
baseline the anomaly part is scaled by `CV_o/CV_m`, which preserves the model's
relative change. Negative results are clipped at `correction.lower_bound`.

### 3.5 `quantile_mapping` (empirical, per grid point)

With `n_quantiles` probabilities `p_k = k/(n−1)`:

```
x_corr = Q_o( F_m(x) )
```

`Q_o`, `Q_m` are empirical quantile tables of the baseline reference/model series;
the mapping between quantiles is linear and constant beyond the sampled range.
The corrected historical distribution equals the observed one **by
construction** — which is why KS drops to the sampling noise.

**Warning.** Applied to future values it maps them onto *today's* distribution
and therefore erases the projected change in the tails. Use it for the historical
period or for a deliberately stationary counterfactual only.

### 3.6 `quantile_delta_mapping` (QDM, Cannon 2015) — additive form

For a future value `x` with future-CDF position `p = F_fut(x)`:

```
x_corr = x + [ Q_o(p) − Q_m(p) ]
```

The model's own quantile change is added to the observed quantile, so the
projected distribution shift survives while the historical distribution is
anchored. Multiplicative form: `x_corr = x · Q_o(p)/Q_m(p)`.

* For the historical period `p = F_m(x)` and QDM reduces to QM.
* QDM preserves the quantile-wise change, **not** the mean trend exactly: a
  slowly varying correction term remains in the mean. Diagnostic runs show
  ~10–15 % deviation in mean-trend slope on synthetic data; state this if trend
  magnitude is the headline result.

### 3.7 `delta_change`

Alias of `mean`: the impact-modelling name for a first-moment correction
(`obs climatology + model change`). Identical arithmetic, kept because papers and
users refer to it by that name.

### 3.8 Trend correction (`correction.trend_correction`)

Applied to the index *series* after correction, relative to
`plots.timeseries.reference_year` `t0`:

```
difference:  x_corr(t) = x(t) − (b_m − b_o)·(t − t0)
ratio:       x_corr(t) = x(t) · b_o/b_m
```

with `b_m` the model's historical trend and `b_o` the observed historical trend.
This forces the projected trend to equal the observed one; it is a strong
assumption and must be reported. Leave it disabled unless the user asks for it.

---

## 4. Ensemble statistics

Matrix `X[m, t]` (models × years, NaN where a model has no data).

* centre: `mean` or `median` over models (`ensemble.statistic`);
* bands: percentiles over models (`ensemble.bands`), typically 10/90 and 20/80;
* `ensemble.weights: skill` uses inverse-RMSE weights derived from the historical
  scores (`w_m ∝ 1/RMSE_m²`, normalised);
* a model contributes only where it has data, so bands narrow where fewer models
  are available — the number of models is drawn in the legend.

Percentile bands describe **model spread**, not a confidence interval, and they
are not a probability distribution: do not attach a probability to "outside the
10–90 % band".

---

## 5. Trend tests

Ordinary least squares (`scipy.stats.linregress`) on the annual series, reported
as `slope` (per year) and `slope_per_decade` (×10) with `p`, `r²`, `stderr`.
Theil–Sen and Mann–Kendall (with tie correction) are available in
`cmipkit.indices` for robust/serial-correlation-aware checks.

Significance of the *difference* between the ensemble-mean future trend and the
historical trend is tested with a one-sample t-test on the per-model slope
differences (`compare_slopes`), which is a model-uncertainty test only.

Field significance for map stippling uses a different criterion: **model
agreement** — a grid point is stippled when at least `plots.maps.agreement`
(default 80 %) of the models agree on the sign of the change. This is a
consensus statement, not a statistical test; say so in the caption.

For multiple testing, `indices.fdr_bh` implements Benjamini–Hochberg FDR control.

---

## 6. Extreme-value analysis

### 6.1 Distribution

`scipy.stats.genextreme` with shape parameter `c = −ξ` (Coles 2001 convention):

```
GEV CDF:  F(x) = exp( −[1 + c (x−μ)/σ]^(−1/c) ),   c ≠ 0
          F(x) = exp( −exp(−(x−μ)/σ) ),             c = 0  (Gumbel)
```

* `c > 0` — Weibull-type, bounded upper tail (common for temperature maxima);
* `c < 0` — Fréchet, heavy tail (common for precipitation);
* `c = 0` — Gumbel (`extremes.distribution: gumbel`).

### 6.2 Return level and return period

```
return level:   z_T = F^{-1}(1 − 1/T)
return period:  T(x) = 1 / (1 − F(x))
```

`T = 100` means "exceeded once per 100 years on average". The return period of a
*record* (e.g. the 2025 value) is the more intuitive number for impact stories.

### 6.3 Fitting

Maximum likelihood by default (`fit_method: mle`), falling back to method of
moments and finally to a Gumbel fit with moment-matched scale. Fits with fewer
than `min_years` (default 20) values, or with zero variance, are skipped. For each
fit the pipeline records the log-likelihood, AIC, and a KS test of the fitted CDF.

### 6.4 Multi-model aggregation

Each model gets its own fit; the plotted curve is the **median of return levels
across models** at each return period, with percentile bands across models. This
mixes only model uncertainty — sampling uncertainty of each fit is separate. Set
`extremes.bootstrap: 200` to add a percentile bootstrap of the return levels for
a single fit (written to `extremes/gev_<scenario>.npz`).

### 6.5 Shape comparison

`kruskal_scenarios` runs a Kruskal–Wallis test on the per-model shape parameters
across scenarios; a small `p` means the scenarios produce systematically
different tail behaviour. With fewer than two scenarios the test is `nan` and the
report omits it.

### 6.6 Assumptions and caveats

* The fit assumes the annual values are (approximately) independent block maxima.
  An annual-mean index is *not* a block maximum — the fit is then a tail model for
  a smoothed quantity, and should be described as such.
* Short future windows (e.g. 25 years) give unstable shape parameters; compare
  medians across models rather than a single fit.
* GEV extrapolation beyond the sampled return period (e.g. T > 5×n) depends
  almost entirely on the fitted shape; report it with the fit's CI or not at all.

---

## 7. Multi-member ensembles

When a model provides several realisations (`r1i1p1f1`, `r2i1p1f1`, …),
`data.models.member_combine` selects the treatment:

* `first` — use one member (the configured order, or the one with the most files).
  Each member is an equally valid realisation; one is enough for a multi-model
  ensemble and it keeps the models independent.
* `mean` / `median` — load every configured (or discovered) member for the same
  model/experiment, align them on the common years and average per grid point:

  ```
  x_model(t,i,j) = mean_m x_m(t,i,j)
  ```

Internal variability is reduced by roughly `1/sqrt(n_members)`, so the model's
forced response is better estimated; the price is that the ensemble spread across
models narrows and part of the "model uncertainty" becomes "member uncertainty".
Fits, trends and GEV results are all computed on the member-averaged series, and
the report records how many members each model used.

## 8. Daily extremes: POT/GPD and daily block maxima

The annual pipeline cannot represent a heatwave or a storm: it has already
averaged the season away. `daily.*` re-reads the daily data for a season, reduces
it to a daily regional index, corrects that series, and fits an extreme-value
model to the daily values.

### 8.1 Peaks over threshold (GPD)

Choose a threshold `u` (observed baseline quantile, a per-series quantile, or an
absolute value). Decluster by runs: exceedances separated by fewer than
`decluster_days` days belong to one event, which contributes its maximum. With
`n` peaks over `N` years the exceedance rate is `λ = n/N`.

The excesses `y = x − u` follow a Generalized Pareto distribution:

```
G(y) = 1 − (1 + ξ y/σ)^(−1/ξ),   y > 0, 1 + ξy/σ > 0
```

and the return level for a return period `T` (years) is

```
z_T = u + (σ/ξ)[(λT)^ξ − 1]        ξ ≠ 0
z_T = u + σ ln(λT)                 ξ = 0
```

The return period of a value `x > u` is

```
T(x) = (1/λ) [1 + ξ (x−u)/σ]^(−1/ξ)
```

The Poisson–GPD combination is the standard way to convert a daily series into
"once in T years" statements (Coles 2001, ch. 4). Threshold choice is the main
subjective step: too low violates the asymptotic justification, too high leaves
too few peaks (the stage requires `min_peaks`, default 30). `threshold_mode`
controls whether the threshold is fixed to the observed climate
(`observed`), recomputed per series (`quantile`), or physically fixed
(`absolute`).

### 8.2 Daily block maxima

`kind: block` takes the maximum of the daily series in each year and fits a GEV —
the classical block-maxima approach, equivalent to `aggregate.stat: max` on daily
input. It is more robust to threshold choice but wastes the within-season
information and needs ≥ 20 years.

### 8.3 Correction and uncertainty

The daily regional series is corrected with the same transfer functions as the
annual pipeline (`correction.fit_series`), by default QDM. The correction is
applied before the extreme fit, so the fitted tail reflects the bias-corrected
climate. Return levels are aggregated across models with the median and
percentile bands, exactly like the annual GEV.

## 9. References

Cannon, A. J. (2015): *Bias correction of GCM precipitation by quantile mapping:
how well do methods preserve changes in quantiles and extremes?* J. Climate 28,
6938–6959. — QDM.

Coles, S. (2001): *An Introduction to Statistical Modeling of Extreme Values.*
Springer. — GEV conventions, return levels.

Eyring, V. et al. (2016): *Overview of the Coupled Model Intercomparison Project
Phase 6 (CMIP6).* GMD 9, 1937–1958. — experiment design.

Maraun, D. (2013): *Bias correction, quantile mapping, and downscaling: revisiting
the inflation issue.* J. Climate 26, 2137–2143. — QM pitfalls.

Maraun, D. et al. (2010): *Precipitation downscaling under climate change.*
Rev. Geophys. 48, RG3003. — delta change vs statistical correction.

Tebaldi, C. & Knutti, R. (2007): *The use of the multi-model ensemble in
probabilistic climate projections.* Phil. Trans. R. Soc. A 365, 2053–2075. —
ensemble interpretation.

Zhang, X., Zwiers, F. W. & Li, G. (2004): *Monte Carlo experiments on the
detection of trends in extreme values.* J. Climate 17, 1945–1952. — threshold
choice and declustering sensitivity.
