# Figures and captions

Use the existing matplotlib functions for supported outputs. Figures are written
to the run's `figures` directory. Set `plots.formats: [png,pdf]` for raster preview
and vector export, then inspect the actual rendering. Plot styling does not
replace scientific checks of the data being displayed.

## Choose useful figures

| Name in `plots.figures` | Product | What to explain in its caption |
|---|---|---|
| `timeseries` | Scenario time series | Index definition, historical baseline, model counts and model-spread bands. |
| `maps` | Period-mean maps per scenario | Absolute/change/anomaly quantity, periods, units, region and agreement rule. |
| `kde` | Observed/raw/adjusted distributions | Whether samples are annual indices or pooled space-time fields; calibration period. |
| `gev` | Annual return-level curves | Sample kind, stationary fit windows, SciPy shape convention and model spread. |
| `combined` | Time series with annual GEV panels | Both sets of assumptions; only request when GEV is scientifically justified. |
| `shape_box` | Annual GEV shape summaries | Estimated shapes, sample sizes and exploratory scenario-comparison limitations. |
| `daily` | Daily POT/block return-level figure | Regional operator, threshold/clustering or block definition, exposure and fitting window. |

Daily plotting is also controlled by `daily.plots`. Annual GEV panels require
annual `extremes.enabled`; enabling daily extremes does not create annual GEV
results. Per-model time series and PDF figures use
`plots.timeseries.per_model` and `plots.kde.per_model`.

## Style controls

| Setting | Use |
|---|---|
| `plots.dpi` | 120–150 for drafts, 300 when raster publication output is requested. |
| `plots.formats` | Requested matplotlib-supported file formats, e.g. png/pdf/svg. |
| `plots.figsize_scale` | Adjust figure/font dimensions for the destination. |
| `plots.font` | Installed font family; do not assume a server-specific CJK font exists. |
| `plots.language` | Automatic label language. Inspect untranslated custom labels. |
| `plots.colors.scenario` | Stable colors for scenarios across all figures. |

Use a sequential color map for nonnegative absolute quantities and a diverging
map for signed differences. Use common color limits for cross-scenario maps
when direct comparison matters. Avoid very pale lines and crowded annotations.

## Time series

The centre is the configured ensemble mean or median. Bands are cross-model
percentiles, not confidence intervals and not prediction probabilities. The mean
is not mathematically required to lie between the 10th and 90th percentiles of a
small or skewed ensemble. Do not clamp it to the band to make a figure look right.

`show_members` draws model values; `show_bands` and `show_trend` control those
layers. `reference_year` is an optional marker, not permission to extrapolate
observations or join missing years. Report the actual historical/projection
coverage, including gaps. Trend units are index units per decade, not necessarily
the input field's units if the index is standardized or dimensionless.

## Maps

Use `plots.maps.variable` to distinguish change, absolute and anomaly products.
Use `levels`, `n_levels`, `cmap` and `extent` for a consistent comparison.
Verify that the numerical field, colorbar and caption describe the same baseline.

Agreement stippling is the fraction of available models sharing the change sign.
Caption it as agreement, such as “at least 80% model sign agreement”; do not call
it a 95% significance test. Model counts and masks can vary across grid cells.
A single model cannot establish robust multi-model agreement. On absolute-value
maps, an agreement pattern about change needs explicit explanation or omission.

Cartopy coastlines need local data or a permitted download. When unavailable,
plain longitude/latitude maps remain useful; do not silently label a different
map projection. Do not delete a working coastline cache to address unrelated
plot errors.

## Distribution diagnostics

`kde.mode: series` plots annual index samples. `field` pools spatial-temporal
field values; `auto` selects the available route. Those are different sampling
units, so their KS statistics can differ. Pooled cells are spatially dependent;
a smooth KDE is not evidence of thousands of independent samples.

Calibration distributions often improve by construction. Do not use a favorable
calibration KS p-value as proof of future skill, or remove models solely to make
curves/p-values look good. Name any independent validation period in the caption.
Degenerate or too-short samples should remain visibly unassessed rather than be
replaced by artificial noise.

## Extreme-value figures

Return-level curves summarize fitted distributions under their windows and event
definitions. Bound the plotted return-period range to the scientific request;
long extrapolations need uncertainty and sensitivity analysis. The plotting code
may provide a broad default axis, which is not validation of its most distant
return levels.

A bounded GEV approaches its endpoint at long T. For a reference value above an
upper endpoint, an infinite return period is a consequence of that fitted model,
not a physical impossibility established from the data. Preserve this distinction.

Across-model bands represent variation among fitted model curves. Sampling
uncertainty for each fit requires separate treatment. Pipeline bootstrap is not
an implemented multi-model confidence-interval method. Shape comparisons across
scenarios use repeated models, so their independence assumptions need scrutiny.

## Finishing and extensions

Inspect label clipping, units, actual scenario names, empty panels, colorbars,
threshold lines and model counts at the final figure size. If the user asks only
for a new figure from saved results, reuse those results with a small explicit
plotting adapter; there is no obligation to recalculate all climate fields.

For a reusable new figure type, add its plotting function, register the config
name, connect it to available stage outputs and verify the result on a small
known dataset. Prefer native matplotlib layouts to stitching raster panels.
