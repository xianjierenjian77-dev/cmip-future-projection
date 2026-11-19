# Figures

All figures are produced by `cmipkit/plotting.py` and written to
`<outdir>/figures/`. They are matplotlib objects saved with
`bbox_inches="tight"` at `plots.dpi`; the default backend is `Agg` when no
display is attached, so runs over SSH never block.

Set `plots.formats: [png, pdf]` to get vector output for journals.

## Catalogue

| Name in `plots.figures` | File | Content |
|---|---|---|
| `timeseries` | `fig_timeseries.png` | one row per scenario: observed history, individual model years, ensemble mean, 10–90 and 20–80 bands, historical and future trend lines with `/10a` values, reference-year marker |
| `gev` | `fig_gev_return_levels.png` | one row per scenario: historical GEV curve (black), future multi-model median + bands, return-period grid at 10/100 a, optional empirical points, reference value line |
| `combined` | `fig_scenarios_combined.png` | n×2 panel: time series on the left, GEV on the right (the 6.25 Fig. 4 layout) |
| `maps` | `fig_maps_<scenario>.png` | period-mean change maps, one panel per `periods.future_windows` entry, model-agreement stipple, region boxes, shared horizontal colour bar |
| `kde` | `fig_diagnostics_kde.png` | per-model distributions: observed vs raw vs corrected, KS statistics in each panel |
| `kde` (per model) | `fig_pdf_<model>.png` | one larger PDF panel per model, with sample counts and KS values (`plots.kde.per_model`) |
| `timeseries` (per model) | `fig_timeseries_<model>.png` | one model's own future projection for every scenario, ensemble mean dashed (`plots.timeseries.per_model`) |
| `daily` | `fig_daily_return_levels.png` | GPD (POT) or GEV (block) return-level curves of the daily extremes, observed history + each scenario with model spreads |
| `shape_box` | `fig_gev_shape.png` | box plot of the GEV shape parameter per scenario + Kruskal–Wallis p |

Any name not in the list is rejected at config validation, so a typo cannot
silently produce nothing.

## Sizing and typography

| Knob | Effect |
|---|---|
| `plots.dpi` | 300 for print, 110–150 for drafts |
| `plots.figsize_scale` | 1.0 default; 0.9 for two-column journals, 1.15 for slides |
| `plots.font` | family name; `null` keeps DejaVu Sans |
| `plots.language: zh` | forces a CJK font even without Chinese in the auto labels |

Chinese labels: WenQuanYi Micro/Zen Hei are installed on this server and are
selected automatically when any label contains CJK characters. If a label renders
as boxes, set `plots.font: WenQuanYi Micro Hei` explicitly.

## Time-series panel

* Individual model years are drawn as small grey points (`plots.timeseries.show_members`).
  With more than ~20 models this becomes a grey cloud; turn it off.
* Bands: the first pair in `ensemble.bands` is the lighter 10–90 %, the second the
  20–80 %. Provide four values, e.g. `[10, 90, 20, 80]`.
* The observed line and the future ensemble are joined through
  `plots.timeseries.reference_year` (dotted grey), which is what makes the
  "history meets projection" visual work. Set it to the last observed year if the
  reference series ends earlier (e.g. 2014).
* Trend annotations are computed on the plotted series (historical period and
  future period separately) and expressed in `variable.units` per decade.
* If a scenario's series is missing, the row is omitted rather than drawn empty.

## GEV panel

* Return period axis is logarithmic from `plots.gev.xlim[0]` to `xlim[1]`.
* An upper-bounded GEV (`c > 0`) yields `NaN` above its bound; matplotlib simply
  stops the curve there — that is correct, not a bug.
* Multi-model bands are percentiles **across models**, not fit uncertainty. For
  fit uncertainty use `extremes.bootstrap` (written to the NPZ; not drawn).
* The purple horizontal line is the reference value (e.g. the observed record);
  the star marks its median return period across models.

## Maps

* Colour limits default to the 98th percentile of |change| (symmetric) — robust to
  one outlier grid point. Fix them with `plots.maps.levels: [-10, 10]` when
  comparing several scenarios side by side.
* `RdBu_r` for signed changes, `BrBG_r` for moisture/flux divergence,
  `viridis`/`YlOrRd` for absolute values. Set `plots.maps.cmap`.
* Stippling is **model agreement**, not a significance test: a point is dotted
  when at least `plots.maps.agreement` of the models agree on the sign. Say
  "stippling: ≥80 % model agreement" in the caption.
* Coastlines/borders use the local cartopy cache in
  `~/.local/share/cartopy` (110 m). If the cache is missing and the machine is
  offline, the figure is still produced without coastlines and a warning is
  logged. Never delete that cache to "fix" anything.
* `plots.maps.extent` defaults to the data extent; set it explicitly for a
  consistent framing across figures.
* Region boxes come from `index.region` and `index.regions`.

## KDE / PDF panel

* **What is plotted.** By default (`plots.kde.mode: auto`) the samples are every
  grid point of the study region over the baseline period, pooled per model —
  the same construction as the 5.19 diagnostic KDE, which makes the density
  smooth and describes the variable rather than a 65-value index. `series` mode
  plots the annual index instead; `field` forces the gridded sampling.
* Because of that, the KS printed on the figure is computed on the *plotted
  samples* and can differ slightly from `tables/model_scores.csv`, whose KS is
  computed on the annual index series. Both are reported; say which one you quote.
* `plots.kde.per_model: true` (default) writes one `fig_pdf_<model>.png` per model
  with the full-size curves and sample counts.
* Grid adapts to the number of models (`kde.grid` is a maximum), so three models
  give one row, not a 4×6 grid of empty panels.
* KS values are computed on the **absolute** historical index series (not on
  anomalies), matching `tables/model_scores.csv`.
* Use this figure to justify excluding a model: a corrected KS far above the
  ensemble median, or a distribution that does not move towards the observed one.

## Shape-parameter panel

* Each box is the distribution of per-model GEV shape parameters in one scenario.
* A positive median shape (`c > 0`) means a bounded upper tail; negative means
  heavy-tailed. Shifts between scenarios indicate changed tail behaviour.
* The Kruskal–Wallis p-value is printed only when at least two scenarios have
  fits.

## Adding or changing a panel

1. Add a function to `cmipkit/plotting.py` that takes `cfg` plus dataclasses and
   returns `_save(fig, outbase, cfg)`.
2. Register the name in `cmipkit/config.py::validate_config` (`valid_figs`).
3. Call it from `Pipeline.figures` with data already available in `self.ens`,
   `self.gev`, `self.series` or `self.map_acc` — do not read files in the plotting
   layer.
4. Add its config block under `plots.<name>` in `DEFAULTS`.
5. Re-run `run.py selftest` and eyeball the new figure.

## Stitching separately generated panels

If a figure must combine panels with different sizes (e.g. two maps and a time
series), the 6.25 workflow stitched PNGs with PIL. Prefer a native matplotlib
`GridSpec` (as `scenarios_combined` does) because it keeps fonts and DPI
consistent. If stitching is unavoidable:

```python
from PIL import Image
imgs = [Image.open(p) for p in paths]
width = max(i.width for i in imgs)
imgs = [i if i.width == width else i.resize((width, round(i.height * width / i.width)), Image.LANCZOS)
        for i in imgs]
canvas = Image.new("RGB", (width, sum(i.height for i in imgs) + 10 * (len(imgs) - 1)), "white")
y = 0
for im in imgs:
    canvas.paste(im, (0, y)); y += im.height + 10
canvas.save(out, dpi=(300, 300))
```
