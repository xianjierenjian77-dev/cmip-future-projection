# Operating procedure

This is the checklist an agent follows for a new projection task. Every step ends
with a **checkpoint** — do not continue while a checkpoint fails.

Conventions: `PY=/disk1/wsr/miniconda3/envs/climate/bin/python`,
`SK=/disk1/wsr/.dsh/skills/cmip-future-projection`.
Work in a project directory (e.g. `/disk1/wsr/proj/<name>/`), never inside the
skill directory.

---

## Step 0 — environment

```bash
$PY $SK/scripts/run.py doctor
```

Checkpoint: `RESULT: ready`. If packages are missing, stop and report; do not
`pip install` into the shared `climate` env without the user's approval.

Optional end-to-end sanity check of the installation (≈10 s, writes to a temp dir):

```bash
$PY $SK/scripts/run.py selftest
```

---

## Step 1 — reconnaissance of the data

Never write a config from assumptions about file layout, variable names, units or
calendars. Inspect first.

**Gridded NetCDF (CMIP6 or any model):**

```bash
$PY - <<'EOF'
import xarray as xr, glob
f = sorted(glob.glob('/disk5/cmip6/monthly/historical/ACCESS-CM2/tas_Amon*.nc'))[0]
ds = xr.open_dataset(f, use_cftime=True)
print(ds)                       # variable names, dims, plev, calendar
v = ds['tas']
print(v.dims, v.shape, v.attrs.get('units'))
print(float(v.min()), float(v.max()))
print(v.time.values[[0, -1]])
EOF
```

**GRIB (ERA5):**

```bash
$PY - <<'EOF'
import xarray as xr
ds = xr.open_dataset('/disk3/era5-1.0/monthly/surface/t2m.2000.grib', engine='cfgrib',
                     backend_kwargs={'indexpath': ''})
print(ds); print(float(ds.t2m.min()), float(ds.t2m.max()))
EOF
```

Record and verify:

| Property | Why it matters | Red flags |
|---|---|---|
| variable name in file | `data.*.variable` | `t2m` vs `tas`, `tp` vs `pr` |
| units | `units_transform` | K vs °C (≈273), m vs mm (1000), kg m⁻² s⁻¹ vs mm/day (86400) |
| time axis | aggregation, calendars | cftime 360-day, monthly vs daily, `valid_time` vs `time` |
| dimensions | plev needed? | `plev` in Pa (CMIP6) vs hPa (ERA5) |
| latitude order | regridding | 90→−90 vs −90→90 (handled, but check bbox) |
| longitude convention | region selection | 0–360 vs −180–180 (handled) |
| grid resolution | runtime & memory | 1° global 75 yr ≈ 65 k points, fine; 0.25° is heavy |
| members available | consistency | `r1i1p1f1` vs `r101i1p1f1` |
| missing years | trends, GEV | gaps are skipped silently — check `n` in the report |

Checkpoint: you can state min/max/mean, units, calendar, and the exact file glob
that matches **historical** and each **scenario** for **one** model. Then confirm
the same glob works for a second model (layouts differ: some archives keep
`<var>/` sub-directories).

---

## Step 2 — analysis design

Write the design down (it becomes the config comments):

1. **Variable & units out.** e.g. `tas` in °C (transform `{scale: 1, offset: -273.15}`),
   or `pr` in mm/day (`{scale: 86400}` if the file is kg m⁻² s⁻¹).
2. **Season.** `[9, 10]` for Sep–Oct, `"JJA"`, `"ANN"`, or `"DJF"`
   (December is assigned to the following winter year).
3. **Aggregation.** `mean` (temperature, height), `sum` (seasonal precipitation),
   `max` (TXx), `min` (TNn), `p95` (hot-day threshold), `count_gt`
   (days above a threshold; requires `threshold`).
4. **Index.** See the table in `SKILL.md`. For a wave-train style analysis use
   `projection_index` with the reference field as template
   (`index.projection.template_source: reference`).
5. **Periods.** Baseline for fitting: usually 1950–2014 (the CMIP6 historical
   end). Projection: 2026–2100 (or the horizon the user cares about). Use
   `periods.future_windows` for the map panels.
6. **Method.** `mean_std` for a general variable; `quantile_delta_mapping` for
   extremes/precipitation; `mean` when the model trend must be preserved
   unchanged; `none` only for a control.
7. **Extremes.** Only meaningful when one value per year is a block maximum or
   an annual extreme index. Say so explicitly in the report ("GEV fitted to the
   annual regional-mean index").
8. **Models/scenarios.** Keep the list explicit in the config when the user named
   models; `models: all` only for a first exploration.

Checkpoint: the design can be explained in five sentences without reference to
code.

---

## Step 3 — config

```bash
$PY $SK/scripts/run.py init /disk1/wsr/proj/<name>
# edit /disk1/wsr/proj/<name>/config.yaml
```

Start from a template when it fits:

| Template | Scenario |
|---|---|
| `templates/config_grid_projection.yaml` | index + 3 scenarios + maps + GEV (6.25-style) |
| `templates/config_regional_extremes.yaml` | one region, extreme index, return period of a record |
| `templates/config_precip_qdm.yaml` | precipitation, multiplicative QDM, non-negative |
| `templates/config_series.yaml` | pre-computed annual series for reference + models |

Rules:

* Absolute paths, or paths relative to the working directory (not to the config).
* `{year}` in `data.reference.paths` is expanded over `data.reference.years`.
* `{root}`, `{model}`, `{experiment}`, `{member}`, `{variable}` are expanded in
  `data.models.pattern(s)`.
* Leave `data.models.patterns` empty to use the built-in CMIP6 layouts
  (flat and `<var>/` sub-directory); add explicit patterns only for other archives.
* Put the scientific justification in YAML comments — the next agent reads them.

Checkpoint: `$PY $SK/scripts/run.py run config.yaml --stages prepare` succeeds
and prints `reference field: N years ... grid AxB` and the model list.

---

## Step 4 — pilot run (cheap)

Restrict to one model and one short projection period before the full run:

```yaml
data: {models: {models: [ACCESS-CM2], scenarios: [ssp245]}}
periods: {future: [2026, 2050]}
plots: {figures: [timeseries], dpi: 110}
```

```bash
$PY $SK/scripts/run.py run config.yaml
```

Checkpoint:

* `report.md` exists and the bias drops from raw to corrected;
* `figures/fig_timeseries.png` shows history joining the future;
* the runtime is acceptable; the corrected files are the size you expect.

Only then widen models/scenarios/period and raise `plots.dpi`.

---

## Step 5 — production run

```bash
cd /disk1/wsr/proj/<name>
nohup $PY $SK/scripts/run.py run config.yaml > run.out 2>&1 &
# then poll: tail -f run.out
```

Long runs (many models × yearly files) take minutes to tens of minutes. While it
runs, do not start a second run into the same `outdir`.

Checkpoint: `=== done in N s ===` in the log and `report.md` updated.

---

## Step 6 — quality control

Read `report.md` and check numerically:

| Metric | Expected | If not |
|---|---|---|
| `bias_corr` | ≈ 0 (by construction) | units/baseline mismatch |
| `ks_corr` | clearly below `ks_raw` | method too weak or index not comparable |
| `sd_ratio_corr` | 0.9–1.1 | `mean_std`/QDM not applied to the same quantity |
| ensemble trend p-value | small for strong scenarios | short period or few models |
| GEV `shape` | \|c\| < 1, stable across models | too few years, or index not a block maximum |
| GEV `ks_p` | > 0.05 | misfit: try Gumbel, longer period, or a different index |
| maps stippling | present where the signal is strong | too few models or threshold too strict |

Then look at every figure (open it, do not just check the file exists):
axis labels and units, overlapping text, colour-bar range, panel letters,
legend readability, and that the historical and future branches join.

---

## Step 7 — iterate

Common adjustments, in order of impact:

1. **Domain/region** — a too-large region dilutes the signal.
2. **Season** — the annual mean often hides the projected change.
3. **Index definition** — zonal anomaly / projection indices sharpen the signal.
4. **Correction method** — QDM vs mean_std changes tails, not means.
5. **Period** — 2026–2050 vs 2081–2100.
6. **Model selection** — exclude models with `ks_corr` far above the ensemble
   median, or with |c| ≫ 1.

Re-run only the affected stages (`--stages correct,ensemble,extremes,figures`)
after a plotting or GEV change; re-run everything after a data/index change.
Cached yearly fields in `<outdir>/cache/` make repeated runs fast; delete the
cache only if the source files changed.

---

## Step 8 — report to the user

State, in this order:

1. variable, units, season, region, index definition;
2. reference dataset and period; model list, member, scenarios;
3. correction method and baseline;
4. headline numbers from `report.md` (bias/KS improvement, ensemble trend with
   p-value, return level or return period of the recorded event);
5. caveats that matter (e.g. single member, stationary-bias assumption, QM vs
   QDM, index is not a block maximum);
6. absolute paths of `report.md`, the figures, and the tables.

Never present a figure without saying which method and period produced it.
