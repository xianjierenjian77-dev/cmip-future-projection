# cmip-future-projection — installation and quick start

A config-driven pipeline that turns raw observational + CMIP6 (or any model
ensemble) files into a complete future-projection study: bias-corrected fields,
index series, ensemble statistics, GEV / POT extremes, publication figures,
CSV tables and a written `report.md`. Any variable, region, season or index.

This archive is a **skill** (a folder with `SKILL.md` plus resources). Read
`SKILL.md` first; it is the operating manual. `references/` has the equations,
the full config dictionary, per-variable recipes and troubleshooting;
`templates/` has ready-to-edit configs.

## What is inside

```
cmip-future-projection/
├── SKILL.md              operating manual (start here)
├── README.md             this file
├── references/           methods, config reference, plotting, variables, workflow, troubleshooting
├── templates/            5 ready-to-edit YAML configs
├── scripts/
│   ├── run.py            CLI entry point
│   └── cmipkit/          the package (~5.5k lines)
├── tests/
│   ├── test_methods.py   unit tests for the correction / index / GEV kernels
│   ├── test_daily.py     POT/GPD + daily block maxima integration test
│   └── preview/          example figures produced from real CMIP6 + ERA5 data
```

## Install

The skill is a folder. Put it in any skill root that your agent tool scans:

| Target | Directory |
|---|---|
| This server (already installed) | `/disk1/wsr/.dsh/skills/cmip-future-projection` |
| DSH, user level | `~/.dsh/skills/cmip-future-projection` (Windows: `%USERPROFILE%\.dsh\skills\`) |
| DSH, project level | `<project>/.dsh/skills/cmip-future-projection` |
| Shared with other agent tools | `~/.agents/skills/` or `<project>/.agents/skills/` |

Requirements: Python ≥ 3.9 with `numpy scipy xarray pandas matplotlib pyyaml`
(required) and `netCDF4 cftime cartopy cfgrib` (optional, needed for NetCDF,
CMIP6 calendars, maps and GRIB). On this server use:

```
/disk1/wsr/miniconda3/envs/climate/bin/python
```

## Quick start

```bash
PY=/disk1/wsr/miniconda3/envs/climate/bin/python
SK=/disk1/wsr/.dsh/skills/cmip-future-projection

$PY $SK/scripts/run.py doctor          # check the environment
$PY $SK/scripts/run.py selftest        # synthetic end-to-end test (~15 s)
$PY $SK/scripts/run.py init myrun      # write myrun/config.yaml
# edit myrun/config.yaml, then
$PY $SK/scripts/run.py run myrun/config.yaml
```

Outputs go to the `project.outdir` of the config: `report.md`, `figures/`,
`tables/`, `corrected/`, `series/`, `ensemble/`, `extremes/`, `params/`.

## Paths in the templates

The templates and `SKILL.md` use this server's data layout (`/disk3/era5-*`,
`/disk5/cmip6`, `/disk1/wsr/...`). On another machine, replace the paths and the
Python executable; nothing else is server-specific.

## Tests

```bash
$PY $SK/scripts/run.py selftest        # full pipeline on synthetic data
$PY $SK/tests/test_methods.py          # correction / index / GEV kernels
$PY $SK/tests/test_daily.py            # daily POT and block-maxima paths
```

All three must pass after any change to the package.
