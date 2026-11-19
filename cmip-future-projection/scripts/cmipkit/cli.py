"""Command line interface.

Commands
--------
``run CONFIG [--stages ...]``
    Execute the pipeline.  ``--stages`` accepts a comma separated subset of
    ``prepare,correct,ensemble,extremes,figures,report``.
``init DIR``
    Write a fully commented ``config.yaml`` and a short README into ``DIR``.
``doctor``
    Check the python environment (packages, fonts, cartopy data, write access).
``selftest [--keep]``
    Build a synthetic data set and run the whole pipeline end to end.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import textwrap
from typing import List, Optional

TEMPLATE = """\
# cmipkit configuration -- every key shown here is optional except the paths.
project:
  name: my_projection
  outdir: runs/my_projection

variable:
  name: tas                 # label used in titles
  units: degC               # units AFTER units_transform
  long_name: Near-surface air temperature

data:
  reference:                # observations / reanalysis
    kind: grid
    paths:
      - /path/to/era5/tas_{year}.nc     # {year} is expanded, globs are allowed
    years: [1950, 2014]
    variable: t2m           # variable inside the file (omit to auto-detect)
    units_in: K
    units_transform: {scale: 1.0, offset: -273.15}
    aggregate: {months: [9, 10], stat: mean}
  models:
    kind: grid
    root: /path/to/cmip6/monthly
    pattern: "{root}/{experiment}/{model}/{variable}_Amon_{model}_{experiment}_{member}_*.nc"
    variable: tas
    units_in: K
    units_transform: {scale: 1.0, offset: -273.15}
    historical_experiment: historical
    scenarios: [ssp245, ssp585]
    models: all             # or an explicit list
    aggregate: {months: [9, 10], stat: mean}

grid:                       # common grid every field is interpolated to
  lat: {start: 20, stop: 60, step: 1}
  lon: {start: 60, stop: 140, step: 1}

periods:
  historical: [1950, 2014]
  future: [2026, 2100]
  reference_climate: [1950, 2014]
  future_windows: [[2026, 2045], [2046, 2065], [2066, 2085], [2086, 2100]]

index:                      # what one number per year means
  type: region_mean         # region_mean | region_mean_anomaly | region_diff | projection_index | series
  region: {lat: [25, 55], lon: [70, 130]}
  anomaly: true

correction:
  method: mean_std          # none | mean | mean_std | quantile_mapping | quantile_delta_mapping
  params_period: [1950, 2014]
  lower_bound: null         # e.g. 0 for precipitation

ensemble:
  statistic: mean
  bands: [10, 90, 20, 80]

extremes:
  enabled: true
  distribution: gev
  return_periods: [2, 5, 10, 20, 50, 100, 200]

plots:
  dpi: 300
  figures: [timeseries, gev, maps, kde, shape_box]
"""


def cmd_init(args) -> int:
    os.makedirs(args.dir, exist_ok=True)
    cfg_path = os.path.join(args.dir, "config.yaml")
    if os.path.exists(cfg_path) and not args.force:
        print(f"refusing to overwrite {cfg_path} (use --force)")
        return 2
    with open(cfg_path, "w", encoding="utf-8") as fh:
        fh.write(TEMPLATE)
    readme = os.path.join(args.dir, "README.md")
    with open(readme, "w", encoding="utf-8") as fh:
        fh.write(
            "# cmipkit run\n\n"
            "1. Edit `config.yaml` (paths, variable, season, region, correction method).\n"
            "2. Check the environment: `python <skill>/scripts/run.py doctor`\n"
            "3. Run everything: `python <skill>/scripts/run.py run config.yaml`\n"
            "4. Read `runs/<name>/report.md` and look at `runs/<name>/figures/`.\n"
        )
    print(f"wrote {cfg_path}")
    print(f"wrote {readme}")
    return 0


def cmd_doctor(args) -> int:
    import importlib

    print(f"python      : {sys.version.split()[0]}  ({sys.executable})")
    print(f"platform    : {platform.platform()}")
    ok = True
    for name in ("numpy", "scipy", "xarray", "matplotlib", "yaml", "pandas", "netCDF4"):
        try:
            mod = importlib.import_module(name)
            print(f"{name:<12}: {getattr(mod, '__version__', 'ok')}")
        except Exception as exc:
            ok = False
            print(f"{name:<12}: MISSING ({exc})")
    optional = {}
    for name in ("cartopy", "cfgrib", "cftime", "dask"):
        try:
            mod = importlib.import_module(name)
            optional[name] = getattr(mod, "__version__", "ok")
            print(f"{name:<12}: {optional[name]} (optional)")
        except Exception:
            print(f"{name:<12}: not installed (optional)")
    # fonts
    import matplotlib.font_manager as fm

    names = {f.name for f in fm.fontManager.ttflist}
    cjk = [n for n in ("WenQuanYi Micro Hei", "WenQuanYi Zen Hei", "Noto Sans CJK SC", "Source Han Sans SC") if n in names]
    print(f"cjk fonts   : {cjk if cjk else 'none found (Chinese labels will fall back to boxes)'}")
    # cartopy data
    if optional.get("cartopy"):
        import cartopy

        data_dir = cartopy.config.get("data_dir", "")
        print(f"cartopy data: {data_dir} ({'present' if os.path.isdir(str(data_dir)) else 'missing'})")
    # writability
    import tempfile

    try:
        with tempfile.NamedTemporaryFile(delete=True) as fh:
            fh.write(b"ok")
        print("write access: ok")
    except Exception as exc:
        ok = False
        print(f"write access: FAILED ({exc})")
    print("\nRESULT:", "ready" if ok else "missing required packages")
    return 0 if ok else 1


def cmd_selftest(args) -> int:
    from .selftest import run_selftest

    result = run_selftest(workdir=args.workdir, keep=args.keep)
    print(json.dumps({k: v for k, v in result.items() if k != "checks"}, indent=2, ensure_ascii=False))
    if not result["ok"]:
        print("\nfailed checks:", ", ".join(result["failed"]))
        return 1
    print("\nselftest passed:", len(result["figures"]), "figures,",
          f"bias {result['checks']['raw_bias']:.3f} -> {result['checks']['corr_bias']:.3f}")
    return 0


def cmd_run(args) -> int:
    from .config import ConfigError, load_config
    from .pipeline import Pipeline

    try:
        cfg = load_config(args.config)
    except ConfigError as exc:
        print(f"configuration error:\n{exc}", file=sys.stderr)
        return 2
    stages = args.stages.split(",") if args.stages else None
    pipe = Pipeline(cfg)
    pipe.run(stages)
    print(f"\noutput: {cfg.outdir}")
    print(f"report: {os.path.join(cfg.outdir, 'report.md')}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="cmipkit",
        description="Bias-corrected CMIP6 future projections, extremes and figures.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run the pipeline from a YAML config")
    r.add_argument("config")
    r.add_argument("--stages", default=None,
                   help="comma separated subset of prepare,correct,ensemble,extremes,figures,report")
    r.set_defaults(func=cmd_run)

    i = sub.add_parser("init", help="write a commented config template")
    i.add_argument("dir")
    i.add_argument("--force", action="store_true")
    i.set_defaults(func=cmd_init)

    d = sub.add_parser("doctor", help="check the environment")
    d.set_defaults(func=cmd_doctor)

    s = sub.add_parser("selftest", help="synthetic end-to-end test")
    s.add_argument("--workdir", default=None)
    s.add_argument("--keep", action="store_true")
    s.set_defaults(func=cmd_selftest)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
