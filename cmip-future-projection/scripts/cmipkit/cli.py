"""Small, portable command line interface; heavy dependencies load on demand."""
from __future__ import annotations

import argparse
import importlib
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
from typing import List, Optional

TEMPLATES = {
    "grid": "config_grid_projection.yaml",
    "series": "config_series.yaml",
    "precip": "config_precip_qdm.yaml",
    "regional-extremes": "config_regional_extremes.yaml",
    "daily": "config_daily_extremes.yaml",
}


def _positive_int(value: str) -> int:
    try:
        number = int(value)
    except (ValueError, TypeError) as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def cmd_init(args) -> int:
    source = Path(__file__).resolve().parents[2] / "templates" / TEMPLATES[args.template]
    destination = Path(args.dir).expanduser().resolve()
    target = destination / "config.yaml"
    if target.exists() and not args.force:
        print(f"refusing to overwrite {target}; choose a new directory or use --force", file=sys.stderr)
        return 2
    destination.mkdir(parents=True, exist_ok=True)
    target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    print(f"wrote {target}")
    print("Edit the paths and scientific choices, then validate CONFIG --check-inputs.")
    print("Relative paths are resolved against the configuration file's directory.")
    return 0


def cmd_doctor(args) -> int:
    print(f"python      : {sys.version.split()[0]} ({sys.executable})")
    print(f"platform    : {platform.platform()}")
    ok = sys.version_info >= (3, 9)
    available = {}
    for name in ("numpy", "scipy", "xarray", "pandas", "matplotlib", "yaml"):
        try:
            mod = importlib.import_module(name)
            available[name] = mod
            print(f"{name:<12}: {getattr(mod, '__version__', 'available')}")
        except Exception as exc:
            ok = False
            print(f"{name:<12}: MISSING ({type(exc).__name__}: {exc})")
    for name in ("netCDF4", "h5netcdf", "cftime", "cartopy", "cfgrib", "dask"):
        try:
            mod = importlib.import_module(name)
            available[name] = mod
            print(f"{name:<12}: {getattr(mod, '__version__', 'available')} (format/feature dependent)")
        except Exception:
            print(f"{name:<12}: not installed (format/feature dependent)")
    try:
        mod = importlib.import_module("threadpoolctl")
        print(f"threadpoolctl: {getattr(mod, '__version__', 'available')} (model multiprocessing)")
    except Exception:
        print("threadpoolctl: not installed; workers > 1 needs: python -m pip install threadpoolctl")
    if not any(name in available for name in ("netCDF4", "h5netcdf")):
        print("NetCDF-4/HDF5 reading needs netCDF4 or h5netcdf; scipy supports NetCDF-3.")
    if "cftime" not in available:
        print("Non-standard CMIP calendars need cftime; do not coerce them to Gregorian dates.")
    if "matplotlib" in available:
        import matplotlib.font_manager as fm
        names = {f.name for f in fm.fontManager.ttflist}
        cjk = [n for n in ("WenQuanYi Micro Hei", "Noto Sans CJK SC", "Source Han Sans SC") if n in names]
        print(f"CJK fonts   : {cjk or 'none detected; use English labels or configure a font'}")
    if "cartopy" in available:
        print(f"cartopy data: {available['cartopy'].config.get('data_dir')}")
        print("Coastlines need an existing Natural Earth cache; offline maps can use plain axes.")
    try:
        with tempfile.TemporaryFile() as fh:
            fh.write(b"ok")
        print("temporary IO: available")
    except OSError as exc:
        ok = False
        print(f"temporary IO: {exc}")
    print("RESULT:", "core dependencies available" if ok else "missing core dependencies")
    return 0 if ok else 1


def _validation_result(args, check_inputs=False):
    from .config import ConfigError, load_config, plan_stages, quality_requires_values
    from .preflight import inspect_config
    try:
        cfg = load_config(args.config)
        for name in ("workers", "worker_threads"):
            override = getattr(args, name, None)
            if override is not None:
                cfg.set_path(f"runtime.{name}", override)
    except ConfigError as exc:
        return None, {"report_type": "cmipkit.input_quality", "ok": False, "config": str(Path(args.config).resolve()),
                      "outdir": None, "inputs_checked": False, "models": [], "datasets": [],
                      "errors": [str(exc)], "warnings": [], "stages": [],
                      "quality": {"level": "config", "data_values_checked": False,
                                  "summary": {"datasets": 0, "files": 0, "errors": 1,
                                              "warnings": 0, "values_checked": 0}, "datasets": []}}
    requested = args.stages.split(",") if getattr(args, "stages", None) else None
    check_data = bool(getattr(args, "check_data", False))
    if getattr(args, "command", None) == "run":
        check_data = check_data or quality_requires_values(cfg)
    result = inspect_config(cfg, check_inputs=check_inputs or check_data, check_data=check_data)
    result["report_type"] = "cmipkit.input_quality"
    result["parallel"] = {
        "requested_workers": cfg.get_path("runtime.workers", 1),
        "worker_threads": cfg.get_path("runtime.worker_threads", 1),
        "effective_workers": (min(cfg.get_path("runtime.workers", 1), len(result["models"]))
                              if result.get("models") else None),
    }
    try:
        result["stages"] = plan_stages(cfg, requested)
    except ConfigError as exc:
        result["stages"] = []
        result["errors"].append(str(exc))
        result["ok"] = False
    return cfg, result


def _show_validation(result, as_json=False):
    if as_json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return
    print("Configuration:", result["config"])
    print("Output:", result["outdir"])
    print("Stages:", " -> ".join(result["stages"]) or "(none enabled)")
    if result.get("parallel"):
        settings = result["parallel"]
        effective = settings["effective_workers"]
        print("Model workers:", settings["requested_workers"], "requested;",
              effective if effective is not None else "undetermined", "planned after model discovery;",
              settings["worker_threads"], "numerical thread(s) per child")
    if result.get("inputs_checked"):
        print("Models found:", ", ".join(result.get("models", [])) or "(none)")
        for item in result.get("datasets", []):
            print(f"  {item['label']}: {len(item['files'])} file(s)")
    for warning in result.get("warnings", [])[:20]:
        print("WARNING:", warning)
    for error in result.get("errors", [])[:20]:
        print("ERROR:", error, file=sys.stderr)
    if len(result.get("warnings", [])) > 20 or len(result.get("errors", [])) > 20:
        print("Additional findings omitted from console; use --json or --report-json for the full report.")
    print("RESULT:", "preflight passed" if result["ok"] else "preflight failed")
    if not result.get("inputs_checked"):
        print("Input existence and file headers were not inspected; add --check-inputs.")
    elif result.get("quality", {}).get("data_values_checked"):
        print("All selected file headers, time coverage and selected data values checked in bounded chunks.")
    elif result.get("quality", {}).get("level") == "full":
        print("Full-value scan requested but not completed for every selected input; see the reported errors.")
    else:
        print("All selected file headers and time coordinates checked; values not scanned. Add --check-data.")


def cmd_validate(args) -> int:
    cfg, result = _validation_result(args, check_inputs=args.check_inputs)
    _save_quality_report(args, cfg, result)
    _show_validation(result, args.json)
    return 0 if result["ok"] else 2


def _save_quality_report(args, cfg, result):
    destination = getattr(args, "report_json", None)
    if not destination:
        return
    from .config import ConfigError, write_json_atomic
    target = Path(destination).expanduser().resolve()
    protected = {Path(args.config).resolve()}
    protected.update(Path(filename).resolve() for item in result.get("datasets", [])
                     for filename in item.get("files", []))
    if cfg is not None:
        for group in ("reference", "models"):
            source = cfg.get_path(f"data.{group}", {})
            if source.get("series_file"):
                protected.add(Path(source["series_file"]).resolve())
            for raw_path in source.get("paths", []):
                # Literal inputs must be protected even without discovery.
                if not any(char in raw_path for char in "*?[]{}"):
                    protected.add(Path(raw_path).resolve())
    input_suffixes = {".nc", ".nc4", ".cdf", ".grib", ".grb", ".grib2", ".grb2",
                      ".csv", ".npz", ".npy", ".yaml", ".yml"}
    if target in protected or target.suffix.lower() in input_suffixes:
        raise ConfigError("quality report must not overwrite a configuration or an input data file")
    if target.exists():
        try:
            previous = json.loads(target.read_text(encoding="utf-8"))
        except (OSError, ValueError, UnicodeError) as exc:
            raise ConfigError("quality report must not overwrite an existing non-report file; choose a new JSON path") from exc
        if not isinstance(previous, dict) or previous.get("report_type") != "cmipkit.input_quality":
            raise ConfigError("quality report must not overwrite an existing non-report file; choose a new JSON path")
    write_json_atomic(str(target), result)


def cmd_selftest(args) -> int:
    from .selftest import run_selftest
    result = run_selftest(workdir=args.workdir, keep=args.keep)
    print(json.dumps({k: v for k, v in result.items() if k != "checks"}, indent=2, ensure_ascii=False))
    if not result["ok"]:
        print("failed checks:", ", ".join(result["failed"]), file=sys.stderr)
        return 1
    print(f"selftest passed: {len(result['figures'])} figure files, synthetic calibration bias "
          f"{result['checks']['raw_bias']:.3f} -> {result['checks']['corr_bias']:.3f}")
    return 0


def cmd_run(args) -> int:
    cfg, result = _validation_result(args, check_inputs=True)
    _save_quality_report(args, cfg, result)
    if args.dry_run:
        _show_validation(result, args.json)
        return 0 if result["ok"] else 2
    if not result["ok"]:
        _show_validation(result, args.json)
        return 2
    if args.overwrite:
        cfg.set_path("runtime.overwrite", True)
    from .pipeline import Pipeline
    pipe = Pipeline(cfg, preflight_result=result)
    pipe.input_inventory = result.get("datasets", [])
    finished = pipe.run(result["stages"])
    print(f"output: {cfg.outdir}")
    if finished.get("report"):
        print(f"report: {finished['report']}")
    else:
        print("Selected stages completed; no new report was requested.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cmipkit", description="Reproducible CMIP projections with optional bias adjustment and extremes."
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="execute a validated configuration")
    run.add_argument("config")
    run.add_argument("--stages", default=None,
                     help="comma-separated stages; required upstream stages are added automatically")
    run.add_argument("--dry-run", action="store_true",
                     help="inspect inputs and execution plan without analysis output (explicit --report-json still writes a report)")
    run.add_argument("--json", action="store_true", help="machine-readable dry-run summary")
    run.add_argument("--check-data", action="store_true", help="scan all selected raw values before analysis")
    run.add_argument("--report-json", default=None, metavar="PATH",
                     help="save the complete preflight report, including failures, to PATH")
    run.add_argument("--workers", type=_positive_int, default=None, metavar="N",
                     help="run up to N models in separate processes; overrides runtime.workers")
    run.add_argument("--worker-threads", type=_positive_int, default=None, metavar="N",
                     help="limit supported BLAS/OpenMP threads per child; overrides runtime.worker_threads")
    run.add_argument("--overwrite", action="store_true", help="allow an existing output directory with a changed configuration")
    run.set_defaults(func=cmd_run)
    init = sub.add_parser("init", help="copy a portable starter configuration")
    init.add_argument("dir")
    init.add_argument("--template", choices=list(TEMPLATES), default="grid")
    init.add_argument("--force", action="store_true", help="replace an existing config.yaml")
    init.set_defaults(func=cmd_init)
    validate = sub.add_parser("validate", help="validate settings without performing an analysis")
    validate.add_argument("config")
    validate.add_argument("--check-inputs", action="store_true", help="discover files and check every selected header/time axis")
    validate.add_argument("--check-data", action="store_true", help="also scan all selected values in bounded chunks")
    validate.add_argument("--report-json", default=None, metavar="PATH", help="save the complete quality report to PATH")
    validate.add_argument("--stages", default=None)
    validate.add_argument("--json", action="store_true")
    validate.set_defaults(func=cmd_validate)
    doctor = sub.add_parser("doctor", help="check dependencies without importing the pipeline")
    doctor.set_defaults(func=cmd_doctor)
    selftest = sub.add_parser("selftest", help="run a small synthetic integration example")
    selftest.add_argument("--workdir", default=None)
    selftest.add_argument("--keep", action="store_true")
    selftest.set_defaults(func=cmd_selftest)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except KeyboardInterrupt:
        print("Interrupted; completed outputs were left in place.", file=sys.stderr)
        return 130
    except ModuleNotFoundError as exc:
        print(f"missing dependency: {exc.name}; run doctor and install the required package", file=sys.stderr)
        return 2
    except Exception as exc:
        # The detailed traceback is retained in run.log for pipeline failures.
        from .config import ConfigError
        if isinstance(exc, (ConfigError, ValueError, OSError, RuntimeError)):
            print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
            return 2
        raise


if __name__ == "__main__":
    raise SystemExit(main())
