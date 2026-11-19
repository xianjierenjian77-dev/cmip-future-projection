"""Read-only config and file-discovery checks; no analysis output is created."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict

from .config import Config, validate_config


def inspect_config(cfg: Config, *, check_inputs: bool = False) -> Dict[str, Any]:
    errors, warnings = validate_config(cfg)
    result = {
        "ok": not errors, "config": cfg.get_path("_config_path"), "outdir": cfg.outdir,
        "inputs_checked": check_inputs, "models": [], "datasets": [],
        "errors": list(errors), "warnings": list(warnings),
    }
    if errors or not check_inputs:
        return result
    from . import io
    source = io.DataSource(cfg, create_dirs=False)

    def add(label, files, spec):
        files = sorted(set(files))
        item = {"label": label, "files": files}
        result["datasets"].append(item)
        if not files:
            result["errors"].append(f"{label}: no files matched the configured paths/member")
            return
        item["files_metadata"] = []
        for filename in files:
            p = Path(filename)
            if not p.is_file():
                result["errors"].append(f"{label}: file not found: {filename}")
            else:
                stat = p.stat()
                item["files_metadata"].append({
                    "path": str(p), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                })
        if len(item["files_metadata"]) != len(files):
            return
        try:
            if spec.get("kind") == "series":
                years, values = io.read_series_file(files[0], spec.get("series_columns"))
                item["sample_header"] = {"n_years": int(years.size),
                                         "first_year": int(years[0]), "last_year": int(years[-1])}
                return
            # Read only one representative header per dataset; time coverage and
            # every file's values are validated in the loader, not claimed here.
            import xarray as xr
            filename = files[0]
            suffix = Path(filename).suffix.lower()
            if suffix in (".tif", ".tiff"):
                raise ValueError("GeoTIFF is not supported: convert with documented calendar/unit metadata first")
            kwargs = {"decode_times": False}
            if suffix in (".grib", ".grb", ".grib2", ".grb2"):
                kwargs.update(engine="cfgrib", backend_kwargs={"indexpath": ""})
            with xr.open_dataset(filename, **kwargs) as ds:
                var = spec.get("variable")
                if var and var not in ds.data_vars:
                    raise ValueError(f"variable {var!r} absent; available: {list(ds.data_vars)}")
                if not var:
                    variables = [name for name in ds.data_vars if ds[name].ndim >= 2]
                    if len(variables) != 1:
                        raise ValueError(f"ambiguous data variable; explicitly choose from {list(ds.data_vars)}")
                    var = variables[0]
                data = ds[var]
                item["sample_header"] = {
                    "variable": var, "dims": dict(data.sizes),
                    "units": str(data.attrs.get("units", "")),
                    "calendar": str(ds["time"].attrs.get("calendar", "standard")) if "time" in ds else None,
                }
                declared = spec.get("units_in")
                found = data.attrs.get("units")
                if declared and found and str(declared).strip() != str(found).strip():
                    result["warnings"].append(
                        f"{label}: configured units_in={declared!r}, header units={found!r}; "
                        "verify aliases and units_transform explicitly"
                    )
                for coord in ("lat", "latitude", "lon", "longitude", "nav_lat", "nav_lon"):
                    if coord in ds.coords and ds[coord].ndim != 1:
                        raise ValueError("curvilinear grids need an external, validated regridding step")
        except Exception as exc:
            result["errors"].append(f"{label}: {type(exc).__name__}: {exc}")

    ref = cfg.get_path("data.reference")
    if ref["kind"] == "series":
        add("reference", [ref["series_file"]], ref)
    else:
        files = []
        for template in ref["paths"]:
            period = ref.get("years") or cfg.period("historical")
            files.extend(io.expand_paths(template, io._input_years(period, cfg.months("reference"))))
        add("reference", files, ref)

    try:
        models = source.discover_models()
    except Exception as exc:
        models = []
        result["errors"].append(f"model discovery: {exc}")
    result["models"] = models
    if not models:
        result["errors"].append("no models found: check root, patterns, models, members and historical experiment")
    spec = cfg.get_path("data.models")
    hist = spec.get("historical_experiment", "historical")
    experiments = [hist] + cfg.scenario_list()
    for model in models:
        locked_members = None
        for experiment in experiments:
            years = cfg.period("historical" if experiment == hist else "future")
            if spec["kind"] == "series":
                root = Path(spec["series_dir"])
                options = ([root / model / f"{experiment}{suffix}" for suffix in (".csv", ".npz")]
                           + [root / f"{model}_{experiment}{suffix}" for suffix in (".csv", ".npz")])
                files = [str(p) for p in options if p.is_file()]
                if len(files) > 1:
                    result["errors"].append(f"{model}/{experiment}: multiple series files match; select one layout unambiguously")
                add(f"{model}/{experiment}", files[:1], spec)
            else:
                try:
                    if spec.get("member_combine") != "first":
                        members = spec.get("members") or ["*"]
                        for member in members:
                            local_spec = dict(spec, members=[member])
                            files = io._files_for_grid_source(local_spec, experiment, model,
                                                              io._input_years(years, cfg.months("models")))
                            add(f"{model}/{experiment}/{member}", files, spec)
                    else:
                        selected = dict(spec)
                        if locked_members:
                            selected["members"] = locked_members
                        files = io._files_for_grid_source(selected, experiment, model,
                                                          io._input_years(years, cfg.months("models")))
                        if experiment == hist:
                            members = {io._member_of(path, model, experiment, spec) for path in files}
                            members.discard(None)
                            locked_members = sorted(members) if members else None
                        add(f"{model}/{experiment}", files, spec)
                except Exception as exc:
                    result["errors"].append(f"{model}/{experiment}: {exc}")
    result["ok"] = not result["errors"]
    return result
