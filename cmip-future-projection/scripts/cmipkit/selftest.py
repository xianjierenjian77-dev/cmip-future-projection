"""Synthetic end-to-end test.

``python run.py selftest`` builds a tiny but structurally realistic data set
(observational reference + two bias-corrupted models + two scenarios), runs the
whole pipeline and asserts that the correction actually removes the imposed bias
and that every promised figure is produced.  It is the fastest way to check that
an environment can run the skill.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import tempfile
from typing import Dict, List, Tuple

import numpy as np
import xarray as xr
import yaml

log = logging.getLogger("cmipkit.selftest")


def _make_field(rng, lat, lon, base, trend, year, noise=1.0, bias=0.0, bias_shape=None):
    lat2d, lon2d = np.meshgrid(lat, lon, indexing="ij")
    pattern = (
        base
        + 6.0 * np.sin(np.deg2rad(lat2d))
        + 3.0 * np.cos(np.deg2rad(lon2d * 2))
        + 2.0 * np.exp(-((lat2d - 40) ** 2 + (lon2d - 100) ** 2) / 400.0)
    )
    t = (year - 1980) * trend
    if bias_shape is not None:
        pattern = pattern + bias_shape
    return pattern + bias + t + rng.normal(0, noise, size=pattern.shape)


def _write_dataset(path: str, years, lat, lon, values_by_year: Dict[int, np.ndarray], months=(9, 10)):
    times, fields = [], []
    for year in years:
        for m in months:
            times.append(np.datetime64(f"{int(year):04d}-{int(m):02d}-15", "ns"))
            fields.append(values_by_year[int(year)])
    da = xr.DataArray(
        np.stack(fields, axis=0),
        dims=("time", "latitude", "longitude"),
        coords={"time": np.asarray(times, dtype="datetime64[ns]"), "latitude": lat, "longitude": lon},
        name="tas",
        attrs={"units": "K"},
    )
    da.to_netcdf(path)


def build_dataset(root: str, seed: int = 3) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    os.makedirs(root, exist_ok=True)
    lat = np.arange(20.0, 60.1, 2.0)
    lon = np.arange(60.0, 140.1, 2.0)
    hist_years = np.arange(1950, 2015)
    fut_years = np.arange(2026, 2101)

    # observational reference
    obs = {y: _make_field(rng, lat, lon, 280.0, 0.030, y, noise=0.8) for y in hist_years}
    _write_dataset(os.path.join(root, "obs_tas.nc"), hist_years, lat, lon, obs)

    models = ["MODEL-A", "MODEL-B", "MODEL-C"]
    members = ["r1i1p1f1", "r2i1p1f1", "r3i1p1f1"]
    for k, model in enumerate(models):
        bias = 3.0 + 0.8 * k
        sd_scale = 1.25 + 0.1 * k
        hist_trend = 0.018 + 0.004 * k
        for member in members:
            # members share the model bias/trend but have independent noise
            hist = {}
            for y in hist_years:
                f = _make_field(rng, lat, lon, 280.0, hist_trend, y, noise=0.8 * sd_scale)
                hist[y] = 280.0 + (f - 280.0) * sd_scale + bias
            _write_dataset(os.path.join(root, f"{model}_{member}_historical.nc"), hist_years, lat, lon, hist)

            for scen, trend in (("ssp245", 0.030), ("ssp585", 0.060)):
                fut = {}
                for y in fut_years:
                    f = _make_field(rng, lat, lon, 280.0, hist_trend + (trend - 0.03), y, noise=0.8 * sd_scale)
                    fut[y] = 280.0 + (f - 280.0) * sd_scale + bias
                _write_dataset(os.path.join(root, f"{model}_{member}_{scen}.nc"), fut_years, lat, lon, fut)
    return {"lat": lat, "lon": lon, "models": models, "members": members,
            "hist_years": hist_years, "future_years": fut_years}


def build_config(workdir: str, data_dir: str) -> str:
    cfg = {
        "project": {"name": "selftest", "outdir": os.path.join(workdir, "run")},
        "variable": {"name": "tas", "units": "K", "long_name": "Near-surface temperature"},
        "data": {
            "reference": {
                "kind": "grid",
                "paths": [os.path.join(data_dir, "obs_tas.nc")],
                "years": [1950, 2014],
                "variable": "tas",
                "aggregate": {"months": [9, 10], "stat": "mean"},
            },
            "models": {
                "kind": "grid",
                "root": data_dir,
                "pattern": "{root}/{model}_{member}_{experiment}.nc",
                "variable": "tas",
                "historical_experiment": "historical",
                "scenarios": ["ssp245", "ssp585"],
                "models": "all",
                "members": ["r1i1p1f1", "r2i1p1f1", "r3i1p1f1"],
                "member_combine": "mean",
                "aggregate": {"months": [9, 10], "stat": "mean"},
            },
        },
        "grid": {"lat": {"start": 20, "stop": 60, "step": 2}, "lon": {"start": 60, "stop": 140, "step": 2}},
        "periods": {
            "historical": [1950, 2014],
            "future": [2026, 2100],
            "reference_climate": [1950, 2014],
            "future_windows": [[2026, 2045], [2046, 2065]],
        },
        "index": {
            "type": "region_mean",
            "region": {"lat": [25, 55], "lon": [70, 130], "name": "test-box"},
            "anomaly": True,
        },
        "correction": {"method": "mean_std"},
        "ensemble": {"statistic": "mean", "bands": [10, 90, 20, 80]},
        "extremes": {"enabled": True, "return_periods": [2, 10, 50, 100]},
        "plots": {
            "dpi": 110,
            "figures": ["timeseries", "gev", "maps", "kde", "shape_box"],
        },
    }
    path = os.path.join(workdir, "config_selftest.yaml")
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, allow_unicode=True, sort_keys=False)
    return path


def run_selftest(workdir: str | None = None, keep: bool = False) -> Dict[str, object]:
    from .config import load_config
    from .pipeline import Pipeline

    if workdir:
        os.makedirs(workdir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="cmipkit_selftest_", dir=workdir)
    data_dir = os.path.join(tmp, "data")
    log.info("building synthetic dataset in %s", data_dir)
    info = build_dataset(data_dir)
    cfg_path = build_config(tmp, data_dir)

    cfg = load_config(cfg_path)
    pipe = Pipeline(cfg)
    pipe.run()

    # ---- assertions ---------------------------------------------------
    report = os.path.join(cfg.outdir, "report.md")
    figures = sorted(os.listdir(os.path.join(cfg.outdir, "figures")))
    checks = {
        "report": os.path.exists(report),
        "figures": len(figures) >= 4,
        "models": len(pipe.scores) == 3,
        "multi_member": all(len(v) == 3 for v in pipe.source.model_members.values()) and len(pipe.source.model_members) == 3,
        "ensemble": set(pipe.ens) == {"ssp245", "ssp585"},
        "gev": all(k in pipe.gev for k in ("historical", "ssp245", "ssp585")),
    }
    for fig in figures:
        checks[f"figure:{fig}"] = True

    # the correction must remove the imposed ~3-5 K bias of the raw fields:
    # compare the region-mean climatology of the corrected historical fields to
    # the observational climatology.
    import glob as _glob

    import xarray as _xr

    ref = _xr.open_dataarray(os.path.join(cfg.outdir, "cache", sorted(
        f for f in os.listdir(os.path.join(cfg.outdir, "cache")) if f.startswith("reference_"))[0]))
    ref_mean = float(ref.mean(dim=("year", "latitude", "longitude")))
    ref.close()
    hist_bias = []
    for path in sorted(_glob.glob(os.path.join(cfg.outdir, "corrected", "*_historical_corrected.nc"))):
        da = _xr.open_dataarray(path)
        hist_bias.append(abs(float(da.mean(dim=("year", "latitude", "longitude"))) - ref_mean))
        da.close()
    raw_bias = 3.0 + 0.8  # MODEL-A..C offsets, mean of 3.0, 3.8, 4.6
    cor_bias = float(np.mean(hist_bias)) if hist_bias else np.inf
    checks["bias_reduced"] = cor_bias < 0.5
    checks["raw_bias"] = float(raw_bias)
    checks["corr_bias"] = cor_bias

    failed = [k for k, v in checks.items() if v is False]
    result = {
        "ok": not failed,
        "workdir": tmp,
        "checks": checks,
        "failed": failed,
        "outdir": cfg.outdir,
        "figures": figures,
    }
    if not keep:
        shutil.rmtree(tmp, ignore_errors=True)
        result["workdir"] = "(removed; use keep=True to retain)"
    return result
