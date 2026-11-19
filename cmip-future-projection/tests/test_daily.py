"""Integration test for the daily-extremes stage (POT/GPD and block maxima).

Builds a small daily data set (3 models, 30 historical years + 25 future years,
JJA, 3x3 grid), runs the pipeline with ``daily.enabled: true`` and checks that
the GPD fits, the return levels and the figure come out.

    python tests/test_daily.py

Exits non-zero on failure.  Takes a few seconds.
"""

import os
import shutil
import sys
import tempfile

import numpy as np
import xarray as xr
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

from cmipkit import daily as daily_mod  # noqa: E402
from cmipkit.config import load_config  # noqa: E402
from cmipkit.pipeline import Pipeline  # noqa: E402


def _daily_field(rng, lat, lon, year, season_len=92, bias=0.0, trend=0.0, tail=1.0):
    lat2d, lon2d = np.meshgrid(lat, lon, indexing="ij")
    base = 28 + 4 * np.sin(np.deg2rad(lat2d))
    days = np.arange(season_len)
    seasonal = 3.0 * np.sin(np.pi * days / season_len)
    noise = rng.normal(0, 1.0, size=(season_len,) + base.shape)
    spikes = rng.gamma(2.0, 2.0 * tail, size=noise.shape) * (rng.random(noise.shape) < 0.12)
    return base[None] + seasonal[:, None, None] + bias + trend * (year - 1980) + noise + spikes


def build(root):
    rng = np.random.default_rng(11)
    lat = np.array([30.0, 32.0, 34.0])
    lon = np.array([100.0, 102.0, 104.0])
    times = []
    for year in list(range(1985, 2015)) + list(range(2026, 2051)):
        start = np.datetime64(f"{year}-06-01")
        times.extend(start + np.arange(92))
    times = np.asarray(times, dtype="datetime64[ns]")
    years = np.array([int(str(t)[:4]) for t in times])
    os.makedirs(root, exist_ok=True)

    def write(path, bias, trend, tail):
        values = np.concatenate([_daily_field(rng, lat, lon, y, bias=bias, trend=trend, tail=tail)
                                 for y in np.unique(years)], axis=0)
        da = xr.DataArray(values, dims=("time", "latitude", "longitude"),
                          coords={"time": times, "latitude": lat, "longitude": lon},
                          name="tas", attrs={"units": "degC"})
        da.to_netcdf(path)

    write(os.path.join(root, "obs.nc"), bias=0.0, trend=0.02, tail=1.0)
    for k, model in enumerate(["MODEL-A", "MODEL-B", "MODEL-C"]):
        write(os.path.join(root, f"{model}_historical.nc"), bias=1.5 + 0.5 * k, trend=0.015 + 0.003 * k,
              tail=1.3 + 0.2 * k)
        for scen, extra in (("ssp245", 0.03), ("ssp585", 0.06)):
            write(os.path.join(root, f"{model}_{scen}.nc"), bias=1.5 + 0.5 * k,
                  trend=0.015 + 0.003 * k + extra, tail=1.3 + 0.2 * k)


def config(root, kind):
    cfg = {
        "project": {"name": f"daily_{kind}", "outdir": os.path.join(root, f"run_{kind}")},
        "variable": {"name": "tas", "units": "degC", "long_name": "daily maximum temperature"},
        "data": {
            "reference": {
                "kind": "grid", "paths": [os.path.join(root, "obs.nc")], "years": [1985, 2014],
                "variable": "tas", "aggregate": {"months": [6, 7, 8], "stat": "max"},
                "bbox": {"lat": [30, 34], "lon": [100, 104]},
            },
            "models": {
                "kind": "grid", "root": root, "pattern": "{root}/{model}_{experiment}.nc",
                "variable": "tas", "historical_experiment": "historical",
                "scenarios": ["ssp245", "ssp585"], "models": "all",
                "bbox": {"lat": [30, 34], "lon": [100, 104]},
                "aggregate": {"months": [6, 7, 8], "stat": "max"},
            },
        },
        "grid": {"lat": {"start": 30, "stop": 34, "step": 2}, "lon": {"start": 100, "stop": 104, "step": 2}},
        "periods": {"historical": [1985, 2014], "future": [2026, 2050],
                    "reference_climate": [1985, 2014], "future_windows": [[2026, 2040], [2041, 2050]]},
        "index": {"type": "region_mean", "region": {"lat": [30, 34], "lon": [100, 104]}, "anomaly": False},
        "correction": {"method": "mean_std", "multiplicative": False},
        "extremes": {"enabled": True, "return_periods": [2, 10, 100]},
        "daily": {
            "enabled": True,
            "months": [6, 7, 8],
            "index": "region_mean",
            "correction": "quantile_delta_mapping",
            "extremes": {
                "kind": kind,
                "threshold_mode": "quantile",
                "threshold_quantile": 0.90,
                "decluster_days": 3,
                "min_peaks": 20,
                "min_years": 15,
                "return_periods": [2, 5, 10, 50, 100],
            },
        },
        "plots": {"dpi": 90, "figures": ["timeseries"]},
    }
    path = os.path.join(root, f"config_{kind}.yaml")
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, allow_unicode=True, sort_keys=False)
    return path


def test_unit():
    """Declustering + GPD sanity, no files."""
    rng = np.random.default_rng(0)
    values = rng.gamma(2.0, 1.0, size=3650)
    threshold = float(np.quantile(values, 0.9))
    idx, peaks = daily_mod.decluster(values, threshold, run_length=3)
    assert peaks.size > 0 and peaks.size < 500, "declustering should reduce the exceedances"
    assert np.all(peaks > threshold)
    fit = daily_mod.fit_pot(values, threshold, n_years=10.0, decluster_days=3, min_peaks=20)
    assert fit is not None and fit.scale > 0
    levels = fit.return_levels([2, 10, 100])
    assert np.all(np.diff(levels) > 0), f"return levels must increase: {levels}"
    assert np.isfinite(fit.return_period(levels[1]))
    dates = np.datetime64("2000-01-01") + np.arange(len(values))
    years, maxima = daily_mod.annual_maxima(dates, values)
    assert years.size == 10 and maxima.size == 10
    print(f"unit: {peaks.size} peaks, shape={fit.shape:+.3f}, scale={fit.scale:.3f}, "
          f"rate={fit.rate:.2f}/yr, RL(100)={levels[2]:.2f}")
    return True


def test_pipeline(kind):
    root = tempfile.mkdtemp(prefix=f"cmipkit_daily_{kind}_")
    try:
        build(root)
        cfg = load_config(config(root, kind))
        pipe = Pipeline(cfg)
        pipe.run()
        checks = {
            "daily stage": bool(pipe.daily),
            "historical fit": pipe.daily.get("historical") is not None,
            "two scenarios": set(pipe.daily.get("scenarios", {})) == {"ssp245", "ssp585"},
            "three models": all(v["n_models"] == 3 for v in pipe.daily.get("scenarios", {}).values()),
            "figure": os.path.exists(os.path.join(cfg.outdir, "figures", "fig_daily_return_levels.png")),
            "table": os.path.exists(os.path.join(cfg.outdir, "tables", "daily_extremes.csv")),
            "report": os.path.exists(os.path.join(cfg.outdir, "report.md")),
        }
        report = open(os.path.join(cfg.outdir, "report.md"), encoding="utf-8").read()
        checks["report section"] = "Daily extremes" in report
        med = pipe.daily["scenarios"]["ssp585"]["stats"]["median"]
        checks["levels increase"] = bool(np.all(np.diff(med[np.isfinite(med)]) > 0))
        failed = [k for k, v in checks.items() if not v]
        print(f"{kind}: " + ", ".join(f"{k}={'ok' if v else 'FAIL'}" for k, v in checks.items()))
        return not failed
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    ok = test_unit()
    ok = test_pipeline("pot") and ok
    ok = test_pipeline("block") and ok
    print("\nRESULT:", "daily tests passed" if ok else "FAILED")
    sys.exit(0 if ok else 1)
