"""Unit tests for the correction / index / extremes kernels.

Run with the climate environment:

    python tests/test_methods.py

No pytest required; the script exits non-zero on the first failure.
"""

import os
import sys

import numpy as np
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

from cmipkit import correction, ensemble, extremes, indices  # noqa: E402


def make_fields(seed=0, n=60, lat=None, lon=None, positive=False):
    """Synthetic (obs, model) fields; ``positive`` mimics precipitation."""
    rng = np.random.default_rng(seed)
    lat = np.arange(25.0, 50.1, 2.5) if lat is None else lat
    lon = np.arange(80.0, 130.1, 2.5) if lon is None else lon
    years = np.arange(1955, 1955 + n)
    if positive:
        base = np.broadcast_to(2.0 + 1.5 * np.sin(np.deg2rad(lat))[:, None], (lat.size, lon.size)).copy()
        obs = np.array([base * np.exp(rng.normal(0, 0.25, base.shape)) * (1 + 0.004 * (y - 1980)) for y in years])
        mod = np.array([1.6 * base * np.exp(rng.normal(0, 0.35, base.shape)) * (1 + 0.002 * (y - 1980)) for y in years])
        coords = {"year": years, "latitude": lat, "longitude": lon}
        dims = ("year", "latitude", "longitude")
        return xr.DataArray(obs, dims=dims, coords=coords), xr.DataArray(mod, dims=dims, coords=coords)
    base = np.broadcast_to(280 + 5 * np.sin(np.deg2rad(lat))[:, None], (lat.size, lon.size)).copy()
    obs = np.array([base + rng.normal(0, 1.0, base.shape) + 0.02 * (y - 1980) for y in years])
    mod = np.array([base + 3.0 + rng.normal(0, 1.4, base.shape) + 0.01 * (y - 1980) for y in years])
    coords = {"year": years, "latitude": lat, "longitude": lon}
    dims = ("year", "latitude", "longitude")
    return xr.DataArray(obs, dims=dims, coords=coords), xr.DataArray(mod, dims=dims, coords=coords)


def test_methods():
    ok = True
    for positive, mult in [(False, False), (True, True)]:
        obs, mod = make_fields(positive=positive)
        methods = ["mean", "mean_std", "quantile_mapping", "quantile_delta_mapping"]
        for method in methods:
            params = correction.fit(obs, mod, method=method, multiplicative=mult, n_quantiles=60)
            corr = correction.apply(params, mod, multiplicative=mult)
            diag = correction.diagnostics(obs.values.mean(axis=(1, 2)), mod.values.mean(axis=(1, 2)),
                                          corr.values.mean(axis=(1, 2)))
            rel = abs(diag["bias_corr"]) / max(abs(diag["bias_raw"]), 1e-9)
            improved = rel < 0.1
            print(f"{'precip' if positive else 'temp':6s} {method:24s} mult={mult!s:5s} "
                  f"bias {diag['bias_raw']:+.3f} -> {diag['bias_corr']:+.3f} "
                  f"KS {diag['ks_raw']:.3f} -> {diag['ks_corr']:.3f}  {'OK' if improved else 'FAIL'}")
            ok &= bool(improved)
    return ok


def test_qdm_preserves_trend():
    obs, mod = make_fields()
    params = correction.fit(obs, mod, method="quantile_delta_mapping", n_quantiles=60)
    years = mod.year.values.astype(float)
    fut = mod.copy()
    fut.values = fut.values + np.linspace(0, 6, fut.shape[0])[:, None, None]
    corr = correction.apply(params, fut)
    raw_slope = np.polyfit(years, fut.values.mean(axis=(1, 2)), 1)[0]
    cor_slope = np.polyfit(years, corr.values.mean(axis=(1, 2)), 1)[0]
    # Cannon's QDM preserves the model's quantile change, not the mean trend
    # exactly, so allow a modest drift here.
    print(f"QDM trend preservation: raw {raw_slope:+.4f} -> corr {cor_slope:+.4f} per year")
    return abs(cor_slope - raw_slope) < 0.25 * abs(raw_slope)


def test_indices():
    obs, mod = make_fields()
    region = {"lat": [28, 45], "lon": [85, 120]}
    spec = {"type": "region_mean", "region": region, "anomaly": False}
    y, v = indices.build_index(obs, spec, [1960, 1980])
    assert y.size == v.size and np.isfinite(v).all()
    trend = indices.linear_trend(y, v)
    print(f"region_mean trend {trend['slope_per_decade']:+.4f} per decade (p={trend['p_value']:.3g})")
    assert abs(trend["slope_per_decade"] - 0.2) < 0.15

    zonal = indices.remove_zonal_mean(obs)
    np.testing.assert_allclose(zonal.mean("longitude").values, 0, atol=1e-10)

    spec2 = {"type": "region_diff", "regions": [{"lat": [30, 40], "lon": [80, 95]},
                                                {"lat": [30, 40], "lon": [110, 125]}]}
    y2, v2 = indices.build_index(obs, spec2, [1960, 1980])
    print(f"region_diff range {np.nanmin(v2):+.2f} .. {np.nanmax(v2):+.2f}")
    return np.isfinite(v2).all()


def test_extremes():
    rng = np.random.default_rng(7)
    sample = 30 + 5 * rng.standard_gamma(2.0, 80)
    fit = extremes.fit_distribution(sample, "gev", "mle", 20)
    assert fit is not None
    rl = fit.return_levels([2, 10, 100])
    print(f"GEV c={fit.shape:+.3f} RL(2,10,100)={np.round(rl, 2)} ks_p={fit.ks_p:.3f}")
    assert rl[0] < rl[1] < rl[2]
    rp = fit.return_period(float(np.max(sample)))
    print(f"return period of the sample maximum: {rp:.1f} years")
    return np.isfinite(rp)


def test_ensemble():
    _, mod = make_fields()
    mat = np.vstack([mod.values.mean(axis=(1, 2)) + k for k in range(5)])
    st = ensemble.stats(mat, "mean", [10, 90])
    assert st["center"].shape == (mat.shape[1],)
    assert np.all(st["p10"] <= st["p90"])
    print(f"ensemble stats ok: center[0]={st['center'][0]:.2f}, p10={st['p10'][0]:.2f}, p90={st['p90'][0]:.2f}")
    return True


if __name__ == "__main__":
    checks = {
        "methods": test_methods(),
        "qdm_trend": test_qdm_preserves_trend(),
        "indices": test_indices(),
        "extremes": test_extremes(),
        "ensemble": test_ensemble(),
    }
    failed = [k for k, v in checks.items() if not v]
    print("\nRESULT:", "all tests passed" if not failed else f"FAILED: {failed}")
    sys.exit(1 if failed else 0)
