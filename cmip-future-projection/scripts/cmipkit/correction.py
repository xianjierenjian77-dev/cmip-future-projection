"""Bias correction of model fields against an observational reference.

Five methods are available.  They all reduce to a per-grid-point transfer
function fitted on a common baseline period and then applied unchanged to the
historical and future fields:

``none``
    no correction (useful as a control run).
``mean`` / ``delta_change``
    first-moment correction.  Additive: ``x + (mu_o - mu_m)``; multiplicative:
    ``x * (mu_o / mu_m)``.  ``delta_change`` is the impact-modelling name for
    exactly the same computation, kept as an alias because papers call it that.
``mean_std``
    first + second moment.  Additive (the method used in the reference 6.25 /
    5.19 workflow): ``(x - mu_m) * sd_o/sd_m + mu_o``.  Multiplicative (for
    positive, ratio-scale variables): anchor the mean and the coefficient of
    variation, ``x' = x * mu_o/mu_m`` then
    ``x'' = mu_o + (x' - mu_o) * (CV_o / CV_m)``.
``quantile_mapping``
    empirical CDF mapping, per grid point: the corrected historical
    distribution matches the observed one by construction.  Best for variables
    whose full distribution matters (precipitation, wind speed).
``quantile_delta_mapping``
    trend-preserving variant: the model's *change* at each quantile is added to
    the observed quantile, so the future distribution keeps the model's own
    climate-change signal instead of being forced onto today's distribution.
    Recommended when the corrected future field feeds an extremes/GEV analysis.

All arithmetic is vectorised over the grid; only the empirical quantile
mapping loops over grid points (one ``numpy.interp`` call per point).
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, Optional, Tuple

import numpy as np
import xarray as xr

from .config import Config

log = logging.getLogger("cmipkit.correction")

METHODS = ("none", "mean", "mean_std", "delta_change", "quantile_mapping", "quantile_delta_mapping")
_PRECIP_KEYWORDS = ("pr", "precip", "rain", "snow", "convective", "shower", "runoff", "discharge")
_POSITIVE_KEYWORDS = ("wind", "speed", "evap", "radiation", "rsds", "rlds", "swdown", "par")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def resolve_multiplicative(cfg: Config, obs: np.ndarray, mod: np.ndarray) -> bool:
    """Decide whether the correction should be multiplicative."""
    mode = cfg.get_path("correction.multiplicative", "auto")
    if mode is True or mode is False:
        return bool(mode)
    label = " ".join(
        str(cfg.get_path(k, "") or "") for k in ("variable.name", "variable.long_name", "data.models.variable")
    ).lower()
    data_min = float(np.nanmin(np.concatenate([obs.ravel(), mod.ravel()]))) if obs.size and mod.size else np.nan
    non_negative = np.isfinite(data_min) and data_min >= 0
    if not non_negative:
        return False
    if any(k in label for k in _PRECIP_KEYWORDS):
        return True
    return False


def _masked_stats(arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Mean and sample std over axis 0, ignoring NaN."""
    mu = np.nanmean(arr, axis=0)
    sd = np.nanstd(arr, axis=0, ddof=1)
    return mu, sd


def _safe_ratio(num: np.ndarray, den: np.ndarray, fallback: float = 1.0) -> np.ndarray:
    out = np.full(np.broadcast(num, den).shape, fallback, dtype=np.float64)
    ok = np.isfinite(num) & np.isfinite(den) & (np.abs(den) > 1e-12)
    out[ok] = num[ok] / den[ok]
    bad = ~np.isfinite(out)
    out[bad] = fallback
    return out


def _empirical_quantiles(arr: np.ndarray, n_quantiles: int) -> np.ndarray:
    """Quantile table of ``arr`` over axis 0 -> shape (n_quantiles, n_points)."""
    probs = np.linspace(0.0, 1.0, int(n_quantiles))
    with np.errstate(invalid="ignore"):
        q = np.nanquantile(arr, probs, axis=0)
    # fill fully-NaN columns with zeros so downstream interp does not explode
    bad = ~np.isfinite(q).any(axis=0)
    if bad.any():
        q[:, bad] = 0.0
    return q


def _interp_columns(values: np.ndarray, x_table: np.ndarray, y_table: np.ndarray) -> np.ndarray:
    """Per-column linear interpolation of ``values`` from x_table -> y_table."""
    nt, npts = values.shape
    out = np.full((nt, npts), np.nan, dtype=np.float64)
    for j in range(npts):
        xs = x_table[:, j]
        ys = y_table[:, j]
        ok = np.isfinite(xs) & np.isfinite(ys)
        if ok.sum() < 2:
            out[:, j] = values[:, j]
            continue
        xs, ys = xs[ok], ys[ok]
        if np.allclose(xs, xs[0]):
            out[:, j] = ys[0]
            continue
        out[:, j] = np.interp(values[:, j], xs, ys)
    return out


def _rank_within_columns(values: np.ndarray) -> np.ndarray:
    """Empirical CDF position (0..1) of every value within its own column."""
    nt, npts = values.shape
    out = np.full((nt, npts), np.nan)
    for j in range(npts):
        col = values[:, j]
        ok = np.isfinite(col)
        n = int(ok.sum())
        if n < 2:
            out[ok, j] = 0.5
            continue
        ref = np.sort(col[ok])
        # average rank for ties, mid-point plotting position
        lo = np.searchsorted(ref, col[ok], side="left")
        hi = np.searchsorted(ref, col[ok], side="right")
        out[ok, j] = (lo + hi + 1) / (2.0 * (n + 1.0))
    return out


# ---------------------------------------------------------------------------
# parameter fitting
# ---------------------------------------------------------------------------
def fit(
    obs_hist: xr.DataArray,
    mod_hist: xr.DataArray,
    method: str = "mean_std",
    multiplicative: bool = False,
    n_quantiles: int = 100,
    lower_bound: Optional[float] = None,
) -> xr.Dataset:
    """Fit the transfer function on a common baseline period."""
    method = str(method).lower()
    if method == "delta_change":
        method = "mean"
    if method not in METHODS:
        raise ValueError(f"unknown correction method {method!r}; choose from {METHODS}")
    if method == "none":
        return xr.Dataset(attrs={"method": "none", "multiplicative": int(multiplicative)})

    if "year" in obs_hist.dims and "year" in mod_hist.dims:
        common = np.intersect1d(obs_hist.year.values, mod_hist.year.values)
        if common.size < 5:
            raise ValueError("fewer than 5 common years between reference and model for bias correction")
        obs_hist = obs_hist.sel(year=common)
        mod_hist = mod_hist.sel(year=common)
    elif obs_hist.sizes.get("time", 0) != mod_hist.sizes.get("time", 0):
        raise ValueError("reference and model fields must share a time/year axis to fit the correction")

    obs = obs_hist.values.astype(np.float64)
    mod = mod_hist.values.astype(np.float64)
    mu_o, sd_o = _masked_stats(obs)
    mu_m, sd_m = _masked_stats(mod)

    coords = {"latitude": obs_hist.latitude, "longitude": obs_hist.longitude}
    ds = xr.Dataset(
        {
            "mu_obs": (("latitude", "longitude"), mu_o),
            "sd_obs": (("latitude", "longitude"), sd_o),
            "mu_mod": (("latitude", "longitude"), mu_m),
            "sd_mod": (("latitude", "longitude"), sd_m),
        },
        coords=coords,
        attrs={
            "method": method,
            "multiplicative": int(multiplicative),
            "n_years": int(obs.shape[0]),
            "lower_bound": -999.0 if lower_bound is None else float(lower_bound),
        },
    )

    if method == "mean_std":
        if multiplicative:
            # anchor both the mean and the coefficient of variation:
            #   x' = x * mu_o/mu_m ; x'' = mu_o + (x' - mu_o) * (CV_o / CV_m)
            # which makes the corrected baseline mean mu_o and sd sd_o exactly.
            ratio = _safe_ratio(mu_o, mu_m)
            cv_ratio = _safe_ratio(sd_o / np.where(np.abs(mu_o) > 1e-12, mu_o, np.nan),
                                   sd_m / np.where(np.abs(mu_m) > 1e-12, mu_m, np.nan))
            ds["factor"] = (("latitude", "longitude"), ratio)
            ds["cv_ratio"] = (("latitude", "longitude"), cv_ratio)
        else:
            ds["factor"] = (("latitude", "longitude"), _safe_ratio(sd_o, sd_m))
            ds["offset"] = (("latitude", "longitude"), mu_o - mu_m * ds["factor"].values)
    elif method == "mean":
        if multiplicative:
            ds["factor"] = (("latitude", "longitude"), _safe_ratio(mu_o, mu_m))
        else:
            ds["offset"] = (("latitude", "longitude"), mu_o - mu_m)
    elif method in ("quantile_mapping", "quantile_delta_mapping"):
        nq = max(10, int(n_quantiles))
        q_mod = _empirical_quantiles(mod, nq)
        q_obs = _empirical_quantiles(obs, nq)
        probs = np.linspace(0.0, 1.0, nq)
        ds["q_mod"] = (("quantile", "latitude", "longitude"), q_mod.reshape(nq, *mu_o.shape))
        ds["q_obs"] = (("quantile", "latitude", "longitude"), q_obs.reshape(nq, *mu_o.shape))
        ds = ds.assign_coords(quantile=probs)
    return ds


# ---------------------------------------------------------------------------
# application
# ---------------------------------------------------------------------------
def apply(
    params: xr.Dataset,
    field: xr.DataArray,
    multiplicative: Optional[bool] = None,
    lower_bound: Optional[float] = None,
    upper_bound: Optional[float] = None,
) -> xr.DataArray:
    """Apply a fitted transfer function to a field of shape (year, lat, lon)."""
    method = str(params.attrs.get("method", "none"))
    if method == "none":
        return field.copy()
    if params.attrs.get("multiplicative") is not None and multiplicative is None:
        multiplicative = bool(params.attrs["multiplicative"])
    multiplicative = bool(multiplicative)
    if lower_bound is None and params.attrs.get("lower_bound", -999.0) != -999.0:
        lower_bound = float(params.attrs["lower_bound"])

    values = field.values.astype(np.float64)
    nt = values.shape[0]
    flat = values.reshape(nt, -1)
    nlat = field.sizes["latitude"]
    nlon = field.sizes["longitude"]

    if method in ("mean", "mean_std"):
        if multiplicative and method == "mean_std":
            factor = params["factor"].values.astype(np.float64).reshape(-1)
            cvr = params["cv_ratio"].values.astype(np.float64).reshape(-1)
            mu_o = params["mu_obs"].values.astype(np.float64).reshape(-1)
            out = mu_o[None, :] + (flat * factor[None, :] - mu_o[None, :]) * cvr[None, :]
        elif multiplicative:
            factor = params["factor"].values.astype(np.float64).reshape(-1)
            out = flat * factor[None, :]
        else:
            offset = params["offset"].values.astype(np.float64).reshape(-1)
            if method == "mean_std":
                factor = params["factor"].values.astype(np.float64).reshape(-1)
                mu_m = params["mu_mod"].values.astype(np.float64).reshape(-1)
                out = (flat - mu_m[None, :]) * factor[None, :] + (offset + mu_m * factor)[None, :]
            else:
                out = flat + offset[None, :]
    elif method == "quantile_mapping":
        q_mod = params["q_mod"].values.astype(np.float64).reshape(params.sizes["quantile"], -1)
        q_obs = params["q_obs"].values.astype(np.float64).reshape(params.sizes["quantile"], -1)
        out = _interp_columns(flat, q_mod, q_obs)
    elif method == "quantile_delta_mapping":
        q_mod = params["q_mod"].values.astype(np.float64).reshape(params.sizes["quantile"], -1)
        q_obs = params["q_obs"].values.astype(np.float64).reshape(params.sizes["quantile"], -1)
        probs = np.asarray(params["quantile"].values, dtype=np.float64)
        rank = _rank_within_columns(flat)
        delta = np.empty_like(flat)
        for j in range(flat.shape[1]):
            if multiplicative:
                base = np.interp(rank[:, j], probs, q_mod[:, j])
                target = np.interp(rank[:, j], probs, q_obs[:, j])
                ratio = np.where(np.abs(base) > 1e-12, target / np.where(np.abs(base) > 1e-12, base, 1.0), 1.0)
                delta[:, j] = flat[:, j] * ratio
            else:
                base = np.interp(rank[:, j], probs, q_mod[:, j])
                target = np.interp(rank[:, j], probs, q_obs[:, j])
                delta[:, j] = flat[:, j] + (target - base)
        out = delta
    else:
        raise ValueError(f"cannot apply unknown method {method!r}")

    if lower_bound is not None:
        out = np.where(out < lower_bound, lower_bound, out)
    if upper_bound is not None:
        out = np.where(out > upper_bound, upper_bound, out)

    arr = out.reshape(nt, nlat, nlon)
    return xr.DataArray(
        arr,
        dims=("year", "latitude", "longitude"),
        coords={"year": field.year.values, "latitude": field.latitude, "longitude": field.longitude},
        name="field",
        attrs={"correction": method, "multiplicative": int(multiplicative)},
    )


# ---------------------------------------------------------------------------
# trend correction (series level, "trend scaling")
# ---------------------------------------------------------------------------
def trend_correct_series(
    years: np.ndarray,
    values: np.ndarray,
    slope_model_hist: float,
    slope_obs_hist: float,
    reference_year: float,
    mode: str = "difference",
) -> np.ndarray:
    """Rescale a future series so that its trend matches the observed trend.

    ``difference`` (default):  ``x - (b_mod - b_obs) * (t - t0)``
    ``ratio``:                 ``x * b_obs / b_mod`` (only sensible for
                               strictly positive series far from zero)
    """
    x = np.asarray(values, dtype=np.float64)
    yr = np.asarray(years, dtype=np.float64)
    dt = yr - float(reference_year)
    if not np.isfinite(slope_model_hist) or abs(slope_model_hist) < 1e-12:
        return x
    if mode == "ratio":
        if not np.isfinite(slope_obs_hist):
            return x
        return x * (slope_obs_hist / slope_model_hist)
    if not np.isfinite(slope_obs_hist):
        return x
    return x - (slope_model_hist - slope_obs_hist) * dt


# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------
def ks_statistic(a: np.ndarray, b: np.ndarray) -> Tuple[float, float]:
    from scipy import stats

    av = np.asarray(a, dtype=float).ravel()
    bv = np.asarray(b, dtype=float).ravel()
    av = av[np.isfinite(av)]
    bv = bv[np.isfinite(bv)]
    if av.size < 2 or bv.size < 2:
        return np.nan, np.nan
    res = stats.ks_2samp(av, bv)
    return float(res.statistic), float(res.pvalue)


def kge(obs: np.ndarray, sim: np.ndarray) -> float:
    """Kling-Gupta efficiency."""
    o = np.asarray(obs, dtype=float).ravel()
    s = np.asarray(sim, dtype=float).ravel()
    ok = np.isfinite(o) & np.isfinite(s)
    o, s = o[ok], s[ok]
    if o.size < 3:
        return np.nan
    r = float(np.corrcoef(o, s)[0, 1])
    sd_o = float(np.std(o, ddof=1))
    alpha = float(np.std(s, ddof=1) / sd_o) if sd_o else np.nan
    # KGE is only defined for series with a non-zero mean; for anomaly-like
    # series the beta ratio explodes, so return NaN instead of a nonsense value
    if abs(float(np.mean(o))) < 0.1 * (sd_o or 1.0):
        return np.nan
    beta = float(np.mean(s) / np.mean(o)) if np.mean(o) else np.nan
    return float(1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def diagnostics(obs: np.ndarray, raw: np.ndarray, corrected: np.ndarray) -> Dict[str, float]:
    """Distribution / bias metrics between observed, raw and corrected series."""
    o = np.asarray(obs, dtype=float)
    r = np.asarray(raw, dtype=float)
    c = np.asarray(corrected, dtype=float)
    ok = np.isfinite(o)
    out: Dict[str, float] = {}
    out["bias_raw"] = float(np.nanmean(r[ok]) - np.nanmean(o[ok])) if ok.any() else np.nan
    out["bias_corr"] = float(np.nanmean(c[ok]) - np.nanmean(o[ok])) if ok.any() else np.nan
    out["rmse_raw"] = float(np.sqrt(np.nanmean((r - o) ** 2)))
    out["rmse_corr"] = float(np.sqrt(np.nanmean((c - o) ** 2)))
    sd_o = np.nanstd(o, ddof=1)
    out["sd_ratio_raw"] = float(np.nanstd(r, ddof=1) / sd_o) if sd_o else np.nan
    out["sd_ratio_corr"] = float(np.nanstd(c, ddof=1) / sd_o) if sd_o else np.nan
    ok2 = np.isfinite(o) & np.isfinite(r)
    out["corr_raw"] = float(np.corrcoef(o[ok2], r[ok2])[0, 1]) if ok2.sum() > 2 else np.nan
    ok3 = np.isfinite(o) & np.isfinite(c)
    out["corr_corr"] = float(np.corrcoef(o[ok3], c[ok3])[0, 1]) if ok3.sum() > 2 else np.nan
    ks_r, p_r = ks_statistic(o, r)
    ks_c, p_c = ks_statistic(o, c)
    out.update({"ks_raw": ks_r, "ks_p_raw": p_r, "ks_corr": ks_c, "ks_p_corr": p_c})
    out["kge_raw"] = kge(o, r)
    out["kge_corr"] = kge(o, c)
    return out


def save_params(params: xr.Dataset, path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    enc = {}
    for name, var in params.data_vars.items():
        if var.ndim >= 2:
            enc[name] = {"zlib": True, "complevel": 4, "dtype": "float32"}
    try:
        params.to_netcdf(path, encoding=enc)
    except Exception:
        params.to_netcdf(path)
    return path


def load_params(path: str) -> xr.Dataset:
    ds = xr.open_dataset(path)
    ds.load()
    ds.close()
    return ds

# ---------------------------------------------------------------------------
# scalar-series variant (data kind == "series")
# ---------------------------------------------------------------------------
def fit_series(
    obs,
    mod,
    method: str = "mean_std",
    multiplicative: bool = False,
    n_quantiles: int = 100,
) -> Dict[str, Any]:
    """Fit the same transfer functions to two 1-D annual series."""
    method = str(method).lower()
    if method == "delta_change":
        method = "mean"
    if method not in METHODS:
        raise ValueError(f"unknown correction method {method!r}")
    o = np.asarray(obs, dtype=float)
    m = np.asarray(mod, dtype=float)
    o = o[np.isfinite(o)]
    m = m[np.isfinite(m)]
    if o.size < 5 or m.size < 5:
        raise ValueError("need at least 5 valid years in both series to fit a correction")
    params: Dict[str, Any] = {
        "method": method,
        "multiplicative": bool(multiplicative),
        "mu_o": float(np.mean(o)),
        "sd_o": float(np.std(o, ddof=1)),
        "mu_m": float(np.mean(m)),
        "sd_m": float(np.std(m, ddof=1)),
    }
    if method == "none":
        return params
    if method == "mean":
        params["factor"] = params["mu_o"] / params["mu_m"] if multiplicative and abs(params["mu_m"]) > 1e-12 else 1.0
        params["offset"] = params["mu_o"] - params["mu_m"]
    elif method == "mean_std":
        if multiplicative:
            ratio = params["mu_o"] / params["mu_m"] if abs(params["mu_m"]) > 1e-12 else 1.0
            cv_o = params["sd_o"] / params["mu_o"] if abs(params["mu_o"]) > 1e-12 else 1.0
            cv_m = params["sd_m"] / params["mu_m"] if abs(params["mu_m"]) > 1e-12 else 1.0
            params["factor"] = ratio
            params["cv_ratio"] = cv_o / cv_m if abs(cv_m) > 1e-12 else 1.0
        else:
            params["factor"] = params["sd_o"] / params["sd_m"] if params["sd_m"] > 1e-12 else 1.0
            params["offset"] = params["mu_o"] - params["mu_m"] * params["factor"]
    else:
        nq = max(10, int(n_quantiles))
        probs = np.linspace(0.0, 1.0, nq)
        params["probs"] = probs
        params["q_mod"] = np.nanquantile(m, probs)
        params["q_obs"] = np.nanquantile(o, probs)
    return params


def apply_series(params: Dict[str, Any], values) -> np.ndarray:
    """Apply a series transfer function (NaN-preserving)."""
    v = np.asarray(values, dtype=float)
    method = str(params.get("method", "none"))
    if method == "none":
        return v.copy()
    mult = bool(params.get("multiplicative", False))
    out = v.copy()
    if method == "mean":
        out = v * params["factor"] if mult else v + params["offset"]
    elif method == "mean_std":
        if mult:
            out = params["mu_o"] + (v * params["factor"] - params["mu_o"]) * params["cv_ratio"]
        else:
            out = (v - params["mu_m"]) * params["factor"] + params["mu_o"]
    elif method == "quantile_mapping":
        probs = np.asarray(params["probs"])
        q_mod = np.asarray(params["q_mod"])
        q_obs = np.asarray(params["q_obs"])
        out = np.full_like(v, np.nan)
        ok = np.isfinite(v)
        out[ok] = np.interp(v[ok], probs, q_obs)
    elif method == "quantile_delta_mapping":
        probs = np.asarray(params["probs"])
        q_mod = np.asarray(params["q_mod"])
        q_obs = np.asarray(params["q_obs"])
        ok = np.isfinite(v)
        out = v.copy()
        if ok.sum() >= 2:
            ranks = _rank_within_columns(v[ok][:, None])[:, 0]
            base = np.interp(ranks, probs, q_mod)
            target = np.interp(ranks, probs, q_obs)
            out[ok] = v[ok] * (target / np.where(np.abs(base) > 1e-12, base, 1.0)) if mult else v[ok] + (target - base)
    lower = params.get("lower_bound")
    if lower is not None:
        out = np.where(out < float(lower), float(lower), out)
    return out
