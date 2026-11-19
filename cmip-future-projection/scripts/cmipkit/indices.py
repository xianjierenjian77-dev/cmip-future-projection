"""Index construction and time-series statistics.

A *field* here is always ``(year, latitude, longitude)``.  This module turns a
field into the one-dimensional index that the rest of the pipeline (ensemble
statistics, trend tests, GEV) actually works with, and provides the trend /
significance tools used everywhere.

All reductions are vectorised (no python loop over grid points) because the
reference workflow, with a 1x1 degree grid and 75 years, would otherwise spend
hours in per-gridpoint python loops.
"""

from __future__ import annotations

import logging
import math
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import xarray as xr

log = logging.getLogger("cmipkit.indices")


# ---------------------------------------------------------------------------
# weights and masks
# ---------------------------------------------------------------------------
def latitude_weights(lat: np.ndarray) -> np.ndarray:
    return np.cos(np.deg2rad(np.asarray(lat, dtype=float)))


def area_weights(da: xr.DataArray) -> xr.DataArray:
    _require_regular_grid(da)
    w = latitude_weights(da.latitude.values)
    return xr.DataArray(w, dims="latitude", coords={"latitude": da.latitude})


def _require_regular_grid(da: xr.DataArray) -> None:
    """Cosine-latitude weights require equal angular grid-cell widths."""
    lat = np.sort(np.asarray(da.latitude.values, dtype=float))
    lon = np.sort(np.unique(np.mod(np.asarray(da.longitude.values, dtype=float), 360)))
    if lat.size > 2 and not np.allclose(np.diff(lat), np.diff(lat)[0], rtol=1e-5, atol=1e-8):
        raise ValueError("cosine-latitude weights require a regular latitude grid; regrid first or implement cell-area weights")
    if lon.size > 2:
        gaps = np.diff(np.r_[lon, lon[0] + 360])
        covered_gaps = np.delete(gaps, int(np.argmax(gaps)))
        if not np.allclose(covered_gaps, covered_gaps[0], rtol=1e-5, atol=1e-8):
            raise ValueError("cosine-latitude weights require a regular longitude grid; regrid first or implement cell-area weights")


def longitude_mask(lon: np.ndarray, bounds: Sequence[float]) -> np.ndarray:
    """Select the eastward interval from west to east, including wrapped boxes."""
    if len(bounds) != 2:
        raise ValueError("longitude bounds must be [west, east]")
    west, east = map(float, bounds)
    lon = np.asarray(lon, dtype=float)
    if abs(east - west) >= 360.0 - 1e-8:
        return np.isfinite(lon)
    width = (east - west) % 360.0
    distance = np.mod(lon - west, 360.0)
    return np.isfinite(lon) & (distance <= width + 1e-8)


def region_mask(lat: np.ndarray, lon: np.ndarray, region: Dict[str, Any]) -> np.ndarray:
    """Boolean (lat, lon) mask for a rectangular region, wrap-around aware."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat_lo, lat_hi = sorted(float(v) for v in region["lat"])
    lat_ok = (lat >= lat_lo) & (lat <= lat_hi)
    lon_ok = longitude_mask(lon, region["lon"])
    return lat_ok[:, None] & lon_ok[None, :]


def _weighted_flat_mean(values: np.ndarray, mask: np.ndarray, w2d: np.ndarray) -> np.ndarray:
    """Weighted mean of (..., lat, lon) arrays over the mask, NaN-safe."""
    shape = values.shape[:-2]
    flat = values.reshape(-1, values.shape[-2], values.shape[-1])
    w = np.where(mask, w2d, 0.0)
    wsum = np.nansum(w)
    if not np.isfinite(wsum) or wsum == 0:
        return np.full(flat.shape[0], np.nan).reshape(shape)
    data = np.where(mask[None, :, :], flat, np.nan)
    num = np.nansum(data * w[None, :, :], axis=(1, 2))
    den = np.nansum(np.where(np.isfinite(data), w[None, :, :], 0.0), axis=(1, 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out = np.where(den > 0, out, np.nan)
    return out.reshape(shape)


# ---------------------------------------------------------------------------
# field transforms
# ---------------------------------------------------------------------------
def remove_zonal_mean(da: xr.DataArray) -> xr.DataArray:
    return da - da.mean(dim="longitude", skipna=True)


def anomaly(da: xr.DataArray, period: Sequence[int], climatology: Optional[xr.DataArray] = None) -> xr.DataArray:
    """Subtract a baseline climatology.

    ``climatology`` must be supplied when ``da`` does not overlap ``period``
    (a future field has no baseline years); otherwise it is computed from ``da``.
    """
    if climatology is None:
        lo, hi = int(period[0]), int(period[1])
        sel = da.sel(year=slice(lo, hi))
        if sel.year.size == 0:
            raise ValueError(
                f"cannot build an anomaly for years {int(da.year[0])}-{int(da.year[-1])} with a "
                f"{period[0]}-{period[1]} baseline: pass the historical climatology explicitly"
            )
        climatology = sel.mean(dim="year", skipna=True)
    if not np.isfinite(climatology.values).any():
        raise ValueError("baseline climatology is all NaN")
    return da - climatology


def detrend_field(da: xr.DataArray) -> xr.DataArray:
    """Remove the least-squares linear trend at every grid point."""
    years = da.year.values.astype(float)
    y = da.values.astype(np.float64)
    x = years.reshape((-1,) + (1,) * (y.ndim - 1))
    valid = np.isfinite(y)
    n = valid.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        xmean = np.sum(np.where(valid, x, 0), axis=0) / n
        ymean = np.nansum(y, axis=0) / n
        xc = x - xmean
        num = np.sum(np.where(valid, xc * (y - ymean), 0), axis=0)
        den = np.sum(np.where(valid, xc ** 2, 0), axis=0)
        slope = np.where((n >= 2) & (den > 0), num / den, np.nan)
    fit = ymean + slope * xc
    out = y - fit
    return xr.DataArray(out, dims=da.dims, coords=da.coords, name=da.name)


# ---------------------------------------------------------------------------
# index builders
# ---------------------------------------------------------------------------
def region_mean_series(field: xr.DataArray, region: Dict[str, Any]) -> np.ndarray:
    _require_regular_grid(field)
    lat = field.latitude.values
    lon = field.longitude.values
    mask = region_mask(lat, lon, region)
    if not mask.any():
        raise ValueError(f"region {region} does not intersect the grid")
    w2d = np.broadcast_to(latitude_weights(lat)[:, None], mask.shape).copy()
    return _weighted_flat_mean(field.values.astype(np.float64), mask, w2d)


def make_template(
    field: xr.DataArray,
    spec: Dict[str, Any],
    fallback_period: Sequence[int],
    model_fields: Optional[Dict[str, xr.DataArray]] = None,
) -> xr.DataArray:
    """Build the anomaly pattern used by ``projection_index``."""
    source = str(spec.get("template_source", "reference")).lower()
    if source == "year":
        yr = int(spec["template_year"])
        if yr not in field.year.values:
            raise ValueError(f"template_year {yr} is not available in the reference field (years {field.year.values[[0, -1]]})")
        tmpl = field.sel(year=yr)
    elif source == "model":
        if not model_fields:
            raise ValueError("template_source='model' but no model fields were supplied")
        name = spec.get("template_model")
        if name is not None and name not in model_fields:
            raise ValueError(f"template_model {name!r} is not available")
        key = name if name is not None else list(model_fields)[0]
        mfield = model_fields[key]
        period = spec.get("template_period") or fallback_period
        tmpl = mfield.sel(year=slice(int(period[0]), int(period[1]))).mean("year", skipna=True)
    else:
        period = spec.get("template_period") or fallback_period
        tmpl = field.sel(year=slice(int(period[0]), int(period[1]))).mean("year", skipna=True)
    return tmpl - tmpl.mean("longitude", skipna=True)


def projection_series(field: xr.DataArray, template: xr.DataArray, region: Optional[Dict[str, Any]] = None) -> np.ndarray:
    """Project a field onto an anomaly template (wave-train style index).

    ``PI(t) = sum(A' T w) / sum(T^2 w)`` with ``A'`` the zonal anomaly of the
    field and ``T`` the zonal-anomaly template.
    """
    _require_regular_grid(field)
    anom = field - field.mean("longitude", skipna=True)
    lat = anom.latitude.values
    lon = anom.longitude.values
    if "latitude" in template.dims:
        from .io import regrid
        tmpl = regrid(template, {"latitude": lat, "longitude": lon})
    else:
        tmpl = template
    tmpl_vals = np.asarray(tmpl.values, dtype=np.float64)
    A = anom.values.astype(np.float64)

    if region:
        mask = region_mask(lat, lon, region)
        A = np.where(mask[None, :, :], A, np.nan)
        tmpl_vals = np.where(mask, tmpl_vals, np.nan)

    w = np.broadcast_to(latitude_weights(lat)[:, None], tmpl_vals.shape)
    template_den = np.nansum(np.where(np.isfinite(tmpl_vals), tmpl_vals ** 2 * w, 0.0))
    if template_den == 0 or not np.isfinite(template_den):
        raise ValueError("projection template has zero variance over the target region")
    valid = np.isfinite(A) & np.isfinite(tmpl_vals)[None, :, :]
    den = np.sum(np.where(valid, (tmpl_vals ** 2 * w)[None, :, :], 0.0), axis=(1, 2))
    numer = np.nansum(A * tmpl_vals[None, :, :] * w[None, :, :], axis=(1, 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, numer / den, np.nan)


def build_index(
    field: xr.DataArray,
    spec: Dict[str, Any],
    reference_period: Sequence[int],
    template: Optional[xr.DataArray] = None,
    climatology: Optional[xr.DataArray] = None,
    *,
    standardization: Optional[Tuple[float, float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(years, values)`` for one field according to ``index`` config."""
    years = np.asarray(field.year.values, dtype=int)
    itype = str(spec.get("type", "region_mean"))

    if itype == "series":
        vals = np.asarray(field.values, dtype=np.float64).reshape(-1)
        return years, vals

    work = field
    if bool(spec.get("remove_zonal_mean", False)):
        work = remove_zonal_mean(work)
        if climatology is not None:
            climatology = remove_zonal_mean(climatology)
    if bool(spec.get("anomaly", False)) or itype == "region_mean_anomaly":
        work = anomaly(work, reference_period, climatology)

    if itype in ("region_mean", "region_mean_anomaly"):
        region = spec.get("region")
        if not region:
            raise ValueError("index.region is required")
        vals = region_mean_series(work, region)
    elif itype == "region_diff":
        regs = spec.get("regions") or []
        if len(regs) != 2:
            raise ValueError("index.regions must be [reference_region, target_region]")
        zonal = remove_zonal_mean(work)
        ref_mean = region_mean_series(zonal, regs[0])
        tgt_mean = region_mean_series(zonal, regs[1])
        vals = tgt_mean - ref_mean
    elif itype == "projection_index":
        if template is None:
            raise ValueError("projection_index requires a template field")
        vals = projection_series(work, template, spec.get("region"))
    else:
        raise ValueError(f"unknown index.type {itype!r}")

    vals = np.asarray(vals, dtype=np.float64)
    if bool(spec.get("detrend", False)):
        vals = detrend_series(years.astype(float), vals)
    if bool(spec.get("standardize", False)):
        if standardization is None:
            baseline = (years >= int(reference_period[0])) & (years <= int(reference_period[1])) & np.isfinite(vals)
            if baseline.sum() < 2:
                raise ValueError("standardization needs at least two historical baseline samples; pass fixed (mean, std) for future data")
            mu = float(np.mean(vals[baseline]))
            sd = float(np.std(vals[baseline], ddof=1))
        else:
            mu, sd = map(float, standardization)
        if not np.isfinite(mu) or not np.isfinite(sd) or sd <= 0:
            raise ValueError("standardization baseline must have finite mean and positive sample standard deviation")
        vals = (vals - mu) / sd
    return years, vals


# ---------------------------------------------------------------------------
# series statistics
# ---------------------------------------------------------------------------
def linear_trend(years: Sequence[float], values: Sequence[float]) -> Dict[str, float]:
    from scipy import stats

    x = np.asarray(years, dtype=float)
    y = np.asarray(values, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return {"slope": np.nan, "slope_per_decade": np.nan, "intercept": np.nan,
                "p_value": np.nan, "r2": np.nan, "stderr": np.nan, "n": int(x.size)}
    res = stats.linregress(x, y)
    return {
        "slope": float(res.slope),
        "slope_per_decade": float(res.slope * 10.0),
        "intercept": float(res.intercept),
        "p_value": float(res.pvalue),
        "r2": float(res.rvalue ** 2),
        "stderr": float(res.stderr),
        "n": int(x.size),
    }


def theil_sen_trend(years: Sequence[float], values: Sequence[float]) -> Dict[str, float]:
    from scipy import stats

    x = np.asarray(years, dtype=float)
    y = np.asarray(values, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return {"slope": np.nan, "slope_per_decade": np.nan, "intercept": np.nan, "p_value": np.nan}
    slope, intercept, lo, hi = stats.theilslopes(y, x, 0.95)
    return {
        "slope": float(slope),
        "slope_per_decade": float(slope * 10.0),
        "intercept": float(intercept),
        "ci_low": float(lo),
        "ci_high": float(hi),
        "p_value": float("nan"),
    }


def detrend_series(years: Sequence[float], values: Sequence[float]) -> np.ndarray:
    x = np.asarray(years, dtype=float)
    y = np.asarray(values, dtype=float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return y
    slope, intercept = np.polyfit(x[ok], y[ok], 1)
    return y - (slope * x + intercept)


def mann_kendall(values: Sequence[float]) -> Dict[str, float]:
    """Two-sided Mann-Kendall trend test (no external dependency)."""
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 4:
        return {"tau": np.nan, "s": np.nan, "z": np.nan, "p_value": np.nan, "trend": "insufficient data"}
    s = 0.0
    for k in range(n - 1):
        s += np.sign(x[k + 1:] - x[k]).sum()
    # variance with tie correction
    _, counts = np.unique(x, return_counts=True)
    tie = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var = (n * (n - 1) * (2 * n + 5) - tie) / 18.0
    if var <= 0:
        return {"tau": np.nan, "s": s, "z": 0.0, "p_value": 1.0, "trend": "no trend"}
    if s > 0:
        z = (s - 1) / math.sqrt(var)
    elif s < 0:
        z = (s + 1) / math.sqrt(var)
    else:
        z = 0.0
    from scipy import stats

    p = 2 * (1 - stats.norm.cdf(abs(z)))
    if p < 0.05:
        trend = "increasing" if s > 0 else "decreasing"
    else:
        trend = "no trend"
    return {"tau": float(s / (0.5 * n * (n - 1))), "s": float(s), "z": float(z), "p_value": float(p), "trend": trend}


def fdr_bh(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg FDR control; returns a boolean rejection mask."""
    p = np.asarray(p_values, dtype=float)
    out = np.zeros(p.shape, dtype=bool)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    n = ranked.size
    thresh = alpha * (np.arange(1, n + 1) / n)
    below = np.where(ranked <= thresh)[0]
    if below.size == 0:
        return out
    kmax = below[-1]
    reject_sorted = np.zeros(n, dtype=bool)
    reject_sorted[: kmax + 1] = True
    rej = np.zeros(n, dtype=bool)
    rej[order] = reject_sorted
    out[ok] = rej
    return out


def significance_mask(field: np.ndarray, std: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Two-sided 1.96-sigma style mask (matches the reference workflow)."""
    from scipy import stats

    z = stats.norm.ppf(1 - alpha / 2.0)
    return np.abs(field) > z * std


def compare_slopes(slopes: np.ndarray, reference: float) -> Dict[str, float]:
    """Is the ensemble mean slope different from a single reference slope?"""
    from scipy import stats

    s = np.asarray(slopes, dtype=float)
    s = s[np.isfinite(s)]
    if s.size < 2:
        return {"mean": float(np.nanmean(s)) if s.size else np.nan, "diff": np.nan,
                "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan, "significant": False}
    diff = s - reference
    n = diff.size
    mean_diff = float(np.mean(diff))
    se = float(np.std(diff, ddof=1) / math.sqrt(n))
    if se == 0:
        return {"mean": float(np.mean(s)), "diff": mean_diff, "ci_low": mean_diff,
                "ci_high": mean_diff, "p_value": 0.0, "significant": mean_diff != 0}
    tcrit = stats.t.ppf(0.975, n - 1)
    p = float(2 * (1 - stats.t.cdf(abs(mean_diff / se), n - 1)))
    return {
        "mean": float(np.mean(s)),
        "diff": mean_diff,
        "ci_low": mean_diff - tcrit * se,
        "ci_high": mean_diff + tcrit * se,
        "p_value": p,
        "significant": (mean_diff - tcrit * se) * (mean_diff + tcrit * se) > 0,
    }
