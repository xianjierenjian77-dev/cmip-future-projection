"""Multi-model ensemble statistics.

The ensemble is stored as a plain ``(n_models, n_years)`` matrix of NaN-padded
series, which keeps the plotting and GEV code trivial while tolerating models
with different coverage.
"""

from __future__ import annotations

import logging
import warnings
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .indices import linear_trend

log = logging.getLogger("cmipkit.ensemble")


def align(series: Dict[str, Tuple[np.ndarray, np.ndarray]], years: Sequence[int]) -> Tuple[np.ndarray, List[str]]:
    """Stack ``{model: (years, values)}`` onto a common year axis."""
    years = np.asarray(years, dtype=int)
    if years.ndim != 1 or np.any(np.diff(years) <= 0):
        raise ValueError("ensemble target years must be one-dimensional, sorted and unique")
    names = sorted(series)
    mat = np.full((len(names), years.size), np.nan)
    for i, name in enumerate(names):
        yrs, vals = series[name]
        yrs = np.asarray(yrs, dtype=int)
        vals = np.asarray(vals, dtype=float)
        if yrs.ndim != 1 or vals.ndim != 1 or yrs.size != vals.size:
            raise ValueError(f"{name}: years and values must be equally long one-dimensional arrays")
        if np.unique(yrs).size != yrs.size:
            raise ValueError(f"{name}: duplicate years cannot be aligned unambiguously")
        if yrs.size == 0 or years.size == 0:
            continue
        idx = np.searchsorted(years, yrs)
        ok = (idx < years.size) & (years[np.clip(idx, 0, years.size - 1)] == yrs)
        mat[i, idx[ok]] = vals[ok]
    return mat, names


def stats(
    matrix: np.ndarray,
    statistic: str = "mean",
    bands: Sequence[int] = (10, 90, 20, 80),
    weights: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """Center and inter-model spread (not a confidence interval).

    Weighted quantiles use the inverse empirical CDF: the smallest value whose
    cumulative normalized weight reaches q. Weighting applies to both the
    center and bands; effective model counts use Kish's definition.
    """
    mat = np.asarray(matrix, dtype=float)
    if mat.ndim != 2:
        raise ValueError("ensemble matrix must have shape (models, years)")
    if statistic not in {"mean", "median"}:
        raise ValueError("ensemble statistic must be mean or median")
    if any(not np.isfinite(q) or q < 0 or q > 100 for q in bands):
        raise ValueError("percentile bands must lie between 0 and 100")
    mat = np.where(np.isfinite(mat), mat, np.nan)
    if weights is None:
        count = np.isfinite(mat).sum(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            center = (np.nanmedian(mat, axis=0) if statistic == "median" else
                      np.divide(np.nansum(mat, axis=0), count,
                                out=np.full(mat.shape[1], np.nan), where=count > 0))
            out = {"center": center, "n_models": count, "n_effective": count.astype(float)}
            for q in bands:
                out[f"p{int(q)}"] = (np.nanpercentile(mat, float(q), axis=0) if mat.shape[0]
                                      else np.full(mat.shape[1], np.nan))
        return out

    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or w.size != mat.shape[0] or not np.isfinite(w).all() or np.any(w < 0) or w.sum() <= 0:
        raise ValueError("weights must be finite nonnegative values, one per model, with positive sum")
    w = w / w.sum()
    out = {"center": np.full(mat.shape[1], np.nan), "n_models": np.zeros(mat.shape[1], dtype=int),
           "n_effective": np.zeros(mat.shape[1], dtype=float)}
    for q in bands:
        out[f"p{int(q)}"] = np.full(mat.shape[1], np.nan)
    for j in range(mat.shape[1]):
        valid = np.isfinite(mat[:, j]) & (w > 0)
        if not valid.any():
            continue
        vals, ww = mat[valid, j], w[valid]
        order = np.argsort(vals)
        vals, ww = vals[order], ww[order]
        ww = ww / ww.sum()
        cumulative = np.cumsum(ww)
        cumulative[-1] = 1.0
        quantile = lambda q: vals[min(np.searchsorted(cumulative, q, side="left"), vals.size - 1)]
        out["center"][j] = np.dot(vals, ww) if statistic == "mean" else quantile(0.5)
        out["n_models"][j] = vals.size
        out["n_effective"][j] = 1.0 / np.sum(ww ** 2)
        for q in bands:
            out[f"p{int(q)}"][j] = quantile(float(q) / 100.0)
    return out


def model_slopes(matrix: np.ndarray, years: Sequence[int]) -> np.ndarray:
    years = np.asarray(years, dtype=float)
    slopes = np.full(matrix.shape[0], np.nan)
    for i in range(matrix.shape[0]):
        t = linear_trend(years, matrix[i])
        slopes[i] = t["slope"]
    return slopes


def skill_weights(scores: Dict[str, float], power: float = 2.0) -> np.ndarray:
    """Inverse-error weights from a ``{model: rmse}`` mapping (sorted by name)."""
    names = sorted(scores)
    err = np.asarray([abs(scores[n]) for n in names], dtype=float)
    ok = np.isfinite(err)
    w = np.zeros_like(err)
    if not ok.any():
        return np.ones_like(err)
    inv = 1.0 / np.maximum(err[ok], 1e-12) ** power
    w[ok] = inv / inv.sum()
    return w


def significance_of_change(hist: np.ndarray, future: np.ndarray, alpha: float = 0.05) -> Dict[str, float]:
    """Welch t-test between the historical and future multi-model means."""
    from scipy import stats as st

    h = np.asarray(hist, dtype=float)
    f = np.asarray(future, dtype=float)
    h = h[np.isfinite(h)]
    f = f[np.isfinite(f)]
    if h.size < 3 or f.size < 3:
        return {"diff": np.nan, "p_value": np.nan, "significant": False}
    res = st.ttest_ind(f, h, equal_var=False)
    diff = float(np.mean(f) - np.mean(h))
    return {"diff": diff, "p_value": float(res.pvalue), "significant": bool(res.pvalue < alpha)}
