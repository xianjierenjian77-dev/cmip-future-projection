"""Multi-model ensemble statistics.

The ensemble is stored as a plain ``(n_models, n_years)`` matrix of NaN-padded
series, which keeps the plotting and GEV code trivial while tolerating models
with different coverage.
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .indices import linear_trend

log = logging.getLogger("cmipkit.ensemble")


def align(series: Dict[str, Tuple[np.ndarray, np.ndarray]], years: Sequence[int]) -> Tuple[np.ndarray, List[str]]:
    """Stack ``{model: (years, values)}`` onto a common year axis."""
    years = np.asarray(years, dtype=int)
    names = sorted(series)
    mat = np.full((len(names), years.size), np.nan)
    for i, name in enumerate(names):
        yrs, vals = series[name]
        yrs = np.asarray(yrs, dtype=int)
        vals = np.asarray(vals, dtype=float)
        if yrs.size == 0:
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
    """Ensemble central tendency and percentile bands along axis 0."""
    mat = np.asarray(matrix, dtype=float)
    with np.errstate(invalid="ignore"):
        if weights is not None:
            w = np.asarray(weights, dtype=float)
            w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
            w = w / w.sum() if w.sum() > 0 else None
        else:
            w = None
        if w is None:
            center = np.nanmedian(mat, axis=0) if statistic == "median" else np.nanmean(mat, axis=0)
        else:
            valid = np.isfinite(mat)
            wmat = np.where(valid, w[:, None], 0.0)
            wsum = wmat.sum(axis=0)
            center = np.where(wsum > 0, np.nansum(np.where(valid, mat, 0.0) * wmat, axis=0) / np.where(wsum > 0, wsum, 1.0), np.nan)
        out: Dict[str, np.ndarray] = {"center": center, "n_models": np.sum(np.isfinite(mat), axis=0)}
        for q in bands:
            out[f"p{int(q)}"] = np.nanpercentile(mat, float(q), axis=0)
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
