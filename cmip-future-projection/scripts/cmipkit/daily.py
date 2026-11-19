"""Daily-scale extreme analysis: peaks-over-threshold (GPD) and daily block maxima.

The annual pipeline works on one value per year.  Many variables need the daily
scale instead: a heatwave is a run of hot days, a storm is a peak above a
meteorological threshold, and a return period is defined on the daily series.
This module provides the two standard routes:

``pot``   peaks over a threshold, declustered by runs, fitted with a Generalized
          Pareto distribution (GPD).  The Poisson rate of exceedances converts
          the GPD into a return level for a return period in years.
``block`` annual maxima of the daily series, fitted with a GEV — the classical
          block-maxima approach, equivalent to ``aggregate.stat: max`` on daily
          input but computed without re-reading the files.

Both return the same small interface (:class:`POTFit` and
:class:`~cmipkit.extremes.Fit`), so the plotting code is shared with the annual
GEV figures.

References: Coles (2001) ch. 4 for the GPD/Poisson formulation; Zhang et al.
(2004) for the threshold-choice caveat.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

from .extremes import Fit, fit_distribution

log = logging.getLogger("cmipkit.daily")


@dataclass
class POTFit:
    """Peaks-over-threshold fit with a Poisson exceedance rate."""

    threshold: float
    shape: float          # GPD xi
    scale: float          # GPD sigma
    rate: float           # exceedances per year
    n_peaks: int
    n_years: float
    decluster_days: int = 3
    ks_stat: float = np.nan
    ks_p: float = np.nan
    method: str = "mle"
    ok: bool = True
    extra: Dict[str, float] = field(default_factory=dict)

    def return_level(self, period: float) -> float:
        """Level exceeded once every ``period`` years (period > 0)."""
        if period <= 0 or not np.isfinite(self.rate) or self.rate <= 0:
            return np.nan
        lam_t = self.rate * float(period)
        if lam_t <= 0:
            return np.nan
        if abs(self.shape) < 1e-8:
            return float(self.threshold + self.scale * np.log(lam_t))
        return float(self.threshold + self.scale / self.shape * (lam_t ** self.shape - 1.0))

    def return_levels(self, periods: Sequence[float]) -> np.ndarray:
        return np.asarray([self.return_level(t) for t in periods], dtype=float)

    def return_period(self, value: float) -> float:
        if value <= self.threshold or not np.isfinite(self.rate) or self.rate <= 0:
            return np.nan
        y = (value - self.threshold) / self.scale
        if abs(self.shape) < 1e-8:
            return float(np.exp(y) / self.rate)
        base = 1.0 + self.shape * y
        if base <= 0:
            return np.inf
        return float(base ** (-1.0 / self.shape) / self.rate)

    def as_dict(self) -> Dict[str, float]:
        return {
            "kind": "pot",
            "threshold": self.threshold,
            "shape": self.shape,
            "scale": self.scale,
            "rate": self.rate,
            "n_peaks": self.n_peaks,
            "n_years": self.n_years,
            "decluster_days": self.decluster_days,
            "ks_stat": self.ks_stat,
            "ks_p": self.ks_p,
        }


# ---------------------------------------------------------------------------
# declustering
# ---------------------------------------------------------------------------
def decluster(
    values: np.ndarray,
    threshold: float,
    run_length: int = 3,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the peak value and index of each cluster above ``threshold``.

    Two exceedances belong to the same cluster when the gap between them is at
    most ``run_length`` days; the cluster contributes its maximum.
    """
    v = np.asarray(values, dtype=float)
    idx = np.where(np.isfinite(v) & (v > float(threshold)))[0]
    if idx.size == 0:
        return np.array([], dtype=int), np.array([], dtype=float)
    clusters = [[idx[0]]]
    for k in idx[1:]:
        if k - clusters[-1][-1] <= int(run_length):
            clusters[-1].append(k)
        else:
            clusters.append([k])
    peak_idx = np.array([max(c, key=lambda i: v[i]) for c in clusters], dtype=int)
    return peak_idx, v[peak_idx]


def fit_pot(
    values: np.ndarray,
    threshold: float,
    n_years: float,
    decluster_days: int = 3,
    min_peaks: int = 30,
) -> Optional[POTFit]:
    """Fit a GPD to declustered excesses and attach the Poisson rate."""
    from scipy import stats

    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    _, peaks = decluster(v, threshold, decluster_days)
    if peaks.size < int(min_peaks):
        log.debug("POT: only %d peaks above %.3f (< %d)", peaks.size, threshold, min_peaks)
        return None
    excess = peaks - float(threshold)
    excess = excess[excess > 0]
    if excess.size < int(min_peaks) or np.allclose(excess, excess[0]):
        return None
    shape, loc, scale = np.nan, 0.0, np.nan
    method = "mle"
    try:
        shape, loc, scale = stats.genpareto.fit(excess, floc=0.0)
    except Exception as exc:
        log.debug("GPD MLE failed (%s); using moments", exc)
    if not (np.isfinite(shape) and np.isfinite(scale) and scale > 0):
        # method of moments starting values, still with loc fixed at zero
        m1 = float(np.mean(excess))
        m2 = float(np.var(excess))
        if m1 <= 0 or m2 <= 0:
            return None
        shape = float(0.5 * (1.0 - m1 ** 2 / m2))
        scale = float(0.5 * m1 * (m1 ** 2 / m2 + 1.0))
        method = "moments"
    fit = POTFit(
        threshold=float(threshold),
        shape=float(shape),
        scale=float(scale),
        rate=float(peaks.size / max(float(n_years), 1e-9)),
        n_peaks=int(peaks.size),
        n_years=float(n_years),
        decluster_days=int(decluster_days),
        method=method,
    )
    try:
        ks = stats.kstest(excess, lambda x: stats.genpareto.cdf(x, fit.shape, loc=0.0, scale=fit.scale))
        fit.ks_stat, fit.ks_p = float(ks.statistic), float(ks.pvalue)
    except Exception:
        pass
    return fit


def annual_maxima(dates: np.ndarray, values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Annual maxima of a daily series -> (years, maxima)."""
    v = np.asarray(values, dtype=float)
    d = np.asarray(dates)
    years = np.asarray([getattr(x, "year", int(str(x)[:4])) for x in d], dtype=int)
    out_years, out_max = [], []
    for year in np.unique(years):
        sel = v[years == year]
        sel = sel[np.isfinite(sel)]
        if sel.size == 0:
            continue
        out_years.append(int(year))
        out_max.append(float(np.max(sel)))
    return np.asarray(out_years, dtype=int), np.asarray(out_max, dtype=float)


def fit_block(
    dates: np.ndarray,
    values: np.ndarray,
    distribution: str = "gev",
    min_years: int = 20,
) -> Optional[Fit]:
    """Annual-maxima GEV/Gumbel fit of a daily series."""
    years, maxima = annual_maxima(dates, values)
    if years.size < int(min_years):
        return None
    return fit_distribution(maxima, distribution=distribution, method="mle", min_years=min_years)


def threshold_value(values: np.ndarray, quantile: Optional[float] = None,
                    absolute: Optional[float] = None) -> float:
    """Threshold for POT: an absolute value or a quantile of the series."""
    if absolute is not None:
        return float(absolute)
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    q = 0.95 if quantile is None else float(quantile)
    return float(np.nanquantile(v, q))


def curve_stats(curves: np.ndarray, bands: Sequence[int] = (10, 90, 20, 80)) -> Dict[str, np.ndarray]:
    if curves.size == 0:
        return {}
    out = {
        "median": np.nanmedian(curves, axis=0),
        "mean": np.nanmean(curves, axis=0),
        "n_models": int(curves.shape[0]),
    }
    for q in bands:
        out[f"p{int(q)}"] = np.nanpercentile(curves, float(q), axis=0)
    return out
