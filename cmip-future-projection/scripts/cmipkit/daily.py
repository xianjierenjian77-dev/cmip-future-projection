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
import calendar
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
    ks_p_naive: float = np.nan
    method: str = "mle"
    ok: bool = True
    extra: Dict[str, float] = field(default_factory=dict)

    def return_level(self, period: float) -> float:
        """Level exceeded once every ``period`` years (period > 0)."""
        if not np.isfinite(period) or period <= 0 or not np.isfinite(self.rate) or self.rate <= 0:
            return np.nan
        lam_t = self.rate * float(period)
        if lam_t < 1:
            return np.nan
        if abs(self.shape) < 1e-8:
            return float(self.threshold + self.scale * np.log(lam_t))
        return float(self.threshold + self.scale / self.shape * (lam_t ** self.shape - 1.0))

    def return_levels(self, periods: Sequence[float]) -> np.ndarray:
        return np.asarray([self.return_level(t) for t in periods], dtype=float)

    def return_period(self, value: float) -> float:
        if not np.isfinite(value) or value < self.threshold or not np.isfinite(self.rate) or self.rate <= 0:
            return np.nan
        y = (value - self.threshold) / self.scale
        if abs(self.shape) < 1e-8:
            return float(np.exp(y) / self.rate)
        base = 1.0 + self.shape * y
        if base <= 0:
            return np.inf
        return float(base ** (1.0 / self.shape) / self.rate)

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
            "ks_p_naive": self.ks_p_naive,
            "method": self.method,
        }


# ---------------------------------------------------------------------------
# declustering
# ---------------------------------------------------------------------------
def decluster(
    values: np.ndarray,
    threshold: float,
    run_length: int = 3,
    *,
    dates: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the peak value and index of each cluster above ``threshold``.

    Two exceedances belong to the same cluster when the gap between them is at
    most ``run_length`` days; the cluster contributes its maximum.
    """
    v = np.asarray(values, dtype=float)
    if v.ndim != 1 or int(run_length) < 1:
        raise ValueError("decluster requires a one-dimensional series and run_length >= 1")
    elapsed = _date_parts(dates, v.size)[2] if dates is not None else np.arange(v.size, dtype=float)
    idx = np.where(np.isfinite(v) & (v > float(threshold)))[0]
    if idx.size == 0:
        return np.array([], dtype=int), np.array([], dtype=float)
    clusters = [[idx[0]]]
    for k in idx[1:]:
        last = clusters[-1][-1]
        gap = elapsed[k] - elapsed[last]
        continuous = gap == k - last and np.isfinite(v[last:k + 1]).all()
        if continuous and gap <= int(run_length):
            clusters[-1].append(k)
        else:
            clusters.append([k])
    peak_idx = np.array([max(c, key=lambda i: v[i]) for c in clusters], dtype=int)
    return peak_idx, v[peak_idx]


def _date_parts(dates, size):
    """Validate actual daily timestamps without converting CF calendars."""
    d = np.asarray(dates)
    if d.ndim != 1 or d.size != size:
        raise ValueError("dates and daily values must be one-dimensional and equally long")
    if not d.size:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
    if np.issubdtype(d.dtype, np.datetime64):
        if np.isnat(d).any():
            raise ValueError("daily dates contain NaT")
        day = d.astype("datetime64[D]")
        elapsed = np.asarray((d - d[0]) / np.timedelta64(1, "D"), dtype=float)
        years = day.astype("datetime64[Y]").astype(int) + 1970
        months = day.astype("datetime64[M]").astype(int) % 12 + 1
    else:
        try:
            years = np.array([x.year for x in d], dtype=int)
            months = np.array([x.month for x in d], dtype=int)
            elapsed = np.array([(x - d[0]).total_seconds() / 86400.0 for x in d], dtype=float)
        except (AttributeError, TypeError) as exc:
            raise ValueError("daily dates must be numpy datetime64 or Python/CF datetime objects") from exc
    if np.any(np.diff(elapsed) <= 0):
        raise ValueError("daily timestamps must be sorted and unique; duplicate days are not independent samples")
    if np.any(np.abs(elapsed - np.round(elapsed)) > 1e-8):
        raise ValueError("daily extrema require a daily, not subdaily, time axis")
    return years, months, elapsed


def _season_layout(dates, size, months=None):
    years, actual_months, elapsed = _date_parts(dates, size)
    wanted = list(range(1, 13)) if months is None else [int(m) for m in months]
    if not wanted or len(set(wanted)) != len(wanted) or any(m < 1 or m > 12 for m in wanted):
        raise ValueError("months must contain unique integers from 1 through 12")
    # Match IO aggregation: DJF/NDJFM end-year labels also work when callers
    # supply sorted month lists. The largest omitted interval separates years.
    early = set()
    if len(wanted) < 12 and 1 in wanted and 12 in wanted:
        ordered = np.array(sorted(wanted), dtype=int)
        start = int(ordered[int(np.argmax(np.diff(ordered))) + 1])
        early = {m for m in wanted if m >= start}
    labels = years + np.isin(actual_months, list(early))
    return labels.astype(int), actual_months, elapsed, wanted, early


def _expected_days(year, months, early, template):
    total = 0
    for month in months:
        calendar_year = int(year) - int(month in early)
        if hasattr(template, "calendar"):
            total += int(type(template)(calendar_year, month, 1).daysinmonth)
        else:
            total += calendar.monthrange(calendar_year, month)[1]
    return total


def exposure_years(dates, values, months=None) -> float:
    """Observed season-years: valid days / expected days in each named season.

    A complete JJA contributes one season-year, including in a 360-day
    calendar. Partial seasons contribute fractional exposure; they are not
    silently counted as complete years. This assumes missing days are not
    preferentially extreme, which the caller must assess.
    """
    v = np.asarray(values, dtype=float)
    labels, actual_months, _, wanted, early = _season_layout(dates, v.size, months)
    if not v.size:
        return 0.0
    active = np.isin(actual_months, wanted)
    total = 0.0
    for year in np.unique(labels[active]):
        observed = int(np.count_nonzero((labels == year) & active & np.isfinite(v)))
        total += observed / _expected_days(year, wanted, early, np.asarray(dates)[0])
    return float(total)


def fit_pot(
    values: np.ndarray,
    threshold: float,
    n_years: float,
    decluster_days: int = 3,
    min_peaks: int = 30,
    *,
    dates: Optional[np.ndarray] = None,
    months: Optional[Sequence[int]] = None,
) -> Optional[POTFit]:
    """Fit a GPD to declustered excesses and attach the Poisson rate."""
    from scipy import stats

    v = np.asarray(values, dtype=float)
    if v.ndim != 1 or not np.isfinite(threshold):
        raise ValueError("POT needs a one-dimensional series and finite threshold")
    if dates is not None:
        n_years = exposure_years(dates, v, months)
        if months is not None:
            _, actual_months, _ = _date_parts(dates, v.size)
            v = np.where(np.isin(actual_months, months), v, np.nan)
    if not np.isfinite(n_years) or n_years <= 0:
        raise ValueError("POT exposure n_years must be positive and finite")
    if v.size == 0:
        return None
    _, peaks = decluster(v, threshold, decluster_days, dates=dates)
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
    if not np.isfinite(stats.genpareto.logpdf(excess, shape, loc=0.0, scale=scale)).all():
        log.warning("GPD fit places observations outside its support; skipping fit")
        return None
    fit = POTFit(
        threshold=float(threshold),
        shape=float(shape),
        scale=float(scale),
        rate=float(peaks.size / float(n_years)),
        n_peaks=int(peaks.size),
        n_years=float(n_years),
        decluster_days=int(decluster_days),
        method=method,
    )
    try:
        ks = stats.kstest(excess, lambda x: stats.genpareto.cdf(x, fit.shape, loc=0.0, scale=fit.scale))
        fit.ks_stat, fit.ks_p_naive = float(ks.statistic), float(ks.pvalue)
    except Exception:
        pass
    return fit


def annual_maxima(dates: np.ndarray, values: np.ndarray, *, months=None,
                  min_coverage: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """Annual maxima of a daily series -> (years, maxima)."""
    v = np.asarray(values, dtype=float)
    d = np.asarray(dates)
    if v.ndim != 1 or not (0 <= float(min_coverage) <= 1):
        raise ValueError("annual_maxima needs a 1-D series and min_coverage in [0,1]")
    years, actual_months, _, wanted, early = _season_layout(d, v.size, months)
    out_years, out_max = [], []
    for year in np.unique(years):
        sel = v[(years == year) & np.isin(actual_months, wanted)]
        sel = sel[np.isfinite(sel)]
        if sel.size == 0:
            continue
        expected = _expected_days(year, wanted, early, d[0])
        if sel.size / expected + 1e-12 < float(min_coverage):
            log.warning("excluding incomplete extreme block %s: %d/%d valid days", year, sel.size, expected)
            continue
        out_years.append(int(year))
        out_max.append(float(np.max(sel)))
    return np.asarray(out_years, dtype=int), np.asarray(out_max, dtype=float)


def fit_block(
    dates: np.ndarray,
    values: np.ndarray,
    distribution: str = "gev",
    min_years: int = 20,
    *,
    months: Optional[Sequence[int]] = None,
    min_coverage: float = 1.0,
) -> Optional[Fit]:
    """Annual-maxima GEV/Gumbel fit of a daily series."""
    years, maxima = annual_maxima(dates, values, months=months, min_coverage=min_coverage)
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
        "n_valid": np.sum(np.isfinite(curves), axis=0),
    }
    for q in bands:
        out[f"p{int(q)}"] = np.nanpercentile(curves, float(q), axis=0)
    return out
