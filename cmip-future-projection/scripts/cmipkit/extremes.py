"""Extreme-value analysis: GEV / Gumbel fits, return levels and uncertainty.

Conventions
-----------
* ``scipy.stats.genextreme`` is used, whose shape parameter ``c`` equals
  ``-xi`` in the classical GEV notation used by Coles (2001): ``c > 0`` is a
  Weibull-type upper-bounded tail, ``c < 0`` is Frechet / heavy-tailed.
* Annual means are not block maxima. A GEV fitted to an annual index is an
  exploratory distribution assumption, not an extreme-value limit theorem.
* KS distance is descriptive. The ordinary KS p-value is uncalibrated after
  fitting on the same sample, so it is stored as ``ks_p_naive`` and ``ks_p``
  remains missing until a calibrated test is actually performed.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

log = logging.getLogger("cmipkit.extremes")

RETURN_POINTS_DEFAULT = (1.1, 2, 5, 10, 20, 30, 50, 100, 200, 500, 1000)


@dataclass
class Fit:
    distribution: str
    shape: float
    loc: float
    scale: float
    n: int = 0
    method: str = "mle"
    loglik: float = np.nan
    aic: float = np.nan
    ks_stat: float = np.nan
    ks_p: float = np.nan
    ks_p_naive: float = np.nan
    ok: bool = True
    extra: Dict[str, float] = field(default_factory=dict)

    # -- distribution helpers -------------------------------------------
    def _dist(self):
        from scipy import stats

        if self.distribution == "gumbel":
            return stats.gumbel_r
        return stats.genextreme

    def _args(self):
        if self.distribution == "gumbel":
            return (self.loc, self.scale)
        return (self.shape, self.loc, self.scale)

    def cdf(self, x):
        return self._dist().cdf(x, *self._args())

    def ppf(self, p):
        return self._dist().ppf(p, *self._args())

    def return_level(self, period: float):
        """Quantile exceeded once every ``period`` years (for period > 1)."""
        if not np.isfinite(period) or period <= 1:
            return np.nan
        return float(self._dist().isf(1.0 / float(period), *self._args()))

    def return_levels(self, periods: Sequence[float]) -> np.ndarray:
        return np.asarray([self.return_level(t) for t in periods], dtype=float)

    def return_period(self, value: float) -> float:
        """Expected waiting time, using the survival function for tail accuracy."""
        if not np.isfinite(value):
            return np.nan
        exceedance = float(self._dist().sf(value, *self._args()))
        if not np.isfinite(exceedance):
            return np.nan
        return np.inf if exceedance <= 0 else 1.0 / exceedance

    def as_dict(self) -> Dict[str, float]:
        return {
            "distribution": self.distribution,
            "shape": self.shape,
            "shape_xi": -self.shape,
            "loc": self.loc,
            "scale": self.scale,
            "n": self.n,
            "method": self.method,
            "loglik": self.loglik,
            "aic": self.aic,
            "ks_stat": self.ks_stat,
            "ks_p": self.ks_p,
            "ks_p_naive": self.ks_p_naive,
        }


def _clean(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values) if not isinstance(values, np.ndarray) else values, dtype=float).ravel()
    return arr[np.isfinite(arr)]


def fit_distribution(
    values: Iterable[float],
    distribution: str = "gev",
    method: str = "mle",
    min_years: int = 20,
) -> Optional[Fit]:
    """Fit a stationary GEV/Gumbel; the caller establishes block validity."""
    distribution = str(distribution).lower()
    method = str(method).lower()
    if distribution not in {"gev", "gumbel"}:
        raise ValueError("distribution must be gev or gumbel")
    if method not in {"mle", "mm"}:
        raise ValueError("fit method must be mle or mm")
    if int(min_years) < 3:
        raise ValueError("min_years must be at least 3")
    x = _clean(values)
    if x.size < int(min_years):
        log.debug("not enough years for a GEV fit (%d < %d)", x.size, min_years)
        return None
    if np.ptp(x) <= np.finfo(float).eps * max(1.0, np.max(np.abs(x))):
        return None
    from scipy import stats

    candidates = []
    if distribution == "gumbel":
        try:
            loc, scale = stats.gumbel_r.fit(x, method=method.upper())
            if np.isfinite([loc, scale]).all() and scale > 0:
                candidates.append(("gumbel", 0.0, loc, scale, method))
        except Exception:
            pass
        if not candidates:
            scale = float(np.std(x, ddof=1) * np.sqrt(6) / np.pi)
            loc = float(np.mean(x) - np.euler_gamma * scale)
            candidates.append(("gumbel", 0.0, loc, scale, "moments"))
    else:
        for meth in ([str(method).lower()] if str(method).lower() != "mle" else []) + ["mle", "mm"]:
            if any(c[4] == meth for c in candidates):
                continue
            try:
                c, loc, scale = stats.genextreme.fit(x, method=meth)
                if (np.isfinite([c, loc, scale]).all() and scale > 0
                        and np.isfinite(stats.genextreme.logpdf(x, c, loc=loc, scale=scale)).all()):
                    candidates.append(("gev", float(c), float(loc), float(scale), meth))
                    break
            except Exception as exc:
                log.debug("GEV %s fit failed: %s", meth, exc)
        if not candidates:
            log.warning("GEV fit failed; no unrelated Gumbel fallback is substituted")
            return None

    name, shape, loc, scale, used = candidates[0]
    if not np.isfinite(scale) or scale <= 0:
        return None
    fit = Fit(distribution=name, shape=shape, loc=loc, scale=scale, n=int(x.size), method=used)
    try:
        logs = fit._dist().logpdf(x, *fit._args())
        if not np.isfinite(logs).all():
            return None
        fit.loglik = float(np.sum(logs))
        k = 2 if name == "gumbel" else 3
        fit.aic = float(-2 * fit.loglik + 2 * k)
    except Exception:
        pass
    try:
        ks = stats.kstest(x, lambda v: _cdf(name, v, shape, loc, scale))
        fit.ks_stat, fit.ks_p_naive = float(ks.statistic), float(ks.pvalue)
    except Exception:
        pass
    return fit


def _cdf(name: str, x, shape, loc, scale):
    from scipy import stats

    if name == "gumbel":
        return stats.gumbel_r.cdf(x, loc=loc, scale=scale)
    return stats.genextreme.cdf(x, shape, loc=loc, scale=scale)


def _pdf(name: str, x, shape, loc, scale):
    from scipy import stats

    if name == "gumbel":
        return stats.gumbel_r.pdf(x, loc=loc, scale=scale)
    return stats.genextreme.pdf(x, shape, loc=loc, scale=scale)


def fit_many(
    series: Dict[str, Iterable[float]],
    distribution: str = "gev",
    method: str = "mle",
    min_years: int = 20,
) -> Dict[str, Fit]:
    fits: Dict[str, Fit] = {}
    for name, values in series.items():
        fit = fit_distribution(values, distribution, method, min_years)
        if fit is not None:
            fits[name] = fit
    return fits


def return_curves(fits: Dict[str, Fit], periods: Sequence[float]) -> np.ndarray:
    """Matrix ``(n_models, n_periods)`` of return levels."""
    if not fits:
        return np.empty((0, len(periods)))
    return np.vstack([fits[name].return_levels(periods) for name in sorted(fits)])


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


def bootstrap_return_levels(
    values: Iterable[float],
    periods: Sequence[float],
    n_boot: int = 200,
    seed: int = 0,
    min_years: int = 20,
) -> np.ndarray:
    """Percentile bootstrap of the return-level curve."""
    x = _clean(values)
    if n_boot <= 0 or x.size < min_years:
        return np.empty((0, len(periods)))
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(int(n_boot)):
        sample = rng.choice(x, size=x.size, replace=True)
        fit = fit_distribution(sample, "gev", "mle", min_years)
        if fit is not None:
            rows.append(fit.return_levels(periods))
    return np.vstack(rows) if rows else np.empty((0, len(periods)))


def empirical_return_levels(values: Iterable[float]) -> Tuple[np.ndarray, np.ndarray]:
    """Empirical return levels with Gringorten plotting positions."""
    x = np.sort(_clean(values))
    n = x.size
    if n < 2:
        return np.array([]), np.array([])
    p = (np.arange(1, n + 1) - 0.44) / (n + 0.12)
    return 1.0 / (1.0 - p), x


def shape_summary(fits) -> Dict[str, float]:
    """Summary of shape parameters; accepts a dict of :class:`Fit` or a sequence of floats."""
    items = list(fits.values()) if isinstance(fits, dict) else list(fits)
    shapes = np.asarray([f.shape if isinstance(f, Fit) else float(f) for f in items], dtype=float)
    shapes = shapes[np.isfinite(shapes)]
    if shapes.size == 0:
        return {"n": 0}
    return {
        "n": int(shapes.size),
        "mean": float(np.mean(shapes)),
        "median": float(np.median(shapes)),
        "std": float(np.std(shapes, ddof=1)) if shapes.size > 1 else np.nan,
        "min": float(np.min(shapes)),
        "max": float(np.max(shapes)),
    }


def kruskal_scenarios(groups: Dict[str, Sequence[float]]) -> Dict[str, float]:
    """Kruskal-Wallis test that the shape parameters differ between scenarios."""
    from scipy import stats

    samples = [np.asarray(v, dtype=float) for v in groups.values() if len(v) >= 2]
    if len(samples) < 2:
        return {"H": np.nan, "p_value": np.nan, "significant": False}
    try:
        res = stats.kruskal(*samples)
        return {"H": float(res.statistic), "p_value": float(res.pvalue), "significant": bool(res.pvalue < 0.05)}
    except Exception:
        return {"H": np.nan, "p_value": np.nan, "significant": False}
