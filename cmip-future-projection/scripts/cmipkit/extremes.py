"""Extreme-value analysis: GEV / Gumbel fits, return levels and uncertainty.

Conventions
-----------
* ``scipy.stats.genextreme`` is used, whose shape parameter ``c`` equals
  ``-xi`` in the classical GEV notation used by Coles (2001): ``c > 0`` is a
  Weibull-type upper-bounded tail, ``c < 0`` is Frechet / heavy-tailed.
* Fits are done on block maxima (one value per year by construction of the
  yearly fields).  If the series is not a block-maxima series (e.g. an annual
  mean index), the fit is still valid as a tail model but should be reported as
  such.
* Every fit records diagnostics (log-likelihood, AIC, KS on the fitted CDF) so
  that bad fits can be filtered instead of silently plotted.
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
        if period <= 1:
            return np.nan
        return float(self.ppf(1.0 - 1.0 / float(period)))

    def return_levels(self, periods: Sequence[float]) -> np.ndarray:
        return np.asarray([self.return_level(t) for t in periods], dtype=float)

    def return_period(self, value: float) -> float:
        """Expected return period of a given value (nan if outside the support)."""
        p = float(self.cdf(value))
        if not np.isfinite(p) or p >= 1.0:
            return np.inf if p >= 1.0 else np.nan
        if p <= 0.0:
            return np.nan
        return 1.0 / (1.0 - p)

    def as_dict(self) -> Dict[str, float]:
        return {
            "distribution": self.distribution,
            "shape": self.shape,
            "loc": self.loc,
            "scale": self.scale,
            "n": self.n,
            "method": self.method,
            "loglik": self.loglik,
            "aic": self.aic,
            "ks_stat": self.ks_stat,
            "ks_p": self.ks_p,
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
    """Fit a GEV (or Gumbel) to a sample of annual values."""
    x = _clean(values)
    if x.size < int(min_years):
        log.debug("not enough years for a GEV fit (%d < %d)", x.size, min_years)
        return None
    if np.allclose(x, x[0]):
        return None
    distribution = str(distribution).lower()
    from scipy import stats

    candidates = []
    if distribution == "gumbel":
        try:
            loc, scale = stats.gumbel_r.fit(x, method="MLE")
            candidates.append(("gumbel", 0.0, loc, scale, "mle"))
        except Exception:
            pass
        if not candidates:
            loc = float(np.mean(x))
            scale = float(np.std(x, ddof=1) * np.sqrt(6) / np.pi)
            candidates.append(("gumbel", 0.0, loc, scale, "moments"))
    else:
        for meth in ([str(method).lower()] if str(method).lower() != "mle" else []) + ["mle", "mm"]:
            if any(c[4] == meth for c in candidates):
                continue
            try:
                c, loc, scale = stats.genextreme.fit(x, method=meth)
                if np.isfinite([c, loc, scale]).all() and scale > 0:
                    candidates.append(("gev", float(c), float(loc), float(scale), meth))
                    break
            except Exception as exc:
                log.debug("GEV %s fit failed: %s", meth, exc)
        if not candidates:
            loc = float(np.mean(x))
            scale = float(np.std(x, ddof=1) * np.sqrt(6) / np.pi)
            candidates.append(("gev", 0.0, loc, scale, "moments"))

    name, shape, loc, scale, used = candidates[0]
    if not np.isfinite(scale) or scale <= 0:
        return None
    fit = Fit(distribution=name, shape=shape, loc=loc, scale=scale, n=int(x.size), method=used)
    try:
        fit.loglik = float(np.sum(np.log(np.maximum(_pdf(name, x, shape, loc, scale), 1e-300))))
        k = 2 if name == "gumbel" else 3
        fit.aic = float(-2 * fit.loglik + 2 * k)
    except Exception:
        pass
    try:
        ks = stats.kstest(x, lambda v: _cdf(name, v, shape, loc, scale))
        fit.ks_stat, fit.ks_p = float(ks.statistic), float(ks.pvalue)
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
