"""Publication-quality figures.

Six figure types are produced by the default suite:

``timeseries``  observed history + multi-model future (members, ensemble mean,
                percentile bands, linear trends in units/decade)
``gev``         GEV return-level curves for history and each scenario
``combined``    3x2 panel: time series and GEV side by side
``maps``        period-mean change maps with stippling and region boxes
``kde``         observed vs raw vs bias-corrected distributions, one panel per
                model, with KS statistics
``shape_box``   GEV shape parameter by scenario + Kruskal-Wallis p

Every function takes a :class:`~cmipkit.config.Config` plus small dataclasses
defined here, writes one file per format and returns the list of paths.  Nothing
in this module reads data from disk: the pipeline hands it ready-made arrays.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib

if not os.environ.get("DISPLAY") and os.environ.get("MPLBACKEND") is None:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.patches import Rectangle

from .config import Config

log = logging.getLogger("cmipkit.plotting")

_HAS_CARTOPY = True
try:  # cartopy is optional: the map figure degrades to a plain axes
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
except Exception:  # pragma: no cover
    _HAS_CARTOPY = False

CJK_FONTS = ["WenQuanYi Micro Hei", "WenQuanYi Zen Hei", "Noto Sans CJK SC", "Source Han Sans SC", "AR PL UMing CN"]


# ---------------------------------------------------------------------------
# data containers
# ---------------------------------------------------------------------------
@dataclass
class HistoricalSeries:
    years: np.ndarray
    values: np.ndarray
    trend: Dict[str, float] = field(default_factory=dict)
    label: str = "Historical"


@dataclass
class ScenarioSeries:
    scenario: str
    label: str
    years: np.ndarray
    matrix: np.ndarray
    center: np.ndarray
    bands: Dict[int, np.ndarray]
    trend: Dict[str, float] = field(default_factory=dict)
    n_models: int = 0


@dataclass
class GEVScenario:
    scenario: str
    label: str
    periods: np.ndarray
    median: np.ndarray
    bands: Dict[int, np.ndarray]
    n_models: int = 0
    reference_period: float = float("nan")
    empirical: Optional[Tuple[np.ndarray, np.ndarray]] = None


@dataclass
class GEVHistorical:
    periods: np.ndarray
    levels: np.ndarray
    reference_period: float = float("nan")
    empirical: Optional[Tuple[np.ndarray, np.ndarray]] = None


@dataclass
class MapPanel:
    label: str
    field: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    stipple: Optional[np.ndarray] = None


# ---------------------------------------------------------------------------
# style
# ---------------------------------------------------------------------------
def _has_cjk(text: str) -> bool:
    return any("\u4e00" <= ch <= "\u9fff" for ch in str(text))


def apply_style(cfg: Config, extra_text: str = "") -> None:
    """Consistent publication style; switches to a CJK font when needed."""
    scale = float(cfg.get_path("plots.figsize_scale", 1.0) or 1.0)
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": int(cfg.get_path("plots.dpi", 300)),
            "savefig.bbox": "tight",
            "font.size": 10 * scale,
            "axes.labelsize": 11 * scale,
            "axes.titlesize": 12 * scale,
            "xtick.labelsize": 9 * scale,
            "ytick.labelsize": 9 * scale,
            "legend.fontsize": 8.5 * scale,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": False,
            "figure.facecolor": "white",
        }
    )
    font = cfg.get_path("plots.font")
    language = str(cfg.get_path("plots.language", "en")).lower()
    if font:
        plt.rcParams["font.sans-serif"] = [font, "DejaVu Sans"]
    elif language.startswith("zh") or _has_cjk(extra_text):
        available = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
        picks = [f for f in CJK_FONTS if f in available]
        plt.rcParams["font.sans-serif"] = picks + ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False


def _units(cfg: Config) -> str:
    label = cfg.get_path("index.units_label")
    if label:
        return str(label).strip()
    return cfg.index_units.strip()


def _var_label(cfg: Config) -> str:
    return str(cfg.get_path("variable.long_name", "") or cfg.varname)


def _panel(ax, letter: str, x: float = -0.02, y: float = 1.02, size: Optional[float] = None) -> None:
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=size or 12,
            fontweight="bold", va="bottom", ha="right")


def _save(fig, outbase: str, cfg: Config) -> List[str]:
    paths = []
    dpi = int(cfg.get_path("plots.dpi", 300))
    formats = cfg.get_path("plots.formats", ["png"]) or ["png"]
    os.makedirs(os.path.dirname(os.path.abspath(outbase)) or ".", exist_ok=True)
    for fmt in formats:
        path = f"{outbase}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def _robust_levels(field: np.ndarray, n: int = 21, symmetric: bool = True) -> np.ndarray:
    data = np.asarray(field, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return np.linspace(-1, 1, n)
    if symmetric:
        vmax = float(np.nanpercentile(np.abs(data), 98))
        vmax = vmax if vmax > 0 else 1.0
        return np.linspace(-vmax, vmax, n)
    lo, hi = np.nanpercentile(data, [2, 98])
    if lo == hi:
        hi = lo + 1.0
    return np.linspace(float(lo), float(hi), n)


def _nice_ticks(ax, axis: str = "x", nbins: int = 6) -> None:
    locator = mticker.MaxNLocator(nbins=nbins, prune=None)
    (ax.xaxis if axis == "x" else ax.yaxis).set_major_locator(locator)


# ---------------------------------------------------------------------------
# time series
# ---------------------------------------------------------------------------
def _plot_series_panel(
    ax,
    cfg: Config,
    hist: Optional[HistoricalSeries],
    fut: ScenarioSeries,
    panel_letter: Optional[str] = None,
    show_legend: bool = True,
) -> None:
    color = cfg.color(fut.scenario)
    ref_year = cfg.get_path("plots.timeseries.reference_year")
    ref_value = cfg.get_path("plots.timeseries.reference_value")

    if hist is not None and hist.values.size:
        ax.plot(hist.years, hist.values, "-", color=cfg.get_path("plots.colors.historical", "#1f4e79"),
                lw=1.8, label=f"Observed ({int(hist.years[0])}–{int(hist.years[-1])})", zorder=6)
        if cfg.get_path("plots.timeseries.show_trend", True):
            tr = hist.trend or {}
            if tr.get("slope") is not None and np.isfinite(tr.get("slope", np.nan)):
                years = np.asarray(hist.years, dtype=float)
                line = tr["slope"] * (years - years.mean()) + np.nanmean(hist.values)
                unit = _units(cfg)
                ax.plot(years, line, "--", color=cfg.get_path("plots.colors.historical", "#1f4e79"),
                        lw=1.2, alpha=0.75, zorder=5,
                        label=f"Hist. trend {tr['slope_per_decade']:+.3f} {unit}/10a")

    if ref_year is not None and ref_value is not None:
        ax.plot([float(ref_year)], [float(ref_value)], marker="s", ms=7, color="#2ca02c",
                mec="k", mew=0.6, zorder=12, label=f"{int(float(ref_year))} observed")

    if cfg.get_path("plots.timeseries.show_members", True) and fut.matrix.size:
        yrs = np.tile(fut.years, fut.matrix.shape[0])
        vals = fut.matrix.reshape(-1)
        counts = np.sum(np.isfinite(fut.matrix), axis=0)
        count_label = str(int(counts[0])) if np.all(counts == counts[0]) else f"{int(counts.min())}–{int(counts.max())}"
        ax.scatter(yrs, vals, s=6, color=cfg.get_path("plots.colors.models", "#9e9e9e"),
                   alpha=0.45, linewidths=0, zorder=1, label=f"Models (available n={count_label})")

    if cfg.get_path("plots.timeseries.show_bands", True):
        need = set(cfg.get_path("ensemble.bands", [10, 90, 20, 80]))
        if {10, 90} <= need and 10 in fut.bands and 90 in fut.bands:
            ax.fill_between(fut.years, fut.bands[10], fut.bands[90], color=color, alpha=0.15,
                            lw=0, label="10–90% range", zorder=2)
        if {20, 80} <= need and 20 in fut.bands and 80 in fut.bands:
            ax.fill_between(fut.years, fut.bands[20], fut.bands[80], color=color, alpha=0.30,
                            lw=0, label="20–80% range", zorder=3)

    ax.plot(fut.years, fut.center, "-", color=color, lw=2.2,
            label=f"Multi-model {'median' if cfg.get_path('ensemble.statistic') == 'median' else 'mean'}",
            zorder=8)

    if cfg.get_path("plots.timeseries.show_trend", True):
        tr = fut.trend or {}
        if tr.get("slope") is not None and np.isfinite(tr.get("slope", np.nan)):
            years = np.asarray(fut.years, dtype=float)
            ok = np.isfinite(fut.center)
            if ok.sum() > 2:
                line = tr["slope"] * (years - years[ok].mean()) + np.nanmean(fut.center[ok])
                ax.plot(years, line, "--", color=color, lw=1.4, alpha=0.9, zorder=9,
                        label=f"Fut. trend {tr['slope_per_decade']:+.3f} {_units(cfg)}/10a")

    if (cfg.get_path("index.anomaly") or cfg.get_path("index.standardize")
            or cfg.get_path("index.type") in ("region_mean_anomaly", "region_diff", "projection_index")):
        ax.axhline(0, color="0.55", ls="--", lw=0.8, zorder=0)
    if ref_value is not None and float(ref_value) != 0:
        ax.axhline(float(ref_value), color="#6a3d9a", ls="-", lw=1.1, alpha=0.7, zorder=4)

    ax.set_ylabel(f"{_var_label(cfg)} ({_units(cfg)})".strip())
    ax.grid(alpha=0.28, lw=0.6)
    counts = np.sum(np.isfinite(fut.matrix), axis=0) if fut.matrix.size else np.array([0])
    available = str(int(counts[0])) if np.all(counts == counts[0]) else f"{int(counts.min())}–{int(counts.max())}"
    ax.set_title(f"{fut.label} · available models: {available}", loc="left", fontweight="bold")
    if panel_letter:
        _panel(ax, panel_letter)
    if show_legend:
        ax.legend(loc="upper left", ncol=1, handlelength=1.6, borderaxespad=0.3)
    xlim = cfg.get_path("plots.timeseries.xlim")
    if xlim:
        ax.set_xlim(*xlim)
    ylim = cfg.get_path("plots.timeseries.ylim")
    if ylim:
        ax.set_ylim(*ylim)


def timeseries(
    cfg: Config,
    hist: Optional[HistoricalSeries],
    scenarios: Sequence[ScenarioSeries],
    outbase: str,
) -> List[str]:
    """One column, one row per scenario, shared x axis."""
    if not scenarios:
        return []
    n = len(scenarios)
    fig, axes = plt.subplots(n, 1, figsize=(8.2, 2.9 * n + 0.6), sharex=True, squeeze=False)
    axes = axes[:, 0]
    letters = "abcdefghij"
    for i, fut in enumerate(scenarios):
        _plot_series_panel(axes[i], cfg, hist, fut, panel_letter=letters[i], show_legend=(i == 0))
        if i < n - 1:
            axes[i].tick_params(labelbottom=False)
    axes[-1].set_xlabel("Year")
    _nice_ticks(axes[-1], "x", 8)
    fig.tight_layout(h_pad=0.8)
    return _save(fig, outbase, cfg)


# ---------------------------------------------------------------------------
# GEV
# ---------------------------------------------------------------------------
def _plot_gev_panel(ax, cfg: Config, hist: Optional[GEVHistorical], fut: GEVScenario,
                    panel_letter: Optional[str] = None, show_legend: bool = True) -> None:
    color = cfg.color(fut.scenario)
    if hist is not None and hist.levels.size:
        ax.plot(hist.periods, hist.levels, "-", color=cfg.get_path("plots.colors.observation", "#111111"),
                lw=2.0, label="Historical fit", zorder=5)
        if hist.empirical is not None and len(hist.empirical[0]):
            ax.plot(*hist.empirical, "o", ms=3, mfc="none", mec="0.35", mew=0.6, label="Empirical", zorder=4)

    need = set(cfg.get_path("ensemble.bands", [10, 90, 20, 80]))
    if {10, 90} <= need and 10 in fut.bands and 90 in fut.bands:
        ax.fill_between(fut.periods, fut.bands[10], fut.bands[90], color=color, alpha=0.15, lw=0,
                        label="10–90% range")
    if {20, 80} <= need and 20 in fut.bands and 80 in fut.bands:
        ax.fill_between(fut.periods, fut.bands[20], fut.bands[80], color=color, alpha=0.30, lw=0,
                        label="20–80% range")
    ax.plot(fut.periods, fut.median, "-", color=color, lw=2.2,
            label=f"Multi-model median ({fut.n_models})")
    if fut.empirical is not None and len(fut.empirical[0]):
        ax.plot(*fut.empirical, "s", ms=3, mfc="none", mec=color, mew=0.6, alpha=0.8)

    ref_value = cfg.get_path("plots.timeseries.reference_value")
    if ref_value is not None:
        ax.axhline(float(ref_value), color="#6a3d9a", ls="-", lw=1.1, alpha=0.7)
        if np.isfinite(fut.reference_period):
            ax.plot([fut.reference_period], [float(ref_value)], "*", ms=9, color="#6a3d9a", zorder=10)
    for t in (10, 100):
        ax.axvline(t, color="0.7", ls=":", lw=0.7, zorder=0)

    ax.set_xscale("log")
    xlim = cfg.get_path("plots.gev.xlim", [1.1, 1100])
    ax.set_xlim(*xlim)
    fmt = mticker.FuncFormatter(lambda v, _: f"{v:g}")
    ax.xaxis.set_major_formatter(fmt)
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xlabel("Return period (years)")
    ax.set_ylabel(f"{_var_label(cfg)} ({_units(cfg)})".strip())
    ax.grid(alpha=0.28, lw=0.6, which="both")
    ax.set_title(fut.label, loc="left", fontweight="bold")
    if panel_letter:
        _panel(ax, panel_letter)
    ylim = cfg.get_path("plots.gev.ylim")
    if ylim:
        ax.set_ylim(*ylim)
    if show_legend:
        ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.3)


def gev(cfg: Config, hist: Optional[GEVHistorical], scenarios: Sequence[GEVScenario], outbase: str) -> List[str]:
    if not scenarios:
        return []
    n = len(scenarios)
    fig, axes = plt.subplots(n, 1, figsize=(6.4, 3.0 * n + 0.5), sharex=True, squeeze=False)
    axes = axes[:, 0]
    letters = "abcdefghij"
    for i, fut in enumerate(scenarios):
        _plot_gev_panel(axes[i], cfg, hist, fut, panel_letter=letters[i], show_legend=(i == 0))
        if i < n - 1:
            axes[i].tick_params(labelbottom=False)
    fig.tight_layout(h_pad=0.8)
    return _save(fig, outbase, cfg)


def scenarios_combined(
    cfg: Config,
    hist: Optional[HistoricalSeries],
    series: Sequence[ScenarioSeries],
    gev_hist: Optional[GEVHistorical],
    gev_scen: Sequence[GEVScenario],
    outbase: str,
) -> List[str]:
    """3x2 (or n x 2) panel: left time series, right GEV."""
    n = max(len(series), len(gev_scen))
    if n == 0:
        return []
    fig, axes = plt.subplots(n, 2, figsize=(14.5, 3.6 * n + 0.8), sharex="col")
    axes = np.atleast_2d(axes)
    letters_l = "acegikmo"
    letters_r = "bdfhjlnp"
    for i in range(n):
        if i < len(series):
            _plot_series_panel(axes[i, 0], cfg, hist, series[i], panel_letter=letters_l[i], show_legend=(i == 0))
        else:
            axes[i, 0].axis("off")
        if i < len(gev_scen):
            _plot_gev_panel(axes[i, 1], cfg, gev_hist, gev_scen[i], panel_letter=letters_r[i], show_legend=(i == 0))
        else:
            axes[i, 1].axis("off")
        if i < n - 1:
            axes[i, 0].tick_params(labelbottom=False)
            axes[i, 1].tick_params(labelbottom=False)
    axes[-1, 0].set_xlabel("Year")
    _nice_ticks(axes[-1, 0], "x", 8)
    fig.tight_layout(h_pad=0.6, w_pad=1.6)
    return _save(fig, outbase, cfg)


# ---------------------------------------------------------------------------
# maps
# ---------------------------------------------------------------------------
def _map_axes(fig, rect, extent, projection: str = "PlateCarree"):
    if _HAS_CARTOPY:
        proj = ccrs.PlateCarree() if projection == "PlateCarree" else ccrs.Robinson()
        ax = fig.add_axes(rect, projection=proj)
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        try:
            ax.add_feature(cfeature.COASTLINE.with_scale("110m"), lw=0.5, edgecolor="0.25")
            ax.add_feature(cfeature.BORDERS.with_scale("110m"), lw=0.35, edgecolor="0.45", alpha=0.7)
        except Exception as exc:  # coastlines may need a download
            log.warning("cartopy features unavailable (%s); drawing data only", exc)
        ax.set_xticks(np.arange(-180, 181, 30), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(-90, 91, 15), crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.tick_params(labelsize=8)
        return ax, ccrs.PlateCarree()
    ax = fig.add_axes(rect)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    # Callers pass this explicitly to meshes, stippling and region boxes.
    # None would select display coordinates rather than longitude/latitude.
    return ax, ax.transData


def _region_boxes(ax, cfg: Config, transform=None) -> None:
    if not cfg.get_path("plots.maps.regions", True):
        return
    regions = []
    if cfg.get_path("index.region"):
        regions.append(cfg.get_path("index.region"))
    regions.extend(cfg.get_path("index.regions") or [])
    for reg in regions:
        if not reg:
            continue
        lat = sorted(reg["lat"])
        lon = sorted(reg["lon"])
        ax.add_patch(
            Rectangle((lon[0], lat[0]), lon[1] - lon[0], lat[1] - lat[0], fill=False,
                      edgecolor="k", lw=1.3, transform=transform, zorder=12)
        )
        if reg.get("name"):
            ax.text(lon[0] + 0.5, lat[1] + 0.5, reg["name"], fontsize=8, transform=transform, zorder=13,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))


def maps(
    cfg: Config,
    panels: Sequence[Sequence[MapPanel]],
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    outbase: str,
    cbar_label: str = "",
    symmetric: bool = True,
    extent: Optional[Sequence[float]] = None,
) -> List[str]:
    """Grid of maps; ``panels[r][c]`` is a :class:`MapPanel`.

    When cartopy is unavailable the axes fall back to plain matplotlib axes so
    the figure is still produced (without coastlines).
    """
    nrows, ncols = len(panels), len(panels[0]) if panels else 0
    if nrows == 0 or ncols == 0:
        return []
    lon = panels[0][0].lon
    lat = panels[0][0].lat
    if extent is None:
        extent = [float(np.min(lon)), float(np.max(lon)), float(np.min(lat)), float(np.max(lat))]

    cmap_name = cfg.get_path("plots.maps.cmap", "RdBu_r")
    levels_cfg = cfg.get_path("plots.maps.levels")
    all_vals = np.concatenate([p.field[np.isfinite(p.field)].ravel() for row in panels for p in row if p is not None])
    if levels_cfg:
        levels = np.asarray(levels_cfg if len(levels_cfg) > 2 else np.linspace(levels_cfg[0], levels_cfg[1],
                                                                              int(cfg.get_path("plots.maps.n_levels", 21))))
    else:
        levels = _robust_levels(all_vals, int(cfg.get_path("plots.maps.n_levels", 21)), symmetric)
    from matplotlib.colors import BoundaryNorm

    norm = BoundaryNorm(levels, ncolors=plt.get_cmap(cmap_name).N, extend="both")

    fig_w = 4.6 * ncols + 0.8
    fig_h = 3.4 * nrows + 1.2
    fig = plt.figure(figsize=(fig_w, fig_h))
    left, right = 0.06, 0.90
    top, bottom = 0.94, 0.10
    w = (right - left) / ncols
    h = (top - bottom) / nrows
    mesh = None
    for r in range(nrows):
        for c in range(ncols):
            rect = [left + c * w + 0.005, top - (r + 1) * h + 0.01, w - 0.02, h - 0.03]
            ax, transform = _map_axes(fig, rect, extent, cfg.get_path("plots.maps.projection", "PlateCarree"))
            panel = panels[r][c]
            if panel is None:
                ax.axis("off")
                continue
            mesh = ax.pcolormesh(panel.lon, panel.lat, panel.field, cmap=cmap_name, norm=norm,
                                 shading="auto", transform=transform, rasterized=True)
            if cfg.get_path("plots.maps.stipple", True) and panel.stipple is not None and panel.stipple.any():
                mlon, mlat = np.meshgrid(panel.lon, panel.lat)
                ax.scatter(mlon[panel.stipple], mlat[panel.stipple], s=2.5, c="k", marker=".",
                           alpha=float(cfg.get_path("plots.maps.stipple_alpha", 0.05)) * 8, linewidths=0,
                           transform=transform, zorder=9)
            _region_boxes(ax, cfg, transform)
            ax.set_title(panel.label, fontsize=10, loc="left")
            if c == 0 and r < len(row_labels):
                ax.text(-0.06, 0.5, row_labels[r], transform=ax.transAxes, rotation=90, va="center",
                        ha="right", fontsize=10, fontweight="bold")
    if mesh is not None:
        cax = fig.add_axes([left + 0.02, 0.045, right - left - 0.06, 0.018])
        cb = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="both", ticks=levels[:: max(1, len(levels) // 7)])
        cb.set_label(cbar_label or f"{_var_label(cfg)} anomaly ({_units(cfg)})")
        cb.ax.tick_params(labelsize=8)
    return _save(fig, outbase, cfg)


# ---------------------------------------------------------------------------
# KDE diagnostics
# ---------------------------------------------------------------------------
def kde(
    cfg: Config,
    records: Sequence[Dict[str, Any]],
    outbase: str,
) -> List[str]:
    """Grid of per-model distributions: observed vs raw vs corrected.

    ``records`` items: ``{"model": str, "obs": arr, "raw": arr, "corr": arr,
    "ks_raw": float, "ks_corr": float}``.
    """
    from scipy.stats import gaussian_kde

    if not records:
        return []
    max_models = int(cfg.get_path("plots.kde.max_models", 24))
    records = records[:max_models]
    grid = cfg.get_path("plots.kde.grid", [4, 6]) or [4, 6]
    max_cols = max(1, int(grid[1]))
    ncols = max(1, min(len(records), max_cols))
    nrows = int(np.ceil(len(records) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 2.6 * nrows), sharex=True, sharey=True,
                             squeeze=False)
    flat = axes.ravel()
    handles = None
    for i, rec in enumerate(records):
        ax = flat[i]
        datasets = [("Obs", rec["obs"], cfg.get_path("plots.colors.observation", "#111111")),
                    ("Raw", rec["raw"], "#d62728"),
                    ("Corrected", rec["corr"], "#1f77b4")]
        allv = np.concatenate([np.asarray(d, float)[np.isfinite(d)] for _, d, _ in datasets])
        if allv.size < 3 or np.nanstd(allv) == 0:
            ax.axis("off")
            continue
        xs = np.linspace(np.nanmin(allv), np.nanmax(allv), 200)
        lines = []
        for label, data, color in datasets:
            arr = np.asarray(data, float)
            arr = arr[np.isfinite(arr)]
            if arr.size < 3:
                lines.append(None)
                continue
            try:
                y = gaussian_kde(arr)(xs)
            except Exception:
                lines.append(None)
                continue
            (line,) = ax.plot(xs, y, color=color, lw=1.4, label=label)
            lines.append(line)
        if handles is None:
            handles = [ln for ln in lines if ln is not None]
        txt = f"KS raw {rec.get('ks_raw', float('nan')):.3f}\nKS corr {rec.get('ks_corr', float('nan')):.3f}"
        ax.text(0.96, 0.95, txt, transform=ax.transAxes, ha="right", va="top", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.85))
        ax.set_title(rec["model"], fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.2, lw=0.5)
        if i == 0:
            _panel(ax, "a", x=-0.12, y=1.05, size=10)
    for j in range(len(records), nrows * ncols):
        flat[j].axis("off")
    for r in range(nrows):
        flat[r * ncols].set_ylabel("Density")
    for c in range(ncols):
        flat[(nrows - 1) * ncols + c].set_xlabel(_units(cfg) or "value")
    if handles:
        fig.legend(handles, [h.get_label() for h in handles], loc="upper center", ncol=3,
                   bbox_to_anchor=(0.5, 1.0), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=0.6, w_pad=0.4)
    return _save(fig, outbase, cfg)


# ---------------------------------------------------------------------------
# shape parameter box plot
# ---------------------------------------------------------------------------
def shape_box(
    cfg: Config,
    shapes: Dict[str, Sequence[float]],
    labels: Dict[str, str],
    outbase: str,
    kruskal: Optional[Dict[str, float]] = None,
) -> List[str]:
    keys = [k for k in labels if len(shapes.get(k, [])) > 0]
    if not keys:
        return []
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    data = [np.asarray(shapes[k], float) for k in keys]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55, showmeans=True,
                    meanprops=dict(marker="^", mfc="white", mec="k", ms=5),
                    medianprops=dict(color="k", lw=1.2))
    for patch, key in zip(bp["boxes"], keys):
        patch.set_facecolor(cfg.color(key))
        patch.set_alpha(0.45)
        patch.set_edgecolor("k")
    ax.axhline(0, color="0.6", ls="--", lw=0.8)
    ax.set_xticks(range(1, len(keys) + 1))
    ax.set_xticklabels([labels[k] for k in keys])
    ax.set_ylabel("GEV shape parameter $c$")
    ax.set_title("GEV shape parameter by scenario", loc="left", fontweight="bold")
    if kruskal and np.isfinite(kruskal.get("p_value", np.nan)):
        ax.text(0.98, 0.03, f"Kruskal–Wallis p = {kruskal['p_value']:.3f}",
                transform=ax.transAxes, ha="right", fontsize=9)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    return _save(fig, outbase, cfg)

# ---------------------------------------------------------------------------
# single-model PDF / projection figures
# ---------------------------------------------------------------------------
def pdf_single(cfg: Config, record: Dict[str, Any], outbase: str, subtitle: str = "") -> List[str]:
    """One probability-density figure for one model: observed vs raw vs corrected."""
    from scipy.stats import gaussian_kde

    datasets = [
        ("Observed", np.asarray(record["obs"], float), cfg.get_path("plots.colors.observation", "#111111")),
        ("Model raw", np.asarray(record["raw"], float), "#d62728"),
        ("Model corrected", np.asarray(record["corr"], float), "#1f77b4"),
    ]
    allv = np.concatenate([d[np.isfinite(d)] for _, d, _ in datasets])
    if allv.size < 5 or np.nanstd(allv) == 0:
        return []
    lo, hi = np.nanpercentile(allv, [0.2, 99.8])
    pad = 0.05 * (hi - lo if hi > lo else 1.0)
    xs = np.linspace(lo - pad, hi + pad, 400)

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    for label, data, color in datasets:
        arr = data[np.isfinite(data)]
        if arr.size < 3:
            continue
        try:
            density = gaussian_kde(arr)(xs)
        except Exception:
            continue
        ax.plot(xs, density, color=color, lw=2.0, label=f"{label} (n={arr.size})")
        ax.fill_between(xs, 0, density, color=color, alpha=0.06, lw=0)

    ks_raw = record.get("ks_raw", float("nan"))
    ks_corr = record.get("ks_corr", float("nan"))
    ax.text(0.97, 0.96, f"KS$_{{raw}}$ = {ks_raw:.3f}\nKS$_{{corr}}$ = {ks_corr:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.9))
    ax.set_xlabel(f"{_var_label(cfg)} ({_units(cfg)})".strip())
    ax.set_ylabel("Probability density")
    title = record.get("model", "")
    if subtitle:
        title = f"{title} — {subtitle}"
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    return _save(fig, outbase, cfg)


def per_model_timeseries(
    cfg: Config,
    hist: Optional[HistoricalSeries],
    series: Dict[str, ScenarioSeries],
    ensemble: Optional[Dict[str, ScenarioSeries]],
    model: str,
    outbase: str,
) -> List[str]:
    """Future projection of one individual model, all scenarios in one panel."""
    from .indices import linear_trend

    if not series:
        return []
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    if hist is not None and hist.values.size:
        ax.plot(hist.years, hist.values, "-", color=cfg.get_path("plots.colors.historical", "#1f4e79"),
                lw=1.8, label=f"Observed ({int(hist.years[0])}–{int(hist.years[-1])})", zorder=6)
    for sc, fut in series.items():
        color = cfg.color(sc)
        years = np.asarray(fut.years, dtype=float)
        trend = linear_trend(years, fut.center)
        ax.plot(years, fut.center, "-", color=color, lw=2.0, zorder=5,
                label=f"{fut.label}  ({trend['slope_per_decade']:+.3f} {_units(cfg)}/10a)")
    for sc, fut in (ensemble or {}).items():
        if sc in series:
            ax.plot(fut.years, fut.center, "--", color="0.45", lw=1.1, zorder=3,
                    label="Ensemble mean" if sc == next(iter(series)) else None)
    if hist is not None and hist.values.size and series:
        first = next(iter(series.values()))
        ax.plot([hist.years[-1], first.years[0]], [hist.values[-1], first.center[0]],
                ":", color="0.4", lw=1.0)
    ax.axhline(0, color="0.55", ls="--", lw=0.8, zorder=0)
    ax.set_xlabel("Year")
    ax.set_ylabel(f"{_var_label(cfg)} ({_units(cfg)})".strip())
    ax.set_title(f"{model}", loc="left", fontweight="bold")
    ax.grid(alpha=0.28, lw=0.6)
    ax.legend(loc="upper left", fontsize=8.5, frameon=False)
    xlim = cfg.get_path("plots.timeseries.xlim")
    if xlim:
        ax.set_xlim(*xlim)
    ylim = cfg.get_path("plots.timeseries.ylim")
    if ylim:
        ax.set_ylim(*ylim)
    fig.tight_layout()
    return _save(fig, outbase, cfg)
