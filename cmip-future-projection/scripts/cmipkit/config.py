"""Configuration handling for cmipkit.

A run is fully described by one YAML file.  The loader deep-merges the user file
on top of :data:`DEFAULTS`, resolves ``{name}``/``{root}`` placeholders inside
string values, and validates the result before any heavy work starts.

Design notes
------------
* Everything is a plain nested ``dict`` (wrapped in :class:`Config`) so that an
  agent can dump / patch any part without fighting dataclass plumbing.
* Unknown keys are *reported* (warnings) but never fatal, so the config format
  can grow without breaking old files.
* ``Config.get_path("a.b.c", default)`` is the only accessor used by the rest of
  the package; it keeps the pipeline code short.
"""

from __future__ import annotations

import copy
import json
import logging
import os
from typing import Any, Dict, Iterable, List, Optional

import yaml

log = logging.getLogger("cmipkit.config")

SEASONS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ANN": list(range(1, 13)),
}

AGG_STATS = (
    "mean",
    "sum",
    "max",
    "min",
    "median",
    "std",
    "p05",
    "p10",
    "p90",
    "p95",
    "p99",
    "count_gt",
    "count_lt",
    "count_ge",
    "count_le",
    "max_abs",
    "mean_abs",
)

DEFAULTS: Dict[str, Any] = {
    "project": {
        "name": "cmip_projection",
        "outdir": "runs/{name}",
        "description": "",
    },
    "variable": {
        "name": "variable",
        "units": "",
        "long_name": "",
    },
    "data": {
        "reference": {
            "kind": "grid",  # grid | series
            "paths": [],
            "years": None,
            "variable": None,
            "units_in": None,
            "units_transform": {"scale": 1.0, "offset": 0.0},
            "plev": None,
            "bbox": None,
            "aggregate": {"months": "ANN", "stat": "mean", "threshold": None},
            "series_file": None,
            "series_columns": {"year": "year", "value": "value"},
        },
        "models": {
            "kind": "grid",  # grid | series
            "root": None,
            "pattern": None,  # single pattern (legacy); prefer `patterns`
            "patterns": [],  # empty -> built-in CMIP6 candidates (flat and per-variable sub-directory)
            "variable": None,
            "units_in": None,
            "units_transform": {"scale": 1.0, "offset": 0.0},
            "plev": None,
            "bbox": None,
            "historical_experiment": "historical",
            "scenarios": ["ssp245"],
            "scenario_labels": {},
            "members": ["r1i1p1f1"],
            "member_combine": "first",  # first | mean | median (multi-member ensemble)
            "models": "all",
            "exclude": [],
            "aggregate": {"months": "ANN", "stat": "mean", "threshold": None},
                "series_dir": None,
        },
    },
    "grid": {
        "lat": None,  # {"start": .., "stop": .., "step": ..} or explicit list
        "lon": None,
        "method": "linear",  # linear | nearest
        "lon_360": True,
    },
    "periods": {
        "historical": [1950, 2014],
        "future": [2026, 2100],
        "reference_climate": [1950, 2014],
        "future_windows": [[2026, 2045], [2046, 2065], [2066, 2085], [2086, 2100]],
    },
    "index": {
        "type": "region_mean",  # region_mean | region_mean_anomaly | region_diff | projection_index | series
        "region": None,  # {"lat": [..], "lon": [..], "name": "..."}
        "regions": None,  # for region_diff: [reference_region, target_region]
        "remove_zonal_mean": False,
        "anomaly": False,
        "detrend": False,
        "standardize": False,
        "projection": {
            "template_source": "reference",  # reference | model | year
            "template_period": None,  # defaults to periods.reference_climate
            "template_year": None,
            "template_model": None,
        },
        "units_label": None,
    },
    "correction": {
        "method": "mean_std",  # none | mean | mean_std | delta_change | quantile_mapping | quantile_delta_mapping
        "params_period": None,  # defaults to periods.reference_climate
        "apply_to": ["historical", "future"],
        "multiplicative": "auto",  # true | false | auto
        "lower_bound": None,
        "upper_bound": None,
        "n_quantiles": 100,
        "trend_correction": {
            "enabled": False,
            "mode": "difference",  # difference | ratio
        },
    },
    "ensemble": {
        "statistic": "mean",  # mean | median
        "bands": [10, 90, 20, 80],
        "weights": "equal",  # equal | skill
        "exclude_models": [],
        "min_models": 3,
    },
    "extremes": {
        "enabled": True,
        "source": "future",  # future | historical | both
        "fit_window": None,  # {"future": [2026, 2100], "historical": [1950, 2024]}
        "distribution": "gev",  # gev | gumbel
        "fit_method": "mle",  # mle | pwm (method of moments fallback is automatic)
        "return_periods": [1.1, 2, 5, 10, 20, 30, 50, 100, 200, 500, 1000],
        "reference_value": None,
        "bootstrap": 0,
        "show_empirical": False,
        "min_years": 20,
        "shape_comparison": True,
    },
    "daily": {
        "enabled": False,
        "months": [6, 7, 8],  # JJA
        "index": "region_mean",  # region_mean | region_mean_anomaly | region_diff
        "baseline": None,  # defaults to periods.reference_climate
        "correction": "quantile_delta_mapping",  # none | mean | mean_std | quantile_mapping | quantile_delta_mapping
        "extremes": {
            "kind": "pot",  # pot | block
            "threshold_mode": "observed",  # observed | quantile | absolute
            "threshold_quantile": 0.95,
            "threshold_abs": None,
            "decluster_days": 3,
            "min_peaks": 30,
            "min_years": 20,
            "distribution": "gev",  # for kind: block
            "return_periods": [2, 5, 10, 20, 50, 100],
        },
        "plots": True,
    },
    "plots": {
        "dpi": 300,
        "formats": ["png"],
        "figsize_scale": 1.0,
        "font": None,
        "language": "en",  # en | zh (only affects automatic titles/labels)
        "colors": {
            "historical": "#1f4e79",
            "observation": "#111111",
            "models": "#9e9e9e",
            "scenario": {"ssp126": "#1b7837", "ssp245": "#e08214", "ssp370": "#b2182b", "ssp585": "#d73027"},
        },
        "figures": ["timeseries", "gev", "combined", "maps", "kde", "shape_box"],
        "timeseries": {
            "show_members": True,
            "show_bands": True,
            "show_trend": True,
            "reference_year": 2025,
            "reference_value": None,
            "ylim": None,
            "xlim": None,
            "per_model": False,   # also write one projection figure per model
        },
        "gev": {"xlim": [1.1, 1100], "ylim": None, "show_empirical": False},
        "maps": {
            "variable": "change",  # change | absolute | anomaly
            "cmap": "RdBu_r",
            "levels": None,  # [vmin, vmax] or list of levels
            "n_levels": 21,
            "stipple": True,
            "stipple_alpha": 0.05,
            "agreement": 0.8,
            "regions": True,
            "extent": None,
            "projection": "PlateCarree",
        },
        "kde": {
            "grid": [4, 6],
            "max_models": 24,
            "mode": "auto",       # auto | field | series
            "per_model": True,    # also write one PDF figure per model
        },
    },
    "runtime": {
        "verbose": True,
        "cache": True,
        "strict": False,
    },
}

_SCENARIO_ALIASES = {
    "ssp126": "SSP1-2.6",
    "ssp245": "SSP2-4.5",
    "ssp370": "SSP3-7.0",
    "ssp585": "SSP5-8.5",
    "ssp119": "SSP1-1.9",
    "ssp434": "SSP4-3.4",
    "ssp460": "SSP4-6.0",
    "historical": "Historical",
}


class ConfigError(ValueError):
    """Raised when a configuration file cannot be used."""


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    out = copy.deepcopy(base)
    for key, val in (override or {}).items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], val)
        else:
            out[key] = copy.deepcopy(val)
    return out


def _iter_leaf_paths(node: Any, prefix: str = "") -> Iterable[str]:
    if isinstance(node, dict):
        if not node and prefix:
            # an empty dict accepts an arbitrary sub-tree (e.g. scenario_labels)
            yield prefix
        for key, val in node.items():
            yield from _iter_leaf_paths(val, f"{prefix}.{key}" if prefix else str(key))
    else:
        yield prefix


class Config(dict):
    """Nested dict + dotted accessor + scenario/period helpers."""

    # -- access ---------------------------------------------------------
    def get_path(self, dotted: str, default: Any = None) -> Any:
        node: Any = self
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def set_path(self, dotted: str, value: Any) -> None:
        parts = dotted.split(".")
        node = self
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = value

    def require(self, dotted: str, message: str = "") -> Any:
        val = self.get_path(dotted)
        if val in (None, "", []):
            raise ConfigError(message or f"missing required configuration key: {dotted}")
        return val

    # -- shortcuts ------------------------------------------------------
    @property
    def name(self) -> str:
        return str(self.get_path("project.name", "cmip_projection"))

    @property
    def outdir(self) -> str:
        return str(self.get_path("project.outdir", "runs/" + self.name))

    @property
    def varname(self) -> str:
        return str(self.get_path("variable.name", "variable"))

    @property
    def units(self) -> str:
        return str(self.get_path("variable.units", "") or "")

    def scenario_list(self) -> List[str]:
        sc = self.get_path("data.models.scenarios", ["ssp245"])
        if isinstance(sc, str):
            sc = [sc]
        return list(sc)

    def scenario_label(self, scenario: str) -> str:
        labels = self.get_path("data.models.scenario_labels", {}) or {}
        if scenario in labels:
            return str(labels[scenario])
        return _SCENARIO_ALIASES.get(scenario, scenario.upper())

    def period(self, key: str) -> List[int]:
        val = self.get_path(f"periods.{key}")
        if val is None and key == "params_period":
            val = self.get_path("periods.reference_climate")
        if val is None:
            raise ConfigError(f"periods.{key} is not configured")
        if len(val) != 2:
            raise ConfigError(f"periods.{key} must be [start, end], got {val!r}")
        return [int(val[0]), int(val[1])]

    def future_windows(self) -> List[List[int]]:
        wins = self.get_path("periods.future_windows")
        if not wins:
            lo, hi = self.period("future")
            wins = [[lo, hi]]
        return [[int(a), int(b)] for a, b in wins]

    def months(self, which: str = "reference") -> List[int]:
        spec = self.get_path(f"data.{which}.aggregate.months", "ANN")
        return resolve_months(spec)

    def color(self, scenario: str) -> str:
        table = self.get_path("plots.colors.scenario", {}) or {}
        if scenario in table:
            return table[scenario]
        fallback = ["#1b7837", "#e08214", "#b2182b", "#d73027", "#4575b4", "#762a83"]
        scs = self.scenario_list()
        idx = scs.index(scenario) if scenario in scs else 0
        return fallback[idx % len(fallback)]

    # -- serialisation --------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return copy.deepcopy(dict(self))

    def dump(self, path: str) -> None:
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(self.to_dict(), fh, allow_unicode=True, sort_keys=False)


def resolve_months(spec: Any) -> List[int]:
    """Turn ``[9, 10]`` / ``"SON"`` / ``"ANN"`` into a list of month numbers."""
    if spec is None:
        return list(range(1, 13))
    if isinstance(spec, int):
        return [spec]
    if isinstance(spec, str):
        key = spec.strip().upper()
        if key in SEASONS:
            return list(SEASONS[key])
        if key in ("ALL", "YEAR", "ANNUAL"):
            return list(range(1, 13))
        # allow "9-10" or "9,10"
        parts = key.replace("-", ",").split(",")
        return sorted({int(p) for p in parts if p.strip()})
    return sorted({int(m) for m in spec})


class _KeepUnknown(dict):
    """Mapping that turns missing ``{placeholders}`` back into their own text."""

    def __missing__(self, key):
        return "{" + key + "}"


def safe_format(template: str, mapping: Dict[str, Any]) -> str:
    """``str.format`` that leaves unknown placeholders untouched."""
    import string

    try:
        return string.Formatter().vformat(template, (), _KeepUnknown(mapping))
    except (IndexError, ValueError):
        return template


def _resolve_placeholders(cfg: Dict[str, Any], env: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve ``{name}`` style placeholders inside strings (paths, pattern)."""

    def walk(node: Any) -> Any:
        if isinstance(node, dict):
            return {k: walk(v) for k, v in node.items()}
        if isinstance(node, list):
            return [walk(v) for v in node]
        if isinstance(node, str):
            return safe_format(node, env)
        return node

    return walk(cfg)


def _env_for_placeholders(cfg: Dict[str, Any]) -> Dict[str, Any]:
    root = cfg.get("data", {}).get("models", {}).get("root") or ""
    name = cfg.get("project", {}).get("name", "cmip_projection")
    # deliberately small: {model}, {experiment}, {member}, {scenario} and {year}
    # must stay unresolved so that the per-experiment file resolution still works
    return {
        "name": name,
        "root": root,
        "outdir": cfg.get("project", {}).get("outdir", "runs/" + name),
        "variable": (cfg.get("data", {}).get("models", {}).get("variable") or cfg.get("variable", {}).get("name", "var")),
    }


def load_config(path: str) -> Config:
    if not os.path.exists(path):
        raise ConfigError(f"config file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        user = yaml.safe_load(fh) or {}
    if not isinstance(user, dict):
        raise ConfigError(f"top level of {path} must be a mapping")

    merged = _deep_merge(DEFAULTS, user)
    if not merged.get("project", {}).get("outdir"):
        merged["project"]["outdir"] = "runs/" + str(merged["project"].get("name", "cmip_projection"))
    merged = _resolve_placeholders(merged, _env_for_placeholders(merged))
    cfg = Config(merged)
    cfg["_config_path"] = os.path.abspath(path)
    cfg["_config_dir"] = os.path.dirname(os.path.abspath(path))
    errors, warnings = validate_config(cfg, user_keys=set(_iter_leaf_paths(user)))
    for msg in warnings:
        log.warning(msg)
    if errors:
        raise ConfigError("invalid configuration:\n  - " + "\n  - ".join(errors))
    return cfg


def validate_config(cfg: Config, user_keys: Optional[set] = None):
    """Return ``(errors, warnings)``; never raise so callers can report all."""
    errors: List[str] = []
    warnings: List[str] = []

    for key in ("project.name", "variable.name"):
        if not cfg.get_path(key):
            errors.append(f"{key} must be set")

    # data.reference
    ref = cfg.get_path("data.reference", {}) or {}
    kind = ref.get("kind", "grid")
    if kind not in ("grid", "series"):
        errors.append(f"data.reference.kind must be 'grid' or 'series', got {kind!r}")
    if kind == "grid":
        if not ref.get("paths"):
            errors.append("data.reference.paths must list at least one file or glob pattern")
        if not ref.get("variable"):
            warnings.append("data.reference.variable is empty; the loader will auto-detect the first data variable")
        years = ref.get("years")
        if years is not None and (not isinstance(years, (list, tuple)) or len(years) != 2):
            errors.append("data.reference.years must be [start, end]")
    else:
        if not ref.get("series_file"):
            errors.append("data.reference.series_file is required when data.reference.kind == 'series'")

    # data.models
    mod = cfg.get_path("data.models", {}) or {}
    mkind = mod.get("kind", "grid")
    if mkind not in ("grid", "series"):
        errors.append(f"data.models.kind must be 'grid' or 'series', got {mkind!r}")
    if mkind == "grid":
        if not mod.get("root"):
            errors.append("data.models.root is required for grid input")
    else:
        if not mod.get("series_dir"):
            errors.append("data.models.series_dir is required when data.models.kind == 'series'")
    if not cfg.scenario_list():
        errors.append("data.models.scenarios must contain at least one scenario")
    for sc in cfg.scenario_list():
        if sc == "historical":
            errors.append("data.models.scenarios must not contain 'historical' (that is the baseline experiment)")

    # periods
    for key in ("historical", "future", "reference_climate"):
        try:
            lo, hi = cfg.period(key)
        except ConfigError as exc:
            errors.append(str(exc))
            continue
        if lo > hi:
            errors.append(f"periods.{key} start {lo} > end {hi}")
    try:
        if cfg.period("historical")[1] >= cfg.period("future")[0]:
            warnings.append("periods.historical overlaps periods.future; check the split year")
    except ConfigError:
        pass

    # index
    itype = cfg.get_path("index.type", "region_mean")
    valid_index = ("region_mean", "region_mean_anomaly", "region_diff", "projection_index", "series")
    if itype not in valid_index:
        errors.append(f"index.type must be one of {valid_index}, got {itype!r}")
    if itype in ("region_mean", "region_mean_anomaly", "projection_index", "region_diff"):
        if itype == "region_diff":
            regs = cfg.get_path("index.regions")
            if not regs or len(regs) != 2:
                errors.append("index.regions must list exactly two regions [reference, target] for region_diff")
        elif not cfg.get_path("index.region") and itype != "projection_index":
            errors.append(f"index.region is required for index.type == {itype!r}")

    # correction
    method = cfg.get_path("correction.method", "mean_std")
    valid_methods = ("none", "mean", "mean_std", "delta_change", "quantile_mapping", "quantile_delta_mapping")
    if method not in valid_methods:
        errors.append(f"correction.method must be one of {valid_methods}, got {method!r}")
    combine = cfg.get_path("data.models.member_combine", "first")
    if combine not in ("first", "mean", "median"):
        errors.append("data.models.member_combine must be 'first', 'mean' or 'median'")
    mult = cfg.get_path("correction.multiplicative", "auto")
    if mult not in ("auto", True, False):
        errors.append("correction.multiplicative must be true, false or 'auto'")
    lb = cfg.get_path("correction.lower_bound")
    if lb is not None and not isinstance(lb, (int, float)):
        errors.append("correction.lower_bound must be numeric or null")

    # daily extremes
    daily = cfg.get_path("daily", {}) or {}
    if daily.get("enabled"):
        if daily.get("index") not in ("region_mean", "region_mean_anomaly", "region_diff"):
            errors.append("daily.index must be region_mean, region_mean_anomaly or region_diff")
        if str(daily.get("correction", "none")) not in valid_methods:
            errors.append(f"daily.correction must be one of {valid_methods}")
        ext = daily.get("extremes") or {}
        kind = ext.get("kind", "pot")
        if kind not in ("pot", "block"):
            errors.append("daily.extremes.kind must be 'pot' or 'block'")
        if ext.get("threshold_mode", "observed") not in ("observed", "quantile", "absolute"):
            errors.append("daily.extremes.threshold_mode must be observed, quantile or absolute")
        if ext.get("threshold_mode") == "absolute" and ext.get("threshold_abs") is None:
            errors.append("daily.extremes.threshold_abs is required when threshold_mode is 'absolute'")
        if cfg.get_path("data.reference.kind") != "grid" or cfg.get_path("data.models.kind") != "grid":
            warnings.append("daily extremes need gridded input; the stage will be skipped")

    # plots
    figs = cfg.get_path("plots.figures", [])
    valid_figs = ("timeseries", "gev", "combined", "maps", "kde", "shape_box")
    for fig in figs:
        if fig not in valid_figs:
            errors.append(f"plots.figures contains unknown figure {fig!r}; valid: {valid_figs}")
    if "maps" in figs and mkind != "grid":
        warnings.append("plots.figures requests 'maps' but model data are series; the map figure will be skipped")

    # unknown user keys (a key is known if it, or one of its ancestors, is a
    # default leaf -- several defaults are ``None`` placeholders that accept an
    # arbitrary sub-dict, e.g. ``grid.lat: {start, stop, step}``)
    if user_keys:
        known = set(_iter_leaf_paths(DEFAULTS)) | {"_config_path", "_config_dir"}

        def is_known(key: str) -> bool:
            parts = key.split(".")
            return any(".".join(parts[:i]) in known for i in range(1, len(parts) + 1))

        for key in sorted(user_keys):
            if not is_known(key):
                warnings.append(f"unknown configuration key (ignored): {key}")
    return errors, warnings


def write_resolved(cfg: Config, path: str) -> None:
    """Dump the fully merged config next to the run for reproducibility."""
    payload = cfg.to_dict()
    payload.pop("_config_path", None)
    payload.pop("_config_dir", None)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(payload, fh, allow_unicode=True, sort_keys=False)


def write_manifest(cfg: Config, path: str, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    import platform
    import sys

    import numpy
    import scipy
    import xarray

    manifest = {
        "cmipkit": __import__("cmipkit").__version__,
        "project": cfg.name,
        "config": cfg.get_path("_config_path"),
        "created": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": {"numpy": numpy.__version__, "scipy": scipy.__version__, "xarray": xarray.__version__},
        "variable": cfg.get_path("variable.name"),
        "units": cfg.get_path("variable.units"),
        "scenarios": cfg.scenario_list(),
        "correction": cfg.get_path("correction.method"),
    }
    if extra:
        manifest.update(extra)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)
    return manifest
