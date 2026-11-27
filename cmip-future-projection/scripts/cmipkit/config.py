"""Configuration handling for cmipkit.

A run is fully described by one YAML file.  The loader deep-merges the user file
on top of :data:`DEFAULTS`, resolves ``{name}``/``{root}`` placeholders inside
string values, and validates the result before any heavy work starts.

Design notes
------------
* Everything is a plain nested ``dict`` (wrapped in :class:`Config`) so that an
  agent can dump / patch any part without fighting dataclass plumbing.
* Unknown keys are errors by default; migration may explicitly retain warnings.
* ``Config.get_path("a.b.c", default)`` is the only accessor used by the rest of
  the package; it keeps the pipeline code short.
"""

from __future__ import annotations

import copy
import json
import logging
import os
import re
import math
import tempfile
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
            "units_transform": {"scale": 1.0, "offset": 0.0, "units": None},
            "plev": None,
            "bbox": None,
            "aggregate": {"months": "ANN", "stat": "mean", "threshold": None,
                          "frequency": "auto", "weight_by_days": True,
                          "require_complete": True, "min_valid_fraction": 1.0},
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
            "units_transform": {"scale": 1.0, "offset": 0.0, "units": None},
            "plev": None,
            "bbox": None,
            "historical_experiment": "historical",
            "scenarios": ["ssp245"],
            "scenario_labels": {},
            "members": ["r1i1p1f1"],
            "member_combine": "first",  # first | mean | median (multi-member ensemble)
            "models": "all",
            "exclude": [],
            "aggregate": {"months": "ANN", "stat": "mean", "threshold": None,
                          "frequency": "auto", "weight_by_days": True,
                          "require_complete": True, "min_valid_fraction": 1.0},
            "series_dir": None,
            "series_columns": {"year": "year", "value": "value"},
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
        "validation": None,
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
        "method": "none",  # opt into a scientifically justified adjustment
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
        "enabled": False,
        "sample_kind": "annual_index",  # annual_index (descriptive) | block_maxima
        "source": "future",  # future | historical | both
        "fit_window": None,  # {"future": [2026, 2100], "historical": [1950, 2024]}
        "distribution": "gev",  # gev | gumbel
        "fit_method": "mle",  # mle | mm (method of moments; no silent fallback)
        "return_periods": [1.1, 2, 5, 10, 20, 30, 50, 100, 200, 500, 1000],
        "reference_value": None,
        "bootstrap": 0,
        "show_empirical": False,
        "min_years": 20,
        "shape_comparison": False,
    },
    "daily": {
        "enabled": False,
        "months": [6, 7, 8],  # JJA
        "index": "region_mean",  # absolute daily index: region_mean | region_diff
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
        "figures": ["timeseries", "maps", "kde"],
        "timeseries": {
            "show_members": True,
            "show_bands": True,
            "show_trend": True,
            "reference_year": None,
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
    "quality": {
        "scan": "headers",  # headers | full; thresholds also require full scanning
        "coverage": "warn",  # warn | error for incomplete requested time coverage
        "require_units": False,
        "max_missing_fraction": None,
        "valid_min": None,  # bounds in explicitly transformed input units
        "valid_max": None,
        "chunk_size": 32,
        "max_chunk_values": 2000000,
    },
    "runtime": {
        "verbose": True,
        "cache": True,
        "workers": 1,  # independent models; 1 preserves in-process execution
        "worker_threads": 1,  # supported BLAS/OpenMP thread pools per child
        "max_input_values": 30000000,
        "strict": True,
        "strict_config": True,
        "overwrite": False,
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

    @property
    def index_units(self) -> str:
        if self.get_path("index.standardize") or self.get_path("index.type") == "projection_index":
            return "dimensionless"
        return self.units

    def scenario_list(self) -> List[str]:
        sc = self.get_path("data.models.scenarios", ["ssp245"])
        if isinstance(sc, str):
            sc = [sc]
        return list(sc) if isinstance(sc, (list, tuple)) else []

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
        if not isinstance(val, (list, tuple)) or len(val) != 2:
            raise ConfigError(f"periods.{key} must be [start, end], got {val!r}")
        if any(isinstance(x, bool) or not isinstance(x, int) for x in val):
            raise ConfigError(f"periods.{key} must contain integer years")
        return list(val)

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
    if isinstance(spec, int) and not isinstance(spec, bool):
        spec = [spec]
    if isinstance(spec, str):
        key = spec.strip().upper()
        if key in SEASONS:
            return list(SEASONS[key])
        if key in ("ALL", "YEAR", "ANNUAL"):
            return list(range(1, 13))
        # allow "9-10" or "9,10"
        if re.fullmatch(r"\d{1,2}-\d{1,2}", key):
            lo, hi = map(int, key.split("-"))
            spec = list(range(lo, hi + 1)) if lo <= hi else list(range(lo, 13)) + list(range(1, hi + 1))
        else:
            spec = [int(p) for p in key.split(",") if p.strip()]
    if not isinstance(spec, (list, tuple)) or not spec:
        raise ConfigError("months must be a non-empty month list or season name")
    if any(isinstance(m, bool) or not isinstance(m, int) or not 1 <= m <= 12 for m in spec):
        raise ConfigError("months must contain integers from 1 to 12")
    if len(set(spec)) != len(spec):
        raise ConfigError("months must not contain duplicates")
    return list(spec)


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
    try:
        with open(path, "r", encoding="utf-8") as fh:
            user = yaml.safe_load(fh) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(f"cannot read YAML config {path}: {exc}") from exc
    if not isinstance(user, dict):
        raise ConfigError(f"top level of {path} must be a mapping")

    _validate_structure(user, DEFAULTS)
    merged = _deep_merge(DEFAULTS, user)
    if "future" in user.get("periods", {}) and "future_windows" not in user.get("periods", {}):
        merged["periods"]["future_windows"] = None
    if not merged.get("project", {}).get("outdir"):
        merged["project"]["outdir"] = "runs/" + str(merged["project"].get("name", "cmip_projection"))
    config_dir = os.path.dirname(os.path.abspath(path))
    root = merged["data"]["models"].get("root")
    if root:
        merged["data"]["models"]["root"] = _resolve_path(root, config_dir)
    merged = _resolve_placeholders(merged, _env_for_placeholders(merged))
    for group in ("reference", "models"):
        spec = merged["data"][group]
        if spec["units_transform"].get("units") is None and merged["variable"].get("units"):
            spec["units_transform"]["units"] = merged["variable"]["units"]
        for key in ("root", "series_file", "series_dir", "pattern"):
            if spec.get(key):
                spec[key] = _resolve_path(spec[key], config_dir)
        for key in ("paths", "patterns"):
            if spec.get(key):
                if not isinstance(spec[key], list):
                    raise ConfigError(f"data.{group}.{key} must be a list")
                spec[key] = [_resolve_path(p, config_dir) for p in spec[key]]
    merged["project"]["outdir"] = _resolve_path(merged["project"]["outdir"], config_dir)
    if isinstance(merged["data"]["models"].get("scenarios"), str):
        merged["data"]["models"]["scenarios"] = [merged["data"]["models"]["scenarios"]]
    cfg = Config(merged)
    cfg["_config_path"] = os.path.abspath(path)
    cfg["_config_dir"] = os.path.dirname(os.path.abspath(path))
    errors, warnings = validate_config(cfg, user_keys=set(_iter_leaf_paths(user)))
    for msg in warnings:
        log.warning(msg)
    if errors:
        raise ConfigError("invalid configuration:\n  - " + "\n  - ".join(errors))
    return cfg


def _validate_structure(user: Dict[str, Any], defaults: Dict[str, Any], prefix: str = "") -> None:
    for key, value in user.items():
        dotted = f"{prefix}.{key}" if prefix else str(key)
        if not isinstance(key, str):
            raise ConfigError("configuration keys must be strings")
        expected = defaults.get(key)
        if isinstance(expected, dict):
            if not isinstance(value, dict):
                raise ConfigError(f"{dotted} must be a mapping")
            _validate_structure(value, expected, dotted)


def _resolve_path(value: str, config_dir: str) -> str:
    if not isinstance(value, str):
        raise ConfigError(f"paths must be strings, got {value!r}")
    expanded = os.path.expanduser(os.path.expandvars(value))
    if re.search(r"\$(?:\{[^}]+\}|[A-Za-z_][A-Za-z0-9_]*)", expanded):
        raise ConfigError(f"undefined environment variable in path: {value}")
    if "://" in expanded:
        raise ConfigError("remote URLs are not local input paths; materialize the files first")
    return os.path.abspath(os.path.join(config_dir, expanded))


def validate_config(cfg: Config, user_keys: Optional[set] = None):
    """Return ``(errors, warnings)``; never raise so callers can report all."""
    errors: List[str] = []
    warnings: List[str] = []

    def number(key: str, lo=None, hi=None, integer=False, nullable=False):
        value = cfg.get_path(key)
        if value is None and nullable:
            return
        valid = (isinstance(value, int) if integer else isinstance(value, (int, float)))
        if not valid or isinstance(value, bool) or not math.isfinite(value):
            errors.append(f"{key} must be a finite {'integer' if integer else 'number'}")
        elif (lo is not None and value < lo) or (hi is not None and value > hi):
            errors.append(f"{key} must be in the range [{lo}, {hi}]")

    def year_pair(value, label):
        if (not isinstance(value, (list, tuple)) or len(value) != 2
                or any(isinstance(x, bool) or not isinstance(x, int) for x in value)
                or value[0] > value[1]):
            errors.append(f"{label} must be [start_year, end_year] in ascending order")
            return False
        return True

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
    if kind != mkind:
        errors.append("reference and models must use the same data kind (grid or series)")
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
        if not isinstance(sc, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", sc):
            errors.append("scenario names must contain only letters, numbers, dot, underscore or hyphen")
        if sc == "historical":
            errors.append("data.models.scenarios must not contain 'historical' (that is the baseline experiment)")
    if len(cfg.scenario_list()) != len(set(map(str, cfg.scenario_list()))):
        errors.append("data.models.scenarios must not contain duplicates")
    models = mod.get("models")
    if models != "all" and (not isinstance(models, list) or not models or
                             any(not isinstance(m, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", m) for m in models)):
        errors.append("data.models.models must be 'all' or a non-empty list of safe model names")
    if isinstance(models, list) and len(models) != len(set(map(str, models))):
        errors.append("data.models.models must not contain duplicates")
    if not isinstance(mod.get("members"), list) or not mod.get("members"):
        errors.append("data.models.members must be a non-empty list")

    for group in ("reference", "models"):
        aggregate = cfg.get_path(f"data.{group}.aggregate")
        try:
            resolve_months(aggregate.get("months"))
        except (ConfigError, TypeError, ValueError) as exc:
            errors.append(f"data.{group}.aggregate.months: {exc}")
        if aggregate.get("stat") not in AGG_STATS:
            errors.append(f"data.{group}.aggregate.stat must be one of {AGG_STATS}")
        if str(aggregate.get("stat", "")).startswith("count_") and aggregate.get("threshold") is None:
            errors.append(f"data.{group}.aggregate.threshold is required for counts")
        if aggregate.get("frequency") not in ("auto", "monthly", "daily", "annual"):
            errors.append(f"data.{group}.aggregate.frequency must be auto, monthly, daily or annual")
        number(f"data.{group}.aggregate.min_valid_fraction", 0, 1)
        number(f"data.{group}.units_transform.scale")
        number(f"data.{group}.units_transform.offset")
        for flag in ("weight_by_days", "require_complete"):
            if not isinstance(aggregate.get(flag), bool):
                errors.append(f"data.{group}.aggregate.{flag} must be boolean")
    if kind == mkind == "grid":
        if ref.get("aggregate", {}).get("stat") != mod.get("aggregate", {}).get("stat"):
            errors.append("reference and models must use the same aggregation statistic")
        try:
            if set(cfg.months("reference")) != set(cfg.months("models")):
                errors.append("reference and models must use the same season/months")
        except (ValueError, TypeError, ConfigError):
            pass

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
    for key in ("correction.params_period", "daily.baseline", "periods.validation"):
        if cfg.get_path(key) is not None:
            year_pair(cfg.get_path(key), key)
    if ref.get("years") is not None:
        year_pair(ref.get("years"), "data.reference.years")
    windows = cfg.get_path("periods.future_windows")
    if windows is not None:
        if not isinstance(windows, list):
            errors.append("periods.future_windows must be a list of year pairs")
        else:
            for win in windows:
                if year_pair(win, "periods.future_windows entry"):
                    try:
                        lo, hi = cfg.period("future")
                        if win[0] < lo or win[1] > hi:
                            errors.append("periods.future_windows must be contained in periods.future")
                    except ConfigError:
                        pass

    # index
    itype = cfg.get_path("index.type", "region_mean")
    valid_index = ("region_mean", "region_mean_anomaly", "region_diff", "projection_index", "series")
    if itype not in valid_index:
        errors.append(f"index.type must be one of {valid_index}, got {itype!r}")
    if (kind == "series") != (itype == "series"):
        errors.append("series input requires index.type: series; grid input requires a gridded index")
    if kind == "series" and any(cfg.get_path(f"index.{key}") for key in
                                 ("anomaly", "detrend", "standardize", "remove_zonal_mean")):
        errors.append("series input is already an index: preprocess anomaly/detrend/standardize upstream")
    if cfg.get_path("index.detrend"):
        warnings.append("index.detrend removes the fitted trend separately from each supplied period; use for variability diagnostics only")
    if cfg.get_path("index.standardize"):
        warnings.append("standardized index values and trend units are dimensionless; set index.units_label accordingly")
    projection = cfg.get_path("index.projection")
    if projection.get("template_source") not in ("reference", "year", "model"):
        errors.append("index.projection.template_source must be reference, year or model")
    if itype == "projection_index":
        if projection.get("template_source") == "model" and not projection.get("template_model"):
            errors.append("index.projection.template_model is required for a model template")
        if projection.get("template_source") == "year" and not isinstance(projection.get("template_year"), int):
            errors.append("index.projection.template_year must be an integer")
        if projection.get("template_period") is not None:
            year_pair(projection["template_period"], "index.projection.template_period")
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
    number("correction.upper_bound", nullable=True)
    number("correction.lower_bound", nullable=True)
    number("correction.n_quantiles", 5, 10000, integer=True)
    ub = cfg.get_path("correction.upper_bound")
    if isinstance(lb, (int, float)) and isinstance(ub, (int, float)) and lb >= ub:
        errors.append("correction.lower_bound must be smaller than upper_bound")
    if not isinstance(cfg.get_path("correction.apply_to"), list) or any(
            x not in ("historical", "future") for x in cfg.get_path("correction.apply_to", [])):
        errors.append("correction.apply_to must be a list drawn from historical, future")
    if method == "quantile_mapping" and "future" in (cfg.get_path("correction.apply_to") or []):
        warnings.append("quantile_mapping can distort future quantile changes; justify its use or choose quantile_delta_mapping")
    if method == "mean_std":
        warnings.append("mean_std rescales modeled variability and generally its trend; historical agreement is calibration, not independent validation")
    if cfg.get_path("correction.trend_correction.enabled"):
        warnings.append("trend_correction imposes an observed-trend adjustment; report it as a sensitivity experiment")
    if cfg.get_path("correction.trend_correction.mode") not in ("difference", "ratio"):
        errors.append("correction.trend_correction.mode must be difference or ratio")

    if cfg.get_path("grid.method") not in ("linear", "nearest"):
        errors.append("grid.method must be linear or nearest; conservative regridding is not implemented")
    for axis in ("lat", "lon"):
        values = cfg.get_path(f"grid.{axis}")
        if isinstance(values, dict):
            if not all(isinstance(values.get(k), (int, float)) for k in ("start", "stop", "step")):
                errors.append(f"grid.{axis} needs numeric start, stop and step")
            elif values["step"] <= 0 or values["stop"] < values["start"]:
                errors.append(f"grid.{axis} must increase with a positive step")
        elif values is not None and (not isinstance(values, list) or len(values) < 1):
            errors.append(f"grid.{axis} must be an axis mapping or a non-empty coordinate list")
    if cfg.get_path("ensemble.statistic") not in ("mean", "median"):
        errors.append("ensemble.statistic must be mean or median")
    if cfg.get_path("ensemble.weights") not in ("equal", "skill"):
        errors.append("ensemble.weights must be equal or skill")
    if cfg.get_path("ensemble.weights") == "skill":
        warnings.append("skill weights use in-sample RMSE and do not adjust for model dependence; treat as sensitivity only")
    number("ensemble.min_models", 1, integer=True)
    bands = cfg.get_path("ensemble.bands")
    if not isinstance(bands, list) or any(isinstance(q, bool) or not isinstance(q, int) or not 0 <= q <= 100 for q in bands):
        errors.append("ensemble.bands must list integer percentiles in [0,100]")

    for prefix in ("extremes", "daily.extremes"):
        number(f"{prefix}.min_years", 3, integer=True)
        rp = cfg.get_path(f"{prefix}.return_periods")
        if not isinstance(rp, list) or not rp or any(not isinstance(t, (int, float)) or isinstance(t, bool) or not math.isfinite(t) or t <= 1 for t in rp):
            errors.append(f"{prefix}.return_periods must list finite numbers greater than 1")
        if cfg.get_path(f"{prefix}.distribution") not in ("gev", "gumbel"):
            errors.append(f"{prefix}.distribution must be gev or gumbel")
    if cfg.get_path("extremes.sample_kind") not in ("annual_index", "block_maxima"):
        errors.append("extremes.sample_kind must be annual_index or block_maxima")
    if cfg.get_path("extremes.enabled") and cfg.get_path("extremes.sample_kind") == "annual_index":
        warnings.append("GEV for an annual index is a descriptive distribution fit, not a block-maxima extreme-value justification")
    if cfg.get_path("extremes.enabled") and combine != "first":
        warnings.append("averaging ensemble members before extremes changes variability and return levels")
    if cfg.get_path("extremes.source") not in ("historical", "future", "both"):
        errors.append("extremes.source must be historical, future or both")
    if cfg.get_path("extremes.fit_method") not in ("mle", "mm"):
        errors.append("extremes.fit_method must be mle or mm; PWM is not implemented")
    number("extremes.bootstrap", 0, integer=True)
    if cfg.get_path("extremes.bootstrap") != 0:
        errors.append("extremes.bootstrap must be 0: ensemble-median resampling is not a valid per-model uncertainty workflow")
    fw = cfg.get_path("extremes.fit_window")
    if fw is not None:
        if not isinstance(fw, dict):
            errors.append("extremes.fit_window must be a mapping of historical/future year pairs")
        else:
            for label, value in fw.items():
                if label not in ("historical", "future"):
                    errors.append(f"unknown extremes.fit_window key: {label}")
                year_pair(value, f"extremes.fit_window.{label}")

    # daily extremes
    daily = cfg.get_path("daily", {}) or {}
    if daily.get("enabled"):
        if daily.get("index") not in ("region_mean", "region_diff"):
            errors.append("daily.index must be region_mean or region_diff; daily anomaly climatology is not implemented")
        if combine != "first":
            errors.append("daily extremes require member_combine: first; compute member-specific extremes before combining results")
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
            errors.append("daily extremes require gridded input")
        for group in ("reference", "models"):
            if not cfg.get_path(f"data.{group}.bbox"):
                errors.append(f"data.{group}.bbox is required for daily loading")
        try:
            resolve_months(daily.get("months"))
        except (ConfigError, TypeError, ValueError) as exc:
            errors.append(f"daily.months: {exc}")
        number("daily.extremes.threshold_quantile", 0.001, 0.999999)
        number("daily.extremes.decluster_days", 1, integer=True)
        number("daily.extremes.min_peaks", 3, integer=True)

    # Raw-input quality checks. Value thresholds are explicit scientific choices,
    # not universal temperature/precipitation defaults.
    if cfg.get_path("quality.scan", "headers") not in ("headers", "full"):
        errors.append("quality.scan must be headers or full")
    if cfg.get_path("quality.coverage", "warn") not in ("warn", "error"):
        errors.append("quality.coverage must be warn or error")
    if not isinstance(cfg.get_path("quality.require_units", False), bool):
        errors.append("quality.require_units must be boolean")
    for name in ("chunk_size", "max_chunk_values"):
        number(f"quality.{name}", 1, integer=True)
    number("quality.max_missing_fraction", 0, 1, nullable=True)
    number("quality.valid_min", nullable=True)
    number("quality.valid_max", nullable=True)
    qlo, qhi = cfg.get_path("quality.valid_min"), cfg.get_path("quality.valid_max")
    if (isinstance(qlo, (int, float)) and not isinstance(qlo, bool)
            and isinstance(qhi, (int, float)) and not isinstance(qhi, bool)
            and qlo > qhi):
        errors.append("quality.valid_min must not exceed quality.valid_max")

    # plots
    figs = cfg.get_path("plots.figures", [])
    if not isinstance(figs, list):
        errors.append("plots.figures must be a list")
        figs = []
    valid_figs = ("timeseries", "gev", "combined", "maps", "kde", "shape_box")
    for fig in figs:
        if fig not in valid_figs:
            errors.append(f"plots.figures contains unknown figure {fig!r}; valid: {valid_figs}")
    number("plots.maps.agreement", 0.5, 1)
    number("plots.dpi", 50, 2400, integer=True)
    for key in ("strict", "strict_config", "overwrite", "cache", "verbose"):
        if not isinstance(cfg.get_path(f"runtime.{key}"), bool):
            errors.append(f"runtime.{key} must be boolean")
    number("runtime.max_input_values", 1, integer=True)
    number("runtime.workers", 1, integer=True)
    number("runtime.worker_threads", 1, integer=True)
    for key in ("extremes.enabled", "daily.enabled", "daily.plots", "grid.lon_360",
                "index.anomaly", "index.standardize", "index.detrend", "index.remove_zonal_mean",
                "correction.trend_correction.enabled"):
        if not isinstance(cfg.get_path(key), bool):
            errors.append(f"{key} must be boolean, not a quoted string")
    validation = cfg.get_path("periods.validation")
    calibration = cfg.get_path("correction.params_period") or cfg.get_path("periods.reference_climate")
    if (isinstance(validation, list) and len(validation) == 2 and isinstance(calibration, list)
            and len(calibration) == 2 and all(isinstance(x, int) for x in validation + calibration)):
        if max(validation[0], calibration[0]) <= min(validation[1], calibration[1]):
            warnings.append("periods.validation overlaps the calibration window; these scores are not held-out validation")
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
                target = errors if cfg.get_path("runtime.strict_config", True) else warnings
                target.append(f"unknown configuration key: {key}")
    return errors, warnings


PIPELINE_STAGES = ("prepare", "correct", "ensemble", "extremes", "daily", "figures", "report")


def plan_stages(cfg: Config, requested=None) -> List[str]:
    """Expand required in-memory dependencies; this is not artifact resume."""
    tokens = [s.strip().lower() for item in (requested or ["all"])
              for s in str(item).split(",") if s.strip()]
    unknown = set(tokens) - set(PIPELINE_STAGES) - {"all"}
    if unknown:
        raise ConfigError(f"unknown stages: {sorted(unknown)}; valid stages: {PIPELINE_STAGES}")
    if not tokens:
        raise ConfigError("at least one stage is required")
    if "all" in tokens:
        tokens = list(PIPELINE_STAGES)
    enabled = {"extremes": bool(cfg.get_path("extremes.enabled")),
               "daily": bool(cfg.get_path("daily.enabled"))}
    figs = set(cfg.get_path("plots.figures") or [])
    fig_deps = []
    if figs & {"timeseries", "combined"}:
        fig_deps.append("ensemble")
    if figs & {"maps", "kde"}:
        fig_deps.append("correct")
    if figs & {"gev", "combined", "shape_box"} and enabled["extremes"]:
        fig_deps.append("extremes")
    if enabled["daily"] and cfg.get_path("daily.plots", True):
        fig_deps.append("daily")
    dependencies = {"prepare": [], "correct": ["prepare"], "ensemble": ["correct"],
                    "extremes": ["ensemble"], "daily": ["prepare"],
                    "figures": fig_deps,
                    "report": ["ensemble", "figures"] + [s for s, yes in enabled.items() if yes]}
    needed = set()

    def add(stage):
        if stage in needed or enabled.get(stage) is False:
            return
        for dependency in dependencies[stage]:
            add(dependency)
        needed.add(stage)

    for stage in tokens:
        add(stage)
    return [s for s in PIPELINE_STAGES if s in needed]


def quality_requires_values(cfg: Config) -> bool:
    """Explicit limits must never silently pass a headers-only production run."""
    return (cfg.get_path("quality.scan", "headers") == "full"
            or any(cfg.get_path(f"quality.{key}") is not None
                   for key in ("max_missing_fraction", "valid_min", "valid_max")))


def _atomic_text(path: str, writer) -> None:
    """Keep the previous complete metadata file if serialization/writing fails."""
    parent = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(parent, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".cmipkit-", suffix=".tmp", dir=parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            writer(stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json_atomic(path: str, payload: Any) -> None:
    _atomic_text(path, lambda stream: json.dump(payload, stream, indent=2, ensure_ascii=False))


def write_resolved(cfg: Config, path: str) -> None:
    """Dump the fully merged config next to the run for reproducibility."""
    payload = cfg.to_dict()
    payload.pop("_config_path", None)
    payload.pop("_config_dir", None)
    _atomic_text(path, lambda stream: yaml.safe_dump(payload, stream, allow_unicode=True, sort_keys=False))


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
    write_json_atomic(path, manifest)
    return manifest
