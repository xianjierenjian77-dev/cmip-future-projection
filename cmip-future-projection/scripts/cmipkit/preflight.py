"""Read-only discovery and quality checks; no analysis output is created."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from .config import Config, resolve_months, validate_config


def inspect_config(cfg: Config, *, check_inputs: bool = False,
                   check_data: bool = False) -> Dict[str, Any]:
    """Check every discovered input; optionally scan all selected field values."""
    check_inputs = bool(check_inputs or check_data)
    errors, warnings = validate_config(cfg)
    result = {
        "ok": not errors, "config": cfg.get_path("_config_path"), "outdir": cfg.outdir,
        "inputs_checked": check_inputs, "models": [], "datasets": [],
        "errors": list(errors), "warnings": list(warnings),
        "quality": {
            "level": "full" if check_data else "headers" if check_inputs else "config",
            "data_values_checked": False, "ok": not errors, "datasets": [],
            "summary": {"datasets": 0, "files": 0, "errors": len(errors),
                        "warnings": len(warnings), "values_checked": 0},
        },
    }
    if errors or not check_inputs:
        return result
    from . import io, quality
    source = io.DataSource(cfg, create_dirs=False)
    options = cfg.get_path("quality", {}) or {}
    audited_scopes = {}

    def scope_key(files, spec, years, months, frequency, expected):
        selected = {key: spec.get(key) for key in
                    ("kind", "variable", "units_in", "units_transform", "plev", "bbox")}
        return (tuple(sorted(set(files))), tuple(years), tuple(sorted(months)), frequency,
                json.dumps(selected, sort_keys=True, default=str), json.dumps(expected, sort_keys=True, default=str))

    def add(label, files, spec, years, months, expected=None, *, use="annual"):
        files = sorted(set(files))
        item = {"label": label, "files": files, "files_metadata": []}
        result["datasets"].append(item)
        if spec.get("kind") == "series":
            audit = quality.audit_series(files, spec, label=label, years=years,
                                         quality=options, check_data=check_data)
        else:
            audit = quality.audit_grid(files, spec, label=label, years=years, months=months,
                                       quality=options, check_data=check_data, expected=expected)
        if not files:
            audit["findings"].insert(0, dict(severity="error", code="no_files",
                                            message="no files matched the configured paths/member"))
            audit["ok"] = False
        for header in audit["files"]:
            if "size_bytes" in header:
                item["files_metadata"].append({key: header[key] for key in ("path", "size_bytes", "mtime_ns")})
        if audit["files"]:
            item["sample_header"] = {k: v for k, v in audit["files"][0].items()
                                      if k not in ("path", "size_bytes", "mtime_ns", "values")}
        result["quality"]["datasets"].append(audit)
        audit["uses"] = [use]
        frequency = audit.get("coverage", {}).get("frequency", spec.get("aggregate", {}).get("frequency", "auto"))
        audited_scopes[scope_key(files, spec, years, months, frequency, expected)] = audit
        return audit

    def add_daily(label, files, spec, years, months, expected=None):
        # An annual-context aggregation of daily data can already have scanned
        # exactly this scope. Reuse it, but require genuinely daily samples.
        spec = dict(spec, aggregate=dict(spec.get("aggregate") or {}, frequency="daily"))
        existing = audited_scopes.get(scope_key(files, spec, years, months, "daily", expected))
        if existing is not None:
            if "daily" not in existing["uses"]:
                existing["uses"].append("daily")
            audit = existing
        else:
            audit = add(label, files, spec, years, months, expected, use="daily")
        if audit.get("input_frequency") != "daily":
            finding = dict(severity="error", code="daily_sampling_required",
                           message=f"daily analysis requires daily input; raw sampling resolved to {audit.get('input_frequency')!r}")
            if finding not in audit["findings"]:
                audit["findings"].append(finding)
            audit["ok"] = False
        return audit

    ref = cfg.get_path("data.reference")
    ref_years = ref.get("years") or cfg.period("historical")
    if ref["kind"] == "series":
        add("reference", [ref["series_file"]], ref, ref_years, cfg.months("reference"))
    else:
        files = []
        for template in ref["paths"]:
            files.extend(io.expand_paths(template, io._input_years(ref_years, cfg.months("reference"))))
        add("reference", files, ref, ref_years, cfg.months("reference"))

    try:
        models = source.discover_models()
    except Exception as exc:
        models = []
        result["errors"].append(f"model discovery: {exc}")
    result["models"] = models
    if not models:
        result["errors"].append("no models found: check root, patterns, models, members and historical experiment")
    elif len(models) < int(cfg.get_path("ensemble.min_models", 3)):
        result["warnings"].append(
            f"[insufficient_models] discovered {len(models)} model(s), below ensemble.min_models="
            f"{cfg.get_path('ensemble.min_models', 3)}; this inventory cannot form an ensemble meeting "
            "that minimum. Add models, or explicitly choose a smaller minimum for a pilot run; "
            "the configured threshold has not been changed."
        )
    spec = cfg.get_path("data.models")
    hist = spec.get("historical_experiment", "historical")
    experiments = [hist] + cfg.scenario_list()
    months = cfg.months("models")
    primary_members = {}
    for model in models:
        locked_members = None
        historical_identity = {}
        for experiment in experiments:
            years = cfg.period("historical" if experiment == hist else "future")
            if spec["kind"] == "series":
                root = Path(spec["series_dir"])
                paths = ([root / model / f"{experiment}{suffix}" for suffix in (".csv", ".npz")]
                         + [root / f"{model}_{experiment}{suffix}" for suffix in (".csv", ".npz")])
                files = [str(p) for p in paths if p.is_file()]
                if len(files) > 1:
                    result["errors"].append(f"{model}/{experiment}: multiple series files match; select one layout unambiguously")
                add(f"{model}/{experiment}", files, spec, years, months)
                continue
            try:
                input_years = io._input_years(years, months)
                expected = {"source_id": model, "experiment_id": experiment}
                if spec.get("member_combine", "first") != "first":
                    configured = list(spec.get("members") or ["*"])
                    members = io.discover_members(spec, experiment, model, input_years) if configured == ["*"] else configured
                    if not members:
                        raise ValueError("cannot identify ensemble members; configure an explicit {member} pattern")
                    if locked_members is not None and members != locked_members:
                        result["errors"].append(f"{model}/{experiment}: member sets differ across experiments; choose common members explicitly")
                    if experiment == hist:
                        locked_members = members
                    for member in locked_members or members:
                        files = io._files_for_member(spec, experiment, model, member, input_years)
                        audit = add(f"{model}/{experiment}/{member}", files, spec, years, months,
                                    dict(expected, variant_label=member))
                        _compare_experiments(audit, historical_identity, member, experiment == hist)
                else:
                    selected = dict(spec)
                    if locked_members:
                        selected["members"] = locked_members
                    files = io._files_for_grid_source(selected, experiment, model, input_years)
                    members = {io._member_of(path, model, experiment, spec) for path in files}
                    members.discard(None)
                    if len(members) > 1:
                        raise ValueError("multiple members selected for a single-member dataset")
                    if experiment == hist:
                        locked_members = sorted(members) or None
                    if locked_members:
                        expected["variant_label"] = locked_members[0]
                    audit = add(f"{model}/{experiment}", files, spec, years, months, expected)
                    # Metadata can establish identity even for an explicit filename
                    # pattern without a parseable member component.
                    if experiment == hist and not locked_members:
                        variant = audit["identity"].get("variant_label")
                        locked_members = [variant] if variant else None
                    if experiment == hist:
                        primary_members[model] = locked_members
                    _compare_experiments(audit, historical_identity, "selected", experiment == hist)
            except Exception as exc:
                result["errors"].append(f"{model}/{experiment}: {type(exc).__name__}: {exc}")

    if cfg.get_path("daily.enabled", False):
        # daily_stage uses its own season and periods.historical for reference,
        # irrespective of data.reference.years or the annual-context season.
        daily_months = resolve_months(cfg.get_path("daily.months"))
        hist_years = cfg.period("historical")
        daily_files = []
        for template in ref["paths"]:
            daily_files.extend(io.expand_paths(template, io._input_years(hist_years, daily_months)))
        add_daily("daily reference", daily_files, ref, hist_years, daily_months)
        for model in models:
            daily_members = primary_members.get(model)
            historical_identity = {}
            for experiment in experiments:
                years = cfg.period("historical" if experiment == hist else "future")
                try:
                    selected = dict(spec)
                    if daily_members:
                        selected["members"] = daily_members
                    files = io._files_for_grid_source(selected, experiment, model,
                                                      io._input_years(years, daily_months))
                    if experiment == hist:
                        found_members = {io._member_of(path, model, experiment, spec) for path in files}
                        found_members.discard(None)
                        if len(found_members) > 1:
                            raise ValueError("daily historical files mix multiple members")
                        if found_members:
                            daily_members = sorted(found_members)
                    expected = {"source_id": model, "experiment_id": experiment}
                    if daily_members:
                        expected["variant_label"] = daily_members[0]
                    audit = add_daily(f"daily {model}/{experiment}", files, spec, years, daily_months, expected)
                    if experiment == hist and not daily_members:
                        variant = audit["identity"].get("variant_label")
                        daily_members = [variant] if variant else None
                    _compare_experiments(audit, historical_identity, "selected", experiment == hist)
                except Exception as exc:
                    result["errors"].append(f"daily {model}/{experiment}: {type(exc).__name__}: {exc}")

    discovery_complete = not result["errors"]
    audits = result["quality"]["datasets"]
    for audit in audits:
        for finding in audit["findings"]:
            context = f" ({finding['file']})" if finding.get("file") else ""
            text = f"{audit['label']}{context}: [{finding['code']}] {finding['message']}"
            result["errors" if finding["severity"] == "error" else "warnings"].append(text)
    result["ok"] = not result["errors"]
    result["quality"]["ok"] = result["ok"]
    result["quality"]["data_values_checked"] = bool(check_data and discovery_complete and audits
                                                    and all(a["values_checked"] for a in audits))
    result["quality"]["summary"].update(
        datasets=len(audits), files=sum(len(a["files"]) for a in audits),
        errors=len(result["errors"]), warnings=len(result["warnings"]),
        values_checked=sum(a.get("values", {}).get("count", 0) for a in audits))
    return result


def _compare_experiments(audit, historical, member, is_historical):
    identity = audit["identity"]
    if is_historical:
        historical[member] = identity
        return
    baseline = historical.get(member, {})
    # Experiments themselves necessarily differ. Calendars are compared only
    # within a model/member, never across unrelated climate models.
    for key in ("source_id", "variant_label", "grid_label", "units", "calendar"):
        if key in identity and key in baseline and identity[key] != baseline[key]:
            finding = dict(severity="error", code="historical_scenario_mismatch",
                           message=f"{key}={identity[key]!r} differs from historical {baseline[key]!r}")
            if finding not in audit["findings"]:
                audit["findings"].append(finding)
            audit["ok"] = False
