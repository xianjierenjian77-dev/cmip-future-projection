"""Pipeline orchestration: from raw files to figures and a written report.

Stages
------
``prepare``   resolve directories, load the observational reference, discover
              models, build the projection template
``correct``   fit the bias-correction transfer function per model, apply it to
              the historical and future fields, derive the index series, save
              corrected fields and accumulate the map fields
``ensemble``  align the multi-model series, optional trend correction, ensemble
              statistics and trend tests
``extremes``  GEV fits, return levels, multi-model spread, shape comparison
``figures``   all requested figures
``report``    CSV tables + ``report.md``

Everything is executed in one process so the stages share the (small) in-memory
series; corrected fields are streamed to disk and only accumulated sums are
kept for the maps.
"""

from __future__ import annotations

import json
import logging
import os
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import xarray as xr

from . import correction, daily as daily_mod, ensemble, extremes, indices, io, plotting
from .config import Config, ConfigError, resolve_months, write_manifest, write_resolved, plan_stages

log = logging.getLogger("cmipkit.pipeline")

STAGES = ("prepare", "correct", "ensemble", "extremes", "daily", "figures", "report")


def _setup_logging(outdir: str, verbose: bool = True) -> logging.Logger:
    logger = logging.getLogger("cmipkit")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(name)s | %(message)s", "%H:%M:%S")
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    try:
        fh = logging.FileHandler(os.path.join(outdir, "logs", "run.log"), encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(name)s | %(message)s"))
        logger.addHandler(fh)
    except Exception:
        pass
    logger.propagate = False
    return logger


def _dirs(outdir: str) -> Dict[str, str]:
    d = {
        "root": outdir,
        "cache": os.path.join(outdir, "cache"),
        "corrected": os.path.join(outdir, "corrected"),
        "params": os.path.join(outdir, "params"),
        "series": os.path.join(outdir, "series"),
        "ensemble": os.path.join(outdir, "ensemble"),
        "extremes": os.path.join(outdir, "extremes"),
        "figures": os.path.join(outdir, "figures"),
        "tables": os.path.join(outdir, "tables"),
        "logs": os.path.join(outdir, "logs"),
    }
    for path in d.values():
        os.makedirs(path, exist_ok=True)
    return d


def _diagnostic_spec(spec: Dict[str, Any]) -> Dict[str, Any]:
    """Index spec for skill metrics: absolute values, no anomaly/standardisation.

    Bias, KS and KGE are meaningless on anomaly series (whose mean is zero by
    construction), so the diagnostics are always computed on the absolute index.
    """
    out = dict(spec)
    out["anomaly"] = False
    out["standardize"] = False
    out["detrend"] = False
    if out.get("type") == "region_mean_anomaly":
        out["type"] = "region_mean"
    return out


def _median_finite(values: Sequence[float]) -> float:
    """median of the finite values, NaN when nothing is finite (no warnings)."""
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if arr.size else float("nan")


def _fmt(x: Any, nd: int = 3) -> str:
    try:
        val = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(val):
        return "nan"
    return f"{val:.{nd}f}"


class Pipeline:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.name = cfg.name
        self.outdir = os.path.abspath(cfg.outdir)
        self.previous_run_directory = None
        previous_path = os.path.join(self.outdir, "config.resolved.yaml")
        if os.path.isfile(previous_path):
            import yaml
            with open(previous_path, encoding="utf-8") as fh:
                previous = yaml.safe_load(fh) or {}
            current = cfg.to_dict()
            for value in (previous, current):
                for key in ("_config_path", "_config_dir", "plots", "runtime"):
                    value.pop(key, None)
            if previous != current:
                if not cfg.get_path("runtime.overwrite", False):
                    raise ConfigError("output directory contains a different scientific configuration; choose a new project.outdir or explicitly use --overwrite")
                self.previous_run_directory = self.outdir + ".previous-" + str(time.time_ns())
                os.rename(self.outdir, self.previous_run_directory)
        os.makedirs(self.outdir, exist_ok=True)
        self.dirs = _dirs(self.outdir)
        self.log = _setup_logging(self.outdir, bool(cfg.get_path("runtime.verbose", True)))
        self.source = io.DataSource(cfg)
        self.template: Optional[xr.DataArray] = None
        self.obs_field: Optional[xr.DataArray] = None
        self.obs_years: Optional[np.ndarray] = None
        self.obs_series: Optional[np.ndarray] = None
        self.obs_diag: Optional[np.ndarray] = None
        self.obs_hist: Optional[Any] = None
        self.models: List[str] = []
        self.series: Dict[str, Dict[str, Any]] = {}
        self.raw_series: Dict[str, Dict[str, Any]] = {}
        self.params: Dict[str, xr.Dataset] = {}
        self.scores: Dict[str, Dict[str, float]] = {}
        self.map_acc: Dict[str, Dict[str, Any]] = {}
        self.pdf_data: Dict[str, Dict[str, np.ndarray]] = {}
        self.ens: Dict[str, Any] = {}
        self.gev: Dict[str, Any] = {}
        self.shapes: Dict[str, List[float]] = {}
        self.kruskal: Dict[str, float] = {}
        self.daily: Dict[str, Any] = {}
        self._finished: Dict[str, Any] = {}
        self.manifest: Dict[str, Any] = {}
        self.input_inventory: List[Dict[str, Any]] = []
        self.daily_input_inventory: List[Dict[str, Any]] = []
        self.failures: List[Dict[str, str]] = []
        self._artifacts_before = self._artifact_snapshot()

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def run(self, stages: Optional[Sequence[str]] = None) -> Dict[str, Any]:
        stages = self._normalise_stages(stages)
        t0 = time.time()
        self.log.info("=== %s | stages: %s ===", self.name, ", ".join(stages))
        actions = {"prepare": self.prepare, "correct": self.correct, "ensemble": self.ensemble_stage,
                   "extremes": self.extremes_stage, "daily": self.daily_stage,
                   "figures": self.figures, "report": self.write_report}
        try:
            for stage in stages:
                actions[stage]()
        except Exception as exc:
            self._save_status(stages, "failed", str(exc))
            raise
        self._save_status(stages, "completed_with_skips" if self.failures else "completed")
        self.log.info("=== done in %.1f s | output: %s ===", time.time() - t0, self.outdir)
        return self._finished

    # ------------------------------------------------------------------
    def _normalise_stages(self, stages: Optional[Sequence[str]]) -> List[str]:
        return plan_stages(self.cfg, stages)

    def _save_status(self, stages, status, error=None):
        after = self._artifact_snapshot()
        current_artifacts = [name for name, stamp in after.items() if self._artifacts_before.get(name) != stamp]
        previous_artifacts = [name for name, stamp in after.items() if self._artifacts_before.get(name) == stamp]
        self.manifest.update({"status": status, "requested_stages": list(stages),
                              "completed_stages": list(self._finished),
                              "failures": self.failures,
                              "preflight_inputs": self.input_inventory,
                              "inputs": list(self.source.actual_inputs.values()) + self.daily_input_inventory,
                              "actual_inputs": list(self.source.actual_inputs.values()) + self.daily_input_inventory,
                              "models_loaded": list(self.series),
                              "model_members": self.source.model_members,
                              "artifacts_current_invocation": sorted(current_artifacts),
                              "artifacts_not_produced_in_this_invocation": sorted(previous_artifacts),
                              "previous_run_directory": self.previous_run_directory})
        if error:
            self.manifest["error"] = error
        with open(os.path.join(self.outdir, "manifest.json"), "w", encoding="utf-8") as fh:
            json.dump(self.manifest, fh, indent=2, ensure_ascii=False)

    def _artifact_snapshot(self):
        snapshot = {}
        for directory, _, filenames in os.walk(self.outdir):
            if os.path.basename(directory) in ("cache", "logs"):
                continue
            for filename in filenames:
                if filename == "manifest.json":
                    continue
                path = os.path.join(directory, filename)
                stat = os.stat(path)
                snapshot[os.path.relpath(path, self.outdir)] = (stat.st_size, stat.st_mtime_ns)
        return snapshot

    def _failed_input(self, label, exc):
        self.failures.append({"dataset": label, "reason": str(exc)})
        self.log.warning("%s: %s", label, exc)
        if self.cfg.get_path("runtime.strict", True):
            raise exc

    def _clip(self, values):
        lower = self.cfg.get_path("correction.lower_bound")
        upper = self.cfg.get_path("correction.upper_bound")
        if lower is None and upper is None:
            return values
        return values.clip(min=lower, max=upper) if isinstance(values, xr.DataArray) else np.clip(values, lower, upper)

    # ------------------------------------------------------------------
    # prepare
    # ------------------------------------------------------------------
    def prepare(self) -> None:
        write_resolved(self.cfg, os.path.join(self.outdir, "config.resolved.yaml"))
        self.manifest = write_manifest(self.cfg, os.path.join(self.outdir, "manifest.json"),
                                       {"status": "running", "preflight_inputs": self.input_inventory,
                                        "actual_inputs": []})

        itype = self.cfg.get_path("index.type", "region_mean")
        ref_kind = self.cfg.get_path("data.reference.kind", "grid")
        ref_period = self.cfg.period("reference_climate")

        if ref_kind == "grid":
            self.obs_field = self.source.reference_field()
            self.obs_years = np.asarray(self.obs_field.year.values, dtype=int)
            self.log.info("reference field: %d years (%s-%s), grid %dx%d",
                          self.obs_years.size, self.obs_years[0], self.obs_years[-1],
                          self.obs_field.sizes["latitude"], self.obs_field.sizes["longitude"])
            if itype == "projection_index":
                projection_spec = self.cfg.get_path("index.projection", {})
                model_fields = None
                if projection_spec.get("template_source") == "model":
                    template_model = projection_spec["template_model"]
                    if template_model not in self.source.discover_models():
                        raise ConfigError(f"projection template model {template_model!r} is not in the discovered model set")
                    model_fields = {template_model: self.source.model_field(
                        template_model, self.cfg.get_path("data.models.historical_experiment"),
                        self.cfg.period("historical"))}
                self.template = indices.make_template(
                    self.obs_field, projection_spec, ref_period, model_fields=model_fields
                )
                self.log.info("projection template built from %s",
                              self.cfg.get_path("index.projection.template_source", "reference"))
            elif itype == "series":
                raise ConfigError("index.type='series' requires data of kind 'series'")
            self.obs_years, self.obs_series = indices.build_index(
                self.obs_field, self.cfg.get_path("index", {}), ref_period, self.template
            )
            self.obs_diag = indices.build_index(
                self.obs_field, _diagnostic_spec(self.cfg.get_path("index", {})), ref_period, self.template
            )[1]
        else:
            self.obs_years, self.obs_series = self.source.reference_series()
            self.obs_diag = self.obs_series.copy()
            self.log.info("reference series: %d years (%s-%s)", self.obs_years.size,
                          int(self.obs_years[0]), int(self.obs_years[-1]))

        self.models = self.source.discover_models()
        if not self.models:
            raise ConfigError("no models discovered; check data.models.pattern / models / series_dir")
        self.log.info("%d model(s): %s", len(self.models), ", ".join(self.models[:8]) + (" ..." if len(self.models) > 8 else ""))

        if (self.cfg.get_path("data.models.kind", "grid") == "grid"
                and str(self.cfg.get_path("data.models.member_combine", "first")) != "first"):
            spec = self.cfg.get_path("data.models", {})
            exp = self.cfg.get_path("data.models.historical_experiment", "historical")
            for model in self.models:
                try:
                    self.source.model_members[model] = io.discover_members(
                        spec, exp, model, self.cfg.period("historical"))
                except Exception:
                    pass
            if self.source.model_members:
                self.log.info("multi-member combination: %s", ", ".join(
                    f"{m} ({len(v)})" for m, v in sorted(self.source.model_members.items())))

        hist_period = self.cfg.period("historical")
        hist_mask = (self.obs_years >= hist_period[0]) & (self.obs_years <= hist_period[1])
        if hist_mask.sum() < 5:
            self.log.warning("only %d reference years inside periods.historical", int(hist_mask.sum()))
        obs_hist_years = self.obs_years[hist_mask]
        obs_hist_vals = self.obs_series[hist_mask]
        self.obs_hist = plotting.HistoricalSeries(
            years=obs_hist_years,
            values=obs_hist_vals,
            trend=indices.linear_trend(obs_hist_years, obs_hist_vals),
        )

        ref_year = self.cfg.get_path("plots.timeseries.reference_year")
        ref_value = self.cfg.get_path("plots.timeseries.reference_value")
        if ref_value is None and ref_year is not None:
            if int(ref_year) in set(self.obs_years.astype(int)):
                val = float(self.obs_series[np.where(self.obs_years == int(ref_year))[0][0]])
                self.cfg.set_path("plots.timeseries.reference_value", val)
                self.log.info("reference value for %s taken from the reference series: %.4f", ref_year, val)
            else:
                self.log.info("reference year %s is not in the reference series; no threshold line is drawn", ref_year)

        self._finished["prepare"] = {
            "models": self.models,
            "obs_years": [int(self.obs_years[0]), int(self.obs_years[-1])],
        }

    # ------------------------------------------------------------------
    # correction + index
    # ------------------------------------------------------------------
    def _template_for(self) -> Optional[xr.DataArray]:
        return self.template

    def correct(self) -> None:
        if self.obs_series is None:
            self.prepare()
        method = str(self.cfg.get_path("correction.method", "mean_std")).lower()
        params_period = self.cfg.get_path("correction.params_period") or self.cfg.period("reference_climate")
        if (self.cfg.get_path("data.models.kind", "grid") == "grid"
                and str(self.cfg.get_path("data.models.member_combine", "first")) != "first"):
            spec = self.cfg.get_path("data.models", {})
            exp = self.cfg.get_path("data.models.historical_experiment", "historical")
            for model in self.models:
                try:
                    self.source.model_members[model] = io.discover_members(
                        spec, exp, model, self.cfg.period("historical"))
                except Exception:
                    pass
            if self.source.model_members:
                self.log.info("multi-member combination: %s", ", ".join(
                    f"{m} ({len(v)})" for m, v in sorted(self.source.model_members.items())))

        hist_period = self.cfg.period("historical")
        fut_period = self.cfg.period("future")
        scenarios = self.cfg.scenario_list()
        index_spec = self.cfg.get_path("index", {})
        lower_bound = self.cfg.get_path("correction.lower_bound")
        apply_to = set(self.cfg.get_path("correction.apply_to", ["historical", "future"]) or [])

        for sc in scenarios:
            self.map_acc[sc] = {"windows": [], "sum": {}, "fut_sum": {}, "n": {}, "pos": {}, "neg": {}}

        n_ok = 0
        for mi, model in enumerate(self.models, 1):
            self.log.info("[%d/%d] %s", mi, len(self.models), model)
            try:
                self._process_model(
                    model, method, params_period, hist_period, fut_period,
                    scenarios, index_spec, lower_bound, apply_to,
                )
                n_ok += 1
            except Exception as exc:
                self.log.exception("model %s failed: %s", model, exc)
                self._failed_input(model, exc)

        self.log.info("corrected %d/%d models", n_ok, len(self.models))
        if not n_ok:
            raise ConfigError("no model completed: see logs/run.log and manifest.json for failures")
        self._save_series()
        self._finished["correct"] = {"models_done": n_ok, "scores": self.scores}

    def _process_model(
        self,
        model: str,
        method: str,
        params_period: Sequence[int],
        hist_period: Sequence[int],
        fut_period: Sequence[int],
        scenarios: Sequence[str],
        index_spec: Dict[str, Any],
        lower_bound: Optional[float],
        apply_to: set,
    ) -> None:
        cfg = self.cfg
        grid_mode = cfg.get_path("data.reference.kind", "grid") == "grid"
        rec: Dict[str, Any] = {"historical": None, "future": {}}
        raw: Dict[str, Any] = {"historical": None, "future": {}}

        if not grid_mode:
            yrs, vals = self.source.model_series(model, cfg.get_path("data.models.historical_experiment", "historical"))
            mask = (yrs >= hist_period[0]) & (yrs <= hist_period[1])
            yrs, vals = np.asarray(yrs)[mask], np.asarray(vals, dtype=float)[mask]
            raw["historical"] = (yrs, vals.copy())

            # fit on the years shared with the observation inside the fit period
            common = np.intersect1d(yrs, np.asarray(self.obs_years, dtype=int))
            common = common[(common >= params_period[0]) & (common <= params_period[1])]
            if method != "none" and common.size < 5:
                raise ValueError(f"only {common.size} common baseline years with the reference")
            obs_fit = self.obs_series[np.isin(self.obs_years, common)]
            mod_fit = vals[np.isin(yrs, common)]
            mult = correction.resolve_multiplicative(cfg, obs_fit, mod_fit)
            sparams = correction.fit_series(
                obs_fit, mod_fit, method=method, multiplicative=mult,
                n_quantiles=int(cfg.get_path("correction.n_quantiles", 100)),
            )
            sparams["lower_bound"] = lower_bound
            self.params[model] = sparams
            io.save_npz(os.path.join(self.dirs["params"], f"{model}_params.npz"),
                        **{k: np.asarray(v) for k, v in sparams.items() if not isinstance(v, str)})
            if method == "quantile_mapping":
                self.log.warning("quantile_mapping can distort future quantile changes; "
                                 "use quantile_delta_mapping for projections")

            corr_vals = self._clip(correction.apply_series(sparams, vals)) if "historical" in apply_to else vals.copy()
            rec["historical"] = (yrs, corr_vals)
            self._score_model(model, yrs, vals, yrs, corr_vals)
            for sc in scenarios:
                try:
                    yrs_f, vals_f = self.source.model_series(model, sc)
                except FileNotFoundError as exc:
                    self._failed_input(f"{model}/{sc}", exc)
                    continue
                mask = (yrs_f >= fut_period[0]) & (yrs_f <= fut_period[1])
                yrs_f = np.asarray(yrs_f)[mask]
                vals_f = np.asarray(vals_f, dtype=float)[mask]
                if not yrs_f.size or not np.isfinite(vals_f).any():
                    self._failed_input(f"{model}/{sc}", ValueError("no finite data in the configured future period"))
                    continue
                raw["future"][sc] = (yrs_f, vals_f.copy())
                rec["future"][sc] = (yrs_f, self._clip(correction.apply_series(sparams, vals_f)) if "future" in apply_to else vals_f)
            self.series[model] = rec
            self.raw_series[model] = raw
            return

        # ---- grid mode ------------------------------------------------
        obs_hist = self.obs_field.sel(year=slice(int(params_period[0]), int(params_period[1])))
        mod_hist_all = self.source.model_field(model, cfg.get_path("data.models.historical_experiment", "historical"), hist_period)
        common = np.intersect1d(obs_hist.year.values, mod_hist_all.year.values)
        common = common[(common >= params_period[0]) & (common <= params_period[1])]
        if method != "none" and common.size < 5:
            raise ValueError(f"only {common.size} common baseline years with the reference")
        obs_fit = obs_hist.sel(year=common)
        mod_fit = mod_hist_all.sel(year=common)

        params = correction.fit(
            obs_fit, mod_fit, method=method,
            multiplicative=correction.resolve_multiplicative(cfg, obs_fit.values, mod_fit.values),
            n_quantiles=int(cfg.get_path("correction.n_quantiles", 100)),
            lower_bound=lower_bound,
        )
        self.params[model] = params
        correction.save_params(params, os.path.join(self.dirs["params"], f"{model}_params.nc"))
        mult = bool(params.attrs.get("multiplicative", 0))

        # historical
        mod_hist = mod_hist_all.sel(year=slice(int(hist_period[0]), int(hist_period[1])))
        corr_hist = self._clip(correction.apply(params, mod_hist, multiplicative=mult, lower_bound=lower_bound)) \
            if "historical" in apply_to else mod_hist.rename("field")
        self._save_corrected(model, "historical", corr_hist)

        ref_clim = cfg.period("reference_climate")
        clim_raw = mod_hist.sel(year=slice(int(ref_clim[0]), int(ref_clim[1]))).mean("year", skipna=True)
        clim_corr = corr_hist.sel(year=slice(int(ref_clim[0]), int(ref_clim[1]))).mean("year", skipna=True)
        samples = self._pdf_samples(model, mod_hist, corr_hist)
        if samples is not None:
            self.pdf_data[model] = samples

        def standardization(field, climatology):
            if not index_spec.get("standardize"):
                return None
            spec = dict(index_spec, standardize=False)
            years, values = indices.build_index(field, spec, ref_clim, self._template_for(), climatology)
            mask = (years >= ref_clim[0]) & (years <= ref_clim[1]) & np.isfinite(values)
            if mask.sum() < 2:
                raise ConfigError("index.standardize needs at least two finite historical baseline values")
            return (float(np.mean(values[mask])), float(np.std(values[mask], ddof=1)))

        raw_standardization = standardization(mod_hist, clim_raw)
        corrected_standardization = standardization(corr_hist, clim_corr)
        y_raw, v_raw = indices.build_index(mod_hist, index_spec, ref_clim, self._template_for(), clim_raw,
                                           standardization=raw_standardization)
        y_cor, v_cor = indices.build_index(corr_hist, index_spec, ref_clim, self._template_for(), clim_corr,
                                           standardization=corrected_standardization)
        rec["historical"] = (y_cor, v_cor)
        raw["historical"] = (y_raw, v_raw)
        if method == "quantile_mapping":
            self.log.warning("quantile_mapping can distort future quantile changes; "
                             "use quantile_delta_mapping for projections")

        # baseline for maps
        base_mean = corr_hist.sel(year=slice(int(cfg.period("reference_climate")[0]),
                                             int(cfg.period("reference_climate")[1]))).mean("year", skipna=True)

        for sc in scenarios:
            try:
                fut = self.source.model_field(model, sc, fut_period)
            except Exception as exc:
                self._failed_input(f"{model}/{sc}", exc)
                continue
            fut = fut.sel(year=slice(int(fut_period[0]), int(fut_period[1])))
            if not fut.year.size or not np.isfinite(fut.values).any():
                self._failed_input(f"{model}/{sc}", ValueError("no finite future data in the selected period"))
                continue
            corr_fut = self._clip(correction.apply(params, fut, multiplicative=mult, lower_bound=lower_bound)) \
                if "future" in apply_to else fut.rename("field")
            self._save_corrected(model, sc, corr_fut)
            yf_raw, vf_raw = indices.build_index(fut, index_spec, ref_clim, self._template_for(), clim_raw,
                                                standardization=raw_standardization)
            yf, vf = indices.build_index(corr_fut, index_spec, ref_clim, self._template_for(), clim_corr,
                                         standardization=corrected_standardization)
            rec["future"][sc] = (yf, vf)
            raw["future"][sc] = (yf_raw, vf_raw)
            self._accumulate_maps(sc, corr_fut, base_mean)

        self.series[model] = rec
        self.raw_series[model] = raw

        # diagnostics on the absolute (non-anomaly) historical index series
        diag = _diagnostic_spec(index_spec)
        _, v_diag_raw = indices.build_index(mod_hist, diag, params_period, self._template_for())
        _, v_diag_cor = indices.build_index(corr_hist, diag, params_period, self._template_for())
        self._score_model(model, mod_hist.year.values, v_diag_raw, corr_hist.year.values, v_diag_cor)

    def _pdf_samples(self, model: str, raw_hist: xr.DataArray, corr_hist: xr.DataArray):
        """Field samples for the per-model PDF diagnostics.

        Pool every grid point of the study
        region over the baseline period, so the density is smooth and reflects the
        variable's distribution rather than that of a 65-value index.
        """
        cfg = self.cfg
        mode = str(cfg.get_path("plots.kde.mode", "auto")).lower()
        if mode == "series" or self.obs_field is None:
            return None
        lo, hi = cfg.period("reference_climate")
        region = cfg.get_path("index.region")
        if not region:
            regions = cfg.get_path("index.regions") or []
            region = regions[-1] if regions else None
        max_points = 20000
        rng = np.random.default_rng(0)

        def pool(da: xr.DataArray) -> np.ndarray:
            sel = da.sel(year=slice(int(lo), int(hi)))
            vals = sel.values.astype(np.float64)
            if region:
                mask = indices.region_mask(sel.latitude.values, sel.longitude.values, region)
                vals = vals[:, mask]
            vals = vals[np.isfinite(vals)].ravel()
            if vals.size > max_points:
                vals = rng.choice(vals, max_points, replace=False)
            return vals

        try:
            return {"obs": pool(self.obs_field), "raw": pool(raw_hist), "corr": pool(corr_hist)}
        except Exception as exc:
            self.log.debug("pdf sampling failed for %s: %s", model, exc)
            return None

    def _score_model(self, model: str, years_raw: np.ndarray, raw: np.ndarray,
                     years_cor: np.ndarray, corr: np.ndarray) -> None:
        inter = np.intersect1d(np.intersect1d(years_raw, years_cor), self.obs_years)
        validation = self.cfg.get_path("periods.validation")
        if validation is not None:
            inter = inter[(inter >= validation[0]) & (inter <= validation[1])]
        if inter.size < 10:
            return
        idx_obs = np.isin(self.obs_years, inter)
        self.scores[model] = correction.diagnostics(
            self.obs_diag[idx_obs],
            np.asarray(raw)[np.isin(years_raw, inter)],
            np.asarray(corr)[np.isin(years_cor, inter)],
        )

    def _save_corrected(self, model: str, tag: str, field: xr.DataArray) -> str:
        name = self.cfg.varname
        path = os.path.join(self.dirs["corrected"], f"{model}_{tag}_corrected.nc")
        out = field.rename(name)
        out.attrs.update({"model": model, "period": tag, "correction": str(self.cfg.get_path("correction.method"))})
        return io.save_field(out, path, name=name)

    def _accumulate_maps(self, scenario: str, future: xr.DataArray, base_mean: xr.DataArray) -> None:
        windows = self.cfg.future_windows()
        acc = self.map_acc.setdefault(scenario, {"windows": [], "sum": {}, "fut_sum": {}, "n": {}, "pos": {}, "neg": {}})
        acc["windows"] = windows
        for win in windows:
            key = f"{win[0]}-{win[1]}"
            sel = future.sel(year=slice(int(win[0]), int(win[1])))
            if sel.year.size == 0:
                continue
            fut_mean = sel.mean("year", skipna=True)
            change = (fut_mean - base_mean).values.astype(np.float64)
            if key not in acc["sum"]:
                acc["sum"][key] = np.zeros_like(change)
                acc["fut_sum"][key] = np.zeros_like(change)
                acc["n"][key] = np.zeros_like(change)
                acc["pos"][key] = np.zeros_like(change)
                acc["neg"][key] = np.zeros_like(change)
            ok = np.isfinite(change)
            acc["sum"][key] += np.where(ok, change, 0.0)
            acc["fut_sum"][key] += np.where(ok, fut_mean.values.astype(np.float64), 0.0)
            acc["n"][key] += ok.astype(float)
            acc["pos"][key] += (np.where(ok, change, 0.0) > 0).astype(float)
            acc["neg"][key] += (np.where(ok, change, 0.0) < 0).astype(float)

    def _save_series(self) -> None:
        payload: Dict[str, np.ndarray] = {
            "obs_years": np.asarray(self.obs_years, dtype=int),
            "obs_series": np.asarray(self.obs_series, dtype=float),
            "model_names": np.asarray(sorted(self.series), dtype=str),
        }
        years_hist = None
        years_fut = {}
        for model, rec in self.series.items():
            y, v = rec["historical"]
            payload[f"{model}__historical__years"] = np.asarray(y, dtype=int)
            payload[f"{model}__historical"] = np.asarray(v, dtype=float)
            if years_hist is None:
                years_hist = np.asarray(y, dtype=int)
            if np.array_equal(np.asarray(y, dtype=int), years_hist):
                payload[f"{model}__historical"] = np.asarray(v, dtype=float)
            for sc, (yf, vf) in rec["future"].items():
                payload[f"{model}__{sc}__years"] = np.asarray(yf, dtype=int)
                years_fut.setdefault(sc, np.asarray(yf, dtype=int))
                payload[f"{model}__{sc}"] = np.asarray(vf, dtype=float)
        if years_hist is not None:
            payload["years_historical"] = years_hist
        for sc, yrs in years_fut.items():
            payload[f"years_{sc}"] = yrs
        path = os.path.join(self.dirs["series"], "index_series.npz")
        io.save_npz(path, **payload)
        self.log.info("series written to %s", path)

    # ------------------------------------------------------------------
    # ensemble
    # ------------------------------------------------------------------
    def ensemble_stage(self) -> None:
        if not self.series:
            raise ConfigError("no corrected series in memory; run the 'correct' stage first")
        stat_name = self.cfg.get_path("ensemble.statistic", "mean")
        bands = [int(b) for b in self.cfg.get_path("ensemble.bands", [10, 90, 20, 80])]
        fut_period = self.cfg.period("future")
        do_trend_corr = bool(self.cfg.get_path("correction.trend_correction.enabled", False))
        tc_mode = self.cfg.get_path("correction.trend_correction.mode", "difference")
        ref_year = self.cfg.get_path("plots.timeseries.reference_year") or fut_period[0] - 1

        obs_trend = self.obs_hist.trend if self.obs_hist is not None else {}
        weights = self.cfg.get_path("ensemble.weights", "equal")

        for sc in self.cfg.scenario_list():
            series = {}
            for model, rec in self.series.items():
                if sc in rec["future"]:
                    series[model] = rec["future"][sc]
            if not series:
                self.log.warning("%s: no model series, skipping ensemble", sc)
                continue
            if len(series) < int(self.cfg.get_path("ensemble.min_models", 3)):
                self.log.warning("%s: only %d model(s) (< ensemble.min_models); bands will be degenerate",
                                 sc, len(series))
            years = np.arange(fut_period[0], fut_period[1] + 1)
            mat, names = ensemble.align(series, years)

            if do_trend_corr:
                for i, model in enumerate(names):
                    hist = self.series[model]["historical"]
                    ht = indices.linear_trend(hist[0], hist[1])
                    b_obs = obs_trend.get("slope", np.nan) if obs_trend else np.nan
                    if np.isfinite(ht["slope"]) and np.isfinite(b_obs):
                        mat[i] = correction.trend_correct_series(
                            years, mat[i], ht["slope"], b_obs, float(ref_year), tc_mode
                        )

            skill_w = None
            if weights == "skill" and self.scores:
                skill_w = ensemble.skill_weights(
                    {m: abs(self.scores[m].get("rmse_raw", np.nan)) if m in self.scores else np.nan
                     for m in names}
                )
            st = ensemble.stats(mat, statistic=stat_name, bands=bands, weights=skill_w)
            slopes = ensemble.model_slopes(mat, years)
            center_trend = indices.linear_trend(years, st["center"])
            self.ens[sc] = {
                "years": years,
                "matrix": mat,
                "names": names,
                "center": st["center"],
                "bands": {b: st[f"p{b}"] for b in bands if f"p{b}" in st},
                "n_models": int(np.sum(np.isfinite(mat).any(axis=1))),
                "n_models_by_year": np.sum(np.isfinite(mat), axis=0),
                "trend": center_trend,
                "model_slopes": slopes,
            }
            minimum = int(self.cfg.get_path("ensemble.min_models", 3))
            eligible = self.ens[sc]["n_models_by_year"] >= minimum
            self.ens[sc]["center"] = np.where(eligible, self.ens[sc]["center"], np.nan)
            self.ens[sc]["bands"] = {q: np.where(eligible, values, np.nan)
                                      for q, values in self.ens[sc]["bands"].items()}
            self.ens[sc]["trend"] = indices.linear_trend(years, self.ens[sc]["center"])
            center_trend = self.ens[sc]["trend"]
            if not eligible.any():
                raise ConfigError(f"{sc}: no year satisfies ensemble.min_models={minimum}; select more models or explicitly lower the minimum")
            self._save_ensemble(sc)
            self.log.info("%s: %d models, trend %+.4f %s/decade (p=%.3g)", sc, self.ens[sc]["n_models"],
                          center_trend["slope_per_decade"], self.cfg.index_units, center_trend["p_value"])

        if not self.ens:
            raise ConfigError("no scenario has usable model series; a successful ensemble cannot be reported")
        if obs_trend:
            self.log.info("observed trend %+.4f %s/decade (p=%.3g)", obs_trend["slope_per_decade"],
                          self.cfg.index_units, obs_trend["p_value"])
        self._finished["ensemble"] = {
            sc: {"n_models": v["n_models"], "trend_per_decade": v["trend"]["slope_per_decade"]}
            for sc, v in self.ens.items()
        }

    def _save_ensemble(self, scenario: str) -> None:
        e = self.ens[scenario]
        payload = {
            "years": e["years"],
            "center": e["center"],
            "matrix": e["matrix"],
            "model_names": np.asarray(e["names"], dtype=str),
            "n_models_by_year": e["n_models_by_year"],
            "trend_slope": e["trend"]["slope"],
            "trend_p": e["trend"]["p_value"],
            "model_slopes": e["model_slopes"],
        }
        for q, arr in e["bands"].items():
            payload[f"p{q}"] = arr
        io.save_npz(os.path.join(self.dirs["ensemble"], f"ensemble_{scenario}.npz"), **payload)

    # ------------------------------------------------------------------
    # extremes
    # ------------------------------------------------------------------
    def extremes_stage(self) -> None:
        conf = self.cfg.get_path("extremes", {})
        if not conf.get("enabled", False):
            self.log.info("extremes disabled")
            return
        if not self.ens:
            raise ConfigError("no ensemble in memory; run the 'ensemble' stage first")
        dist = str(conf.get("distribution", "gev"))
        method = str(conf.get("fit_method", "mle"))
        min_years = int(conf.get("min_years", 20))
        fit_window = conf.get("fit_window") or {}
        hist_lo, hist_hi = fit_window.get("historical", self.cfg.period("historical"))
        fut_lo, fut_hi = fit_window.get("future", self.cfg.period("future"))
        source = str(conf.get("source", "future"))
        table_periods = [float(t) for t in conf.get("return_periods", [])]
        periods = np.unique(np.concatenate([np.logspace(np.log10(1.1), np.log10(max(table_periods + [100.])), 120),
                                            np.asarray(table_periods + [100.])]))
        ref_value = conf.get("reference_value")
        if ref_value is None:
            ref_value = self.cfg.get_path("plots.timeseries.reference_value")

        # historical observed fit
        obs_mask = (self.obs_years >= int(hist_lo)) & (self.obs_years <= int(hist_hi))
        hist_fit = extremes.fit_distribution(self.obs_series[obs_mask], dist, method, min_years)
        hist_curve = hist_fit.return_levels(periods) if hist_fit else np.array([])
        hist_ref_rp = hist_fit.return_period(float(ref_value)) if (hist_fit and ref_value is not None) else np.nan
        self.gev["historical"] = {
            "fit": hist_fit, "periods": periods, "levels": hist_curve, "reference_period": hist_ref_rp,
        }
        if hist_fit:
            self.log.info("historical GEV: c=%.3f loc=%.3f scale=%.3f (n=%d)", hist_fit.shape, hist_fit.loc,
                          hist_fit.scale, hist_fit.n)
        elif source in ("historical", "both"):
            raise ConfigError("requested historical extreme fit could not be estimated; check sample size and variability")

        shapes: Dict[str, List[float]] = {}
        for sc in self.cfg.scenario_list():
            if source == "historical":
                break
            if sc not in self.ens:
                continue
            years = self.ens[sc]["years"]
            mat = self.ens[sc]["matrix"]
            win = (years >= int(fut_lo)) & (years <= int(fut_hi))
            mat = mat[:, win]
            years = years[win]
            fits = {}
            for i, model in enumerate(self.ens[sc]["names"]):
                y = mat[i]
                ok = np.isfinite(y)
                fit = extremes.fit_distribution(y[ok], dist, method, min_years)
                if fit is not None:
                    fits[model] = fit
                    fits[model].extra["return_period_ref"] = fit.return_period(float(ref_value)) if ref_value is not None else np.nan
                else:
                    self.failures.append({"dataset": f"extremes {model}/{sc}",
                                          "reason": "insufficient finite years or unsuccessful distribution fit"})
            if not fits:
                self._failed_input(f"extremes {sc}", ConfigError("no valid extreme-value fits for the requested scenario"))
                continue
            if len(fits) < int(self.cfg.get_path("ensemble.min_models", 3)):
                self._failed_input(f"extremes {sc}", ConfigError("fewer valid fits than ensemble.min_models"))
            if source == "historical":
                self.log.info("extremes.source=historical: scenario curves are not reported")
                continue
            curves = extremes.return_curves(fits, periods)
            st = extremes.curve_stats(curves)
            shapes[sc] = [f.shape for f in fits.values()]
            y_med = np.nanmedian(mat, axis=0)
            emp = extremes.empirical_return_levels(y_med) if conf.get("show_empirical", False) else None
            self.gev[sc] = {
                "fits": fits,
                "periods": periods,
                "curves": curves,
                "stats": st,
                "table_periods": table_periods,
                "reference_period": _median_finite([f.return_period(float(ref_value)) for f in fits.values()])
                if ref_value is not None else np.nan,
                "empirical": emp,
            }
            io.save_npz(
                os.path.join(self.dirs["extremes"], f"gev_{sc}.npz"),
                periods=periods,
                curves=curves,
                median=st["median"],
                model_names=np.asarray(sorted(fits), dtype=str),
                shapes=np.asarray([fits[m].shape for m in sorted(fits)]),
                reference_value=np.asarray([ref_value if ref_value is not None else np.nan]),
            )
            self.log.info("%s GEV: %d fits, median return level(100a)=%.3f", sc, len(fits),
                          _median_finite(f.return_level(100.) for f in fits.values()))

        if source != "historical" and not any("fits" in value for value in self.gev.values()):
            raise ConfigError("extreme analysis requested but no scenario produced a valid fit")
        kruskal = extremes.kruskal_scenarios(shapes) if conf.get("shape_comparison", False) else {}
        self.shapes = shapes
        self.kruskal = kruskal
        self._save_gev_tables()
        self._finished["extremes"] = {
            "historical_shape": hist_fit.shape if hist_fit else None,
            "kruskal_p": kruskal.get("p_value", np.nan) if kruskal else None,
        }

    def _save_gev_tables(self) -> None:
        rows = []
        for sc, data in self.gev.items():
            if sc == "historical":
                fit = data.get("fit")
                if fit:
                    rows.append({"scenario": "historical", "model": "OBS", **fit.as_dict(),
                                 "return_period_ref": data.get("reference_period", np.nan)})
                continue
            for model, fit in data.get("fits", {}).items():
                rows.append({"scenario": sc, "model": model, **fit.as_dict(),
                             "return_period_ref": fit.extra.get("return_period_ref", np.nan)})
        if not rows:
            return
        keys = ["scenario", "model", "distribution", "shape", "shape_xi", "loc", "scale", "n", "method",
                "loglik", "aic", "ks_stat", "ks_p", "ks_p_naive", "return_period_ref"]
        path = os.path.join(self.dirs["tables"], "gev_fits.csv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(",".join(keys) + "\n")
            for row in rows:
                fh.write(",".join(_fmt(row.get(k, ""), 6) for k in keys) + "\n")
        rows_rl = []
        for sc, data in self.gev.items():
            if sc == "historical" or "fits" not in data:
                continue
            tp = data.get("table_periods") or []
            for model, fit in data["fits"].items():
                row = {"scenario": sc, "model": model}
                for period in tp:
                    row[f"T{period:g}"] = fit.return_level(period)
                rows_rl.append(row)
        if rows_rl:
            keys = list(rows_rl[0])
            path = os.path.join(self.dirs["tables"], "gev_return_levels.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(",".join(keys) + "\n")
                for row in rows_rl:
                    cells = [str(row["scenario"]), str(row["model"])]
                    cells += [_fmt(row.get(k, np.nan), 4) for k in keys[2:]]
                    fh.write(",".join(cells) + "\n")

        shapes = self.shapes
        if shapes:
            path = os.path.join(self.dirs["tables"], "gev_shape_summary.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("scenario,n,mean,median,std,min,max\n")
                for sc, vals in shapes.items():
                    summ = extremes.shape_summary({str(i): f for i, f in enumerate(vals)})
                    fh.write(",".join([sc] + [_fmt(summ.get(k, np.nan), 4) for k in ("n", "mean", "median", "std", "min", "max")]) + "\n")


    # ------------------------------------------------------------------
    # daily-scale extremes (POT/GPD or daily block maxima)
    # ------------------------------------------------------------------
    def daily_stage(self) -> None:
        conf = self.cfg.get_path("daily", {}) or {}
        if not conf.get("enabled", False):
            return
        if (self.cfg.get_path("data.reference.kind") != "grid"
                or self.cfg.get_path("data.models.kind") != "grid"):
            self.log.warning("daily extremes skipped: gridded input is required")
            return
        ext = conf.get("extremes", {}) or {}
        kind = str(ext.get("kind", "pot")).lower()
        months = resolve_months(conf.get("months", [6, 7, 8]))
        hist_years = self.cfg.period("historical")
        fut_years = self.cfg.period("future")
        base = conf.get("baseline") or self.cfg.period("reference_climate")
        method = str(conf.get("correction", "quantile_delta_mapping"))
        index_spec = dict(self.cfg.get_path("index", {}) or {})
        index_spec["type"] = conf.get("index", "region_mean")
        ref_spec = self.cfg.get_path("data.reference", {})
        mod_spec = self.cfg.get_path("data.models", {})
        hist_exp = self.cfg.get_path("data.models.historical_experiment", "historical")
        min_years = int(ext.get("min_years", 20))
        min_peaks = int(ext.get("min_peaks", 30))
        decluster_days = int(ext.get("decluster_days", 3))
        threshold_mode = str(ext.get("threshold_mode", "observed")).lower()
        threshold_q = ext.get("threshold_quantile", 0.95)
        table_periods = [float(t) for t in ext.get("return_periods", [2, 5, 10, 20, 50, 100])]
        periods = np.unique(np.concatenate([np.logspace(np.log10(1.1), np.log10(max(table_periods + [100.])), 120),
                                            np.asarray(table_periods + [100.])]))

        def years_of(dates):
            return np.asarray([getattr(t, "year", int(str(t)[:4])) for t in dates], dtype=int)

        def mask_between(dates, lo, hi):
            yrs = io._season_year(dates, months)
            return (yrs >= int(lo)) & (yrs <= int(hi))

        def series_threshold(values) -> float:
            if threshold_mode == "absolute":
                return float(ext.get("threshold_abs"))
            if threshold_mode == "quantile":
                return daily_mod.threshold_value(values, threshold_q, None)
            return float(threshold)

        def fit_series_extreme(dates, values):
            if kind == "block":
                return daily_mod.fit_block(dates, values, str(ext.get("distribution", "gev")), min_years,
                                           months=months)
            yrs = years_of(dates)
            n_years = float(len(np.unique(yrs))) if yrs.size else 1.0
            exposure = daily_mod.exposure_years(dates, values, months=months)
            if exposure < min_years:
                self.log.warning("daily exposure %.2f seasons is below min_years=%d", exposure, min_years)
                return None
            return daily_mod.fit_pot(values, series_threshold(values), n_years, decluster_days, min_peaks,
                                     dates=dates, months=months)

        def record_inputs(label, field):
            filenames = json.loads(field.attrs.get("source_files", "[]"))
            records = []
            for filename in filenames:
                stat = os.stat(filename)
                records.append({"path": os.path.abspath(filename), "size_bytes": stat.st_size,
                                "mtime_ns": stat.st_mtime_ns})
            self.daily_input_inventory.append({"label": label, "files_metadata": records})

        self.log.info("daily extremes: %s, months=%s, %s-%s", kind.upper(), months, hist_years[0], hist_years[-1])
        obs_field = io.load_daily_field(ref_spec, hist_years, months)
        record_inputs("daily reference", obs_field)
        obs_dates = np.asarray(obs_field.time.values)
        obs_values = io.daily_index_series(obs_field, index_spec)
        del obs_field
        base_mask = mask_between(obs_dates, base[0], base[1])
        threshold = daily_mod.threshold_value(
            obs_values[base_mask],
            threshold_q,
            ext.get("threshold_abs") if threshold_mode == "absolute" else None,
        )
        if not np.isfinite(threshold):
            raise ConfigError("daily extremes: could not determine a threshold from the reference series")
        self.log.info("daily threshold (%s): %.4f %s", threshold_mode, threshold, self.cfg.units)

        obs_fit = fit_series_extreme(obs_dates[base_mask], obs_values[base_mask])
        record_value = float(np.nanmax(obs_values[base_mask])) if base_mask.any() else np.nan
        obs_rp = obs_fit.return_period(record_value) if obs_fit else np.nan
        if obs_fit:
            self.log.info("historical %s fit: shape=%.3f scale=%.3f%s", kind.upper(), obs_fit.shape,
                          obs_fit.scale, f" rate={obs_fit.rate:.2f}/yr" if kind == "pot" else "")

        scenarios: Dict[str, Any] = {}
        rows: List[Dict[str, Any]] = []
        for model in self.models:
            member_spec = dict(mod_spec)
            if self.source.model_members.get(model):
                member_spec["members"] = list(self.source.model_members[model])
            try:
                historical_files = io._files_for_grid_source(member_spec, hist_exp, model,
                                                               io._input_years(hist_years, months))
                actual_members = {io._member_of(path, model, hist_exp, member_spec) for path in historical_files}
                actual_members.discard(None)
                if len(actual_members) > 1:
                    raise ConfigError("daily historical files mix multiple members")
                if actual_members:
                    member_spec["members"] = sorted(actual_members)
                    self.source.model_members[model] = sorted(actual_members)
                hist_field = io.load_daily_field(member_spec, hist_years, months, experiment=hist_exp, model=model)
                record_inputs(f"daily {model}/{hist_exp}", hist_field)
                hist_dates = np.asarray(hist_field.time.values)
                hist_values = io.daily_index_series(hist_field, index_spec)
                del hist_field
            except Exception as exc:
                self._failed_input(f"daily {model}/historical", exc)
                continue
            mod_base = mask_between(hist_dates, base[0], base[1])
            if mod_base.sum() < 365:
                self.failures.append({"dataset": f"daily {model}/historical", "reason": "too few baseline days"})
                self.log.warning("daily: %s has too few baseline days", model)
                continue
            try:
                params = correction.fit_series(
                    obs_values[base_mask], hist_values[mod_base], method=method,
                    multiplicative=correction.resolve_multiplicative(
                        self.cfg, obs_values[base_mask], hist_values[mod_base]),
                    n_quantiles=int(self.cfg.get_path("correction.n_quantiles", 100)),
                )
            except Exception as exc:
                self._failed_input(f"daily correction {model}", exc)
                continue

            for sc in self.cfg.scenario_list():
                try:
                    fut_field = io.load_daily_field(member_spec, fut_years, months, experiment=sc, model=model)
                    record_inputs(f"daily {model}/{sc}", fut_field)
                    fut_dates = np.asarray(fut_field.time.values)
                    fut_values = io.daily_index_series(fut_field, index_spec)
                    del fut_field
                except Exception as exc:
                    self._failed_input(f"daily {model}/{sc}", exc)
                    continue
                fut_values = self._clip(correction.apply_series(params, fut_values))
                fit = fit_series_extreme(fut_dates, fut_values)
                if fit is None:
                    self.failures.append({"dataset": f"daily {model}/{sc}",
                                          "reason": "insufficient exposure/peaks/complete blocks or unsuccessful fit"})
                    self.log.warning("daily: %s/%s has no valid %s fit", model, sc, kind.upper())
                    continue
                scenarios.setdefault(sc, {})[model] = fit
                rows.append({
                    "scenario": sc, "model": model, **fit.as_dict(),
                    "return_level_100": fit.return_level(100.0),
                    "return_period_of_record": fit.return_period(record_value) if np.isfinite(record_value) else np.nan,
                })

        if not scenarios:
            raise ConfigError("daily extreme analysis requested but no scenario produced a valid fit")
        for scenario in self.cfg.scenario_list():
            if len(scenarios.get(scenario, {})) < int(self.cfg.get_path("ensemble.min_models", 3)):
                self._failed_input(f"daily {scenario}", ConfigError("fewer valid fits than ensemble.min_models"))

        self.daily = {
            "kind": kind,
            "threshold": threshold,
            "periods": periods,
            "historical": obs_fit,
            "record_value": record_value,
            "historical_return_period": obs_rp,
            "scenarios": {},
            "table_periods": table_periods,
        }
        for sc, fits in scenarios.items():
            curves = np.vstack([fits[m].return_levels(periods) for m in sorted(fits)])
            stats = daily_mod.curve_stats(curves)
            rp = _median_finite([f.return_period(record_value) for f in fits.values()]) \
                if np.isfinite(record_value) else np.nan
            self.daily["scenarios"][sc] = {
                "fits": fits, "curves": curves, "stats": stats, "reference_period": rp,
                "n_models": len(fits),
            }
            io.save_npz(
                os.path.join(self.dirs["extremes"], f"daily_{sc}.npz"),
                periods=periods, curves=curves, median=stats["median"],
                model_names=np.asarray(sorted(fits), dtype=str),
                threshold=np.asarray([threshold]),
                kind=np.asarray([kind]),
                record_value=np.asarray([record_value]),
            )
        if rows:
            keys = list(rows[0])
            path = os.path.join(self.dirs["tables"], "daily_extremes.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(",".join(keys) + "\n")
                for row in rows:
                    fh.write(",".join(
                        str(row[k]) if k in ("scenario", "model", "kind") else _fmt(row[k], 6) for k in keys
                    ) + "\n")
        self._finished["daily"] = {
            "kind": kind,
            "threshold": threshold,
            "scenarios": {sc: v["n_models"] for sc, v in self.daily["scenarios"].items()},
        }
        self.log.info("daily extremes: %d scenario(s) fitted", len(self.daily["scenarios"]))

    def _daily_figure(self) -> List[str]:
        if not self.daily:
            return []
        periods = self.daily["periods"]
        hist_fit = self.daily["historical"]
        hist_obj = None
        if hist_fit is not None:
            hist_obj = plotting.GEVHistorical(
                periods=periods,
                levels=hist_fit.return_levels(periods),
                reference_period=self.daily.get("historical_return_period", np.nan),
            )
        objs = []
        for sc, data in self.daily["scenarios"].items():
            st = data["stats"]
            objs.append(plotting.GEVScenario(
                scenario=sc,
                label=f"{self.cfg.scenario_label(sc)} — daily {self.daily['kind'].upper()}",
                periods=periods, median=st["median"],
                bands={q: st[f"p{q}"] for q in (10, 20, 80, 90) if f"p{q}" in st},
                n_models=data["n_models"], reference_period=data.get("reference_period", np.nan),
            ))
        daily_plot_cfg = Config(self.cfg.to_dict())
        daily_plot_cfg.set_path("index.standardize", False)
        daily_plot_cfg.set_path("index.type", self.cfg.get_path("daily.index"))
        daily_plot_cfg.set_path("index.units_label", None)
        daily_plot_cfg.set_path("plots.timeseries.reference_value", self.daily.get("record_value"))
        return plotting.gev(daily_plot_cfg, hist_obj, objs,
                            os.path.join(self.dirs["figures"], "fig_daily_return_levels"))

    # ------------------------------------------------------------------
    # figures
    # ------------------------------------------------------------------
    def figures(self) -> None:
        figs = self.cfg.get_path("plots.figures", []) or []
        extra_text = " ".join(str(v) for v in self.cfg.get_path("variable", {}).values())
        plotting.apply_style(self.cfg, extra_text)
        outdir = self.dirs["figures"]
        produced: List[str] = []

        def series_objects():
            out = []
            for sc, e in self.ens.items():
                out.append(plotting.ScenarioSeries(
                    scenario=sc, label=self.cfg.scenario_label(sc), years=e["years"],
                    matrix=e["matrix"], center=e["center"], bands=e["bands"],
                    trend=e["trend"], n_models=e["n_models"],
                ))
            return out

        def gev_objects():
            out = []
            for sc, g in self.gev.items():
                if sc == "historical" or "stats" not in g:
                    continue
                st = g["stats"]
                out.append(plotting.GEVScenario(
                    scenario=sc, label=self.cfg.scenario_label(sc), periods=g["periods"],
                    median=st["median"], bands={q: st[f"p{q}"] for q in (10, 20, 80, 90) if f"p{q}" in st},
                    n_models=st.get("n_models", 0), reference_period=g.get("reference_period", np.nan),
                    empirical=g.get("empirical"),
                ))
            return out

        def gev_hist_object():
            h = self.gev.get("historical")
            if not h or h.get("fit") is None:
                return None
            return plotting.GEVHistorical(periods=h["periods"], levels=h["levels"],
                                          reference_period=h.get("reference_period", np.nan))

        if "timeseries" in figs and self.ens:
            produced += plotting.timeseries(self.cfg, self.obs_hist, series_objects(),
                                            os.path.join(outdir, "fig_timeseries"))
            if self.cfg.get_path("plots.timeseries.per_model", False):
                produced += self._per_model_series()
        if "gev" in figs and any(k != "historical" and "stats" in self.gev[k] for k in self.gev):
            produced += plotting.gev(self.cfg, gev_hist_object(), gev_objects(),
                                     os.path.join(outdir, "fig_gev_return_levels"))
        if "combined" in figs and self.ens and any(k != "historical" and "stats" in self.gev[k] for k in self.gev):
            produced += plotting.scenarios_combined(self.cfg, self.obs_hist, series_objects(),
                                                    gev_hist_object(), gev_objects(),
                                                    os.path.join(outdir, "fig_scenarios_combined"))
        if "maps" in figs and self.map_acc:
            produced += self._map_figures()
        if "kde" in figs and self.scores:
            produced += self._kde_figure()
        if self.daily and self.cfg.get_path("daily.plots", True):
            produced += self._daily_figure()
        if "shape_box" in figs and self.shapes:
            labels = {sc: self.cfg.scenario_label(sc) for sc in self.shapes}
            produced += plotting.shape_box(self.cfg, self.shapes, labels,
                                           os.path.join(outdir, "fig_gev_shape"),
                                           kruskal=self.kruskal)
        self._finished["figures"] = produced
        self.log.info("%d figure file(s) written to %s", len(produced), outdir)

    def _per_model_series(self) -> List[str]:
        produced = []
        ens_by_sc = {}
        for sc, e in self.ens.items():
            ens_by_sc[sc] = plotting.ScenarioSeries(
                scenario=sc, label=self.cfg.scenario_label(sc), years=e["years"], matrix=e["matrix"],
                center=e["center"], bands=e["bands"], trend=e["trend"], n_models=e["n_models"],
            )
        for model, rec in sorted(self.series.items()):
            if not rec["future"]:
                continue
            sers = {}
            for sc, (years, values) in rec["future"].items():
                values = np.asarray(values, dtype=float)
                sers[sc] = plotting.ScenarioSeries(
                    scenario=sc, label=self.cfg.scenario_label(sc), years=np.asarray(years, dtype=int),
                    matrix=values[None, :], center=values, bands={},
                    trend=indices.linear_trend(years, values), n_models=1,
                )
            produced += plotting.per_model_timeseries(
                self.cfg, self.obs_hist, sers, ens_by_sc, model,
                os.path.join(self.dirs["figures"], f"fig_timeseries_{model}"),
            )
        return produced

    def _map_figures(self) -> List[str]:
        produced = []
        acc0 = next(iter(self.map_acc.values()))
        windows = acc0.get("windows") or self.cfg.future_windows()
        grid = io.resolve_target_grid(self.cfg, self.obs_field)
        lon, lat = grid["longitude"], grid["latitude"]
        agree_thr = float(self.cfg.get_path("plots.maps.agreement", 0.8))
        mode = str(self.cfg.get_path("plots.maps.variable", "change")).lower()
        for sc, acc in self.map_acc.items():
            panels, labels = [], []
            for win in windows:
                key = f"{win[0]}-{win[1]}"
                if key not in acc["sum"]:
                    continue
                n = np.where(acc["n"][key] >= int(self.cfg.get_path("ensemble.min_models", 3)), acc["n"][key], np.nan)
                mean = acc["fut_sum"][key] / n if mode == "absolute" else acc["sum"][key] / n
                agree = np.maximum(acc["pos"][key], acc["neg"][key]) / n
                stipple = np.isfinite(mean) & (agree >= agree_thr)
                panels.append(plotting.MapPanel(label=f"{key}", field=mean, lon=lon, lat=lat, stipple=stipple))
                labels.append(key)
            if not panels:
                continue
            # one row per 2 panels, single scenario per figure
            rows = [panels[i:i + 2] for i in range(0, len(panels), 2)]
            row_labels = [f"{r[0].label}" for r in rows]
            produced += plotting.maps(
                self.cfg, rows, row_labels, [sc], os.path.join(self.dirs["figures"], f"fig_maps_{sc}"),
                cbar_label=(f"{self.cfg.varname} ({self.cfg.units})" if mode == "absolute"
                            else f"{self.cfg.varname} change ({self.cfg.units})"),
                symmetric=(mode != "absolute"),
                extent=self.cfg.get_path("plots.maps.extent"),
            )
        return produced

    def _kde_figure(self) -> List[str]:
        records = []
        for model in sorted(self.scores):
            data = self.pdf_data.get(model)
            if data is not None and data["obs"].size > 5:
                ks_r, _ = correction.ks_statistic(data["obs"], data["raw"])
                ks_c, _ = correction.ks_statistic(data["obs"], data["corr"])
                records.append({"model": model, "obs": data["obs"], "raw": data["raw"],
                                "corr": data["corr"], "ks_raw": ks_r, "ks_corr": ks_c})
                continue
            rec = self.raw_series.get(model, {})
            cor = self.series.get(model, {})
            if not rec.get("historical") or not cor.get("historical"):
                continue
            y_r, v_r = rec["historical"]
            y_c, v_c = cor["historical"]
            common = np.intersect1d(np.intersect1d(y_r, y_c), self.obs_years)
            if common.size < 10:
                continue
            records.append({
                "model": model,
                "obs": self.obs_series[np.isin(self.obs_years, common)],
                "raw": v_r[np.isin(y_r, common)],
                "corr": v_c[np.isin(y_c, common)],
                "ks_raw": self.scores[model].get("ks_raw", np.nan),
                "ks_corr": self.scores[model].get("ks_corr", np.nan),
            })
        if not records:
            return []
        produced = plotting.kde(self.cfg, records, os.path.join(self.dirs["figures"], "fig_diagnostics_kde"))
        if self.cfg.get_path("plots.kde.per_model", False):
            lo, hi = self.cfg.period("reference_climate")
            for rec in records:
                produced += plotting.pdf_single(
                    self.cfg, rec, os.path.join(self.dirs["figures"], f"fig_pdf_{rec['model']}"),
                    subtitle=f"{lo}-{hi}",
                )
        return produced

    # ------------------------------------------------------------------
    # report
    # ------------------------------------------------------------------
    def write_report(self) -> None:
        tables = self._write_tables()
        path = os.path.join(self.outdir, "report.md")
        cfg = self.cfg
        lines: List[str] = []
        lines.append(f"# {self.name}\n")
        lines.append(f"Generated by cmipkit on {time.strftime('%Y-%m-%d %H:%M')}  ")
        lines.append(f"Output directory: `{self.outdir}`\n")
        lines.append("## Configuration\n")
        lines.append(f"- Variable: **{cfg.varname}** ({cfg.units})")
        lines.append(f"- Season / months: {cfg.months('reference')}")
        lines.append(f"- Index: `{cfg.get_path('index.type')}`")
        lines.append(f"- Bias correction: `{cfg.get_path('correction.method')}` fitted on "
                     f"{cfg.get_path('correction.params_period') or cfg.period('reference_climate')}")
        lines.append(f"- Periods: historical {cfg.period('historical')}, future {cfg.period('future')}")
        lines.append(f"- Scenarios: {', '.join(cfg.scenario_label(s) for s in cfg.scenario_list())}")
        lines.append(f"- Stages completed before this report: {', '.join(self._finished)}")
        lines.append(f"- Validation window: {cfg.get_path('periods.validation') or 'not specified; historical diagnostics may overlap calibration'}")
        members = self.source.model_members
        if members:
            members_txt = ", ".join(
                f"{m} ({len(members[m])} members)" if len(members.get(m, [])) > 1 else m
                for m in self.models
            )
            lines.append(f"- Models ({len(self.models)}, member combination = "
                         f"`{cfg.get_path('data.models.member_combine')}`): {members_txt}\n")
        else:
            lines.append(f"- Models ({len(self.models)}): {', '.join(self.models)}\n")

        if self.scores:
            lines.append("## Historical diagnostics against the reference\n")
            lines.append("Calibration agreement is not independent validation. Year-by-year RMSE, correlation and KGE "
                         "are descriptive for free-running models whose internal variability is not synchronized to observations.\n")
            lines.append("| model | KS raw | KS corr | RMSE raw | RMSE corr | bias raw | bias corr | KGE corr |")
            lines.append("|---|---|---|---|---|---|---|---|")
            for m in sorted(self.scores):
                s = self.scores[m]
                lines.append(f"| {m} | {_fmt(s['ks_raw'])} | {_fmt(s['ks_corr'])} | {_fmt(s['rmse_raw'])} | "
                             f"{_fmt(s['rmse_corr'])} | {_fmt(s['bias_raw'])} | {_fmt(s['bias_corr'])} | "
                             f"{_fmt(s['kge_corr'])} |")
            lines.append("")

        if self.ens:
            lines.append("## Interpretation of ensemble results\n")
            lines.append("Percentile bands describe the available model distribution; they are not automatically "
                         "confidence intervals or calibrated probabilities. Model counts can vary by year; "
                         "see `tables/ensemble_coverage.csv`. Years below `ensemble.min_models` are masked. "
                         "OLS trend p-values do not account for serial correlation.\n")
            lines.append("## Trends (units per decade)\n")
            lines.append(f"Observed: {_fmt(self.obs_hist.trend.get('slope_per_decade'))} {cfg.index_units}/10a "
                         f"(p={_fmt(self.obs_hist.trend.get('p_value'))})\n")
            lines.append("| scenario | models | ensemble trend | p | model slope range |")
            lines.append("|---|---|---|---|---|")
            for sc, e in self.ens.items():
                tr = e["trend"]
                sl = e["model_slopes"][np.isfinite(e["model_slopes"])] * 10.0
                rng = f"{_fmt(np.min(sl))} .. {_fmt(np.max(sl))}" if sl.size else "nan"
                lines.append(f"| {cfg.scenario_label(sc)} | {e['n_models']} | {_fmt(tr['slope_per_decade'])} | "
                             f"{_fmt(tr['p_value'])} | {rng} |")
            lines.append("")

        if any(k != "historical" and "stats" in self.gev[k] for k in self.gev):
            lines.append("## GEV fits\n")
            lines.append(f"Sample kind: `{cfg.get_path('extremes.sample_kind')}`. An annual index fit is a "
                         "descriptive distribution assumption unless the samples are justified independent block maxima. "
                         "The fitted window is treated as stationary; warming within it can violate that assumption. "
                         "SciPy shape c = -xi. `ks_p_naive` is uncalibrated after parameter fitting; "
                         "`ks_p` remains missing until a suitable refitting test is implemented.\n")
            h = self.gev.get("historical", {}).get("fit")
            if h:
                lines.append(f"Historical ({h.n} years): c={_fmt(h.shape)}, loc={_fmt(h.loc)}, scale={_fmt(h.scale)}\n")
            lines.append("| scenario | fits | median c | return level 100a | median return period of reference value |")
            lines.append("|---|---|---|---|---|")
            for sc, g in self.gev.items():
                if sc == "historical" or "fits" not in g:
                    continue
                shapes = [f.shape for f in g["fits"].values()]
                idx = int(np.argmin(np.abs(g["periods"] - 100)))
                rp = g.get("reference_period")
                rp_txt = _fmt(rp) if rp is not None and np.isfinite(rp) else "n/a"
                lines.append(f"| {cfg.scenario_label(sc)} | {len(g['fits'])} | {_fmt(np.median(shapes))} | "
                             f"{_fmt(g['stats']['median'][idx])} | {rp_txt} |")
            if self.kruskal and len(self.shapes) >= 2:
                k = self.kruskal
                lines.append(f"\nKruskal–Wallis on the shape parameter across scenarios: "
                             f"H={_fmt(k.get('H'))}, p={_fmt(k.get('p_value'))}\n")

        if self.daily:
            kind = self.daily["kind"].upper()
            lines.append(f"## Daily extremes ({kind})\n")
            lines.append("POT return periods are mean waiting times between independent threshold events; "
                         "block-maxima return periods invert annual block exceedance probability. "
                         "Daily exposure uses actual dates and the selected season's calendar.\n")
            lines.append(f"- Threshold: {_fmt(self.daily['threshold'])} {cfg.units}"
                         + (f", runs declustered over {cfg.get_path('daily.extremes.decluster_days')} days"
                            if self.daily["kind"] == "pot" else ""))
            if np.isfinite(self.daily.get("record_value", np.nan)):
                lines.append(f"- Reference record (observed maximum of the baseline): "
                             f"{_fmt(self.daily['record_value'])} {cfg.units}, return period "
                             f"{_fmt(self.daily.get('historical_return_period'))} a")
            lines.append("")
            lines.append("| scenario | models | median shape | return level 100a | median return period of the record |")
            lines.append("|---|---|---|---|---|")
            idx = int(np.argmin(np.abs(np.asarray(self.daily["periods"]) - 100.0)))
            for sc, data in self.daily["scenarios"].items():
                shapes = [f.shape for f in data["fits"].values()]
                lines.append(f"| {cfg.scenario_label(sc)} | {data['n_models']} | {_fmt(np.median(shapes))} | "
                             f"{_fmt(data['stats']['median'][idx])} | {_fmt(data.get('reference_period'))} |")
            lines.append("")

        if self.failures:
            lines.append("## Skipped or failed datasets\n")
            for failure in self.failures:
                lines.append(f"- {failure['dataset']}: {failure['reason']}")
            lines.append("")
        if self._finished.get("figures"):
            lines.append("")
            lines.append("## Figures\n")
            for p in self._finished["figures"]:
                lines.append(f"- `{os.path.relpath(p, self.outdir)}`")
            lines.append("")

        if tables:
            lines.append("## Tables\n")
            for name, p in tables.items():
                lines.append(f"- [{name}]({os.path.relpath(p, self.outdir)})")
            lines.append("")

        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        self._finished["report"] = path
        self.log.info("report written to %s", path)

    def _write_tables(self) -> Dict[str, str]:
        out: Dict[str, str] = {}
        if self.scores:
            keys = ["ks_raw", "ks_corr", "rmse_raw", "rmse_corr", "bias_raw", "bias_corr",
                    "sd_ratio_raw", "sd_ratio_corr", "corr_raw", "corr_corr", "kge_raw", "kge_corr"]
            path = os.path.join(self.dirs["tables"], "model_scores.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("model," + ",".join(keys) + "\n")
                for m in sorted(self.scores):
                    fh.write(m + "," + ",".join(_fmt(self.scores[m].get(k, np.nan), 5) for k in keys) + "\n")
            out["model_scores"] = path

        if self.ens:
            coverage_path = os.path.join(self.dirs["tables"], "ensemble_coverage.csv")
            with open(coverage_path, "w", encoding="utf-8") as fh:
                fh.write("scenario,year,n_models,meets_minimum\n")
                for sc, ensemble_result in self.ens.items():
                    for year, count in zip(ensemble_result["years"], ensemble_result["n_models_by_year"]):
                        fh.write(f"{sc},{year},{count},{int(count >= self.cfg.get_path('ensemble.min_models', 3))}\n")
            out["ensemble_coverage"] = coverage_path
            path = os.path.join(self.dirs["tables"], "trends.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("scenario,series,slope_per_year,slope_per_decade,p_value,r2\n")
                tr = self.obs_hist.trend
                fh.write(f"historical,obs,{_fmt(tr['slope'],6)},{_fmt(tr['slope_per_decade'],6)},{_fmt(tr['p_value'],6)},{_fmt(tr['r2'],4)}\n")
                for sc, e in self.ens.items():
                    for model, slope in zip(e["names"], e["model_slopes"]):
                        fh.write(f"{sc},{model},{_fmt(slope,6)},{_fmt(slope * 10,6)},,\n")
                    t = e["trend"]
                    fh.write(f"{sc},ensemble_{self.cfg.get_path('ensemble.statistic')},{_fmt(t['slope'],6)},{_fmt(t['slope_per_decade'],6)},{_fmt(t['p_value'],6)},{_fmt(t['r2'],4)}\n")
            out["trends"] = path
        return out
