"""Regression checks for the portable CLI and stage/output contracts (unittest)."""
from __future__ import annotations
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import yaml
import xarray as xr

SKILL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL / "scripts"))
from cmipkit.config import ConfigError, load_config, plan_stages, resolve_months
from cmipkit.pipeline import Pipeline
from cmipkit.preflight import inspect_config

CLI = SKILL / "scripts" / "run.py"


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cmip_workflow_")
        self.root = Path(self.tmp.name)
        (self.root / "data").mkdir()
        (self.root / "elsewhere").mkdir()
        self.hist = np.arange(1981, 2001)
        self.future = np.arange(2031, 2051)
        self.observed = 10.0 + np.arange(20) * 0.2
        self.write_series("observations.csv", self.hist, self.observed)
        for i, model in enumerate(("A", "B", "C")):
            self.write_series(f"{model}_historical.csv", self.hist, self.observed + i + 1)
            years = self.future[self.future != 2035] if model == "B" else self.future
            self.write_series(f"{model}_ssp245.csv", years, 15 + 0.2 * (years - 2031) + i + 1)
        self.config = {
            "project": {"name": "workflow", "outdir": "output"},
            "variable": {"name": "index", "units": "K"},
            "data": {
                "reference": {"kind": "series", "series_file": "data/observations.csv"},
                "models": {"kind": "series", "series_dir": "data", "scenarios": ["ssp245"],
                           "models": ["A", "B", "C"]},
            },
            "periods": {"historical": [1981, 2000], "reference_climate": [1981, 2000],
                        "future": [2031, 2050], "future_windows": [[2031, 2050]]},
            "index": {"type": "series"},
            "correction": {"method": "mean"},
            "ensemble": {"min_models": 2},
            "extremes": {"enabled": False},
            "plots": {"figures": [], "dpi": 80},
            "runtime": {"verbose": False},
        }

    def tearDown(self):
        self.tmp.cleanup()

    def write_series(self, name, years, values):
        with (self.root / "data" / name).open("w") as f:
            f.write("year,value\n")
            for year, value in zip(years, values):
                f.write(f"{year},{value}\n")

    def save_config(self, config=None):
        path = self.root / "config.yaml"
        path.write_text(yaml.safe_dump(self.config if config is None else config))
        return path

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(CLI), *map(str, args)],
                              cwd=self.root / "elsewhere", capture_output=True, text=True)

    def test_relative_paths_and_environment(self):
        os.environ["CMIP_WORKFLOW_DATA"] = str(self.root / "data")
        self.config["data"]["models"]["series_dir"] = "$CMIP_WORKFLOW_DATA"
        cfg = load_config(str(self.save_config()))
        self.assertEqual(cfg.outdir, str(self.root / "output"))
        self.assertEqual(cfg.get_path("data.reference.series_file"), str(self.root / "data/observations.csv"))
        os.environ.pop("CMIP_WORKFLOW_DATA", None)

    def test_undefined_environment_is_error(self):
        self.config["data"]["models"]["series_dir"] = "$CMIP_VARIABLE_THAT_IS_NOT_DEFINED"
        with self.assertRaises(ConfigError):
            load_config(str(self.save_config()))

    def test_unknown_keys_fail_by_default(self):
        self.config["correction"]["methood"] = "mean"
        with self.assertRaisesRegex(ConfigError, "methood"):
            load_config(str(self.save_config()))

    def test_unknown_keys_migration_mode(self):
        self.config["runtime"]["strict_config"] = False
        self.config["correction"]["methood"] = "mean"
        self.assertEqual(load_config(str(self.save_config())).get_path("correction.method"), "mean")

    def test_wrong_section_type_is_actionable(self):
        self.config["data"]["reference"] = None
        with self.assertRaisesRegex(ConfigError, "mapping"):
            load_config(str(self.save_config()))

    def test_quoted_boolean_is_rejected(self):
        self.config["extremes"]["enabled"] = "false"
        with self.assertRaisesRegex(ConfigError, "boolean"):
            load_config(str(self.save_config()))

    def test_month_order_and_ranges(self):
        self.assertEqual(resolve_months("11-2"), [11, 12, 1, 2])
        self.assertEqual(resolve_months([12, 1, 2]), [12, 1, 2])
        with self.assertRaises(ConfigError):
            resolve_months([0, 13])

    def test_incompatible_data_kind_is_rejected(self):
        self.config["data"]["models"].update(kind="grid", root="data")
        with self.assertRaisesRegex(ConfigError, "same data kind"):
            load_config(str(self.save_config()))

    def test_stage_dependencies_are_automatic(self):
        self.config["plots"]["figures"] = ["timeseries"]
        cfg = load_config(str(self.save_config()))
        self.assertEqual(plan_stages(cfg, ["figures"]), ["prepare", "correct", "ensemble", "figures"])
        with self.assertRaises(ConfigError):
            plan_stages(cfg, ["ensembl"])

    def test_preflight_does_not_create_output(self):
        cfg = load_config(str(self.save_config()))
        result = inspect_config(cfg, check_inputs=True)
        self.assertTrue(result["ok"], result)
        self.assertEqual(len(result["datasets"]), 7)
        self.assertFalse((self.root / "output").exists())

    def test_missing_inputs_fail_preflight(self):
        (self.root / "data/C_ssp245.csv").unlink()
        result = inspect_config(load_config(str(self.save_config())), check_inputs=True)
        self.assertFalse(result["ok"])
        self.assertFalse((self.root / "output").exists())

    def test_dry_run_works_from_other_directory(self):
        result = self.run_cli("run", self.save_config(), "--dry-run", "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertTrue(payload["ok"])
        self.assertFalse((self.root / "output").exists())

    def test_partial_stage_does_not_claim_report(self):
        result = self.run_cli("run", self.save_config(), "--stages", "correct")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("no new report", result.stdout)
        self.assertFalse((self.root / "output/report.md").exists())
        manifest = json.loads((self.root / "output/manifest.json").read_text())
        self.assertEqual(manifest["completed_stages"], ["prepare", "correct"])

    def test_model_counts_and_individual_year_axes(self):
        pipe = Pipeline(load_config(str(self.save_config())))
        pipe.run(["report"])
        item = pipe.ens["ssp245"]
        self.assertEqual(item["n_models"], 3)
        j = np.where(item["years"] == 2035)[0][0]
        self.assertEqual(item["n_models_by_year"][j], 2)
        self.assertTrue(np.isfinite(item["center"][j]))
        with np.load(self.root / "output/series/index_series.npz", allow_pickle=False) as saved:
            self.assertEqual(saved["B__ssp245"].size, saved["B__ssp245__years"].size)
            self.assertNotIn(2035, saved["B__ssp245__years"])
            self.assertEqual(saved["model_names"].tolist(), ["A", "B", "C"])
        self.assertTrue((self.root / "output/tables/ensemble_coverage.csv").is_file())

    def test_minimum_masks_missing_year(self):
        self.config["ensemble"]["min_models"] = 3
        pipe = Pipeline(load_config(str(self.save_config())))
        pipe.run(["ensemble"])
        item = pipe.ens["ssp245"]
        self.assertTrue(np.isnan(item["center"][item["years"] == 2035]).all())

    def test_no_eligible_ensemble_is_error(self):
        self.config["ensemble"]["min_models"] = 4
        pipe = Pipeline(load_config(str(self.save_config())))
        with self.assertRaisesRegex(ConfigError, "no year satisfies"):
            pipe.run(["ensemble"])
        manifest = json.loads((self.root / "output/manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")

    def test_changed_scientific_config_cannot_overwrite(self):
        cfg = load_config(str(self.save_config()))
        Pipeline(cfg).run(["prepare"])
        self.config["correction"]["method"] = "none"
        with self.assertRaisesRegex(ConfigError, "different scientific configuration"):
            Pipeline(load_config(str(self.save_config())))

    def test_init_templates_are_portable_and_do_not_overwrite(self):
        for name in ("grid", "series", "precip", "regional-extremes", "daily"):
            destination = self.root / name
            first = self.run_cli("init", destination, "--template", name)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(self.run_cli("validate", destination / "config.yaml").returncode, 0)
            self.assertEqual(self.run_cli("init", destination, "--template", name).returncode, 2)

    def test_fixed_historical_standardization_reaches_future(self):
        lat, lon = [30., 32.], [100., 102.]
        def field(path, years, values):
            da = xr.DataArray(np.broadcast_to(np.array(values)[:, None, None], (len(years), 2, 2)).copy(),
                              dims=("year", "lat", "lon"), coords={"year": years, "lat": lat, "lon": lon},
                              name="tas", attrs={"units": "K"})
            da.to_netcdf(self.root / "data" / path)
        field("obs.nc", self.hist, self.observed)
        field("A_historical.nc", self.hist, self.observed + 2)
        field("A_ssp245.nc", self.future, self.observed + 8)
        config = copy.deepcopy(self.config)
        config["data"] = {
            "reference": {"kind": "grid", "paths": ["data/obs.nc"], "variable": "tas",
                          "aggregate": {"frequency": "annual"}},
            "models": {"kind": "grid", "root": "data", "pattern": "{root}/{model}_{experiment}.nc",
                       "variable": "tas", "models": ["A"], "scenarios": ["ssp245"],
                       "aggregate": {"frequency": "annual"}},
        }
        config["index"] = {"type": "region_mean", "region": {"lat": lat, "lon": lon},
                           "anomaly": True, "standardize": True, "units_label": "Standardized index"}
        config["ensemble"]["min_models"] = 1
        pipe = Pipeline(load_config(str(self.save_config(config))))
        pipe.run(["ensemble"])
        self.assertGreater(float(np.nanmean(pipe.ens["ssp245"]["center"])), 4.0)



    def test_overwrite_preserves_old_run_in_archive(self):
        Pipeline(load_config(str(self.save_config()))).run(["ensemble"])
        self.config["correction"]["method"] = "none"
        self.config["runtime"]["overwrite"] = True
        pipe = Pipeline(load_config(str(self.save_config())))
        pipe.run(["prepare"])
        self.assertFalse((self.root / "output/ensemble/ensemble_ssp245.npz").exists())
        archived = Path(pipe.previous_run_directory)
        self.assertTrue((archived / "ensemble/ensemble_ssp245.npz").is_file())
        manifest = json.loads((self.root / "output/manifest.json").read_text())
        self.assertEqual(manifest["previous_run_directory"], str(archived))

    def test_one_missing_extreme_scenario_fails_strict_run(self):
        self.config["data"]["models"]["scenarios"].append("ssp585")
        for model in ("A", "B", "C"):
            self.write_series(f"{model}_ssp585.csv", self.future, np.ones(self.future.size) * 30)
        self.config["extremes"].update(enabled=True, min_years=3, distribution="gumbel")
        pipe = Pipeline(load_config(str(self.save_config())))
        with self.assertRaisesRegex(ConfigError, "no valid extreme-value fits"):
            pipe.run(["extremes"])
        self.assertEqual(json.loads((self.root / "output/manifest.json").read_text())["status"], "failed")

    def test_all_insufficient_extreme_samples_fail(self):
        self.config["extremes"].update(enabled=True, min_years=30)
        pipe = Pipeline(load_config(str(self.save_config())))
        with self.assertRaises(ConfigError):
            pipe.run(["extremes"])

    def test_actual_inputs_are_separate_from_preflight_candidates(self):
        pipe = Pipeline(load_config(str(self.save_config())))
        pipe.input_inventory = [{"label": "intentionally only a preflight marker", "files": []}]
        pipe.run(["correct"])
        manifest = json.loads((self.root / "output/manifest.json").read_text())
        self.assertEqual(len(manifest["actual_inputs"]), 7)
        self.assertEqual(manifest["preflight_inputs"], pipe.input_inventory)
        self.assertTrue(all(item["files_metadata"] for item in manifest["actual_inputs"]))

    def test_daily_unsupported_anomaly_and_member_average_are_rejected(self):
        for field, value in (("index", "region_mean_anomaly"), ("member_combine", "mean")):
            config = copy.deepcopy(self.config)
            config["daily"] = {"enabled": True}
            if field == "index":
                config["daily"][field] = value
            else:
                config["data"]["models"][field] = value
            with self.assertRaises(ConfigError):
                load_config(str(self.save_config(config)))

    def test_plot_does_not_bridge_history_and_future_and_reports_counts(self):
        import matplotlib.pyplot as plt
        from cmipkit import plotting
        cfg = load_config(str(self.save_config()))
        hist = plotting.HistoricalSeries(self.hist, self.observed)
        matrix = np.ones((3, 20))
        matrix[1, 5] = np.nan
        scenario = plotting.ScenarioSeries("ssp245", "Scenario", self.future, matrix,
                                           np.ones(20), {}, n_models=3)
        fig, ax = plt.subplots()
        plotting._plot_series_panel(ax, cfg, hist, scenario, show_legend=False)
        self.assertIn("2–3", ax.get_title(loc="left"))
        for line in ax.lines:
            x = np.asarray(line.get_xdata(), dtype=float)
            self.assertFalse(np.any(x <= 2000) and np.any(x >= 2031))
        plt.close(fig)

    def test_djf_daily_baseline_and_member_identity(self):
        def dataset(filename, start, stop, level):
            dates = xr.date_range(start, stop, freq="D", use_cftime=False)
            year = np.asarray(dates.year)
            month = np.asarray(dates.month)
            values = level + 0.1 * (year - year.min()) + np.sin(np.arange(len(dates)) * .27)
            values = values[:, None, None] + np.zeros((len(dates), 2, 2))
            xr.DataArray(values, dims=("time", "lat", "lon"),
                         coords={"time": dates, "lat": [30., 32.], "lon": [100., 102.]},
                         name="tas", attrs={"units": "K"}).to_netcdf(self.root / "data" / filename)
        dataset("obs.nc", "1989-12-01", "1994-02-28", 10)
        dataset("A_r2_historical.nc", "1989-12-01", "1994-02-28", 12)
        dataset("A_r2_ssp245.nc", "2029-12-01", "2034-02-28", 14)
        dataset("A_r1_ssp245.nc", "2029-12-01", "2034-02-28", 114)
        bbox = {"lat": [30., 32.], "lon": [100., 102.]}
        config = copy.deepcopy(self.config)
        config["data"] = {
            "reference": {"kind": "grid", "paths": ["data/obs.nc"], "variable": "tas", "bbox": bbox,
                          "aggregate": {"months": "DJF", "frequency": "daily"}},
            "models": {"kind": "grid", "root": "data",
                       "pattern": "{root}/{model}_{member}_{experiment}.nc",
                       "variable": "tas", "models": ["A"], "members": ["*"], "bbox": bbox,
                       "scenarios": ["ssp245"], "aggregate": {"months": "DJF", "frequency": "daily"}},
        }
        config["periods"] = {"historical": [1990, 1994], "reference_climate": [1990, 1994],
                             "future": [2030, 2034], "future_windows": [[2030, 2034]]}
        config["index"] = {"type": "region_mean", "region": bbox}
        config["correction"] = {"method": "none"}
        config["daily"] = {"enabled": True, "months": "DJF", "correction": "none",
                           "plots": False, "extremes": {"kind": "block", "distribution": "gumbel",
                                                       "min_years": 3}}
        config["ensemble"]["min_models"] = 1
        pipe = Pipeline(load_config(str(self.save_config(config))))
        pipe.run(["daily"])
        self.assertEqual(pipe.daily["historical"].n, 5)
        fit = pipe.daily["scenarios"]["ssp245"]["fits"]["A"]
        self.assertEqual(fit.n, 5)
        self.assertLess(fit.loc, 30)
        self.assertEqual(pipe.source.model_members["A"], ["r2"])
        manifest = json.loads((self.root / "output/manifest.json").read_text())
        actual = [entry["path"] for item in manifest["actual_inputs"] for entry in item["files_metadata"]]
        self.assertIn(str(self.root / "data/A_r2_ssp245.nc"), actual)
        self.assertNotIn(str(self.root / "data/A_r1_ssp245.nc"), actual)

if __name__ == "__main__":
    unittest.main(verbosity=2)
