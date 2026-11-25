"""Production-entry regressions for raw-input quality reports and safe failure."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import xarray as xr
import yaml

SKILL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL / "scripts"))
from cmipkit.config import ConfigError, load_config, quality_requires_values, write_json_atomic
from cmipkit.pipeline import Pipeline


class QualityWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cmip_quality_workflow_")
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        self.member = "r1i1p1f1"
        for name, start, offset, experiment in (
            ("observed.nc", 2000, 0, None),
            (f"A_historical_{self.member}.nc", 2000, 1, "historical"),
            (f"A_ssp245_{self.member}.nc", 2030, 3, "ssp245"),
        ):
            dates = np.arange(f"{start}-01", f"{start+3}-01", dtype="datetime64[M]")
            values = 280 + offset + np.arange(36)[:, None, None] * 0.05 + np.zeros((36, 2, 3))
            ds = xr.Dataset({"tas": (("time", "lat", "lon"), values)},
                            coords={"time": dates.astype("datetime64[ns]"), "lat": [10, 20], "lon": [100, 110, 120]})
            ds.tas.attrs["units"] = "K"
            if experiment:
                ds.attrs.update(source_id="A", experiment_id=experiment,
                                variant_label=self.member, grid_label="gn", table_id="Amon")
            ds.to_netcdf(self.root / name)
        self.config = {
            "project": {"name": "quality_workflow", "outdir": "output"},
            "variable": {"name": "tas", "units": "K"},
            "data": {
                "reference": {"paths": ["observed.nc"], "variable": "tas", "units_in": "K"},
                "models": {"root": ".", "patterns": ["{root}/{model}_{experiment}_{member}.nc"],
                           "models": ["A"], "members": [self.member], "scenarios": ["ssp245"],
                           "variable": "tas", "units_in": "K"},
            },
            "periods": {"historical": [2000, 2002], "future": [2030, 2032],
                        "reference_climate": [2000, 2002]},
            "index": {"type": "region_mean", "region": {"lat": [10, 20], "lon": [100, 120]}},
            "correction": {"method": "none"},
            "plots": {"figures": []},
            "runtime": {"cache": False, "verbose": False},
        }

    def save(self):
        path = self.root / "config.yaml"
        path.write_text(yaml.safe_dump(self.config), encoding="utf-8")
        return path

    def cli(self, *arguments):
        return subprocess.run([sys.executable, str(SKILL / "scripts/run.py"), *map(str, arguments)],
                              cwd=self.root, capture_output=True, text=True)

    def test_deep_validation_writes_report_without_analysis_or_input_changes(self):
        hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in self.root.glob("*.nc")}
        result = self.cli("validate", self.save(), "--check-data", "--report-json", "qc.json")
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        report = json.loads((self.root / "qc.json").read_text())
        self.assertTrue(report["ok"], report)
        self.assertEqual(report["quality"]["level"], "full")
        self.assertTrue(report["quality"]["data_values_checked"])
        self.assertFalse((self.root / "output").exists())
        self.assertEqual(hashes, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in hashes})

    def test_header_validation_does_not_claim_value_checks(self):
        result = self.cli("validate", self.save(), "--check-inputs", "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads(result.stdout)
        self.assertFalse(report["quality"]["data_values_checked"])
        self.assertEqual(report["quality"]["level"], "headers")

    def test_explicit_limits_trigger_full_run_preflight_and_save_failure(self):
        self.config["quality"] = {"valid_max": 100}
        cfg_path = self.save()
        self.assertTrue(quality_requires_values(load_config(str(cfg_path))))
        result = self.cli("run", cfg_path, "--report-json", "failed-qc.json")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        report = json.loads((self.root / "failed-qc.json").read_text())
        self.assertFalse(report["ok"])
        self.assertEqual(report["quality"]["level"], "full")
        self.assertFalse((self.root / "output").exists())

    def test_programmatic_failure_preserves_quality_report_and_status(self):
        self.config["quality"] = {"valid_max": 100}
        pipe = Pipeline(load_config(str(self.save())))
        with self.assertRaisesRegex(ConfigError, "input quality checks failed"):
            pipe.run(["prepare"])
        status = json.loads((self.root / "output/manifest.json").read_text())
        self.assertEqual(status["status"], "failed")
        self.assertFalse(status["quality"]["ok"])
        self.assertEqual(status["completed_stages"], [])
        self.assertTrue((self.root / "output/quality_report.json").is_file())

    def test_successful_run_records_quality_scope(self):
        self.config["quality"] = {"scan": "full", "chunk_size": 2, "max_chunk_values": 2}
        pipe = Pipeline(load_config(str(self.save())))
        pipe.run(["prepare"])
        status = json.loads((self.root / "output/manifest.json").read_text())
        self.assertEqual(status["status"], "completed")
        self.assertEqual(status["quality"]["level"], "full")
        self.assertTrue(status["quality"]["data_values_checked"])
        self.assertIn("quality_report.json", status["artifacts_current_invocation"])

    def test_interrupt_records_interrupted_status(self):
        pipe = Pipeline(load_config(str(self.save())))
        with mock.patch.object(pipe, "prepare", side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                pipe.run(["prepare"])
        status = json.loads((self.root / "output/manifest.json").read_text())
        self.assertEqual(status["status"], "interrupted")
        self.assertTrue(status["quality"]["ok"])

    def test_report_cannot_overwrite_input_or_configuration(self):
        cfg_path = self.save()
        for protected in (cfg_path, self.root / "observed.nc"):
            before = protected.read_bytes()
            result = self.cli("validate", cfg_path, "--check-inputs", "--report-json", protected)
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertIn("must not overwrite", result.stderr)
            self.assertEqual(protected.read_bytes(), before)

    def test_atomic_json_failure_keeps_previous_complete_file(self):
        target = self.root / "result.json"
        write_json_atomic(str(target), {"previous": True})
        with self.assertRaises(TypeError):
            write_json_atomic(str(target), {"bad": object()})
        self.assertEqual(json.loads(target.read_text()), {"previous": True})
        self.assertEqual(list(self.root.glob(".cmipkit-*.tmp")), [])

    def test_configuration_failure_can_be_saved_as_json(self):
        self.config["quality"] = {"valid_mxa": 100}
        result = self.cli("validate", self.save(), "--report-json", "invalid-qc.json")
        self.assertEqual(result.returncode, 2)
        report = json.loads((self.root / "invalid-qc.json").read_text())
        self.assertFalse(report["ok"])
        self.assertEqual(report["quality"]["level"], "config")
        self.assertIn("valid_mxa", " ".join(report["errors"]))

    def test_config_only_report_cannot_overwrite_netcdf(self):
        target = self.root / "observed.nc"
        before = target.read_bytes()
        result = self.cli("validate", self.save(), "--report-json", target)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(target.read_bytes(), before)

    def test_report_does_not_overwrite_existing_text_or_json_data(self):
        cfg_path = self.save()
        for name, contents in (("observations.txt", "year,value\n2000,280\n2001,281\n"),
                               ("data.json", '{"observations": [280, 281]}')):
            target = self.root / name
            target.write_text(contents)
            result = self.cli("validate", cfg_path, "--report-json", target)
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertEqual(target.read_text(), contents)

    def test_existing_quality_report_can_be_updated(self):
        cfg_path = self.save()
        for check in ("--check-inputs", "--check-data"):
            result = self.cli("validate", cfg_path, check, "--report-json", "qc.json")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads((self.root / "qc.json").read_text())
        self.assertEqual(report["report_type"], "cmipkit.input_quality")
        self.assertEqual(report["quality"]["level"], "full")

    def test_old_resolved_configuration_accepts_new_quality_options(self):
        cfg = load_config(str(self.save()))
        previous = cfg.to_dict()
        previous.pop("quality")
        previous["runtime"].pop("max_input_values")
        (self.root / "output").mkdir()
        (self.root / "output/config.resolved.yaml").write_text(yaml.safe_dump(previous))
        self.config["quality"] = {"scan": "full"}
        Pipeline(load_config(str(self.save())))

    def test_quality_settings_reject_ambiguous_or_invalid_limits(self):
        for key, value in (("scan", "sample"), ("coverage", "ignore"),
                           ("require_units", "false"), ("chunk_size", 0),
                           ("max_chunk_values", 1.5), ("max_missing_fraction", 1.1),
                           ("valid_min", float("nan")), ("valid_mxa", 100)):
            with self.subTest(key=key):
                self.config["quality"] = {key: value}
                with self.assertRaises(ConfigError):
                    load_config(str(self.save()))
        self.config["quality"] = {"valid_min": 10, "valid_max": 0}
        with self.assertRaises(ConfigError):
            load_config(str(self.save()))


if __name__ == "__main__":
    unittest.main(verbosity=2)
