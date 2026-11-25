"""Synthetic failures and bounded-read checks for raw input quality reports."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import cftime
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from cmipkit import quality
from cmipkit.config import Config, DEFAULTS
from cmipkit.preflight import inspect_config


class QualityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cmip_quality_")
        self.root = Path(self.temp.name)
        self.times = np.array(["2000-01-15", "2000-02-15", "2000-03-15"], dtype="datetime64[ns]")
        self.spec = {"kind": "grid", "variable": "tas", "units_in": "K",
                     "aggregate": {"frequency": "monthly"}, "units_transform": {"scale": 1, "offset": 0, "units": "K"}}
        self.options = {"coverage": "error", "require_units": False, "chunk_size": 2, "max_chunk_values": 6}

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, times=None, values=None, units="K", metadata=None, lat=(10., 20.), lon=(100., 110., 120.)):
        times = self.times if times is None else times
        if values is None:
            values = np.arange(len(times), dtype=float) + 280
        values = np.asarray(values, dtype=float)
        if values.ndim == 1:
            values = np.broadcast_to(values[:, None, None], (len(times), len(lat), len(lon))).copy()
        data = xr.DataArray(values, dims=("time", "lat", "lon"),
                            coords={"time": times, "lat": list(lat), "lon": list(lon)}, name="tas")
        if units is not None:
            data.attrs["units"] = units
        ds = data.to_dataset()
        ds.attrs.update(metadata or {})
        path = self.root / name
        ds.to_netcdf(path)
        return str(path)

    def audit(self, files, *, deep=False, spec=None, options=None, years=(2000, 2000), months=(1, 2, 3), expected=None):
        return quality.audit_grid(files, spec or self.spec, label="sample", years=years, months=months,
                                  quality=options or self.options, check_data=deep, expected=expected)

    def codes(self, report):
        return {f["code"] for f in report["findings"]}

    def test_all_headers_including_second_corrupt_file(self):
        good = self.write("a.nc")
        bad = self.root / "b.nc"
        bad.write_bytes(b"not a netcdf")
        result = self.audit([good, str(bad)], deep=True)
        self.assertFalse(result["ok"])
        self.assertFalse(result["values_checked"])
        self.assertEqual(result["files_values_scanned"], 1)
        self.assertIn("unreadable_input", self.codes(result))

    def test_empty_file_is_explicit_error(self):
        empty = self.root / "empty.nc"; empty.touch()
        result = self.audit([str(empty)])
        self.assertFalse(result["ok"])
        self.assertIn("empty input file", str(result["findings"]))

    def test_headers_do_not_access_field_values(self):
        path = self.write("data.nc")
        original = xr.DataArray.values.fget
        def guarded(data):
            if data.name == "tas":
                raise AssertionError("header scan attempted a data read")
            return original(data)
        with patch.object(xr.DataArray, "values", property(guarded)):
            result = self.audit([path])
        self.assertTrue(result["ok"], result)
        self.assertFalse(result["values_checked"])

    def test_full_scan_uses_spatial_and_time_tiles(self):
        path = self.write("data.nc")
        original = xr.DataArray.values.fget
        reads = []
        def guarded(data):
            if data.name == "tas":
                reads.append(data.size)
                self.assertLessEqual(data.size, 6)
            return original(data)
        with patch.object(xr.DataArray, "values", property(guarded)):
            result = self.audit([path], deep=True)
        self.assertTrue(result["ok"], result)
        self.assertTrue(result["values_checked"])
        self.assertEqual(result["values"]["count"], 18)
        self.assertTrue(reads)
        self.assertLessEqual(result["values"]["max_chunk_values_observed"], 6)

    def test_second_file_units_and_identity_conflicts(self):
        first = self.write("a.nc", self.times[:1], metadata={"source_id": "A", "variant_label": "r1i1p1f1"})
        second = self.write("b.nc", self.times[1:], units="degC", metadata={"source_id": "A", "variant_label": "r2i1p1f1"})
        result = self.audit([first, second])
        self.assertFalse(result["ok"])
        self.assertTrue({"mixed_units", "mixed_identity", "units_mismatch"}.issubset(self.codes(result)))

    def test_unit_aliases_are_not_unit_conversions(self):
        self.assertEqual(quality.canonical_units("degrees_Celsius"), quality.canonical_units("°C"))
        self.assertEqual(quality.canonical_units("kg m**-2 s**-1"), quality.canonical_units("kg/m2/s"))
        self.assertNotEqual(quality.canonical_units("K"), quality.canonical_units("degC"))
        self.assertNotEqual(quality.canonical_units("kg m-2 s-1"), quality.canonical_units("mm/day"))
        self.assertIsNone(quality.canonical_units(" "))
        path = self.write("a.nc", units="kelvin")
        self.assertTrue(self.audit([path])["ok"])

    def test_missing_units_optional_or_required(self):
        path = self.write("a.nc", units=None)
        self.assertTrue(self.audit([path])["ok"])
        self.assertFalse(self.audit([path], options=dict(self.options, require_units=True))["ok"])

    def test_unit_relabel_without_conversion_is_visible(self):
        path = self.write("a.nc")
        spec = copy.deepcopy(self.spec); spec["units_transform"]["units"] = "degC"
        self.assertIn("units_relabelled_without_conversion", self.codes(self.audit([path], spec=spec)))

    def test_missing_month_and_unsorted_dates(self):
        path = self.write("a.nc", self.times[[2, 0]])
        result = self.audit([path])
        self.assertFalse(result["ok"])
        self.assertEqual(result["coverage"]["missing_samples"], 1)
        self.assertIn("unsorted_time", self.codes(result))
        self.assertTrue(self.audit([path], options=dict(self.options, coverage="warn"))["ok"])

    def test_same_month_different_timestamp_is_error(self):
        times = np.array(["2000-01-01", "2000-01-15"], dtype="datetime64[ns]")
        path = self.write("a.nc", times)
        result = self.audit([path], months=(1,))
        self.assertIn("multiple_samples_per_period", self.codes(result))
        self.assertFalse(result["ok"])

    def test_identical_and_conflicting_duplicate_time(self):
        first = self.write("a.nc", self.times[:2], [280, 281])
        second = self.write("b.nc", self.times[1:], [281, 282])
        result = self.audit([first, second], deep=True)
        self.assertTrue(result["ok"], result)
        self.assertIn("duplicate_time", self.codes(result))
        second = self.write("b.nc", self.times[1:], [299, 282])
        result = self.audit([first, second], deep=True)
        self.assertFalse(result["ok"])
        self.assertIn("conflicting_duplicate_time", self.codes(result))

    def test_zero_scale_does_not_hide_raw_duplicate_conflict(self):
        first = self.write("a.nc", self.times, [1, 2, 3])
        second = self.write("b.nc", self.times, [4, 5, 6])
        spec = copy.deepcopy(self.spec); spec["units_transform"]["scale"] = 0
        result = self.audit([first, second], deep=True, spec=spec)
        self.assertIn("conflicting_duplicate_time", self.codes(result))

    def test_cross_year_monthly_coverage(self):
        dates = np.array(["1999-12-15", "2000-01-15", "2000-02-15"], dtype="datetime64[ns]")
        path = self.write("winter.nc", dates)
        result = self.audit([path], months=(12, 1, 2))
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["coverage"]["expected_samples"], 3)

    def test_daily_360_day_and_noleap_coverage(self):
        spec = copy.deepcopy(self.spec); spec["aggregate"]["frequency"] = "daily"
        for calendar, count in (("360_day", 30), ("noleap", 28), ("all_leap", 29)):
            times = xr.date_range("2000-02-01", periods=count, freq="D", calendar=calendar, use_cftime=True)
            path = self.write(calendar + ".nc", times)
            result = self.audit([path], spec=spec, months=(2,))
            self.assertTrue(result["ok"], result)
            self.assertEqual(result["coverage"]["expected_samples"], count)
        path = self.write("missing.nc", times[1:])
        self.assertEqual(self.audit([path], spec=spec, months=(2,))["coverage"]["missing_samples"], 1)

    def test_daily_djf_noleap_counts_prior_december(self):
        times = xr.date_range("1999-12-01", periods=90, freq="D", calendar="noleap", use_cftime=True)
        spec = copy.deepcopy(self.spec); spec["aggregate"]["frequency"] = "daily"
        result = self.audit([self.write("djf.nc", times)], spec=spec, months=(12, 1, 2))
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["coverage"]["expected_samples"], 90)

    def test_mixed_calendars_in_one_dataset_fail(self):
        first = self.write("gregorian.nc", self.times[:1])
        second = self.write("noleap.nc", [cftime.DatetimeNoLeap(2000, 2, 15), cftime.DatetimeNoLeap(2000, 3, 15)])
        self.assertIn("mixed_calendars", self.codes(self.audit([first, second])))

    def test_annual_input_ignores_month_filter(self):
        times = np.array(["2000-06-15", "2001-06-15"], dtype="datetime64[ns]")
        spec = copy.deepcopy(self.spec); spec["aggregate"]["frequency"] = "auto"
        result = self.audit([self.write("annual.nc", times)], spec=spec, years=(2000, 2001), months=(12, 1, 2))
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["coverage"]["unique_samples"], 2)

    def test_static_nan_mask_allowed_until_threshold_set(self):
        values = np.full((3, 2, 3), 280.); values[:, 0, 0] = np.nan
        path = self.write("masked.nc", values=values)
        result = self.audit([path], deep=True)
        self.assertTrue(result["ok"], result)
        self.assertAlmostEqual(result["values"]["missing_fraction"], 1 / 6)
        result = self.audit([path], deep=True, options=dict(self.options, max_missing_fraction=0.1))
        self.assertIn("missing_fraction_exceeded", self.codes(result))

    def test_all_missing_and_infinite_selected_data(self):
        nanfile = self.write("nan.nc", values=[np.nan] * 3)
        self.assertIn("all_values_missing", self.codes(self.audit([nanfile], deep=True)))
        inf = self.write("inf.nc", values=[280, np.inf, 282])
        spec = copy.deepcopy(self.spec); spec["units_transform"]["scale"] = 0
        self.assertIn("infinite_values", self.codes(self.audit([inf], deep=True, spec=spec)))

    def test_limits_apply_after_explicit_unit_conversion(self):
        path = self.write("a.nc", values=[273.15, 283.15, 293.15])
        spec = copy.deepcopy(self.spec)
        spec["units_transform"].update(offset=-273.15, units="degC")
        result = self.audit([path], deep=True, spec=spec, options=dict(self.options, valid_min=-1, valid_max=21))
        self.assertTrue(result["ok"], result)
        self.assertAlmostEqual(result["values"]["min"], 0)
        self.assertAlmostEqual(result["values"]["max"], 20)
        result = self.audit([path], deep=True, spec=spec, options=dict(self.options, valid_max=15))
        self.assertIn("value_range_exceeded", self.codes(result))

    def test_scan_respects_time_and_bbox(self):
        values = np.full((3, 2, 3), np.inf); values[1, 0, 1] = 280
        path = self.write("a.nc", values=values)
        spec = dict(self.spec, bbox={"lat": [9, 11], "lon": [109, 111]})
        result = self.audit([path], deep=True, spec=spec, months=(2,))
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["values"]["count"], 1)

    def test_metadata_identifies_wrong_file(self):
        path = self.write("a.nc", metadata={"source_id": "Other", "experiment_id": "historical"})
        result = self.audit([path], expected={"source_id": "A", "experiment_id": "ssp245"})
        self.assertFalse(result["ok"])
        self.assertIn("identity_mismatch", self.codes(result))

    def test_report_is_json_serializable(self):
        result = self.audit([self.write("a.nc")], deep=True)
        json.dumps(result, allow_nan=False)

    def config(self):
        cfg = Config(copy.deepcopy(DEFAULTS))
        cfg["project"].update(name="qc", outdir=str(self.root / "output"))
        cfg["variable"].update(name="tas", units="K")
        cfg["index"]["region"] = {"lat": [9, 21], "lon": [99, 121]}
        cfg["periods"].update(historical=[2000, 2000], future=[2001, 2001], reference_climate=[2000, 2000], future_windows=[[2001, 2001]])
        cfg["data"]["reference"].update(paths=[self.write("reference.nc")], variable="tas")
        cfg["data"]["reference"]["aggregate"].update(months=[1, 2, 3], frequency="monthly")
        cfg["data"]["models"].update(root=str(self.root), pattern="{root}/{model}_{experiment}_{member}.nc",
                                     variable="tas", models=["A"], members=["*"], member_combine="mean", scenarios=["ssp245"])
        cfg["data"]["models"]["aggregate"].update(months=[1, 2, 3], frequency="monthly")
        for experiment in ("historical", "ssp245"):
            dates = self.times if experiment == "historical" else self.times + np.timedelta64(366, "D")
            for member in ("r1i1p1f1", "r2i1p1f1"):
                self.write(f"A_{experiment}_{member}.nc", dates,
                           metadata={"source_id": "A", "experiment_id": experiment, "variant_label": member, "grid_label": "gn"})
        return cfg

    def test_wildcard_ensemble_inspects_each_member(self):
        cfg = self.config()
        result = inspect_config(cfg, check_inputs=True)
        self.assertTrue(result["ok"], result)
        self.assertEqual(len(result["datasets"]), 5)
        self.assertFalse((self.root / "output").exists())
        (self.root / "A_ssp245_r2i1p1f1.nc").write_bytes(b"bad second member")
        result = inspect_config(cfg, check_data=True)
        self.assertFalse(result["ok"])
        self.assertFalse(result["quality"]["data_values_checked"])

    def test_cross_experiment_grid_identity_mismatch(self):
        cfg = self.config()
        dates = self.times + np.timedelta64(366, "D")
        self.write("A_ssp245_r1i1p1f1.nc", dates,
                   metadata={"source_id": "A", "experiment_id": "ssp245", "variant_label": "r1i1p1f1", "grid_label": "gr"})
        result = inspect_config(cfg, check_inputs=True)
        self.assertFalse(result["ok"])
        self.assertIn("historical_scenario_mismatch", str(result["errors"]))

    def daily_config(self):
        cfg = Config(copy.deepcopy(DEFAULTS))
        cfg["project"].update(name="daily_qc", outdir=str(self.root / "output"))
        cfg["variable"].update(name="tas", units="K")
        region = {"lat": [9, 21], "lon": [99, 121]}
        cfg["index"]["region"] = region
        cfg["periods"].update(historical=[2000, 2000], future=[2001, 2001], reference_climate=[2000, 2000],
                              future_windows=[[2001, 2001]])
        cfg["daily"].update(enabled=True, months=[12, 1, 2])
        cfg["quality"]["coverage"] = "error"
        cfg["ensemble"]["min_models"] = 1
        cfg["data"]["reference"].update(paths=[str(self.root / "reference.nc")], variable="tas", bbox=region)
        cfg["data"]["models"].update(root=str(self.root), pattern="{root}/{model}_{experiment}_{member}.nc",
                                     variable="tas", bbox=region, models=["A"], members=["r1i1p1f1"],
                                     member_combine="first", scenarios=["ssp245"])
        for group in ("reference", "models"):
            cfg["data"][group]["aggregate"].update(months=[6, 7, 8], frequency="daily")
        for label, start, stop in (("reference.nc", "1999-12-01", "2001-09-01"),
                                   ("A_historical_r1i1p1f1.nc", "1999-12-01", "2000-09-01"),
                                   ("A_ssp245_r1i1p1f1.nc", "2000-12-01", "2001-09-01")):
            dates = np.arange(np.datetime64(start), np.datetime64(stop))
            self.write(label, dates, values=np.full(len(dates), 280.))
        return cfg

    def test_daily_scope_checks_winter_values_outside_annual_jja(self):
        cfg = self.daily_config()
        dates = np.arange(np.datetime64("1999-12-01"), np.datetime64("2001-09-01"))
        values = np.full(len(dates), 280.)
        values[dates == np.datetime64("1999-12-15")] = np.inf
        self.write("reference.nc", dates, values=values)
        result = inspect_config(cfg, check_data=True)
        reports = {a["label"]: a for a in result["quality"]["datasets"]}
        self.assertTrue(reports["reference"]["ok"], reports["reference"])
        self.assertIn("infinite_values", self.codes(reports["daily reference"]))
        self.assertFalse(result["ok"])
        self.assertEqual(reports["daily reference"]["coverage"]["months"], [12, 1, 2])

    def test_daily_scope_reports_missing_prior_december(self):
        cfg = self.daily_config()
        dates = np.arange(np.datetime64("2000-01-01"), np.datetime64("2000-09-01"))
        self.write("reference.nc", dates, values=np.full(len(dates), 280.))
        result = inspect_config(cfg, check_inputs=True)
        reports = {a["label"]: a for a in result["quality"]["datasets"]}
        self.assertTrue(reports["reference"]["ok"])
        self.assertFalse(reports["daily reference"]["ok"])
        self.assertEqual(reports["daily reference"]["coverage"]["missing_samples"], 31)
        self.assertIn("1999, 12", str(reports["daily reference"]["coverage"]["missing_examples"]))

    def test_daily_reference_uses_historical_period_not_reference_years(self):
        cfg = self.daily_config()
        cfg["data"]["reference"]["years"] = [2001, 2001]
        dates = np.arange(np.datetime64("1999-12-01"), np.datetime64("2001-09-01"))
        values = np.full(len(dates), 280.)
        values[dates == np.datetime64("1999-12-15")] = np.inf
        self.write("reference.nc", dates, values=values)
        result = inspect_config(cfg, check_data=True)
        reports = {a["label"]: a for a in result["quality"]["datasets"]}
        self.assertEqual(reports["reference"]["coverage"]["requested_years"], [2001, 2001])
        self.assertTrue(reports["reference"]["ok"])
        self.assertEqual(reports["daily reference"]["coverage"]["requested_years"], [2000, 2000])
        self.assertIn("infinite_values", self.codes(reports["daily reference"]))

    def test_identical_daily_and_annual_scopes_reuse_value_scan(self):
        cfg = self.daily_config()
        cfg["daily"]["months"] = [6, 7, 8]
        with patch.object(quality, "_scan_values", wraps=quality._scan_values) as scanner:
            result = inspect_config(cfg, check_data=True)
        self.assertTrue(result["ok"], result)
        self.assertEqual(scanner.call_count, 3)
        self.assertEqual(len(result["datasets"]), 3)
        self.assertTrue(all(a["uses"] == ["annual", "daily"] for a in result["quality"]["datasets"]))

    def test_daily_rejects_monthly_samples_even_when_gap_policy_warns(self):
        cfg = self.daily_config()
        cfg["daily"]["months"] = [6, 7, 8]
        cfg["quality"]["coverage"] = "warn"
        for group in ("reference", "models"):
            cfg["data"][group]["aggregate"]["frequency"] = "monthly"
        for name, year in (("reference.nc", 2000), ("A_historical_r1i1p1f1.nc", 2000), ("A_ssp245_r1i1p1f1.nc", 2001)):
            dates = np.array([f"{year}-06-15", f"{year}-07-15", f"{year}-08-15"], dtype="datetime64[ns]")
            self.write(name, dates)
        result = inspect_config(cfg, check_inputs=True)
        self.assertFalse(result["ok"], result)
        self.assertIn("daily_sampling_required", str(result["errors"]))

    def test_too_few_models_is_nonblocking_and_does_not_change_threshold(self):
        cfg = self.config()
        cfg["ensemble"]["min_models"] = 3
        result = inspect_config(cfg, check_inputs=True)
        self.assertTrue(result["ok"], result)
        self.assertEqual(cfg["ensemble"]["min_models"], 3)
        self.assertIn("insufficient_models", str(result["warnings"]))

    def test_sparse_amon_auto_uses_monthly_season_labels(self):
        dates = np.array(["1999-12-15", "2000-12-15", "2001-12-15"], dtype="datetime64[ns]")
        path = self.write("sparse_monthly.nc", dates, values=[280., 281., 10000.], metadata={"table_id": "Amon"})
        spec = copy.deepcopy(self.spec)
        spec["aggregate"]["frequency"] = "auto"
        result = self.audit([path], deep=True, spec=spec, years=(2000, 2001), months=(12, 1, 2),
                            options=dict(self.options, coverage="warn", valid_max=300))
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["input_frequency"], "monthly")
        self.assertEqual(result["coverage"]["frequency"], "monthly")
        self.assertEqual(result["values"]["max"], 281.)

    def test_shared_frequency_helpers_keep_explicit_choice_and_reject_conflicts(self):
        dates = np.array(["2000-12-15", "2001-12-15"], dtype="datetime64[ns]")
        self.assertEqual(quality.resolve_frequency(dates, "auto", declared_frequencies=["monthly"]), "monthly")
        self.assertEqual(quality.resolve_frequency(dates, "annual", declared_frequencies=["monthly"]), "annual")
        self.assertEqual(quality.resolve_frequency(dates), "annual")
        self.assertEqual(quality.declared_frequency({"_cmip_table_id": "Amon"}), "monthly")
        with self.assertRaises(ValueError):
            quality.declared_frequency({"frequency": "day", "table_id": "Amon"})
        with self.assertRaises(ValueError):
            quality.resolve_frequency(dates, declared_frequencies=["monthly", "daily"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
