"""Regression cases for bounded CMIP reads, source identity and safe writes."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from cmipkit import io
from cmipkit.config import Config, ConfigError


def write_field(path, dates, *, lat=(0., 1.), lon=(0., 1.), values=None,
                units="K", metadata=None):
    shape = (len(dates), len(lat), len(lon))
    values = np.full(shape, 280.) if values is None else np.broadcast_to(values, shape).copy()
    data = xr.DataArray(values, dims=("time", "lat", "lon"),
                        coords={"time": dates, "lat": list(lat), "lon": list(lon)},
                        name="tas", attrs={"units": units})
    data.to_dataset().assign_attrs(metadata or {}).to_netcdf(path)
    return str(path)


def monthly(year):
    return np.array([f"{year}-{month:02d}-15" for month in range(1, 13)], dtype="datetime64[ns]")


class RawIOContracts(unittest.TestCase):
    def test_subset_loaded_while_file_open(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_field(Path(tmp) / "raw.nc", monthly(2000), lat=range(10), lon=range(10))
            opened, closed, loaded = [], [], []
            real_open, real_load = xr.open_dataset, xr.DataArray.load

            def traced_open(*args, **kwargs):
                ds = real_open(*args, **kwargs)
                close = ds._close
                marker = object()
                opened.append(marker)
                def traced_close():
                    closed.append(marker)
                    close()
                ds.set_close(traced_close)
                return ds

            def traced_load(data, *args, **kwargs):
                if data.name == "tas":
                    self.assertGreater(len(opened), len(closed))
                    loaded.append(tuple(data.shape))
                return real_load(data, *args, **kwargs)

            with mock.patch.object(xr, "open_dataset", side_effect=traced_open), \
                    mock.patch.object(xr.DataArray, "load", traced_load):
                result = io.open_dataarray(path, "tas", bbox={"lat": [2, 3], "lon": [4, 5]},
                                           time_indices=[0, 1], max_values=8)
            self.assertEqual(loaded, [(2, 2, 2)])
            self.assertEqual(len(closed), len(opened))
            Path(path).unlink()
            np.testing.assert_array_equal(result.values, np.full((2, 2, 2), 280.))

    def test_yearly_loads_only_requested_year_and_region(self):
        with tempfile.TemporaryDirectory() as tmp:
            dates = np.concatenate([monthly(year) for year in (1999, 2000, 2001)])
            path = write_field(Path(tmp) / "raw.nc", dates, lat=range(10), lon=range(10))
            loaded = []
            real_load = xr.DataArray.load
            def traced_load(data, *args, **kwargs):
                if data.name == "tas":
                    loaded.append(tuple(data.shape))
                return real_load(data, *args, **kwargs)
            with mock.patch.object(xr.DataArray, "load", traced_load):
                result = io.load_yearly_field([path], "tas", [1, 2], "mean", None,
                                              bbox={"lat": [1, 2], "lon": [2, 3]},
                                              year_range=[2000, 2000], frequency="monthly")
            self.assertEqual(loaded, [(2, 2, 2)])
            np.testing.assert_array_equal(result.year, [2000])
            np.testing.assert_allclose(result, 280.)

    def test_daily_memory_guard_precedes_value_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            dates = np.arange(np.datetime64("2000-01-01"), np.datetime64("2000-02-01"))
            path = write_field(Path(tmp) / "raw.nc", dates)
            spec = {"kind": "grid", "paths": [path], "variable": "tas"}
            with mock.patch.object(xr.DataArray, "load", side_effect=AssertionError("values were loaded")):
                with self.assertRaisesRegex(ConfigError, "memory guard"):
                    io.load_daily_field(spec, [2000, 2000], [1], max_values=20)

    def test_year_group_guard_checks_combined_file_size_before_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = [write_field(Path(tmp) / "a.nc", monthly(2000)[:6]),
                     write_field(Path(tmp) / "b.nc", monthly(2000)[6:])]
            with mock.patch.object(xr.DataArray, "load", side_effect=AssertionError("values were loaded")):
                with self.assertRaisesRegex(ConfigError, "memory guard"):
                    io.load_yearly_field(files, "tas", list(range(1, 13)), "mean", None,
                                          frequency="monthly", max_values=30)

    def test_daily_guard_counts_cropped_requested_dates(self):
        with tempfile.TemporaryDirectory() as tmp:
            dates = np.arange(np.datetime64("1999-12-01"), np.datetime64("2000-03-01"))
            path = write_field(Path(tmp) / "raw.nc", dates, lat=range(10), lon=range(10))
            spec = {"kind": "grid", "paths": [path], "variable": "tas",
                    "bbox": {"lat": [0, 0], "lon": [0, 0]}}
            result = io.load_daily_field(spec, [2000, 2000], [2], max_values=29)
            self.assertEqual(result.shape, (29, 1, 1))

    def test_disjoint_years_with_different_units_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = [write_field(Path(tmp) / "a.nc", monthly(2000), units="K"),
                     write_field(Path(tmp) / "b.nc", monthly(2001), units="degC")]
            with self.assertRaisesRegex(ConfigError, "units"):
                io.load_yearly_field(files, "tas", [1], "mean", None, frequency="monthly")

    def test_equivalent_unit_aliases_can_be_combined(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = [write_field(Path(tmp) / "a.nc", monthly(2000), units="K"),
                     write_field(Path(tmp) / "b.nc", monthly(2001), units="kelvin")]
            result = io.load_yearly_field(files, "tas", [1], "mean", None, frequency="monthly")
            np.testing.assert_array_equal(result.year, [2000, 2001])

    def test_declared_source_units_must_match_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_field(Path(tmp) / "raw.nc", monthly(2000))
            with self.assertRaisesRegex(ConfigError, "units_in"):
                io.open_dataarray(path, "tas", units_in="degC")

    def test_disjoint_file_identity_conflicts_fail(self):
        for key, first, second in [("source_id", "M1", "M2"),
                                   ("variant_label", "r1i1p1f1", "r2i1p1f1"),
                                   ("experiment_id", "historical", "ssp245"),
                                   ("grid_label", "gn", "gr")]:
            with self.subTest(key=key), tempfile.TemporaryDirectory() as tmp:
                files = [write_field(Path(tmp) / "a.nc", monthly(2000), metadata={key: first}),
                         write_field(Path(tmp) / "b.nc", monthly(2001), metadata={key: second})]
                with self.assertRaisesRegex(ConfigError, key):
                    io.load_yearly_field(files, "tas", [1], "mean", None, frequency="monthly")

    def test_disjoint_calendar_conflict_fails(self):
        import cftime
        with tempfile.TemporaryDirectory() as tmp:
            noleap = [cftime.DatetimeNoLeap(2001, month, 15) for month in range(1, 13)]
            files = [write_field(Path(tmp) / "a.nc", monthly(2000)),
                     write_field(Path(tmp) / "b.nc", noleap)]
            with self.assertRaisesRegex(ConfigError, "calendar"):
                io.load_yearly_field(files, "tas", [1], "mean", None, frequency="monthly")

    def test_disjoint_grid_conflict_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = [write_field(Path(tmp) / "a.nc", monthly(2000)),
                     write_field(Path(tmp) / "b.nc", monthly(2001), lat=(0., 2.))]
            with self.assertRaisesRegex(ConfigError, "grids"):
                io.load_yearly_field(files, "tas", [1], "mean", None, frequency="monthly")

    def test_history_future_member_from_headers_cannot_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            for experiment, year, member in [("historical", 2000, "r1i1p1f1"), ("ssp245", 2050, "r2i1p1f1")]:
                write_field(Path(tmp) / f"M_{experiment}.nc", monthly(year),
                            metadata={"source_id": "M", "experiment_id": experiment, "variant_label": member})
            cfg = Config({"project": {"outdir": str(Path(tmp) / "out")}, "runtime": {"cache": False},
                          "data": {"reference": {"kind": "series"},
                                   "models": {"kind": "grid", "root": tmp,
                                              "pattern": "{root}/{model}_{experiment}.nc",
                                              "members": ["*"], "variable": "tas",
                                              "aggregate": {"months": [1], "frequency": "monthly"}}}})
            source = io.DataSource(cfg)
            source.model_field("M", "historical", [2000, 2000])
            with self.assertRaisesRegex(ConfigError, "variant_label"):
                source.model_field("M", "ssp245", [2050, 2050])

    def test_bad_time_fails_before_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_field(Path(tmp) / "raw.nc", [2000., 2000.5])
            with self.assertRaisesRegex(ConfigError, "time"):
                io.open_dataarray(path, "tas")

    def test_monthly_duplicate_month_is_rejected_even_with_incomplete_allowed(self):
        data = xr.DataArray([[[1.]], [[2.]]], dims=("time", "latitude", "longitude"),
                            coords={"time": np.array(["2000-01-01", "2000-01-15"], dtype="datetime64[ns]"),
                                    "latitude": [0.], "longitude": [0.]})
        with self.assertRaisesRegex(ConfigError, "one sample per month"):
            io.aggregate_yearly(data, [1], frequency="monthly", require_complete=False)

    def test_sparse_monthly_metadata_preserves_djf_selection(self):
        dates = np.array(["1999-12-15", "2000-12-15", "2001-12-15"], dtype="datetime64[ns]")
        values = np.array([280., 281., 10000.])[:, None, None]
        for metadata in ({"frequency": "mon", "table_id": "Amon"},
                         {"table_id": "Amon"}, {"frequency": "monthly"}):
            with self.subTest(metadata=metadata), tempfile.TemporaryDirectory() as tmp:
                path = write_field(Path(tmp) / "raw.nc", dates, values=values, metadata=metadata)
                result = io.load_yearly_field([path], "tas", [12, 1, 2], "mean", None,
                                              year_range=[2000, 2001], frequency="auto",
                                              require_complete=False)
                np.testing.assert_array_equal(result.year, [2000, 2001])
                np.testing.assert_allclose(result.values[:, 0, 0], [280., 281.])
                with self.assertRaisesRegex(ConfigError, "no complete season-years"):
                    io.load_yearly_field([path], "tas", [12, 1, 2], "mean", None,
                                          year_range=[2000, 2001], frequency="auto")

    def test_frequency_metadata_conflicts_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_field(Path(tmp) / "raw.nc", monthly(2000),
                               metadata={"frequency": "day", "table_id": "Amon"})
            with self.assertRaises(ConfigError):
                io.load_yearly_field([path], "tas", [1], "mean", None)
            path = write_field(Path(tmp) / "raw.nc", monthly(2000), metadata={"frequency": "mon"})
            with self.assertRaisesRegex(ConfigError, "frequency"):
                io.load_yearly_field([path], "tas", [1], "mean", None, frequency="annual")

    def test_sampling_metadata_retained_and_checked_across_chunks(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = write_field(Path(tmp) / "a.nc", monthly(2000), metadata={"frequency": "mon"})
            second = write_field(Path(tmp) / "b.nc", monthly(2001), metadata={"frequency": "monthly"})
            result = io.load_yearly_field([first, second], "tas", [1], "mean", None, frequency="auto")
            self.assertEqual(result.attrs["_cmip_sampling_frequency"], "monthly")
            second = write_field(Path(tmp) / "b.nc", monthly(2001), metadata={"frequency": "day"})
            with self.assertRaisesRegex(ConfigError, "frequency"):
                io.load_yearly_field([first, second], "tas", [1], "mean", None)

    def test_daily_metadata_is_not_reinferred_from_one_selected_day_per_year(self):
        with tempfile.TemporaryDirectory() as tmp:
            dates = np.array(["2000-01-01", "2001-01-01"], dtype="datetime64[ns]")
            path = write_field(Path(tmp) / "raw.nc", dates, metadata={"frequency": "day", "table_id": "day"})
            result = io.load_daily_field({"kind": "grid", "paths": [path], "variable": "tas"},
                                          [2000, 2001], [1])
            self.assertEqual(result.time.size, 2)
            self.assertEqual(result.attrs["_cmip_sampling_frequency"], "daily")

    def test_cf_packed_fill_value_stays_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "packed.nc"
            data = xr.DataArray([[[280.]], [[np.nan]]], dims=("time", "lat", "lon"),
                                coords={"time": monthly(2000)[:2], "lat": [0.], "lon": [0.]},
                                name="tas", attrs={"units": "K"})
            data.to_netcdf(path, encoding={"tas": {"dtype": "int16", "scale_factor": .1,
                                                  "add_offset": 273.15, "_FillValue": -32768}})
            result = io.open_dataarray(str(path), "tas")
            self.assertTrue(np.isnan(result.values[1, 0, 0]))
            self.assertAlmostEqual(float(result.values[0, 0, 0]), 280., delta=.051)

    def test_broken_cache_is_rebuilt_with_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = Config({"project": {"outdir": tmp}, "data": {"reference": {"kind": "series"}}})
            source = io.DataSource(cfg)
            key = {"tag": "test", "loader_version": 3}
            cache = Path(source._cache_path(key)); cache.write_bytes(b"incomplete NetCDF")
            field = xr.DataArray([[[2.]]], dims=("year", "latitude", "longitude"),
                                 coords={"year": [2000], "latitude": [0.], "longitude": [0.]})
            builder = mock.Mock(return_value=field)
            with self.assertLogs("cmipkit.io", level="WARNING") as logged:
                first = source._cached(key, builder)
            second = source._cached(key, builder)
            self.assertEqual(builder.call_count, 1)
            self.assertIn("rebuilding", " ".join(logged.output))
            xr.testing.assert_equal(first, second)

    def test_failed_output_write_preserves_previous_file_and_cleans_temp(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "result.nc"; path.write_bytes(b"old result")
            field = xr.DataArray([1.], dims="x", name="field")
            def interrupted(data, temporary, **kwargs):
                Path(temporary).write_bytes(b"partial result")
                raise OSError("disk full")
            with mock.patch.object(xr.DataArray, "to_netcdf", interrupted):
                with self.assertRaisesRegex(OSError, "disk full"):
                    io.save_field(field, str(path))
            self.assertEqual(path.read_bytes(), b"old result")
            self.assertEqual(list(Path(tmp).iterdir()), [path])

    def test_cached_field_respects_new_memory_limit_before_loading(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = Config({"project": {"outdir": tmp}, "data": {"reference": {"kind": "series"}}})
            source = io.DataSource(cfg)
            key = {"tag": "test", "loader_version": 3}
            field = xr.DataArray([[[2.]], [[3.]]], dims=("year", "latitude", "longitude"),
                                 coords={"year": [2000, 2001], "latitude": [0.], "longitude": [0.]})
            source._cached(key, lambda: field)
            cfg["runtime"] = {"max_input_values": 1}
            with mock.patch.object(xr.DataArray, "load", side_effect=AssertionError("values were loaded")):
                with self.assertRaisesRegex(ConfigError, "memory guard"):
                    source._cached(key, lambda: field)

    def test_failed_npz_write_is_atomic(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "result.npz"; path.write_bytes(b"old result")
            with mock.patch.object(np, "savez_compressed", side_effect=OSError("disk full")):
                with self.assertRaises(OSError):
                    io.save_npz(str(path), values=[1.])
            self.assertEqual(path.read_bytes(), b"old result")
            self.assertEqual(list(Path(tmp).iterdir()), [path])


if __name__ == "__main__":
    unittest.main(verbosity=2)
