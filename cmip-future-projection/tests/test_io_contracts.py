"""Synthetic regression checks for time, coordinates, identity and missingness.

Run: python tests/test_io_contracts.py
"""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from cmipkit import indices, io
from cmipkit.config import Config, ConfigError


def field(times, values, lat=(0.,), lon=(0.,)):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = np.broadcast_to(values[:, None, None], (len(times), len(lat), len(lon))).copy()
    return xr.DataArray(values, dims=("time", "latitude", "longitude"),
                        coords={"time": times, "latitude": list(lat), "longitude": list(lon)}, name="tas")


class IOContracts(unittest.TestCase):
    def test_no_artificial_regional_seam(self):
        data = field(np.arange(2), [1, 2], lon=(100., 110., 120.))
        self.assertEqual(io.align_longitude(data).longitude.size, 3)

    def test_longitude_boxes_both_conventions(self):
        lon = np.array([0., 10., 20., 170., 180., 190., 350.])
        for coords in (lon, (lon + 180) % 360 - 180):
            np.testing.assert_array_equal(indices.longitude_mask(coords, [350, 20]),
                                          [True, True, True, False, False, False, True])
            np.testing.assert_array_equal(indices.longitude_mask(coords, [170, -170]),
                                          [False, False, False, True, True, True, False])
            self.assertTrue(indices.longitude_mask(coords, [-180, 180]).all())

    def test_bbox_wrap(self):
        data = field(np.arange(2), [1, 2], lon=(0., 10., 100., 350.))
        out = io.select_bbox(data, {"lat": [-10, 10], "lon": [350, 10]})
        np.testing.assert_array_equal(out.longitude, [0, 10, 350])

    def test_global_regrid_negative_target(self):
        data = field([2000], np.array([[[1., 0., -1., 0.]]]), lon=(0., 90., 180., 270.))
        out = io.regrid(data, {"latitude": np.array([0.]), "longitude": np.array([-180., -90., 0., 90.])})
        np.testing.assert_allclose(out.values[0, 0], [-1, 0, 1, 0])
        np.testing.assert_array_equal(out.longitude, [-180, -90, 0, 90])

    def test_regional_regrid_does_not_bridge_unsampled_globe(self):
        data = field([2000], np.array([[[0., 10., -10.]]]), lon=(0., 10., 350.))
        out = io.regrid(data, {"latitude": np.array([0.]), "longitude": np.array([5., 180., 355.])})
        np.testing.assert_allclose(out.values[0, 0], [5., np.nan, -5.], equal_nan=True)

    def test_multifile_winter_month_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/'a.nc', Path(tmp)/'b.nc'
            field(np.array(['1999-12-15'], dtype='datetime64[ns]'), [12]).to_netcdf(a)
            field(np.array(['2000-01-15', '2000-02-15'], dtype='datetime64[ns]'), [1, 2]).to_netcdf(b)
            out = io.load_yearly_field([str(a), str(b)], 'tas', [12, 1, 2], 'mean', None,
                                      frequency='monthly')
            np.testing.assert_array_equal(out.year, [2000])
            self.assertAlmostEqual(float(out.values[0, 0, 0]), (12*31 + 1*31 + 2*29)/91)

    def test_extended_cross_year_season(self):
        times = np.array(['1999-11-15', '1999-12-15', '2000-01-15', '2000-02-15', '2000-03-15'], dtype='datetime64[ns]')
        np.testing.assert_array_equal(io._season_year(times, [11, 12, 1, 2, 3]), [2000]*5)

    def test_incomplete_season_is_not_annual_sample(self):
        times = np.array(['1999-12-15', '2000-01-15'], dtype='datetime64[ns]')
        data = field(times, [1, 3])
        self.assertEqual(io.aggregate_yearly(data, [12, 1, 2], frequency='monthly').year.size, 0)
        out = io.aggregate_yearly(data, [12, 1, 2], frequency='monthly', require_complete=False)
        self.assertEqual(out.year.size, 1)

    def test_missing_values_do_not_become_zero_events(self):
        times = np.array(['2000-01-15', '2000-02-15'], dtype='datetime64[ns]')
        data = field(times, [np.nan, np.nan])
        for stat, threshold in [('sum', None), ('count_gt', 0)]:
            out = io.aggregate_yearly(data, [1, 2], stat, threshold, min_valid_fraction=0)
            self.assertTrue(np.isnan(out.values).all())

    def test_valid_fraction(self):
        times = np.array(['2000-01-15', '2000-02-15'], dtype='datetime64[ns]')
        data = field(times, [3, np.nan])
        self.assertTrue(np.isnan(io.aggregate_yearly(data, [1, 2]).values).all())
        self.assertEqual(float(io.aggregate_yearly(data, [1, 2], min_valid_fraction=0.5).values[0, 0, 0]), 3)

    def test_360_day_month_weights(self):
        import cftime
        times = [cftime.Datetime360Day(2000, m, 15) for m in [1, 2, 3]]
        data = field(times, [1, 2, 9])
        self.assertEqual(float(io.aggregate_yearly(data, [1, 2, 3]).values[0, 0, 0]), 4.)

    def test_daily_completeness_and_leap_calendar(self):
        times = np.arange(np.datetime64('2000-02-01'), np.datetime64('2000-03-01'))
        data = field(times, np.ones(29))
        self.assertEqual(float(io.aggregate_yearly(data, [2], 'sum').values[0, 0, 0]), 29.)
        self.assertEqual(io.aggregate_yearly(data.isel(time=slice(1,None)), [2], 'sum').year.size, 0)

    def test_annual_input(self):
        data = field([2000, 2001], [4, 5])
        out = io.aggregate_yearly(data, list(range(1, 13)))
        np.testing.assert_array_equal(out.year, [2000, 2001])
        np.testing.assert_array_equal(out.values[:, 0, 0], [4, 5])

    def test_conflicting_duplicate_times_fail(self):
        data = field(np.array(['2000-01-15', '2000-01-15'], dtype='datetime64[ns]'), [1, 2])
        with self.assertRaises(ConfigError):
            io._deduplicate_time(data)
        data.values[1] = 1
        self.assertEqual(io._deduplicate_time(data).time.size, 1)

    def test_explicit_member_never_falls_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, 'Model_A_r2i1p1f1_historical.nc').touch()
            spec = {'root': tmp, 'kind': 'grid', 'pattern': '{root}/{model}_{member}_{experiment}.nc',
                    'members': ['r1i1p1f1'], 'variable': 'tas'}
            self.assertEqual(io._files_for_grid_source(spec, 'historical', 'Model_A', [2000, 2001]), [])
            self.assertEqual(io.discover_members(spec, 'historical', 'Model_A'), ['r2i1p1f1'])

    def test_explicit_pattern_not_replaced(self):
        spec = {'pattern': '{root}/custom/{model}.nc'}
        self.assertEqual(io.grid_patterns(spec), ['{root}/custom/{model}.nc'])

    def test_absolute_path_model_discovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, 'Model_A_r1i1p1f1_historical.nc').touch()
            cfg = Config({'project': {'outdir': str(Path(tmp)/'out')}, 'runtime': {'cache': True},
                          'data': {'reference': {'kind': 'series'}, 'models': {'kind': 'series', 'series_dir': tmp}}})
            source = io.DataSource(cfg, create_dirs=False)
            spec = {'root': tmp, 'pattern': '{root}/{model}_{member}_{experiment}.nc',
                    'members': ['r1i1p1f1'], 'variable': 'tas'}
            self.assertEqual(source._discover_from_pattern(spec), ['Model_A'])
            self.assertFalse((Path(tmp)/'out').exists())

    def test_cache_includes_config_and_file_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.nc'; path.write_bytes(b'a')
            cfg = Config({'project': {'outdir': tmp}, 'data': {'reference': {'kind': 'series'}, 'models': {'kind': 'series'}}})
            source = io.DataSource(cfg)
            spec = {'units_transform': {'offset': 0}}
            first = source._cache_path(source._source_key('t', [str(path)], spec, [2000, 2001]))
            spec['units_transform']['offset'] = -273.15
            second = source._cache_path(source._source_key('t', [str(path)], spec, [2000, 2001]))
            path.write_bytes(b'changed')
            third = source._cache_path(source._source_key('t', [str(path)], spec, [2000, 2001]))
            self.assertEqual(len({first, second, third}), 3)

    def test_series_retains_missingness_and_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'x.csv'
            path.write_text('y,t\n2001,2\n2000,\n')
            years, values = io.read_series_file(str(path), {'year': 'y', 'value': 't'})
            np.testing.assert_array_equal(years, [2000, 2001])
            self.assertTrue(np.isnan(values[0]))
            path.write_text('year,value\n2000,1\n2000,2\n')
            with self.assertRaises(ConfigError):
                io.read_series_file(str(path))

    def test_npz_requires_named_numeric_arrays(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'x.npz'
            np.savez(path, years=[2001, 2000], values=[2., 1.])
            years, values = io.read_series_file(str(path))
            np.testing.assert_array_equal(years, [2000, 2001])
            np.testing.assert_array_equal(values, [1, 2])
            np.savez(path, x=[2000], y=[1])
            with self.assertRaises(ConfigError):
                io.read_series_file(str(path))

    def test_zero_scale_is_not_replaced_by_one(self):
        out = io.convert_units(xr.DataArray([1., 2.]), {'scale': 0, 'offset': 3})
        np.testing.assert_array_equal(out, [3, 3])

    def test_projection_all_missing_returns_nan(self):
        data = field([2000, 2001], np.array([[[-1., 1.]], [[np.nan, np.nan]]]), lon=(0., 1.)).rename(time='year')
        template = data.isel(year=0)
        result = indices.projection_series(data, template)
        self.assertAlmostEqual(result[0], 1.)
        self.assertTrue(np.isnan(result[1]))

    def test_fixed_baseline_standardization(self):
        data = field([2000, 2001, 2002], [1., 2., 3.]).rename(time='year')
        spec = {'type': 'region_mean', 'region': {'lat': [-1, 1], 'lon': [-1, 1]}, 'standardize': True}
        _, historical = indices.build_index(data, spec, [2000, 2002])
        np.testing.assert_allclose(historical, [-1, 0, 1])
        future = field([2100, 2101], [11., 12.]).rename(time='year')
        with self.assertRaises(ValueError):
            indices.build_index(future, spec, [2000, 2002])
        _, values = indices.build_index(future, spec, [2000, 2002], standardization=(2., 1.))
        np.testing.assert_allclose(values, [9, 10])

    def test_detrending_uses_each_cells_available_years(self):
        data = field(np.arange(2000, 2005), [np.nan, 3., 5., 7., 9.]).rename(time='year')
        values = indices.detrend_field(data).values[:, 0, 0]
        self.assertTrue(np.isnan(values[0]))
        np.testing.assert_allclose(values[1:], 0., atol=1e-12)

    def test_pressure_selection_is_exact_and_units_aware(self):
        data = xr.DataArray([1., 2.], dims='plev', coords={'plev': [50000, 85000]})
        self.assertEqual(float(io.select_plev(data, 500)), 1.)
        with self.assertRaises(ConfigError):
            io.select_plev(data, 700)

    def test_irregular_grid_does_not_claim_cosine_area_weights(self):
        data = field([2000], [1.], lat=(0., 10., 30.)).rename(time='year')
        with self.assertRaises(ValueError):
            indices.region_mean_series(data, {'lat': [-90, 90], 'lon': [-180, 180]})

    def test_inferred_reference_grid_and_cache_reload(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.nc'
            times = np.array(['2000-01-15', '2000-02-15'], dtype='datetime64[ns]')
            field(times, [1., 2.], lon=(0., 1., 2.)).to_netcdf(path)
            cfg = Config({'project': {'outdir': str(Path(tmp)/'output')},
                          'data': {'reference': {'kind': 'grid', 'paths': [str(path)], 'variable': 'tas',
                                                'aggregate': {'months': [1, 2]}}, 'models': {'kind': 'grid'}},
                          'periods': {'historical': [2000, 2000]}})
            first = io.DataSource(cfg); a = first.reference_field()
            second = io.DataSource(cfg); b = second.reference_field()
            xr.testing.assert_equal(a, b)
            np.testing.assert_array_equal(second.grid['longitude'], [0, 1, 2])

    def test_daily_loader_cross_year_paths_preserves_missing_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            for year, start, end in [(1999, '1999-12-01', '2000-01-01'), (2000, '2000-01-01', '2000-03-01')]:
                dates = np.arange(np.datetime64(start), np.datetime64(end))
                values = np.ones(dates.size); values[0] = np.nan
                field(dates, values).to_netcdf(Path(tmp)/f'{year}.nc')
            spec = {'kind': 'grid', 'paths': [str(Path(tmp)/'{year}.nc')], 'variable': 'tas', 'bbox': {'lat': [-1,1], 'lon': [-1,1]}}
            loaded = io.load_daily_field(spec, [2000, 2000], [12, 1, 2])
            self.assertEqual(loaded.time.size, 91)
            self.assertEqual(int(np.isnan(loaded.values).sum()), 2)

    def test_identical_monthly_overlap_not_misclassified_as_daily(self):
        with tempfile.TemporaryDirectory() as tmp:
            dates = np.array(['2000-01-15', '2000-02-15'], dtype='datetime64[ns]')
            files = []
            for name in ('a', 'b'):
                path = Path(tmp)/f'{name}.nc'; field(dates, [1., 2.]).to_netcdf(path); files.append(str(path))
            out = io.load_yearly_field(files, 'tas', [1, 2], 'mean', None)
            self.assertEqual(out.year.size, 1)

    def test_series_units_and_model_name_with_underscore(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, 'Model_A_historical.csv').write_text('Y,T\n2000,273.15\n2001,274.15\n')
            cfg = Config({'project': {'outdir': tmp}, 'data': {'reference': {'kind': 'series'},
                          'models': {'kind': 'series', 'series_dir': tmp, 'series_columns': {'year':'Y', 'value':'T'},
                                     'units_transform': {'offset': -273.15}}}})
            source = io.DataSource(cfg)
            self.assertEqual(source.discover_models(), ['Model_A'])
            _, values = source.model_series('Model_A', 'historical')
            np.testing.assert_allclose(values, [0., 1.])


if __name__ == '__main__':
    unittest.main(verbosity=2)
