"""Regression checks for statistical contracts; synthetic inputs only.

Run: python tests/test_statistical_contracts.py
"""

from pathlib import Path
import sys
import unittest
import warnings

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from cmipkit import correction, daily, ensemble, extremes
from cmipkit.config import Config


def field(values, sample="year"):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None, None]
    dates = np.arange(1980, 1980 + values.shape[0])
    if sample == "time":
        dates = np.datetime64("2000-01-01") + np.arange(values.shape[0])
    return xr.DataArray(values, dims=(sample, "latitude", "longitude"),
                        coords={sample: dates, "latitude": np.arange(values.shape[1]),
                                "longitude": np.arange(values.shape[2])},
                        attrs={"units": "degC"})


class CorrectionContracts(unittest.TestCase):
    def test_scalar_qm_maps_values_not_probabilities(self):
        model = np.arange(280.0, 300.0)
        observed = model - 4
        p = correction.fit_series(observed, model, "quantile_mapping", n_quantiles=20)
        actual = correction.apply_series(p, model)
        np.testing.assert_allclose(actual, observed)
        self.assertGreater(np.ptp(actual), 10)

    def test_scalar_and_grid_qm_agree_with_missing_and_ties(self):
        model = np.array([0, 0, 0, 0, 1, 2, 4, 8, 16, 20], dtype=float)
        observed = model * 2
        values = np.array([np.nan, 0, 2, 20, np.nan])
        scalar = correction.apply_series(correction.fit_series(observed, model, "quantile_mapping"), values)
        grid = correction.apply(correction.fit(field(observed), field(model), "quantile_mapping"), field(values))
        np.testing.assert_allclose(scalar, grid.values.ravel(), equal_nan=True)

    def test_qdm_additive_change_and_multiplicative_change(self):
        model = np.linspace(1, 20, 20)
        for mult, observed, future, expected in [
            (False, model - 3, model + 5, model + 2),
            (True, model * 2, model * 3, model * 6),
        ]:
            params = correction.fit_series(observed, model, "quantile_delta_mapping", mult, 20)
            np.testing.assert_allclose(correction.apply_series(params, future), expected, atol=1e-10)
            gp = correction.fit(field(observed), field(model), "quantile_delta_mapping", mult, 20)
            np.testing.assert_allclose(correction.apply(gp, field(future)).values.ravel(), expected, atol=1e-10)

    def test_missing_grid_baseline_stays_missing(self):
        values = np.tile(np.arange(10.0)[:, None, None], (1, 1, 2))
        values[:, :, 1] = np.nan
        for method in ["mean", "mean_std", "quantile_mapping", "quantile_delta_mapping"]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                params = correction.fit(field(values), field(values + 2), method)
                result = correction.apply(params, field(np.ones_like(values)))
            self.assertTrue(np.isnan(result.values[:, :, 1]).all(), method)

    def test_constant_quantiles_do_not_fill_nan(self):
        params = correction.fit(field(np.ones(10) * 3), field(np.ones(10) * 5), "quantile_mapping")
        result = correction.apply(params, field([5, np.nan, 5]))
        np.testing.assert_allclose(result.values.ravel(), [3, np.nan, 3], equal_nan=True)

    def test_transposed_time_arrays_keep_dimensions_and_units(self):
        obs = field(np.arange(10.0), "time").transpose("longitude", "time", "latitude")
        mod = obs + 3
        mod.attrs["units"] = "degC"
        result = correction.apply(correction.fit(obs, mod, "mean"), mod)
        self.assertEqual(result.dims, mod.dims)
        self.assertEqual(result.attrs["units"], "degC")
        np.testing.assert_allclose(result, obs)

    def test_grid_mismatch_requires_explicit_regridding(self):
        obs = field(np.arange(10.0))
        mod = field(np.arange(10.0)).assign_coords(longitude=[10])
        with self.assertRaisesRegex(ValueError, "coordinates differ"):
            correction.fit(obs, mod, "mean")

    def test_pressure_does_not_auto_select_precipitation_ratios(self):
        cfg = Config({"correction": {"multiplicative": "auto"},
                      "variable": {"name": "ps", "long_name": "surface pressure"}})
        self.assertFalse(correction.resolve_multiplicative(cfg, np.ones(5), np.ones(5)))
        cfg["variable"] = {"name": "pr", "long_name": "precipitation"}
        self.assertTrue(correction.resolve_multiplicative(cfg, np.ones(5), np.ones(5)))

    def test_zero_model_baseline_is_not_a_wet_reference_ratio(self):
        with self.assertRaisesRegex(ValueError, "entirely zero"):
            correction.fit_series(np.ones(10), np.zeros(10), "mean", True)
        with self.assertRaisesRegex(ValueError, "wet-day occurrence"):
            correction.fit_series(np.r_[np.ones(5), np.arange(2, 7)],
                                  np.r_[np.zeros(5), np.arange(1, 6)],
                                  "quantile_delta_mapping", True, 10)
        model = np.r_[np.zeros(5), np.arange(1, 6)]
        params = correction.fit_series(model * 2, model, "quantile_delta_mapping", True, 10)
        np.testing.assert_allclose(correction.apply_series(params, model * 3), model * 6)

    def test_none_does_not_require_baseline_and_bounds_work(self):
        params = correction.fit_series([], [], "none")
        np.testing.assert_equal(correction.apply_series(params, [2, 3]), [2, 3])
        params = correction.fit_series(np.arange(10), np.arange(10), "mean")
        params.update(lower_bound=0, upper_bound=1)
        np.testing.assert_equal(correction.apply_series(params, [-2, 0.5, 2]), [0, 0.5, 1])

    def test_difference_trend_accepts_zero_model_slope(self):
        result = correction.trend_correct_series(np.array([2000, 2001]), np.zeros(2), 0, 2, 2000)
        np.testing.assert_equal(result, [0, 2])


class ExtremeContracts(unittest.TestCase):
    def test_pot_return_period_is_inverse_for_both_shape_signs(self):
        for shape in [-0.2, 0, 0.2]:
            fit = daily.POTFit(threshold=30, shape=shape, scale=2, rate=3, n_peaks=90, n_years=30)
            for period in [1, 2, 10, 100]:
                self.assertAlmostEqual(fit.return_period(fit.return_level(period)), period, places=8)

    def test_pot_does_not_extrapolate_below_threshold(self):
        fit = daily.POTFit(threshold=10, shape=0.2, scale=2, rate=0.1, n_peaks=3, n_years=30)
        self.assertTrue(np.isnan(fit.return_level(2)))
        self.assertAlmostEqual(fit.return_level(10), 10)

    def test_gev_return_period_uses_stable_survival_function(self):
        fit = extremes.Fit("gumbel", 0, 0, 1)
        level = fit.return_level(1e18)
        self.assertTrue(np.isfinite(level))
        self.assertAlmostEqual(fit.return_period(level) / 1e18, 1, places=12)

    def test_fitted_ks_p_is_not_presented_as_calibrated(self):
        rng = np.random.default_rng(42)
        fit = extremes.fit_distribution(rng.gumbel(size=100), "gumbel")
        self.assertIsNotNone(fit)
        self.assertTrue(np.isnan(fit.ks_p))
        self.assertTrue(np.isfinite(fit.ks_p_naive))
        self.assertEqual(fit.as_dict()["shape_xi"], -fit.shape)

    def test_declustering_respects_missing_days_and_season_breaks(self):
        dates = np.array(["2001-08-30", "2001-08-31", "2002-06-01", "2002-06-02"], dtype="datetime64[D]")
        indices, values = daily.decluster(np.array([0, 10, 11, 0]), 5, dates=dates)
        np.testing.assert_equal(indices, [1, 2])
        indices, values = daily.decluster(np.array([10, np.nan, 11]), 5)
        np.testing.assert_equal(indices, [0, 2])

    def test_daily_timestamp_validation_rejects_duplicates(self):
        with self.assertRaisesRegex(ValueError, "sorted and unique"):
            daily.decluster([10, 11], 5, dates=np.array(["2001-01-01", "2001-01-01"], dtype="datetime64[D]"))

    def test_pot_keeps_nan_separation_when_fitting(self):
        # Each observed peak is separated by an unknown day; finite filtering
        # before declustering would collapse all peaks into one cluster.
        values = np.full(80, np.nan)
        values[::2] = np.arange(1, 41)
        fit = daily.fit_pot(values, 0, 10, min_peaks=20)
        self.assertIsNotNone(fit)
        self.assertEqual(fit.n_peaks, 40)

    def test_fractional_seasonal_exposure(self):
        dates = np.datetime64("2001-06-01") + np.arange(92)
        values = np.ones(92)
        self.assertEqual(daily.exposure_years(dates, values, [6, 7, 8]), 1)
        values[0] = np.nan
        self.assertAlmostEqual(daily.exposure_years(dates, values, [6, 7, 8]), 91 / 92)

    def test_djf_and_360_day_calendar_exposure(self):
        import cftime
        import datetime
        dates = np.array([cftime.Datetime360Day(2000, 12, 1) + datetime.timedelta(days=i)
                          for i in range(90)])
        values = np.arange(90.0)
        self.assertEqual(daily.exposure_years(dates, values, [1, 2, 12]), 1)
        years, maxima = daily.annual_maxima(dates, values, months=[12, 1, 2], min_coverage=1)
        np.testing.assert_equal(years, [2001])
        np.testing.assert_equal(maxima, [89])

    def test_block_drops_incomplete_seasons(self):
        dates = np.datetime64("2001-06-01") + np.arange(92)
        values = np.arange(92.0)
        values[0] = np.nan
        years, maxima = daily.annual_maxima(dates, values, months=[6, 7, 8], min_coverage=1)
        self.assertEqual(years.size, 0)


class EnsembleContracts(unittest.TestCase):
    def test_weighted_median_and_bands_share_weights(self):
        result = ensemble.stats(np.array([[0], [10], [100]]), "median", [10, 90], np.array([0.1, 0.8, 0.1]))
        np.testing.assert_equal(result["center"], [10])
        np.testing.assert_equal(result["p90"], [10])
        self.assertAlmostEqual(result["n_effective"][0], 1 / 0.66)

    def test_missing_models_renormalize_weights(self):
        result = ensemble.stats(np.array([[2, np.nan], [10, 10]]), "mean", [10, 90], np.array([3, 1]))
        np.testing.assert_allclose(result["center"], [4, 10])
        np.testing.assert_equal(result["n_models"], [2, 1])

    def test_empty_alignment_and_duplicate_years(self):
        matrix, names = ensemble.align({"m": (np.array([2000]), np.array([1]))}, [])
        self.assertEqual(matrix.shape, (1, 0))
        with self.assertRaisesRegex(ValueError, "duplicate years"):
            ensemble.align({"m": ([2000, 2000], [1, 2])}, [2000])

    def test_invalid_weights_are_not_silently_equal_weights(self):
        for weights in ([1], [-1, 2], [0, 0], [np.nan, 1]):
            with self.assertRaises(ValueError):
                ensemble.stats(np.ones((2, 3)), weights=weights)


if __name__ == "__main__":
    unittest.main(verbosity=2)
