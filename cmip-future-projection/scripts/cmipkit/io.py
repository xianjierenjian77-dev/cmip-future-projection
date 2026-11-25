"""Input / output layer: open gridded or tabular data, harmonise it, aggregate it
to one value per year and put it on a common grid.

Everything downstream works with two very small data structures:

* ``Field``  - ``xarray.DataArray`` with dims ``(year, latitude, longitude)``
* ``Series`` - ``numpy.ndarray`` of shape ``(n_years,)`` + ``years`` vector

The functions here are the only place that knows about file formats, variable
names, GRIB vs NetCDF, 0-360 vs -180-180 longitudes or unit conversions, so a
new data source only has to be taught to produce one of the two structures.
"""

from __future__ import annotations

import glob
import hashlib
import json
import logging
import os
import re
import tempfile
from contextlib import contextmanager
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import xarray as xr

from .config import Config, ConfigError, resolve_months

log = logging.getLogger("cmipkit.io")

GRID_ENGINES = (".grib", ".grb", ".grib2", ".grb2")

PLEV_NAMES = ("plev", "isobaricInhPa", "level", "lev", "pressure")
LAT_NAMES = ("latitude", "lat", "rlat", "y")
LON_NAMES = ("longitude", "lon", "rlon", "x")
TIME_NAMES = ("time", "valid_time", "date", "year")


# ---------------------------------------------------------------------------
# generic helpers
# ---------------------------------------------------------------------------
def _pick_name(da: xr.DataArray, candidates: Sequence[str]) -> Optional[str]:
    for name in candidates:
        if name in da.dims or name in da.coords:
            return name
    return None


def normalize_dims(da: xr.DataArray) -> xr.DataArray:
    """Rename dimensions/coords to latitude/longitude/time/plev and fix order."""
    for name in ("latitude", "lat", "longitude", "lon"):
        if name in da.coords and da[name].ndim != 1:
            raise ConfigError("curvilinear coordinates need a dedicated regridder before this pipeline")
    if (("rlat" in da.dims or "rlon" in da.dims)
            and not any(n in da.coords for n in ("latitude", "lat", "longitude", "lon"))):
        raise ConfigError("rotated grid coordinates are not geographic latitude/longitude; regrid first")
    rename: Dict[str, str] = {}
    lat = _pick_name(da, LAT_NAMES)
    lon = _pick_name(da, LON_NAMES)
    tim = _pick_name(da, TIME_NAMES)
    plev = _pick_name(da, PLEV_NAMES)
    if lat and lat != "latitude":
        rename[lat] = "latitude"
    if lon and lon != "longitude":
        rename[lon] = "longitude"
    if tim and tim != "time":
        rename[tim] = "time"
    if plev and plev != "plev":
        rename[plev] = "plev"
    if rename:
        da = da.rename({k: v for k, v in rename.items() if k in da.dims or k in da.coords})

    # drop auxiliary coordinates that would block concatenation
    drop = [c for c in da.coords if c not in da.dims and c not in ("latitude", "longitude", "time", "plev")]
    if drop:
        da = da.drop_vars([c for c in drop if c in da.coords], errors="ignore")

    for name in ("latitude", "longitude"):
        if name in da.coords and da[name].dims != (name,):
            original = da[name].dims[0]
            da = da.swap_dims({original: name})
        if name in da.dims:
            vals = np.asarray(da[name].values, dtype=float)
            if not np.isfinite(vals).all():
                raise ConfigError(f"{name} coordinates must be finite")
    if "latitude" in da.dims:
        if np.any(np.abs(da.latitude.values) > 90):
            raise ConfigError("latitude must contain geographic degrees in [-90, 90]")
        if np.unique(da.latitude.values).size != da.latitude.size:
            raise ConfigError("duplicate latitude coordinates")
        da = da.sortby("latitude")
    expected = [n for n in ("time", "latitude", "longitude") if n in da.dims]
    if set(da.dims) == set(expected):
        da = da.transpose(*expected)
    return da


def align_longitude(da: xr.DataArray, lon_360: bool = True) -> xr.DataArray:
    """Return a copy whose longitude is sorted and spans 0-360 (or -180-180)."""
    if "longitude" not in da.dims:
        return da
    lon = da.longitude.values.astype(float)
    if lon_360:
        new = np.mod(lon, 360.0)
    else:
        new = ((lon + 180.0) % 360.0) - 180.0
    da = da.assign_coords(longitude=("longitude", new)).sortby("longitude")
    # a 0..360 inclusive grid produces a duplicate 0 after the modulo
    lon_sorted = da.longitude.values
    if lon_sorted.size and np.unique(lon_sorted).size != lon_sorted.size:
        duplicates, counts = np.unique(lon_sorted, return_counts=True)
        for value in duplicates[counts > 1]:
            columns = da.isel(longitude=np.flatnonzero(lon_sorted == value)).values
            axis = da.get_axis_num("longitude")
            copies = np.moveaxis(columns, axis, 0)
            if not np.allclose(copies, copies[0], equal_nan=True, rtol=1e-7, atol=1e-10):
                raise ConfigError(f"duplicate longitude {value} has conflicting values")
        _, keep = np.unique(lon_sorted, return_index=True)
        da = da.isel(longitude=np.sort(keep))
    # Do not append a cyclic column here: it would be counted as a real cell in
    # regional means, and a regional grid has no cyclic seam. Regridding adds
    # temporary halo cells only after identifying a global grid.
    return da


def convert_units(da: xr.DataArray, spec: Optional[Dict[str, Any]]) -> xr.DataArray:
    spec = spec or {}
    scale = float(spec.get("scale", 1.0))
    offset = float(spec.get("offset", 0.0))
    if not np.isfinite(scale) or not np.isfinite(offset):
        raise ConfigError("units_transform scale and offset must be finite")
    if scale == 1.0 and offset == 0.0 and not spec.get("units"):
        return da
    out = da * scale + offset
    out.attrs = dict(da.attrs)
    if spec.get("units"):
        out.attrs["units"] = str(spec["units"])
    else:
        # A stale K label on Celsius values is worse than no label. The
        # pipeline can set variable.units on the final product.
        out.attrs.pop("units", None)
    return out


def select_plev(da: xr.DataArray, plev: Any) -> xr.DataArray:
    """Select a pressure level given in hPa (ERA5) or Pa (CMIP6)."""
    if plev is None:
        return da
    if "plev" not in da.dims and "plev" not in da.coords:
        raise ConfigError("a pressure level was requested but the input has no pressure coordinate")
    vals = np.atleast_1d(np.asarray(da.plev.values, dtype=float))
    target = float(plev)
    if target <= 2000 and vals.max() > 5000:  # user gave hPa, data in Pa
        target *= 100.0
    elif target > 5000 and vals.max() <= 2000:
        target /= 100.0
    candidates = np.abs(vals - target)
    idx = int(np.argmin(candidates))
    if candidates[idx] > max(1e-6, 1e-4 * abs(target)):
        raise ConfigError(f"pressure level {plev} not found; available levels are {vals.tolist()}")
    out = da.isel(plev=idx) if "plev" in da.dims else da
    out = out.drop_vars("plev", errors="ignore")
    return out


def select_bbox(da: xr.DataArray, bbox: Optional[Dict[str, Sequence[float]]]) -> xr.DataArray:
    if not bbox:
        return da
    lat_rng = bbox.get("lat")
    lon_rng = bbox.get("lon")
    if lat_rng is not None and "latitude" in da.dims:
        lo, hi = sorted(float(v) for v in lat_rng)
        da = da.sel(latitude=slice(lo, hi))
    if lon_rng is not None and "longitude" in da.dims:
        from .indices import longitude_mask
        da = da.isel(longitude=longitude_mask(da.longitude.values, lon_rng))
    if any(da.sizes.get(n, 1) == 0 for n in ("latitude", "longitude")):
        raise ConfigError(f"bbox {bbox} does not intersect the input grid")
    return da


# ---------------------------------------------------------------------------
# year aggregation
# ---------------------------------------------------------------------------
def _season_year(time_values, months: Sequence[int]) -> np.ndarray:
    """Map each timestep to the year the *season* belongs to.

    Label a wrapping season by its ending year (DJF, NDJFM, etc.).
    """
    months = list(months)
    times = np.atleast_1d(time_values)
    if times.size and hasattr(times[0], "year"):  # cftime objects
        years = np.asarray([t.year for t in times], dtype=int)
        mon = np.asarray([t.month for t in times], dtype=int)
    else:  # numpy datetime64 (int() would give nanoseconds, not years)
        # Preserve dates outside the nanosecond datetime range (e.g. year 2300).
        dt = np.asarray(times).astype("datetime64[D]")
        years = dt.astype("datetime64[Y]").astype(int) + 1970
        mon = dt.astype("datetime64[M]").astype(int) % 12 + 1
    if months and len(months) < 12 and 1 in months and 12 in months:
        ordered = np.asarray(sorted(set(months)), dtype=int)
        gaps = np.diff(ordered)
        start = int(ordered[int(np.argmax(gaps)) + 1])
        years = np.where(mon >= start, years + 1, years)
    return years


def _frequency(da: xr.DataArray, frequency: str, *, declared_frequencies=()) -> str:
    """Use the same sampling rules as QC, retaining authoritative file hints."""
    from .quality import declared_frequency, resolve_frequency
    declared = [hint for hint in declared_frequencies if hint]
    hint = declared_frequency(da.attrs)
    if hint:
        declared.append(hint)
    if da.attrs.get("_cmip_sampling_frequency"):
        declared.append(da.attrs["_cmip_sampling_frequency"])
    resolved = resolve_frequency(da.time.values, requested=frequency, declared_frequencies=declared)
    if any(hint != resolved for hint in declared):
        raise ConfigError(f"aggregate.frequency={frequency!r} conflicts with source sampling metadata "
                          f"{sorted(set(declared))}")
    if (str(frequency).lower() == "auto" and not declared and resolved == "annual"
            and not np.issubdtype(da.time.dtype, np.number)):
        log.warning("one timestamp per year: treating input as preaggregated annual/seasonal data; "
                    "set aggregate.frequency explicitly if this is sparse monthly input")
    return resolved


def _time_parts(da: xr.DataArray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return tuple(np.asarray(getattr(da.time.dt, n).values, dtype=int)
                 for n in ("year", "month", "day", "days_in_month"))


def _complete_season(da: xr.DataArray, months: Sequence[int], frequency: str) -> bool:
    _, mon, day, month_length = _time_parts(da)
    if set(mon.tolist()) != set(months):
        return False
    for m in months:
        selected = mon == m
        if frequency == "monthly":
            if selected.sum() != 1:
                return False
        elif (selected.sum() != month_length[selected][0]
              or np.unique(day[selected]).size != month_length[selected][0]):
            return False
    return True


def _reduce(values: np.ndarray, stat: str, threshold: Optional[float]) -> np.ndarray:
    """Reduce along axis 0 (time) according to ``stat``."""
    stat = str(stat).lower()
    with np.errstate(invalid="ignore", divide="ignore"):
        if stat == "mean":
            return np.nanmean(values, axis=0)
        if stat == "sum":
            return np.where(np.isfinite(values).any(axis=0), np.nansum(values, axis=0), np.nan)
        if stat == "max":
            return np.nanmax(values, axis=0)
        if stat == "min":
            return np.nanmin(values, axis=0)
        if stat == "median":
            return np.nanmedian(values, axis=0)
        if stat == "std":
            return np.nanstd(values, axis=0, ddof=1)
        if stat == "max_abs":
            return np.nanmax(np.abs(values), axis=0)
        if stat == "mean_abs":
            return np.nanmean(np.abs(values), axis=0)
        if stat.startswith("p") and stat[1:].isdigit():
            q = float(stat[1:])
            return np.nanpercentile(values, q, axis=0)
        if stat.startswith("count_"):
            if threshold is None:
                raise ConfigError(f"aggregate.stat={stat} requires aggregate.threshold")
            op = stat.split("_", 1)[1]
            thr = float(threshold)
            table = {
                "gt": values > thr,
                "ge": values >= thr,
                "lt": values < thr,
                "le": values <= thr,
            }
            if op not in table:
                raise ConfigError(f"unsupported count operator in {stat!r}")
            return np.where(np.isfinite(values).any(axis=0), np.sum(table[op], axis=0).astype(float), np.nan)
    raise ConfigError(f"unsupported aggregate.stat {stat!r}")


def aggregate_yearly(
    da: xr.DataArray,
    months: Sequence[int],
    stat: str = "mean",
    threshold: Optional[float] = None,
    *,
    weight_by_days: bool = True,
    require_complete: bool = True,
    min_valid_fraction: float = 1.0,
    frequency: str = "auto",
) -> xr.DataArray:
    """Aggregate monthly/daily data to one value per (season-)year on (lat, lon)."""
    if "time" not in da.dims:
        raise ConfigError("aggregate_yearly needs a 'time' dimension")
    if not 0 <= float(min_valid_fraction) <= 1:
        raise ConfigError("aggregate.min_valid_fraction must be in [0, 1]")
    if set(da.dims) != {"time", "latitude", "longitude"}:
        raise ConfigError(f"expected time/latitude/longitude dimensions, found {da.dims}; select extra levels first")
    da = da.transpose("time", "latitude", "longitude")
    months = list(months) or list(range(1, 13))
    frequency = _frequency(da, frequency)
    if frequency == "annual":
        vals = np.asarray(da.time.values)
        labels = vals.astype(int) if np.issubdtype(vals.dtype, np.number) else np.asarray(da.time.dt.year.values)
        if np.unique(labels).size != labels.size:
            raise ConfigError("annual input must have exactly one sample per year")
        log.info("annual input is already aggregated; month/stat choices are not re-applied")
        return da.rename(time="year").assign_coords(year=labels).sortby("year")
    if months and sorted(months) != list(range(1, 13)):
        mon = da.time.dt.month
        da = da.isel(time=np.isin(mon.values, months))
    if frequency == "monthly" and da.time.size:
        month_keys = da.time.dt.year.values * 12 + da.time.dt.month.values
        if np.unique(month_keys).size != month_keys.size:
            raise ConfigError("monthly aggregation requires one sample per month; "
                              "resolve duplicate months or aggregate finer-frequency data explicitly")
    if frequency == "daily" and da.time.size:
        yr, mon, day, _ = _time_parts(da)
        day_keys = yr * 10000 + mon * 100 + day
        if np.unique(day_keys).size != day_keys.size:
            raise ConfigError("daily aggregation requires one sample per date; preaggregate subdaily data")
    if da.time.size == 0:
        return da.rename(time="year").assign_coords(year=np.array([], dtype=int))
    labels = _season_year(da.time.values, months)
    out_years = np.unique(labels)
    rows = []
    valid_years = []
    for yr in out_years:
        sel = da.isel(time=np.where(labels == yr)[0])
        if require_complete and not _complete_season(sel, months, frequency):
            log.warning("dropping incomplete %s season-year %s for months %s", frequency, yr, months)
            continue
        values = sel.values.astype(np.float64)
        if frequency == "monthly" and weight_by_days and stat in ("mean", "mean_abs"):
            weights = np.asarray(sel.time.dt.days_in_month.values, dtype=float)[:, None, None]
            valid = np.isfinite(values)
            work = np.abs(values) if stat == "mean_abs" else values
            den = np.sum(np.where(valid, weights, 0.0), axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                red = np.sum(np.where(valid, work * weights, 0.0), axis=0) / den
        else:
            red = _reduce(values, stat, threshold)
        coverage = np.mean(np.isfinite(values), axis=0)
        red = np.where(coverage >= float(min_valid_fraction), red, np.nan)
        rows.append(red)
        valid_years.append(yr)
    if not rows:
        return da.isel(time=slice(0, 0)).rename(time="year").assign_coords(year=np.array([], dtype=int))
    arr = np.stack(rows, axis=0)
    da_out = xr.DataArray(
        arr,
        dims=("year", "latitude", "longitude"),
        coords={
            "year": np.asarray(valid_years, dtype=int),
            "latitude": da.latitude.values if "latitude" in da.coords else np.arange(arr.shape[1]),
            "longitude": da.longitude.values if "longitude" in da.coords else np.arange(arr.shape[2]),
        },
        name=da.name or "field",
        attrs=dict(da.attrs),
    )
    if "plev" in da.coords:
        da_out = da_out.assign_coords(plev=da.plev)
    return da_out


# ---------------------------------------------------------------------------
# file opening
# ---------------------------------------------------------------------------
def _is_grib(path: str) -> bool:
    return path.lower().endswith(GRID_ENGINES)


def _calendar_name(da: xr.DataArray) -> str:
    if "time" not in da.coords:
        return ""
    calendar = da.time.encoding.get("calendar", da.time.attrs.get("calendar"))
    if calendar is None:
        times = np.asarray(da.time.values)
        if times.size and hasattr(times[0], "calendar"):
            calendar = times[0].calendar
        elif np.issubdtype(times.dtype, np.number):
            return "annual_numeric"
        else:
            calendar = "proleptic_gregorian"
    name = str(calendar).strip().lower()
    return {"gregorian": "standard", "365_day": "noleap", "366_day": "all_leap"}.get(name, name)


def _check_time(da: xr.DataArray, path: str) -> None:
    if "time" not in da.dims or da.time.dims != ("time",):
        raise ConfigError(f"{path}: expected a one-dimensional time coordinate")
    times = np.asarray(da.time.values)
    if not times.size:
        raise ConfigError(f"{path}: time coordinate is empty")
    if np.issubdtype(times.dtype, np.datetime64):
        invalid = np.isnat(times).any()
    elif np.issubdtype(times.dtype, np.number):
        invalid = not np.isfinite(times).all() or np.any(times != np.floor(times))
    else:
        invalid = any(t is None or not hasattr(t, "year") for t in times)
    if invalid:
        raise ConfigError(f"{path}: time contains invalid dates or non-integral annual labels")


def _source_metadata(da: xr.DataArray, ds: xr.Dataset, path: str,
                     expected: Optional[Dict[str, str]], units_in: Optional[str]) -> xr.DataArray:
    """Retain input identity after the Dataset owning its metadata is closed."""
    from .quality import canonical_units, declared_frequency
    attrs = dict(da.attrs)
    for name in ("source_id", "model_id", "experiment_id", "variant_label", "grid_label", "table_id", "frequency"):
        value = ds.attrs.get(name, da.attrs.get(name))
        if value is not None and str(value).strip():
            attrs[f"_cmip_{name}"] = str(value).strip()
    hint = declared_frequency(attrs)
    variable_hint = declared_frequency(da.attrs)
    dataset_hint = declared_frequency(ds.attrs)
    if len({value for value in (hint, variable_hint, dataset_hint) if value}) > 1:
        raise ConfigError(f"{path}: variable and dataset sampling metadata conflict "
                          f"({variable_hint!r} versus {dataset_hint or hint!r})")
    if hint or variable_hint:
        attrs["_cmip_sampling_frequency"] = hint or variable_hint
    attrs["_cmip_calendar"] = _calendar_name(da)
    attrs["_cmip_raw_units"] = str(da.attrs.get("units", "")).strip()
    for key, value in (expected or {}).items():
        found = attrs.get(f"_cmip_{key}")
        if found and value and found != str(value):
            raise ConfigError(f"{path}: {key}={found!r} conflicts with expected {value!r}")
    found_units = attrs["_cmip_raw_units"]
    if units_in and found_units and canonical_units(units_in) != canonical_units(found_units):
        raise ConfigError(f"{path}: units={found_units!r} conflict with units_in={units_in!r}; "
                          "set the correct source units and explicit conversion")
    da.attrs = attrs
    return da


@contextmanager
def _open_dataarray_context(path: str, variable: Optional[str] = None, plev: Any = None,
                            *, expected_metadata=None, units_in=None):
    """Own the file handle while callers inspect coordinates or load a subset.

    Longitude normalization is intentionally deferred until after selection:
    checking duplicate cyclic columns can otherwise read an entire source file.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if _is_grib(path):
        try:
            ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
        except Exception as exc:  # pragma: no cover - depends on optional cfgrib
            raise ConfigError(
                f"cannot open GRIB file {path}: {exc}. Install cfgrib/eccodes or convert to NetCDF first."
            )
    else:
        try:
            coder = getattr(getattr(xr, "coders", None), "CFDatetimeCoder", None)
            if coder is not None:
                ds = xr.open_dataset(path, decode_times=coder(use_cftime=None))
            else:
                ds = xr.open_dataset(path, decode_times=True)
        except Exception as exc:
            raise ConfigError(f"cannot decode {path}: {exc}; nonstandard calendars require cftime") from exc
    try:
        if variable:
            if variable not in ds.variables:
                raise ConfigError(f"variable {variable!r} not in {path}; available: {list(ds.data_vars)}")
            da = ds[variable]
        else:
            candidates = [name for name in ds.data_vars if ds[name].ndim >= 2]
            if len(candidates) != 1:
                raise ConfigError(f"choose variable explicitly in {path}; candidates are {list(ds.data_vars)}")
            da = ds[candidates[0]]
        da = da.copy(deep=False)
        da = normalize_dims(da)
        da = select_plev(da, plev)
        _check_time(da, path)
        if set(da.dims) != {"time", "latitude", "longitude"}:
            raise ConfigError(f"{path}: expected time/latitude/longitude dimensions, found {da.dims}; "
                              "select extra levels explicitly")
        da = _source_metadata(da, ds, path, expected_metadata, units_in)
        yield da
    finally:
        ds.close()


def open_dataarray(path: str, variable: Optional[str] = None, plev: Any = None, *,
                   bbox=None, time_indices=None, max_values=None,
                   expected_metadata=None, units_in=None) -> xr.DataArray:
    """Return a detached, normalized array, loading only the requested subset.

    Select time/space before loading values, enforce the optional value-count
    guard before data reads, and close the underlying file on every exit path.
    Callers needing only headers should use the internal context manager.
    """
    with _open_dataarray_context(path, variable, plev, expected_metadata=expected_metadata,
                                 units_in=units_in) as da:
        if time_indices is not None:
            da = da.isel(time=time_indices)
        da = select_bbox(da, bbox)
        if max_values is not None and da.size > float(max_values):
            raise ConfigError(f"{path}: selected field exceeds the memory guard "
                              f"({da.size} values > {float(max_values):g}); crop bbox or shorten the period")
        da = align_longitude(da, lon_360=True)
        return da.load()


def _check_headers(headers: Sequence[Dict[str, Any]]) -> None:
    """Reject incompatible chunks even when they occupy disjoint years."""
    from .quality import canonical_units
    known: Dict[str, Tuple[Any, str]] = {}
    for header in headers:
        path = header["path"]
        attrs = header["attrs"]
        for key in ("_cmip_source_id", "_cmip_model_id", "_cmip_experiment_id",
                    "_cmip_variant_label", "_cmip_grid_label", "_cmip_table_id",
                    "_cmip_calendar", "_cmip_raw_units", "_cmip_sampling_frequency"):
            value = attrs.get(key)
            if not value:
                continue
            if key == "_cmip_raw_units":
                value = canonical_units(value)
            if key in known and known[key][0] != value:
                raise ConfigError(f"incompatible source {key.removeprefix('_cmip_')}: "
                                  f"{known[key][1]} ({known[key][0]!r}) versus {path} ({value!r})")
            known.setdefault(key, (value, path))
        if headers and header is not headers[0]:
            for axis in ("latitude", "longitude"):
                if not np.array_equal(headers[0][axis], header[axis]):
                    raise ConfigError(f"source files have incompatible {axis} grids: "
                                      f"{headers[0]['path']} versus {path}")


def _field_headers(files, variable, plev, bbox, expected_metadata=None, units_in=None):
    headers = []
    for path in dict.fromkeys(files):
        with _open_dataarray_context(path, variable, plev, expected_metadata=expected_metadata,
                                     units_in=units_in) as da:
            da = select_bbox(da, bbox)
            headers.append({"path": path, "time": np.asarray(da.time.values).copy(),
                            "latitude": np.asarray(da.latitude.values).copy(),
                            "longitude": np.unique(np.mod(da.longitude.values, 360.0)),
                            "raw_spatial_size": da.sizes["latitude"] * da.sizes["longitude"],
                            "attrs": dict(da.attrs)})
    _check_headers(headers)
    return headers


def _headers_frequency(headers, requested="auto"):
    times = np.concatenate([header["time"] for header in headers])
    hints = [header["attrs"].get("_cmip_sampling_frequency") for header in headers]
    hints = [hint for hint in hints if hint]
    # Keep the provenance on this synthetic time-only object too. In particular,
    # sparse Amon records remain monthly; their December labels must shift to
    # the following DJF year instead of being mistaken for annual samples.
    attrs = {"_cmip_sampling_frequency": hints[0]} if hints else {}
    time_only = xr.DataArray(np.zeros(times.size), dims="time", coords={"time": times}, attrs=attrs)
    return _frequency(time_only, requested, declared_frequencies=hints)


def expand_paths(template: str, years: Optional[Sequence[int]] = None) -> List[str]:
    """Expand a glob template, optionally replacing ``{year}`` with each year."""
    paths: List[str] = []
    if "{year}" in template:
        years = list(range(1850, 2101)) if not years else list(years)
        for yr in years:
            paths.extend(sorted(glob.glob(template.format(year=yr))))
    else:
        paths = sorted(glob.glob(template))
    # de-duplicate while keeping order
    seen = set()
    out = []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _year_list(years: Optional[Sequence[int]]) -> Optional[List[int]]:
    """Turn a ``[start, end]`` period into the explicit list of years."""
    if years is None:
        return None
    if len(years) == 2 and int(years[1]) >= int(years[0]):
        return list(range(int(years[0]), int(years[1]) + 1))
    return [int(y) for y in years]


def _input_years(years: Sequence[int], months: Sequence[int]) -> List[int]:
    out = _year_list(years) or []
    if out and 1 in months and 12 in months and len(months) < 12:
        out = [min(out) - 1] + out
    return out


def _aggregate_options(spec):
    aggregate = spec.get("aggregate") or {}
    defaults = {"weight_by_days": True, "require_complete": True,
                "min_valid_fraction": 1.0, "frequency": "auto"}
    return {key: aggregate.get(key, default) for key, default in defaults.items()}


BUILTIN_PATTERNS = [
    "{root}/{experiment}/{model}/{variable}_Amon_{model}_{experiment}_{member}_*.nc",
    "{root}/{experiment}/{model}/{variable}/{variable}_Amon_{model}_{experiment}_{member}_*.nc",
]


def grid_patterns(spec: Dict[str, Any]) -> List[str]:
    """Candidate file patterns, in priority order.

    Explicit patterns are tried first; the built-in CMIP6 layouts (flat and
    per-variable sub-directory) are appended as fallbacks so that a single
    archive with a mixed layout still works.
    """
    user: List[str] = []
    if spec.get("pattern"):
        user.append(spec["pattern"])
    user.extend(p for p in (spec.get("patterns") or []) if p)
    if user:
        return list(dict.fromkeys(user))
    out = []
    for pattern in BUILTIN_PATTERNS:
        if pattern not in out:
            out.append(pattern)
    return out


def parse_placeholders(pattern: str, path: str) -> Dict[str, str]:
    """Extract ``{placeholder}`` values from a path using the pattern itself."""
    parts = re.split(r"(\{[A-Za-z_]+\}|\*)", pattern)
    regex, names = [], []
    for part in parts:
        if not part:
            continue
        if part == "*":
            regex.append(".*")
        elif re.fullmatch(r"\{[A-Za-z_]+\}", part):
            name = part[1:-1]
            if name in names:
                regex.append(f"(?P={name})")
            else:
                names.append(name)
                regex.append(f"(?P<{name}>[^/]*)")
        else:
            regex.append(re.escape(part))
    match = re.match("^" + "".join(regex) + "$", path)
    return match.groupdict() if match else {}


def _member_of(path: str, model: str, experiment: str, spec: Optional[Dict[str, Any]] = None) -> Optional[str]:
    """Extract the ensemble member from a file name.

    Tries the configured patterns first (works for any custom layout), then two
    common orders: ``<model>_<experiment>_<member>_...`` (CMIP6 DRS) and
    ``<model>_<member>_<experiment>``.
    """
    if spec is not None:
        for pattern in grid_patterns(spec):
            values = parse_placeholders(pattern.format(
                root=spec.get("root", ""), experiment=experiment, model=model,
                scenario=experiment, year="{year}",
                member="{member}", variable=spec.get("variable") or "*"), path)
            if values.get("member"):
                return values["member"]
    base = os.path.basename(path)
    for pattern in (
        rf"(?:^|_){re.escape(model)}_{re.escape(experiment)}_([A-Za-z0-9]+?)_",
        rf"(?:^|_){re.escape(model)}_([A-Za-z0-9]+?)_?{re.escape(experiment)}",
    ):
        match = re.search(pattern, base)
        if match:
            return match.group(1)
    return None


def _pick_member(files: Sequence[str], model: str, experiment: str,
                 spec: Optional[Dict[str, Any]] = None) -> List[str]:
    """When the member is a wildcard, keep the most common single member."""
    groups: Dict[str, List[str]] = {}
    for path in files:
        member = _member_of(path, model, experiment, spec) or "?"
        groups.setdefault(member, []).append(path)
    if len(groups) <= 1:
        return list(files)
    if "?" in groups:
        raise ConfigError("cannot identify ensemble members; use an explicit {member} pattern")
    best = sorted(groups)[0]
    log.info("  multiple members found for %s (%s); using %s", model, ", ".join(sorted(groups)), best)
    return groups[best]


def _files_for_member(spec: Dict[str, Any], experiment: str, model: Optional[str],
                      member: str, years: Optional[Sequence[int]]) -> List[str]:
    """Files of one experiment for one model and one member (``*`` = wildcard)."""
    yl = _year_list(years)
    for pattern in grid_patterns(spec):
        try:
            template = pattern.format(
                root=spec.get("root", ""),
                experiment=experiment,
                model=model or "*",
                member=member,
                variable=spec.get("variable") or "*",
                scenario=experiment,
                year="{year}",
            )
        except (KeyError, IndexError):
            continue
        paths = expand_paths(template, yl)
        if paths:
            return paths
    return []


def discover_members(spec: Dict[str, Any], experiment: str, model: str,
                     years: Optional[Sequence[int]] = None) -> List[str]:
    """All ensemble members available for one model/experiment."""
    files = _files_for_member(spec, experiment, model, "*", years)
    members = {_member_of(path, model, experiment, spec) for path in files}
    return sorted(m for m in members if m)


def _files_for_grid_source(spec: Dict[str, Any], experiment: str, model: Optional[str], years: Optional[Sequence[int]]):
    """First usable member (used when member_combine == 'first')."""
    if spec.get("kind") == "series":
        return []
    yl = _year_list(years)
    for member in list(spec.get("members") or ["*"]):
        paths = _files_for_member(spec, experiment, model, member, years)
        if paths:
            if member == "*":
                paths = _pick_member(paths, model or "", experiment, spec)
            return paths
    paths = []
    for template in spec.get("paths", []) or []:
        paths.extend(expand_paths(template, yl))
    return paths


def load_yearly_field(
    files: Sequence[str],
    variable: Optional[str],
    months: Sequence[int],
    stat: str,
    threshold: Optional[float],
    units_transform: Optional[Dict[str, Any]] = None,
    plev: Any = None,
    bbox: Optional[Dict[str, Sequence[float]]] = None,
    year_range: Optional[Sequence[int]] = None,
    *,
    weight_by_days: bool = True,
    require_complete: bool = True,
    min_valid_fraction: float = 1.0,
    frequency: str = "auto",
    units_in: Optional[str] = None,
    expected_metadata: Optional[Dict[str, str]] = None,
    max_values: Optional[int] = None,
) -> xr.DataArray:
    """Combine all samples for a season before reducing, regardless of file splits.

    Only one season-year is materialized at a time. Calendars and spatial grids
    must agree across the source files; conflicting duplicate timestamps fail.
    """
    if not files:
        raise ConfigError("no input files matched the configured pattern")
    if max_values is not None and (not np.isfinite(max_values) or max_values <= 0):
        raise ConfigError("max_values must be finite and positive")
    headers = _field_headers(files, variable, plev, bbox, expected_metadata, units_in)
    frequency = _headers_frequency(headers, frequency)
    labels_by_block = []
    for header in headers:
        times = header["time"]
        if frequency == "annual":
            t = np.asarray(times)
            labels = t.astype(int) if np.issubdtype(t.dtype, np.number) else xr.DataArray(t, dims="time").dt.year.values
        else:
            labels = _season_year(times, months)
        labels_by_block.append(np.asarray(labels, dtype=int))
    years_available = np.unique(np.concatenate(labels_by_block))
    if year_range is not None:
        years_available = years_available[(years_available >= int(year_range[0])) & (years_available <= int(year_range[1]))]
    yearly_blocks = []
    for year in years_available:
        selected_parts = []
        selected_size = 0
        for header, labels in zip(headers, labels_by_block):
            selected = labels == year
            if frequency != "annual":
                time_coord = xr.DataArray(header["time"], dims="time")
                selected &= np.isin(time_coord.dt.month.values, list(months) or list(range(1, 13)))
            if selected.any():
                selected_parts.append((header, np.flatnonzero(selected)))
                selected_size += int(selected.sum()) * header["raw_spatial_size"]
        if max_values is not None and selected_size > max_values:
            raise ConfigError(f"season-year {year} exceeds the memory guard ({selected_size} values); "
                              "crop bbox or explicitly increase runtime.max_input_values")
        parts = [open_dataarray(header["path"], variable, plev, bbox=bbox,
                                time_indices=selected, units_in=units_in, max_values=max_values,
                                expected_metadata=expected_metadata) for header, selected in selected_parts]
        if not parts:
            continue
        try:
            group = xr.concat(parts, dim="time", join="exact").sortby("time")
        except (ValueError, TypeError) as exc:
            raise ConfigError(f"source files have incompatible calendars or grids: {exc}") from exc
        group = _deduplicate_time(group)
        group = convert_units(group, units_transform)
        reduced = aggregate_yearly(
            group, months, stat, threshold, weight_by_days=weight_by_days,
            require_complete=require_complete, min_valid_fraction=min_valid_fraction,
            frequency=frequency,
        )
        if reduced.year.size:
            if max_values is not None and sum(b.size for b in yearly_blocks) + reduced.size > max_values:
                raise ConfigError("combined yearly output exceeds the memory guard; crop bbox, "
                                  "shorten the period or explicitly increase runtime.max_input_values")
            yearly_blocks.append(reduced)
    if not yearly_blocks:
        raise ConfigError("no complete season-years remain; inspect input coverage and aggregate settings")
    field = xr.concat(yearly_blocks, dim="year", join="exact").sortby("year")
    if year_range is not None:
        field = field.sel(year=slice(int(year_range[0]), int(year_range[1])))
    field.name = "field"
    return field


def _deduplicate_time(da: xr.DataArray) -> xr.DataArray:
    times = np.asarray(da.time.values)
    unique, first, counts = np.unique(times, return_index=True, return_counts=True)
    for value in unique[counts > 1]:
        copies = da.isel(time=np.flatnonzero(times == value)).values
        if not np.allclose(copies, copies[0], equal_nan=True, rtol=0, atol=0):
            raise ConfigError(f"conflicting duplicate timestamp {value}; select one version/member explicitly")
    if first.size != times.size:
        log.warning("removed %d identical duplicate timestamps", times.size - first.size)
        da = da.isel(time=np.sort(first))
    return da


# ---------------------------------------------------------------------------
# grids
# ---------------------------------------------------------------------------
def make_axis(spec: Any, default_step: float = 1.0) -> Optional[np.ndarray]:
    if spec is None:
        return None
    if isinstance(spec, dict):
        start = float(spec["start"])
        stop = float(spec["stop"])
        step = float(spec.get("step", default_step))
        n = int(round((stop - start) / step)) + 1
        return np.round(start + step * np.arange(n), 6)
    if isinstance(spec, (list, tuple)) and spec and isinstance(spec[0], (list, tuple)):
        start, stop = float(spec[0]), float(spec[1])
        step = float(spec[2]) if len(spec) > 2 else default_step
        n = int(round((stop - start) / step)) + 1
        return np.round(start + step * np.arange(n), 6)
    return np.asarray(spec, dtype=float)


def resolve_target_grid(cfg: Config, sample: Optional[xr.DataArray] = None) -> Dict[str, np.ndarray]:
    lat = make_axis(cfg.get_path("grid.lat"))
    lon = make_axis(cfg.get_path("grid.lon"))
    lon_360 = bool(cfg.get_path("grid.lon_360", True))
    if lat is None and sample is not None:
        lat = np.asarray(sample.latitude.values, dtype=float)
    if lon is None and sample is not None:
        lon = np.asarray(sample.longitude.values, dtype=float)
    if lat is None or lon is None:
        raise ConfigError("grid.lat/grid.lon are not configured and no sample field is available to copy the grid from")
    lon = np.asarray(lon, dtype=float)
    if lon_360:
        lon = np.mod(lon, 360.0)
    else:
        lon = ((lon + 180.0) % 360.0) - 180.0
    return {"latitude": np.sort(np.unique(np.asarray(lat, dtype=float))), "longitude": np.sort(np.unique(lon))}


def regrid(da: xr.DataArray, grid: Dict[str, np.ndarray], method: str = "linear") -> xr.DataArray:
    """Interpolate a field to the target lat/lon (longitude wrap handled)."""
    if "latitude" not in da.dims or "longitude" not in da.dims:
        return da
    da = align_longitude(da, lon_360=True)
    lat_t = np.asarray(grid["latitude"], dtype=float)
    lon_t = np.asarray(grid["longitude"], dtype=float)
    lon_query = np.mod(lon_t, 360.0)
    da = da.sortby("latitude")
    same = (
        da.latitude.size == lat_t.size
        and da.longitude.size == lon_t.size
        and np.allclose(da.latitude.values, lat_t)
        and np.allclose(da.longitude.values, lon_query)
    )
    if same:
        out = da
    else:
        lon = np.asarray(da.longitude.values)
        if lon.size > 1:
            gaps = np.diff(np.r_[lon, lon[0] + 360.0])
            typical = float(np.median(gaps))
            global_grid = np.allclose(gaps, typical, rtol=0.05, atol=1e-6)
            if global_grid:
                left = da.isel(longitude=slice(-1, None)).assign_coords(longitude=[lon[-1] - 360.0])
                right = da.isel(longitude=slice(0, 1)).assign_coords(longitude=[lon[0] + 360.0])
                da = xr.concat([left, da, right], dim="longitude")
            else:
                # Cut at the largest unsampled gap, e.g. [350, 0, 10] becomes
                # [350, 360, 370], never an interpolation across 340 degrees.
                start = (int(np.argmax(gaps)) + 1) % lon.size
                order = np.r_[np.arange(start, lon.size), np.arange(start)]
                unwrapped = np.mod(lon[order] - lon[start], 360.0) + lon[start]
                da = da.isel(longitude=order).assign_coords(longitude=unwrapped)
                lon_query = np.mod(lon_query - lon[start], 360.0) + lon[start]
        out = da
        for axis, target in (("latitude", lat_t), ("longitude", lon_query)):
            if out.sizes[axis] == 1:
                out = out.reindex({axis: target})
            else:
                out = out.interp({axis: target}, method=method)
    out = out.assign_coords(latitude=("latitude", lat_t), longitude=("longitude", lon_t))
    return out


# ---------------------------------------------------------------------------
# source loaders (grid + series)
# ---------------------------------------------------------------------------
def _pattern_to_regex(pattern: str) -> str:
    """Turn a path pattern with ``{placeholders}`` and ``*`` into a regex.

    ``{model}`` becomes a capturing group; every other placeholder matches one
    path component; ``*`` behaves like a shell wildcard.
    """
    parts = re.split(r"(\{[A-Za-z_]+\}|\*)", pattern)
    out = []
    for part in parts:
        if not part:
            continue
        if part == "{model}":
            out.append("(.+?)")
        elif re.fullmatch(r"\{[A-Za-z_]+\}", part):
            out.append("[^/]*")
        elif part == "*":
            out.append(".*")
        else:
            out.append(re.escape(part))
    return "^" + "".join(out) + "$"


class DataSource:
    """Loads reference and model fields/series, with an optional on-disk cache."""

    def __init__(self, cfg: Config, *, create_dirs: bool = True):
        self.cfg = cfg
        self.cache_dir = os.path.join(cfg.outdir, "cache")
        self.use_cache = bool(cfg.get_path("runtime.cache", True))
        if self.use_cache and create_dirs:
            os.makedirs(self.cache_dir, exist_ok=True)
        self.grid = None
        self.model_members: Dict[str, List[str]] = {}
        self.model_metadata: Dict[str, Dict[str, Any]] = {}
        self.actual_inputs: Dict[str, Dict[str, Any]] = {}
        if (create_dirs and cfg.get_path("grid.lat") is not None and cfg.get_path("grid.lon") is not None
                and (cfg.get_path("data.reference.kind", "grid") == "grid"
                     or cfg.get_path("data.models.kind", "grid") == "grid")):
            self.grid = resolve_target_grid(cfg)

    # -- cache ----------------------------------------------------------
    def _cache_path(self, key: Dict[str, Any]) -> str:
        digest = hashlib.md5(json.dumps(key, sort_keys=True, default=str).encode()).hexdigest()[:12]
        return os.path.join(self.cache_dir, f"{key['tag']}_{digest}.nc")

    def _cached(self, key: Dict[str, Any], builder):
        path = self._cache_path(key)
        fingerprint = hashlib.sha256(json.dumps(key, sort_keys=True, default=str).encode()).hexdigest()
        if self.use_cache and os.path.exists(path):
            try:
                with xr.open_dataarray(path) as opened:
                    if opened.attrs.get("_cmip_cache_key") != fingerprint:
                        raise ValueError("cache fingerprint does not match this input/configuration")
                    if set(opened.dims) != {"year", "latitude", "longitude"} or not opened.year.size:
                        raise ValueError("cache does not contain a nonempty yearly field")
                    years = np.asarray(opened.year.values)
                    if not np.isfinite(years).all() or np.any(np.diff(years) <= 0):
                        raise ValueError("cache years are invalid or not strictly increasing")
                    if opened.size > self.cfg.get_path("runtime.max_input_values", 30000000):
                        raise ConfigError("cached field exceeds the memory guard; crop bbox or "
                                          "explicitly increase runtime.max_input_values")
                    da = opened.load()
                log.info("cache hit %s", os.path.basename(path))
                return da
            except ConfigError:
                raise
            except Exception as exc:
                log.warning("discarding unreadable or incompatible cache %s; rebuilding: %s", path, exc)
        da = builder()
        if self.use_cache:
            try:
                da = da.copy(deep=False)
                da.attrs = dict(da.attrs, _cmip_cache_key=fingerprint)
                _atomic_netcdf(da, path)
            except Exception as exc:  # pragma: no cover
                log.warning("cache write failed; continuing with computed data: %s", exc)
        return da

    # -- grid -----------------------------------------------------------
    def reference_field(self) -> xr.DataArray:
        spec = self.cfg.get_path("data.reference", {})
        if spec.get("kind") == "series":
            raise ConfigError("reference_field() called for a series source")
        years = spec.get("years") or self.cfg.period("historical")
        months = self.cfg.months("reference")
        files = []
        input_years = _input_years(years, months)
        for template in spec.get("paths", []):
            files.extend(expand_paths(template, input_years))

        def build():
            log.info("reference: %d file(s)", len(files))
            field = load_yearly_field(
                files,
                spec.get("variable"),
                self.cfg.months("reference"),
                spec.get("aggregate", {}).get("stat", "mean"),
                spec.get("aggregate", {}).get("threshold"),
                spec.get("units_transform"),
                spec.get("plev"),
                spec.get("bbox"),
                years,
                units_in=spec.get("units_in"),
                max_values=self.cfg.get_path("runtime.max_input_values", 30000000),
                **_aggregate_options(spec),
            )
            if self.grid is None:
                self.grid = resolve_target_grid(self.cfg, sample=field)
            field = regrid(field, self.grid, self.cfg.get_path("grid.method", "linear"))
            return field

        field = self._cached(self._source_key("reference", files, spec, years), build)
        if self.grid is None:
            self.grid = resolve_target_grid(self.cfg, sample=field)
        return field

    def model_field(self, model: str, experiment: str, years: Sequence[int]) -> xr.DataArray:
        spec = self.cfg.get_path("data.models", {})
        if spec.get("kind") == "series":
            raise ConfigError("model_field() called for a series source")
        combine = str(spec.get("member_combine", "first")).lower()
        if combine not in ("first", "mean", "median"):
            raise ConfigError(f"data.models.member_combine must be first|mean|median, got {combine!r}")

        if combine == "first":
            requested = dict(spec)
            if model in self.model_members:
                requested["members"] = self.model_members[model]
            files = _files_for_grid_source(requested, experiment, model, _input_years(years, self.cfg.months("models")))
            if not files:
                raise ConfigError(f"no input files matched the configured pattern for {model}/{experiment}")
            actual = {_member_of(f, model, experiment, spec) for f in files}
            actual.discard(None)
            if len(actual) > 1:
                raise ConfigError(f"multiple members selected for {model}/{experiment}: {sorted(actual)}")
            if actual:
                self.model_members[model] = sorted(actual)
            log.info("%s / %s: %d file(s)", model, experiment, len(files))
            expected = {"source_id": model, "experiment_id": experiment}
            if actual:
                expected["variant_label"] = next(iter(actual))
            field = self._load_and_regrid(files, spec, years, cache_key=f"{model}_{experiment}",
                                         expected_metadata=expected)
            self._check_model_metadata(model, field)
            return field

        # ---- multi-member ensemble -----------------------------------
        configured = spec.get("members") or ["*"]
        input_years = _input_years(years, self.cfg.months("models"))
        if "*" in configured or configured == ["*"]:
            members = discover_members(spec, experiment, model, input_years)
        else:
            members = [m for m in configured if _files_for_member(spec, experiment, model, m, input_years)]
            missing = set(configured) - set(members)
            if missing:
                raise ConfigError(f"requested ensemble members missing for {model}/{experiment}: {sorted(missing)}")
        if not members:
            raise ConfigError(f"no ensemble members found for {model}/{experiment}")
        if model in self.model_members and members != self.model_members[model]:
            raise ConfigError(f"member sets differ across experiments for {model}; choose their common members explicitly")
        fields = []
        for member in members:
            files = _files_for_member(spec, experiment, model, member, input_years)
            if not files:
                log.warning("  %s/%s: member %s has no files, skipped", model, experiment, member)
                continue
            member_field = self._load_and_regrid(files, spec, years, cache_key=f"{model}_{experiment}_{member}",
                                                 expected_metadata={"source_id": model,
                                                                    "experiment_id": experiment,
                                                                    "variant_label": member})
            self._check_model_metadata(f"{model}/{member}", member_field)
            fields.append(member_field)
        if not fields:
            raise ConfigError(f"no usable member for {model}/{experiment}")
        self.model_members[model] = members
        if len(fields) == 1:
            return fields[0]
        common = xr.align(*fields, join="inner")
        if not common[0].year.size:
            raise ConfigError(f"no common years across members for {model}/{experiment}")
        stacked = xr.concat(common, dim="member")
        field = stacked.mean("member", skipna=False) if combine == "mean" else stacked.median("member", skipna=False)
        log.info("%s / %s: %d members combined (%s): %s", model, experiment, len(fields), combine, ", ".join(members))
        return field.rename("field")

    def _check_model_metadata(self, model, field):
        """Keep one model's historical/scenario identity stable across calls."""
        from .quality import canonical_units
        previous = self.model_metadata.setdefault(model, {})
        for key in ("_cmip_source_id", "_cmip_variant_label", "_cmip_grid_label",
                    "_cmip_calendar", "_cmip_raw_units"):
            value = field.attrs.get(key)
            if not value:
                continue
            if key == "_cmip_raw_units":
                value = canonical_units(value)
            if key in previous and previous[key] != value:
                raise ConfigError(f"{model}: historical/scenario {key.removeprefix('_cmip_')} differs "
                                  f"({previous[key]!r} versus {value!r})")
            previous.setdefault(key, value)

    def _load_and_regrid(self, files: Sequence[str], spec: Dict[str, Any], years: Sequence[int],
                         cache_key: str, expected_metadata=None) -> xr.DataArray:
        def build():
            field = load_yearly_field(
                files,
                spec.get("variable"),
                self.cfg.months("models"),
                spec.get("aggregate", {}).get("stat", "mean"),
                spec.get("aggregate", {}).get("threshold"),
                spec.get("units_transform"),
                spec.get("plev"),
                spec.get("bbox"),
                years,
                units_in=spec.get("units_in"),
                expected_metadata=expected_metadata,
                max_values=self.cfg.get_path("runtime.max_input_values", 30000000),
                **_aggregate_options(spec),
            )
            if self.grid is None:
                self.grid = resolve_target_grid(self.cfg, sample=field)
            return regrid(field, self.grid, self.cfg.get_path("grid.method", "linear"))

        return self._cached(
            self._source_key(cache_key, files, spec, years), build
        )

    def _source_key(self, tag, files, spec, years):
        from . import __version__
        fingerprints = []
        for path in sorted(set(files)):
            status = os.stat(path)
            fingerprints.append([os.path.abspath(path), status.st_size, status.st_mtime_ns])
        self.actual_inputs[tag] = {"label": tag, "files_metadata": [
            {"path": path, "size_bytes": size, "mtime_ns": modified}
            for path, size, modified in fingerprints], "years": list(years)}
        effective_grid = ({k: np.asarray(v).tolist() for k, v in self.grid.items()}
                          if self.grid is not None and tag != "reference" else None)
        return {"tag": tag, "loader_version": 3, "cmipkit_version": __version__, "years": list(years),
                "spec": spec, "files": fingerprints, "grid": self.cfg.get_path("grid"),
                "effective_grid": effective_grid}

    # -- discovery ------------------------------------------------------
    def discover_models(self) -> List[str]:
        spec = self.cfg.get_path("data.models", {})
        configured = spec.get("models", "all")
        exclude = set(spec.get("exclude", []) or []) | set(self.cfg.get_path("ensemble.exclude_models", []) or [])
        if isinstance(configured, (list, tuple)) and configured and configured != ["all"]:
            return [n for n in configured if n not in exclude]
        if spec.get("kind") == "series":
            root = spec.get("series_dir")
            if not root or not os.path.isdir(root):
                raise ConfigError(f"data.models.series_dir does not exist: {root}")
            experiments = [spec.get("historical_experiment", "historical")] + list(spec.get("scenarios") or [])
            names = []
            for entry in sorted(os.listdir(root)):
                path = os.path.join(root, entry)
                if os.path.isdir(path):
                    historical = spec.get("historical_experiment", "historical")
                    if any(os.path.isfile(os.path.join(path, historical + ext)) for ext in (".csv", ".npz")):
                        names.append(entry)
                elif entry.lower().endswith((".csv", ".npz")):
                    stem = os.path.splitext(entry)[0]
                    for experiment in experiments:
                        suffix = "_" + experiment
                        if stem.endswith(suffix):
                            names.append(stem[:-len(suffix)])
                            break
            return sorted(set(names) - exclude)
        if isinstance(configured, (list, tuple)) and configured and configured != ["all"]:
            names = list(configured)
        else:
            names = self._discover_from_pattern(spec)
        exclude = set(spec.get("exclude", []) or []) | set(self.cfg.get_path("ensemble.exclude_models", []) or [])
        return [n for n in names if n not in exclude]

    def _discover_from_pattern(self, spec: Dict[str, Any]) -> List[str]:
        """Infer model names from the configured file pattern.

        Works whether ``{model}`` is a directory level (CMIP6 archive layout) or
        part of the file name (flat layout), by turning the pattern into a
        regular expression and capturing the ``{model}`` group.
        """
        fmt = {
            "root": spec.get("root", ""),
            "experiment": spec.get("historical_experiment", "historical"),
            "member": (spec.get("members") or ["*"])[0],
            "variable": spec.get("variable") or "*",
            "scenario": spec.get("historical_experiment", "historical"),
            "year": "*",
        }
        names: List[str] = []
        for pattern in grid_patterns(spec):
            if "{model}" not in pattern:
                continue
            glob_pat = pattern.format(model="*", **fmt)
            files = sorted(glob.glob(glob_pat))
            if not files:
                continue
            capture = pattern.format(model="{model}", **fmt)
            for path in files:
                values = parse_placeholders(capture, path)
                if values.get("model"):
                    names.append(values["model"])
        return sorted(set(names))

    # -- series ---------------------------------------------------------
    def reference_series(self) -> Tuple[np.ndarray, np.ndarray]:
        spec = self.cfg.get_path("data.reference", {})
        path = spec.get("series_file")
        self._source_key("reference", [path], spec, self.cfg.period("historical"))
        years, values = read_series_file(path, spec.get("series_columns"))
        values = convert_units(xr.DataArray(values), spec.get("units_transform")).values
        lo, hi = self.cfg.period("historical")
        mask = (years >= lo) & (years <= hi)
        return years[mask], values[mask]

    def model_series(self, model: str, experiment: str) -> Tuple[np.ndarray, np.ndarray]:
        spec = self.cfg.get_path("data.models", {})
        root = spec.get("series_dir")
        cands = [
            os.path.join(root, model, f"{experiment}.csv"),
            os.path.join(root, model, f"{experiment}.npz"),
            os.path.join(root, f"{model}_{experiment}.csv"),
            os.path.join(root, f"{model}_{experiment}.npz"),
        ]
        for path in cands:
            if os.path.exists(path):
                years, values = read_series_file(path, spec.get("series_columns"))
                self._source_key(f"{model}/{experiment}", [path], spec, [int(years[0]), int(years[-1])])
                values = convert_units(xr.DataArray(values), spec.get("units_transform")).values
                return years, values
        raise FileNotFoundError(f"no series file for {model}/{experiment} under {root}")


def read_series_file(path: str, columns: Optional[Dict[str, str]] = None) -> Tuple[np.ndarray, np.ndarray]:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if path.lower().endswith(".npz"):
        with np.load(path, allow_pickle=False) as data:
            if "years" not in data or "values" not in data:
                raise ConfigError(f"{path} must contain numeric 'years' and 'values' arrays")
            return _validate_series(data["years"], data["values"], path)
    import pandas as pd

    df = pd.read_csv(path)
    cols = columns or {"year": "year", "value": "value"}
    if cols["year"] not in df.columns or cols["value"] not in df.columns:
        raise ConfigError(f"{path} must have columns {cols['year']!r} and {cols['value']!r}; found {list(df.columns)}")
    return _validate_series(df[cols["year"]].to_numpy(dtype=float),
                            df[cols["value"]].to_numpy(dtype=float), path)


def _validate_series(years, values, path):
    years, values = np.asarray(years, dtype=float), np.asarray(values, dtype=float)
    if years.ndim != 1 or values.ndim != 1 or years.shape != values.shape or not years.size:
        raise ConfigError(f"{path}: years/values must be nonempty 1-D arrays of equal length")
    if not np.isfinite(years).all() or not np.equal(years, np.floor(years)).all():
        raise ConfigError(f"{path}: year coordinates must be finite integers")
    if np.unique(years).size != years.size:
        raise ConfigError(f"{path}: duplicate year coordinates; aggregate or resolve duplicates explicitly")
    if np.isinf(values).any():
        raise ConfigError(f"{path}: series values cannot contain infinity")
    order = np.argsort(years)
    return years[order].astype(int), values[order]


# ---------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------
def _atomic_netcdf(da: xr.DataArray, path: str, **kwargs) -> None:
    directory = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(directory, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + os.path.basename(path) + ".", suffix=".tmp.nc", dir=directory)
    os.close(fd)
    try:
        da.to_netcdf(temporary, **kwargs)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def save_field(da: xr.DataArray, path: str, name: str = "field") -> str:
    """Write a complete new field atomically, preserving an old file on failure."""
    import importlib.util
    out = da.rename(name) if da.name != name else da
    if importlib.util.find_spec("netCDF4") is not None:
        _atomic_netcdf(out, path, encoding={name: {"zlib": True, "complevel": 4}}, engine="netcdf4")
    else:
        _atomic_netcdf(out, path)
    return path


def save_npz(path: str, **arrays) -> str:
    directory = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(directory, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + os.path.basename(path) + ".", suffix=".tmp.npz", dir=directory)
    try:
        with os.fdopen(fd, "wb") as handle:
            np.savez_compressed(handle, **arrays)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return path

# ---------------------------------------------------------------------------
# daily fields (used by the daily-extremes stage)
# ---------------------------------------------------------------------------
def load_daily_field(
    spec: Dict[str, Any],
    years: Sequence[int],
    months: Sequence[int],
    experiment: Optional[str] = None,
    model: Optional[str] = None,
    bbox: Optional[Dict[str, Sequence[float]]] = None,
    max_values: float = 3.0e7,
) -> xr.DataArray:
    """Load the raw daily field for a season (no yearly aggregation).

    ``experiment=None`` means the observational reference (``spec['paths']``),
    otherwise the model archive is searched for ``model``/``experiment``.

    A guard rail stops a run that would hold more than ``max_values`` floats in
    memory — set ``data.*.bbox`` to the study region for continental daily data.
    """
    if spec.get("kind") == "series":
        raise ConfigError("daily extremes need gridded input; series input has no sub-daily values")
    yl = _input_years(years, months)
    if experiment is None:
        files: List[str] = []
        for template in spec.get("paths", []) or []:
            files.extend(expand_paths(template, yl))
    else:
        files = _files_for_grid_source(spec, experiment, model, yl)
    if not files:
        raise ConfigError("no daily input files matched the configured pattern")
    if not np.isfinite(max_values) or max_values <= 0:
        raise ConfigError("daily max_values must be finite and positive")
    expected = {"source_id": model, "experiment_id": experiment} if experiment else {}
    members = spec.get("members") or []
    if len(members) == 1 and members[0] != "*" and experiment:
        expected["variant_label"] = members[0]
    headers = _field_headers(files, spec.get("variable"), spec.get("plev"),
                             bbox or spec.get("bbox"), expected, spec.get("units_in"))
    if _headers_frequency(headers) != "daily":
        raise ConfigError("daily extremes require daily input, not monthly/annual samples")
    selected_headers = []
    total_values = 0
    for header in headers:
        times = xr.DataArray(header["time"], dims="time")
        if np.issubdtype(times.dtype, np.number):
            raise ConfigError("daily extremes require decoded daily dates, not annual labels")
        months_arr = times.dt.month.values
        season_years = _season_year(header["time"], months)
        sel = (np.isin(months_arr, list(months))
               & (season_years >= int(years[0])) & (season_years <= int(years[1])))
        if not sel.any():
            continue
        indices = np.flatnonzero(sel)
        total_values += indices.size * header["raw_spatial_size"]
        selected_headers.append((header, indices))
        if total_values > max_values:
            raise ConfigError("daily field exceeds the memory guard; crop data.*.bbox or shorten the period")
    blocks = [open_dataarray(header["path"], spec.get("variable"), spec.get("plev"),
                              bbox=bbox or spec.get("bbox"), time_indices=selected,
                              max_values=max_values, expected_metadata=expected,
                              units_in=spec.get("units_in")) for header, selected in selected_headers]
    if not blocks:
        raise ConfigError("no daily timesteps inside the requested months")
    try:
        field = xr.concat(blocks, dim="time", join="exact")
    except (ValueError, TypeError) as exc:
        raise ConfigError(f"daily files have incompatible grids or calendars: {exc}") from exc
    field = field.sortby("time")
    field = _deduplicate_time(field)
    if _frequency(field, "daily") != "daily":
        raise ConfigError("daily extremes require daily input, not monthly/annual samples")
    parts = _time_parts(field)
    day_keys = parts[0] * 10000 + parts[1] * 100 + parts[2]
    if np.unique(day_keys).size != day_keys.size:
        raise ConfigError("daily input contains multiple samples per day; aggregate subdaily data explicitly")
    field = convert_units(field, spec.get("units_transform"))
    field = field.transpose("time", "latitude", "longitude")
    n_values = field.sizes["time"] * field.sizes.get("latitude", 1) * field.sizes.get("longitude", 1)
    if n_values > max_values:
        raise ConfigError(
            f"daily field would hold {n_values/1e6:.1f}M values; "
            f"set data.*.bbox to the study region or shorten the period"
        )
    field.name = "field"
    field.attrs["source_files"] = json.dumps(list(dict.fromkeys(files)))
    log.info("daily field: %d timesteps (%d-%d), %dx%d grid",
             field.sizes["time"], int(years[0]), int(years[1]),
             field.sizes.get("latitude", 0), field.sizes.get("longitude", 0))
    return field


def daily_index_series(field: xr.DataArray, spec: Dict[str, Any]) -> np.ndarray:
    """Reduce a daily field to a 1-D daily index (region mean / region diff)."""
    from . import indices

    itype = str(spec.get("type", "region_mean"))
    work = field
    if spec.get("remove_zonal_mean"):
        work = indices.remove_zonal_mean(work)
    if itype == "region_mean":
        region = spec.get("region")
        if not region:
            raise ConfigError("daily extremes need index.region")
        return indices.region_mean_series(work, region)
    if itype == "region_diff":
        regions = spec.get("regions") or []
        if len(regions) != 2:
            raise ConfigError("index.regions must contain [reference, target] for region_diff")
        zonal = indices.remove_zonal_mean(work)
        return indices.region_mean_series(zonal, regions[1]) - indices.region_mean_series(zonal, regions[0])
    raise ConfigError(f"daily extremes do not support daily.index={itype!r}; use region_mean or region_diff. Daily anomaly climatology must be defined separately.")
