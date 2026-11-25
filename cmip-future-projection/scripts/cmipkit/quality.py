"""Read-only, bounded-memory checks of the input actually requested by a run.

Header scans read coordinates and metadata, never the field values. Full scans
decode CF masks/scales and inspect all selected values in bounded tiles. These
checks do not fill gaps, clip values, or certify every CF/CMIP convention.
"""
from __future__ import annotations

import calendar as pycalendar
import hashlib
import itertools
import re
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np
import xarray as xr

from .config import ConfigError


def canonical_units(value: Any) -> Optional[str]:
    """Recognise spelling aliases only, never perform an implicit conversion."""
    if value is None or not str(value).strip():
        return None
    text = str(value).strip().lower().replace("−", "-").replace("**", "^")
    compact = re.sub(r"[\s_]", "", text).replace("^", "")
    aliases = {
        "k": "K", "kelvin": "K", "kelvins": "K", "degk": "K",
        "degreecelsius": "degC", "degreescelsius": "degC", "celsius": "degC",
        "degc": "degC", "°c": "degC", "degreescelcius": "degC",
        "pa": "Pa", "pascal": "Pa", "pascals": "Pa",
        "hpa": "hPa", "hectopascal": "hPa", "mbar": "hPa", "millibar": "hPa",
        "mm/day": "mm/day", "mm/d": "mm/day", "mmday-1": "mm/day", "mmd-1": "mm/day",
        "kgm-2s-1": "kg m-2 s-1", "kg/m2/s": "kg m-2 s-1", "kgm-2/s": "kg m-2 s-1",
        "m/s": "m s-1", "ms-1": "m s-1", "wm-2": "W m-2", "w/m2": "W m-2",
        "%": "%", "percent": "%", "percentage": "%",
        "1": "1", "dimensionless": "1", "unitless": "1",
    }
    return aliases.get(compact, text)


def canonical_calendar(value: Any) -> str:
    name = str(value or "standard").strip().lower()
    return {"gregorian": "standard", "365_day": "noleap", "366_day": "all_leap"}.get(name, name)


def declared_frequency(metadata: Dict[str, Any]) -> Optional[str]:
    """Read recognised raw frequency/table labels without guessing other tables."""
    value = str(metadata.get("frequency") or metadata.get("_cmip_frequency") or "").strip().lower()
    table = str(metadata.get("table_id") or metadata.get("_cmip_table_id") or "").strip().lower()
    aliases = {"day": "daily", "daily": "daily", "mon": "monthly", "monthly": "monthly",
               "yr": "annual", "year": "annual", "annual": "annual", "yearly": "annual"}
    frequency = aliases.get(value)
    table_frequency = "monthly" if table.endswith("mon") else "daily" if table.endswith("day") else None
    if frequency and table_frequency and frequency != table_frequency:
        raise ConfigError(f"frequency={value!r} conflicts with table_id={table!r}")
    return frequency or table_frequency


def resolve_frequency(time_values, requested="auto", *, declared_frequencies=()) -> str:
    """Share frequency resolution between QC and the runtime loader.

    Explicit requests retain their meaning. For auto, recognised metadata takes
    precedence over time density: sparse Amon input remains monthly input.
    """
    requested = str(requested).lower()
    if requested not in ("auto", "monthly", "daily", "annual"):
        raise ConfigError("aggregate.frequency must be auto|monthly|daily|annual")
    if requested != "auto":
        return requested
    declarations = {str(value) for value in declared_frequencies if value}
    if not declarations.issubset({"monthly", "daily", "annual"}):
        raise ConfigError(f"unrecognised sampling frequency declarations: {sorted(declarations)}")
    if len(declarations) > 1:
        raise ConfigError(f"files declare different sampling frequencies: {sorted(declarations)}")
    values = np.asarray(time_values)
    if not values.size:
        raise ConfigError("cannot infer frequency from an empty time coordinate")
    if np.issubdtype(values.dtype, np.number):
        if declarations - {"annual"}:
            raise ConfigError("numeric annual labels disagree with the declared file frequency")
        return "annual"
    if declarations:
        return next(iter(declarations))
    # Exact duplicate stamps must not masquerade as finer sampling.
    try:
        values = np.unique(values)
        coordinate = xr.DataArray(values, dims="time", coords={"time": values})
        years = np.asarray(coordinate.time.dt.year.values, dtype=int)
        months = np.asarray(coordinate.time.dt.month.values, dtype=int)
    except (TypeError, ValueError, AttributeError) as exc:
        raise ConfigError("cannot resolve frequency from incompatible/undecoded dates") from exc
    if np.unique(years * 12 + months).size < values.size:
        return "daily"
    if values.size > 1 and np.unique(years).size == values.size:
        return "annual"
    return "monthly"


def _finding(report, severity, code, message, **details):
    report["findings"].append(dict(severity=severity, code=code, message=message, **details))


def _open(path):
    suffix = Path(path).suffix.lower()
    if suffix in (".tif", ".tiff"):
        raise ConfigError("GeoTIFF needs an explicit calendar/unit adapter before this workflow")
    kwargs = {}
    if suffix in (".grib", ".grb", ".grib2", ".grb2"):
        kwargs.update(engine="cfgrib", backend_kwargs={"indexpath": ""})
    return xr.open_dataset(path, **kwargs)


def _prepare(ds, spec):
    from . import io
    variable = spec.get("variable")
    if variable and variable not in ds.data_vars:
        raise ConfigError(f"variable {variable!r} absent; available: {list(ds.data_vars)}")
    if not variable:
        candidates = [name for name in ds.data_vars if ds[name].ndim >= 2]
        if len(candidates) != 1:
            raise ConfigError(f"ambiguous data variable; choose explicitly from {list(ds.data_vars)}")
        variable = candidates[0]
    da = io.select_plev(io.normalize_dims(ds[variable]), spec.get("plev"))
    for dim in list(da.dims):
        if dim not in ("time", "latitude", "longitude") and da.sizes[dim] == 1:
            da = da.isel({dim: 0}, drop=True)
    if set(da.dims) != {"time", "latitude", "longitude"}:
        raise ConfigError(f"expected time and one-dimensional latitude/longitude; found {dict(da.sizes)}")
    if da.time.dims != ("time",):
        raise ConfigError("time must be one-dimensional")
    # Sorting/normalising coordinates is lazy. Do not call align_longitude here:
    # its cyclic-column equality test reads data, even during a header-only scan.
    lon = np.mod(np.asarray(da.longitude.values, dtype=float), 360.0)
    da = da.assign_coords(longitude=("longitude", lon)).sortby("longitude")
    da = io.select_bbox(da, spec.get("bbox"))
    if any(da.sizes[d] == 0 for d in da.dims):
        raise ConfigError("empty time or spatial dimension")
    return da.transpose("time", "latitude", "longitude"), variable


def _stamp(value):
    if isinstance(value, np.datetime64):
        return np.datetime_as_string(value, unit="us")
    if hasattr(value, "year"):
        return (f"{value.year:04d}-{value.month:02d}-{value.day:02d}T"
                f"{getattr(value, 'hour', 0):02d}:{getattr(value, 'minute', 0):02d}:"
                f"{getattr(value, 'second', 0):02d}.{getattr(value, 'microsecond', 0):06d}")
    return str(int(value))


def _time_info(da, spec):
    values = np.asarray(da.time.values)
    if not values.size:
        raise ConfigError("empty time coordinate")
    numeric = np.issubdtype(values.dtype, np.number)
    if numeric:
        if not np.isfinite(values).all() or np.any(values != np.floor(values)):
            raise ConfigError("numeric annual labels must be finite integers")
        if spec.get("aggregate", {}).get("frequency", "auto") not in ("auto", "annual"):
            raise ConfigError("numeric time is not decoded CF dates; check time units/calendar")
        if np.any((values < 1) | (values > 9999)):
            raise ConfigError("numeric time does not look like annual year labels; decode CF time explicitly")
        parts = [(int(y), 1, 1) for y in values]
        calendar = "annual_numeric"
    else:
        if np.issubdtype(values.dtype, np.datetime64) and np.isnat(values).any():
            raise ConfigError("time contains missing/invalid dates")
        try:
            years, months, days = (np.asarray(getattr(da.time.dt, x).values, dtype=int)
                                    for x in ("year", "month", "day"))
        except (AttributeError, TypeError, ValueError) as exc:
            raise ConfigError("time is not a decoded calendar coordinate") from exc
        parts = list(zip(years.tolist(), months.tolist(), days.tolist()))
        calendar = da.time.encoding.get("calendar", da.time.attrs.get("calendar"))
        if not calendar:
            calendar = getattr(values[0], "calendar", "proleptic_gregorian")
        calendar = canonical_calendar(calendar)
    return values, parts, calendar


def _wrap_start(months):
    months = sorted(set(months))
    if len(months) < 12 and 1 in months and 12 in months:
        return int(months[int(np.argmax(np.diff(months))) + 1])
    return None


def _selection(parts, years, months, annual=False):
    start = _wrap_start(months)
    lo, hi = map(int, years)
    return np.asarray([
        lo <= y <= hi if annual else
        (m in months and lo <= y + int(start is not None and m >= start) <= hi)
        for y, m, _ in parts], dtype=bool)


def _month_days(year, month, calendar):
    if calendar == "360_day":
        return 30
    if calendar == "noleap":
        return 28 if month == 2 else pycalendar.monthrange(2001, month)[1]
    if calendar == "all_leap":
        return 29 if month == 2 else pycalendar.monthrange(2000, month)[1]
    if calendar == "julian":
        return 29 if month == 2 and year % 4 == 0 else pycalendar.monthrange(2001, month)[1]
    return pycalendar.monthrange(year, month)[1]


def _coverage(report, records, years, months, frequency, calendar, policy):
    unique_stamps = {stamp: part for stamp, part in records}
    if frequency == "annual":
        keys = [part[0] for part in unique_stamps.values()]
    elif frequency == "monthly":
        keys = [part[:2] for part in unique_stamps.values()]
    else:
        keys = [part for part in unique_stamps.values()]
    observed = set(keys)
    if len(keys) != len(observed):
        _finding(report, "error", "multiple_samples_per_period",
                 f"distinct timestamps map to the same {frequency} sample; aggregate sub-frequency data explicitly",
                 count=len(keys) - len(observed))
    start = _wrap_start(months)
    expected, missing, examples = 0, 0, []
    for season_year in range(int(years[0]), int(years[1]) + 1):
        if frequency == "annual":
            candidates = [season_year]
        else:
            candidates = []
            for month in months:
                year = season_year - int(start is not None and month >= start)
                if frequency == "monthly":
                    candidates.append((year, month))
                else:
                    candidates.extend((year, month, day) for day in range(1, _month_days(year, month, calendar) + 1))
        for key in candidates:
            expected += 1
            if key not in observed:
                missing += 1
                if len(examples) < 5:
                    examples.append(str(key))
    report["coverage"] = {"frequency": frequency, "calendar": calendar,
                          "requested_years": list(years), "months": list(months),
                          "expected_samples": expected, "unique_samples": len(observed),
                          "missing_samples": missing, "missing_examples": examples}
    if missing:
        _finding(report, "error" if policy == "error" else "warning", "missing_time_samples",
                 f"{missing} requested {frequency} samples are absent; gaps are not filled",
                 count=missing, examples=examples)


def _new_stats():
    return dict(count=0, finite_count=0, nan_count=0, inf_count=0, raw_inf_count=0, min=None, max=None,
                all_missing_timesteps=0, all_missing_examples=[], below_min=0, above_max=0,
                max_chunk_values_observed=0)


def _scan_values(da, stamps, quality, transform, hashes):
    """Scan all selected data with both temporal and spatial tiling."""
    stats = _new_stats()
    limit = int(quality.get("max_chunk_values", 2000000))
    tstep = min(int(quality.get("chunk_size", 32)), limit)
    spatial_limit = max(1, limit // tstep)
    ntime, nlat, nlon = (da.sizes[n] for n in ("time", "latitude", "longitude"))
    lonstep = min(nlon, spatial_limit)
    latstep = min(nlat, max(1, spatial_limit // lonstep))
    scale, offset = float(transform.get("scale", 1)), float(transform.get("offset", 0))
    lower, upper = quality.get("valid_min"), quality.get("valid_max")
    conflicts = []
    for t0 in range(0, ntime, tstep):
        t1 = min(ntime, t0 + tstep)
        digests = [hashlib.sha256() for _ in range(t1 - t0)]
        finite_per_time = np.zeros(t1 - t0, dtype=np.int64)
        for lat0, lon0 in itertools.product(range(0, nlat, latstep), range(0, nlon, lonstep)):
            block = da.isel(time=slice(t0, t1), latitude=slice(lat0, lat0 + latstep),
                            longitude=slice(lon0, lon0 + lonstep))
            values = np.asarray(block.values, dtype=np.float64)
            raw_values = values
            stats["max_chunk_values_observed"] = max(stats["max_chunk_values_observed"], int(values.size))
            stats["raw_inf_count"] += int(np.isinf(values).sum())
            with np.errstate(over="ignore", invalid="ignore"):
                values = values * scale + offset
            finite = np.isfinite(values)
            finite_per_time += finite.sum(axis=(1, 2))
            stats["count"] += int(values.size)
            stats["finite_count"] += int(finite.sum())
            stats["nan_count"] += int(np.isnan(values).sum())
            stats["inf_count"] += int(np.isinf(values).sum())
            if finite.any():
                lo, hi = float(values[finite].min()), float(values[finite].max())
                stats["min"] = lo if stats["min"] is None else min(stats["min"], lo)
                stats["max"] = hi if stats["max"] is None else max(stats["max"], hi)
            if lower is not None:
                stats["below_min"] += int((finite & (values < lower)).sum())
            if upper is not None:
                stats["above_max"] += int((finite & (values > upper)).sum())
            # Normalise NaN payloads and signed zero for exact equality matching.
            canonical = np.array(raw_values, copy=True)
            canonical[np.isnan(canonical)] = np.nan
            canonical[canonical == 0] = 0
            for j, row in enumerate(canonical):
                digests[j].update(np.ascontiguousarray(row).tobytes())
        for j, digest in enumerate(digests):
            stamp = stamps[t0 + j]
            digest = digest.hexdigest()
            if stamp in hashes and hashes[stamp] != digest:
                conflicts.append(stamp)
            hashes[stamp] = digest
            if finite_per_time[j] == 0:
                stats["all_missing_timesteps"] += 1
                if len(stats["all_missing_examples"]) < 5:
                    stats["all_missing_examples"].append(stamp)
    return stats, conflicts


def _merge_stats(target, source):
    for key in ("count", "finite_count", "nan_count", "inf_count", "raw_inf_count", "all_missing_timesteps", "below_min", "above_max"):
        target[key] += source[key]
    for key, op in (("min", min), ("max", max)):
        if source[key] is not None:
            target[key] = source[key] if target[key] is None else op(target[key], source[key])
    target["max_chunk_values_observed"] = max(target["max_chunk_values_observed"], source["max_chunk_values_observed"])
    target["all_missing_examples"] = (target["all_missing_examples"] + source["all_missing_examples"])[:5]


def _finish_values(report, quality):
    stats = report["values"]
    count = stats["count"]
    stats["missing_fraction"] = ((count - stats["finite_count"]) / count) if count else None
    stats["missing_fraction_denominator"] = "all selected raw grid values, including static masked cells and overlapping file records"
    if count and not stats["finite_count"]:
        _finding(report, "error", "all_values_missing", "all selected dataset values are missing/non-finite")
    if stats["inf_count"] or stats["raw_inf_count"]:
        _finding(report, "error", "infinite_values", "infinite values found in decoded input or after configured unit conversion",
                 raw_count=stats["raw_inf_count"], transformed_count=stats["inf_count"])
    threshold = quality.get("max_missing_fraction")
    if count and threshold is not None and stats["missing_fraction"] > threshold:
        _finding(report, "error", "missing_fraction_exceeded",
                 f"missing fraction {stats['missing_fraction']:.6g} exceeds {threshold}; denominator includes static masks")
    if stats["below_min"] or stats["above_max"]:
        _finding(report, "error", "value_range_exceeded", "values outside configured bounds in transformed output units; no clipping applied",
                 below_min=stats["below_min"], above_max=stats["above_max"])
    if stats["all_missing_timesteps"]:
        _finding(report, "warning", "all_missing_timesteps", "some selected timestamps have no finite cells",
                 count=stats["all_missing_timesteps"], examples=stats["all_missing_examples"])


def audit_grid(files: Sequence[str], spec: Dict[str, Any], *, label: str,
               years: Sequence[int], months: Sequence[int], quality: Dict[str, Any],
               check_data: bool = False, expected: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """Inspect every file of one model/experiment/member or reference dataset."""
    report = dict(label=label, kind="grid", level="full" if check_data else "headers",
                  findings=[], files=[], identity={}, values_checked=False)
    if check_data:
        report["values"] = _new_stats()
    records, calendars, known_units, identities, grids, hashes = [], set(), set(), {}, set(), {}
    pending = []
    requested_frequency = spec.get("aggregate", {}).get("frequency", "auto")
    frequency = requested_frequency
    time_blocks = []
    metadata_frequencies = set()
    missing_units = []
    for path in files:
        header = dict(path=str(path))
        report["files"].append(header)
        try:
            stat = Path(path).stat()
            header.update(size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
            if not stat.st_size:
                raise ConfigError("empty input file")
            with _open(path) as ds:
                da, variable = _prepare(ds, spec)
                values, parts, calendar = _time_info(da, spec)
                numeric = np.issubdtype(values.dtype, np.number)
                units = canonical_units(da.attrs.get("units"))
                declared = canonical_units(spec.get("units_in"))
                if units:
                    known_units.add(units)
                else:
                    missing_units.append(str(path))
                if declared and units and declared != units:
                    _finding(report, "error", "units_mismatch", f"units_in={spec.get('units_in')!r} disagrees with header units={da.attrs.get('units')!r}", file=str(path))
                transform = spec.get("units_transform") or {}
                output_units = canonical_units(transform.get("units"))
                if (output_units and (units or declared) and output_units != (units or declared)
                        and float(transform.get("scale", 1)) == 1 and float(transform.get("offset", 0)) == 0):
                    _finding(report, "warning", "units_relabelled_without_conversion",
                             "unit label changes while scale=1 and offset=0; verify the intended physical conversion", file=str(path))
                calendars.add(calendar)
                stamps = [_stamp(v) for v in values]
                if np.any(values[1:] < values[:-1]):
                    _finding(report, "warning", "unsorted_time", "timestamps are out of order; loader sorts them before aggregation", file=str(path))
                metadata = {}
                for key, aliases in {"source_id": ("source_id", "model_id"),
                                     "experiment_id": ("experiment_id",),
                                     "variant_label": ("variant_label", "member_id"),
                                     "grid_label": ("grid_label",), "table_id": ("table_id",)}.items():
                    found = next((ds.attrs.get(alias) or da.attrs.get(alias) for alias in aliases
                                  if ds.attrs.get(alias) or da.attrs.get(alias)), None)
                    if found:
                        value = str(found).strip()
                        metadata[key] = value
                        identities.setdefault(key, set()).add(value)
                        if expected and expected.get(key) and value != str(expected[key]):
                            _finding(report, "error", "identity_mismatch", f"{key}={value!r} disagrees with selected {expected[key]!r}", file=str(path))
                declared_freq = declared_frequency({"frequency": ds.attrs.get("frequency") or da.attrs.get("frequency"),
                                                     "table_id": metadata.get("table_id")})
                if declared_freq:
                    metadata_frequencies.add(declared_freq)
                grid = hashlib.sha256()
                for coord in ("latitude", "longitude"):
                    array = np.asarray(da[coord].values, dtype=np.float64)
                    grid.update(coord.encode()); grid.update(array.shape.__repr__().encode()); grid.update(array.tobytes())
                grids.add(grid.hexdigest())
                if np.unique(da.longitude.values).size != da.longitude.size:
                    _finding(report, "warning", "cyclic_longitude", "duplicate cyclic longitude coordinates require value equality in the loader", file=str(path))
                header.update(variable=variable, dims=dict(da.sizes), units=str(da.attrs.get("units", "")),
                              calendar=calendar, first_time=stamps[0], last_time=stamps[-1],
                              identity=metadata, grid_fingerprint=grid.hexdigest())
                time_blocks.append(values)
                pending.append((path, header, parts, stamps))
        except Exception as exc:
            _finding(report, "error", "unreadable_input", f"{type(exc).__name__}: {exc}", file=str(path))
    if missing_units:
        _finding(report, "error" if quality.get("require_units", False) else "warning", "missing_units",
                 "input units absent; configured conversions cannot be checked against the file header",
                 count=len(missing_units), examples=missing_units[:5])
    for code, values in (("mixed_units", known_units), ("mixed_calendars", calendars), ("mixed_grids", grids)):
        if len(values) > 1:
            _finding(report, "error", code, f"one dataset mixes incompatible {code[6:]} across files", examples=sorted(values)[:5])
    for key, values in identities.items():
        if len(values) > 1:
            _finding(report, "error", "mixed_identity", f"one dataset mixes {key}: {sorted(values)}")
        else:
            report["identity"][key] = next(iter(values))
    if known_units:
        report["identity"]["units"] = sorted(known_units)[0]
    if len(calendars) == 1:
        report["identity"]["calendar"] = next(iter(calendars))
    if len(metadata_frequencies) > 1:
        _finding(report, "error", "mixed_frequency", "file headers declare different sampling frequencies")
    report["input_frequency"] = None
    if time_blocks:
        try:
            all_times = np.concatenate(time_blocks)
            report["input_frequency"] = resolve_frequency(all_times, "auto", declared_frequencies=metadata_frequencies)
            frequency = resolve_frequency(all_times, requested_frequency, declared_frequencies=metadata_frequencies)
            if (requested_frequency == "auto" and frequency == "annual" and not metadata_frequencies
                    and not np.issubdtype(all_times.dtype, np.number)):
                _finding(report, "warning", "ambiguous_frequency", "one date per year was inferred as annual; set aggregate.frequency explicitly if input is sparse monthly data")
        except (ConfigError, TypeError, ValueError) as exc:
            _finding(report, "error", "frequency_resolution", str(exc))
    if frequency == "auto":
        # Retain an explicit diagnostic shape for unreadable/mixed files; errors
        # above prevent this fallback from being accepted as valid input.
        frequency = "monthly"
    if metadata_frequencies and requested_frequency != "auto" and any(f != frequency for f in metadata_frequencies):
        _finding(report, "error", "frequency_mismatch", f"configured frequency={frequency} disagrees with declared metadata {sorted(metadata_frequencies)}")
    # Determine sampling from all files before selecting months; yearly inputs
    # carry one label per season and must not be filtered by that label's month.
    for path, header, parts, stamps in pending:
        indices = np.flatnonzero(_selection(parts, years, months, annual=frequency == "annual"))
        header["selected_timesteps"] = int(indices.size)
        records.extend((stamps[i], parts[i]) for i in indices)
        if check_data and indices.size:
            try:
                with _open(path) as ds:
                    da, _ = _prepare(ds, spec)
                    stats, conflicts = _scan_values(da.isel(time=indices), [stamps[i] for i in indices], quality,
                                                   spec.get("units_transform") or {}, hashes)
                header["values"] = stats
                _merge_stats(report["values"], stats)
                if conflicts:
                    _finding(report, "error", "conflicting_duplicate_time", "overlapping timestamps contain different selected values",
                             file=str(path), count=len(conflicts), examples=conflicts[:5])
            except Exception as exc:
                _finding(report, "error", "unreadable_values", f"{type(exc).__name__}: {exc}", file=str(path))
    duplicates = len(records) - len({stamp for stamp, _ in records})
    if duplicates:
        _finding(report, "warning", "duplicate_time", "overlapping timestamps found; identical selected values are allowed, conflicting values are rejected" if check_data else
                 "overlapping timestamps found; value consistency is not checked by a header scan",
                 count=duplicates)
    if not records:
        _finding(report, "error", "no_selected_time", "no input samples lie inside the requested years/months")
    elif len(calendars) == 1:
        _coverage(report, records, years, months, frequency, next(iter(calendars)), quality.get("coverage", "warn"))
    if check_data:
        _finish_values(report, quality)
        report["values_checked"] = bool(files and report["values"]["count"] and
            all("selected_timesteps" in h and (h["selected_timesteps"] == 0 or "values" in h)
                for h in report["files"]))
        report["files_values_scanned"] = sum("values" in h for h in report["files"])
    report["ok"] = not any(f["severity"] == "error" for f in report["findings"])
    return report


def audit_series(files, spec, *, label, years, quality, check_data=False):
    from . import io
    report = dict(label=label, kind="series", level="full" if check_data else "headers",
                  findings=[], files=[], identity={}, values_checked=False)
    if check_data:
        report["values"] = _new_stats()
    records = []
    for path in files:
        header = {"path": str(path)}
        report["files"].append(header)
        try:
            stamp = Path(path).stat()
            header.update(size_bytes=stamp.st_size, mtime_ns=stamp.st_mtime_ns)
            y, values = io.read_series_file(path, spec.get("series_columns"))
            selected = (y >= years[0]) & (y <= years[1])
            records.extend((str(int(v)), (int(v), 1, 1)) for v in y[selected])
            header.update(n_years=int(y.size), first_year=int(y[0]), last_year=int(y[-1]), selected_timesteps=int(selected.sum()))
            if check_data and selected.any():
                # Reuse the same bounded scanner with a one-cell spatial grid.
                da = xr.DataArray(values[selected, None, None], dims=("time", "latitude", "longitude"))
                stats, _ = _scan_values(da, [str(int(v)) for v in y[selected]], quality,
                                        spec.get("units_transform") or {}, {})
                header["values"] = stats
                _merge_stats(report["values"], stats)
        except Exception as exc:
            _finding(report, "error", "unreadable_input", f"{type(exc).__name__}: {exc}", file=str(path))
    if not records:
        _finding(report, "error", "no_selected_time", "no annual samples inside requested years")
    else:
        _coverage(report, records, years, list(range(1, 13)), "annual", "annual_numeric", quality.get("coverage", "warn"))
    if check_data:
        _finish_values(report, quality)
        report["values_checked"] = bool(files and report["values"]["count"] and
            all("selected_timesteps" in h and (h["selected_timesteps"] == 0 or "values" in h)
                for h in report["files"]))
        report["files_values_scanned"] = sum("values" in h for h in report["files"])
    report["ok"] = not any(f["severity"] == "error" for f in report["findings"])
    return report
