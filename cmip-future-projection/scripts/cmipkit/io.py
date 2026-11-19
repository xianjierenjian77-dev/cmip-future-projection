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

    if "latitude" in da.dims and da.latitude.size > 1:
        if float(da.latitude[0]) > float(da.latitude[-1]):
            da = da.sortby("latitude")
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
        _, keep = np.unique(lon_sorted, return_index=True)
        da = da.isel(longitude=np.sort(keep))
    # close the seam so linear interpolation works at the wrap-around point
    if da.longitude.size > 1:
        first, last = float(da.longitude[0]), float(da.longitude[-1])
        span = last - first
        if span < (360.0 - 0.5):
            seam = da.isel(longitude=slice(0, 1)).assign_coords(longitude=[last + (float(da.longitude[1]) - first)])
            da = xr.concat([da, seam], dim="longitude")
    return da


def convert_units(da: xr.DataArray, spec: Dict[str, Any]) -> xr.DataArray:
    scale = float(spec.get("scale", 1.0) or 1.0)
    offset = float(spec.get("offset", 0.0) or 0.0)
    if scale == 1.0 and offset == 0.0:
        return da
    return da * scale + offset


def select_plev(da: xr.DataArray, plev: Any) -> xr.DataArray:
    """Select a pressure level given in hPa (ERA5) or Pa (CMIP6)."""
    if plev is None:
        return da
    if "plev" not in da.dims and "plev" not in da.coords:
        return da
    vals = np.asarray(da.plev.values, dtype=float)
    target = float(plev)
    candidates = [abs(v - target) for v in vals]
    if target <= 2000 and vals.max() > 5000:  # user gave hPa, data in Pa
        converted = [abs(v - target * 100.0) for v in vals]
        if min(converted) < min(candidates):
            candidates = converted
    idx = int(np.argmin(candidates))
    if candidates[idx] > max(1.0, 0.02 * target):
        log.warning("pressure level %s not found in %s; nearest is %s", target, vals.tolist(), vals[idx])
    out = da.isel(plev=idx).drop_vars("plev", errors="ignore")
    return out


def select_bbox(da: xr.DataArray, bbox: Optional[Dict[str, Sequence[float]]]) -> xr.DataArray:
    if not bbox:
        return da
    lat_rng = bbox.get("lat")
    lon_rng = bbox.get("lon")
    if lat_rng and "latitude" in da.dims:
        lo, hi = sorted(float(v) for v in lat_rng)
        da = da.sel(latitude=slice(lo, hi))
    if lon_rng and "longitude" in da.dims:
        lo, hi = sorted(float(v) for v in lon_rng)
        lon = da.longitude.values
        if lon.max() > 180 and lo < 0:  # data 0-360, request -180-180
            lo, hi = lo % 360, hi % 360
        da = da.sel(longitude=slice(lo, hi))
    return da


# ---------------------------------------------------------------------------
# year aggregation
# ---------------------------------------------------------------------------
def _season_year(time_values, months: Sequence[int]) -> np.ndarray:
    """Map each timestep to the year the *season* belongs to.

    Only DJF needs special care: December belongs to the following year's winter.
    """
    months = list(months)
    times = np.atleast_1d(time_values)
    if times.size and hasattr(times[0], "year"):  # cftime objects
        years = np.asarray([t.year for t in times], dtype=int)
        mon = np.asarray([t.month for t in times], dtype=int)
    else:  # numpy datetime64 (int() would give nanoseconds, not years)
        dt = np.asarray(times).astype("datetime64[ns]")
        years = dt.astype("datetime64[Y]").astype(int) + 1970
        mon = dt.astype("datetime64[M]").astype(int) % 12 + 1
    if sorted(months) == [1, 2, 12]:
        years = np.where(mon == 12, years + 1, years)
    return years


def _reduce(values: np.ndarray, stat: str, threshold: Optional[float]) -> np.ndarray:
    """Reduce along axis 0 (time) according to ``stat``."""
    stat = str(stat).lower()
    with np.errstate(invalid="ignore", divide="ignore"):
        if stat == "mean":
            return np.nanmean(values, axis=0)
        if stat == "sum":
            return np.nansum(values, axis=0)
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
            return np.sum(table[op], axis=0).astype(float)
    raise ConfigError(f"unsupported aggregate.stat {stat!r}")


def aggregate_yearly(
    da: xr.DataArray,
    months: Sequence[int],
    stat: str = "mean",
    threshold: Optional[float] = None,
) -> xr.DataArray:
    """Aggregate monthly/daily data to one value per (season-)year on (lat, lon)."""
    if "time" not in da.dims:
        raise ConfigError("aggregate_yearly needs a 'time' dimension")
    months = list(months)
    if months and sorted(months) != list(range(1, 13)):
        mon = da.time.dt.month
        da = da.isel(time=np.isin(mon.values, months))
    if da.time.size == 0:
        return xr.DataArray(
            np.empty((0, da.sizes.get("latitude", 0), da.sizes.get("longitude", 0))),
            dims=("year", "latitude", "longitude"),
            coords={"year": [], "latitude": da.latitude, "longitude": da.longitude},
        )
    labels = _season_year(da.time.values, months)
    out_years = np.unique(labels)
    rows = []
    for yr in out_years:
        sel = da.isel(time=np.where(labels == yr)[0])
        red = _reduce(sel.values.astype(np.float64), stat, threshold)
        rows.append(red)
    arr = np.stack(rows, axis=0)
    da_out = xr.DataArray(
        arr,
        dims=("year", "latitude", "longitude"),
        coords={
            "year": np.asarray(out_years, dtype=int),
            "latitude": da.latitude.values if "latitude" in da.coords else np.arange(arr.shape[1]),
            "longitude": da.longitude.values if "longitude" in da.coords else np.arange(arr.shape[2]),
        },
        name=da.name or "field",
    )
    if "plev" in da.coords:
        da_out = da_out.assign_coords(plev=da.plev)
    return da_out


# ---------------------------------------------------------------------------
# file opening
# ---------------------------------------------------------------------------
def _is_grib(path: str) -> bool:
    return path.lower().endswith(GRID_ENGINES)


def open_dataarray(path: str, variable: Optional[str] = None, plev: Any = None) -> xr.DataArray:
    """Open one NetCDF/GRIB file and return a normalised DataArray."""
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
            ds = xr.open_dataset(path, decode_times=True, use_cftime=True)
        except Exception:
            ds = xr.open_dataset(path)
    try:
        if variable:
            if variable not in ds.variables:
                raise ConfigError(f"variable {variable!r} not in {path}; available: {list(ds.data_vars)}")
            da = ds[variable]
        else:
            if not ds.data_vars:
                raise ConfigError(f"no data variables in {path}")
            da = ds[list(ds.data_vars)[0]]
        da = da.copy(deep=False)
        da = normalize_dims(da)
        da = select_plev(da, plev)
        da = align_longitude(da, lon_360=True)
        return da
    finally:
        ds.close()


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
    out = list(user)
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
                regex.append(".*")
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
                member="*", variable=spec.get("variable") or "*"), path)
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


def _pick_member(files: Sequence[str], model: str, experiment: str) -> List[str]:
    """When the member is a wildcard, keep the most common single member."""
    groups: Dict[str, List[str]] = {}
    for path in files:
        member = _member_of(path, model, experiment) or "?"
        groups.setdefault(member, []).append(path)
    if len(groups) <= 1:
        return list(files)
    best = max(groups, key=lambda k: (len(groups[k]), k))
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
    for member in list(spec.get("members") or ["*"]) + ["*"]:
        paths = _files_for_member(spec, experiment, model, member, years)
        if paths:
            if member == "*":
                paths = _pick_member(paths, model or "", experiment)
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
) -> xr.DataArray:
    """Open many files, aggregate each to yearly (lat, lon) and concatenate."""
    if not files:
        raise ConfigError("no input files matched the configured pattern")
    blocks: List[xr.DataArray] = []
    for path in files:
        try:
            da = open_dataarray(path, variable, plev)
        except ConfigError as exc:
            log.warning("skipping %s (%s)", path, exc)
            continue
        if "time" not in da.dims:
            log.warning("skipping %s (no time dimension)", path)
            continue
        da = select_bbox(da, bbox)
        da = convert_units(da, units_transform or {})
        yearly = aggregate_yearly(da, months, stat, threshold)
        blocks.append(yearly)
        del da
    if not blocks:
        raise ConfigError(f"none of the {len(files)} files could be read (variable/plev/time missing)")
    field = xr.concat(blocks, dim="year") if len(blocks) > 1 else blocks[0]
    field = field.sortby("year")
    _, unique_idx = np.unique(field.year.values, return_index=True)
    if len(unique_idx) != field.year.size:
        log.warning("duplicate years found in input; keeping the first occurrence")
        field = field.isel(year=np.sort(unique_idx))
    if year_range is not None:
        field = field.sel(year=slice(int(year_range[0]), int(year_range[1])))
    field.name = "field"
    return field


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
    return {"latitude": np.asarray(lat, dtype=float), "longitude": np.sort(lon)}


def regrid(da: xr.DataArray, grid: Dict[str, np.ndarray], method: str = "linear") -> xr.DataArray:
    """Interpolate a field to the target lat/lon (longitude wrap handled)."""
    if "latitude" not in da.dims or "longitude" not in da.dims:
        return da
    da = align_longitude(da, lon_360=True)
    lat_t = grid["latitude"]
    lon_t = grid["longitude"]
    same = (
        da.latitude.size == lat_t.size
        and da.longitude.size == lon_t.size
        and np.allclose(np.sort(da.latitude.values), np.sort(lat_t))
        and np.allclose(np.sort(da.longitude.values), np.sort(lon_t))
    )
    if same:
        out = da
    else:
        out = da.interp(latitude=lat_t, longitude=lon_t, method=method)
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

    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.cache_dir = os.path.join(cfg.outdir, "cache")
        self.use_cache = bool(cfg.get_path("runtime.cache", True))
        if self.use_cache:
            os.makedirs(self.cache_dir, exist_ok=True)
        self.grid = None
        self.model_members: Dict[str, List[str]] = {}
        if cfg.get_path("data.reference.kind", "grid") == "grid" or cfg.get_path("data.models.kind", "grid") == "grid":
            self.grid = resolve_target_grid(cfg)

    # -- cache ----------------------------------------------------------
    def _cache_path(self, key: Dict[str, Any]) -> str:
        digest = hashlib.md5(json.dumps(key, sort_keys=True, default=str).encode()).hexdigest()[:12]
        return os.path.join(self.cache_dir, f"{key['tag']}_{digest}.nc")

    def _cached(self, key: Dict[str, Any], builder):
        path = self._cache_path(key)
        if self.use_cache and os.path.exists(path):
            try:
                da = xr.open_dataarray(path)
                da.load()
                log.info("cache hit %s", os.path.basename(path))
                return da
            except Exception:
                pass
        da = builder()
        if self.use_cache:
            try:
                da.to_netcdf(path)
            except Exception as exc:  # pragma: no cover
                log.debug("cache write failed: %s", exc)
        return da

    # -- grid -----------------------------------------------------------
    def reference_field(self) -> xr.DataArray:
        spec = self.cfg.get_path("data.reference", {})
        if spec.get("kind") == "series":
            raise ConfigError("reference_field() called for a series source")
        years = spec.get("years") or self.cfg.period("historical")

        def build():
            files = []
            for template in spec.get("paths", []):
                files.extend(expand_paths(template, _year_list(years)))
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
            )
            field = regrid(field, self.grid, self.cfg.get_path("grid.method", "linear"))
            return field

        return self._cached({"tag": "reference", "years": ["%s" % y for y in years], "var": spec.get("variable")}, build)

    def model_field(self, model: str, experiment: str, years: Sequence[int]) -> xr.DataArray:
        spec = self.cfg.get_path("data.models", {})
        if spec.get("kind") == "series":
            raise ConfigError("model_field() called for a series source")
        combine = str(spec.get("member_combine", "first")).lower()
        if combine not in ("first", "mean", "median"):
            raise ConfigError(f"data.models.member_combine must be first|mean|median, got {combine!r}")

        if combine == "first":
            files = _files_for_grid_source(spec, experiment, model, years)
            if not files:
                raise ConfigError(f"no input files matched the configured pattern for {model}/{experiment}")
            log.info("%s / %s: %d file(s)", model, experiment, len(files))
            return self._load_and_regrid(files, spec, years, cache_key=f"{model}_{experiment}")

        # ---- multi-member ensemble -----------------------------------
        configured = spec.get("members") or ["*"]
        if "*" in configured or configured == ["*"]:
            members = discover_members(spec, experiment, model, years)
        else:
            members = [m for m in configured if _files_for_member(spec, experiment, model, m, years)]
            if not members:
                members = discover_members(spec, experiment, model, years)
        if not members:
            raise ConfigError(f"no ensemble members found for {model}/{experiment}")
        fields = []
        for member in members:
            files = _files_for_member(spec, experiment, model, member, years)
            if not files:
                log.warning("  %s/%s: member %s has no files, skipped", model, experiment, member)
                continue
            try:
                fields.append(self._load_and_regrid(files, spec, years, cache_key=f"{model}_{experiment}_{member}"))
            except ConfigError as exc:
                log.warning("  %s/%s: member %s failed (%s), skipped", model, experiment, member, exc)
        if not fields:
            raise ConfigError(f"no usable member for {model}/{experiment}")
        self.model_members[model] = members
        if len(fields) == 1:
            return fields[0]
        stacked = xr.concat(fields, dim="member")
        field = stacked.mean("member", skipna=True) if combine == "mean" else stacked.median("member", skipna=True)
        log.info("%s / %s: %d members combined (%s): %s", model, experiment, len(fields), combine, ", ".join(members))
        return field.rename("field")

    def _load_and_regrid(self, files: Sequence[str], spec: Dict[str, Any], years: Sequence[int],
                         cache_key: str) -> xr.DataArray:
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
            )
            return regrid(field, self.grid, self.cfg.get_path("grid.method", "linear"))

        return self._cached(
            {"tag": cache_key, "years": ["%s" % y for y in years], "var": spec.get("variable")}, build
        )

    # -- discovery ------------------------------------------------------
    def discover_models(self) -> List[str]:
        spec = self.cfg.get_path("data.models", {})
        if spec.get("kind") == "series":
            root = spec.get("series_dir")
            names = []
            for entry in sorted(os.listdir(root)):
                path = os.path.join(root, entry)
                if os.path.isdir(path):
                    names.append(entry)
                elif entry.lower().endswith((".csv", ".npz")):
                    names.append(re.sub(r"[_.].*$", "", entry))
            return sorted(set(names))
        configured = spec.get("models", "all")
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
        }
        names: List[str] = []
        for pattern in grid_patterns(spec):
            if "{model}" not in pattern:
                continue
            glob_pat = pattern.format(model="*", **fmt)
            files = sorted(glob.glob(glob_pat))
            if not files:
                continue
            compiled = re.compile(_pattern_to_regex(pattern))
            for path in files:
                match = compiled.match(path)
                if match:
                    names.append(match.group(1))
        if not names:
            # last resort: directories that contain the expected files
            for pattern in grid_patterns(spec):
                if "{model}" not in pattern:
                    continue
                glob_pat = pattern.format(model="*", **fmt)
                base = glob_pat.split("*")[0]
                base_dir = base if os.path.isdir(base) else os.path.dirname(base.rstrip("/"))
                names.extend(
                    os.path.basename(p) for p in sorted(glob.glob(os.path.join(base_dir, "*")))
                    if os.path.isdir(p)
                )
        return sorted(set(names))

    # -- series ---------------------------------------------------------
    def reference_series(self) -> Tuple[np.ndarray, np.ndarray]:
        spec = self.cfg.get_path("data.reference", {})
        path = spec.get("series_file")
        years, values = read_series_file(path, spec.get("series_columns"))
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
                return read_series_file(path, None)
        raise FileNotFoundError(f"no series file for {model}/{experiment} under {root}")


def read_series_file(path: str, columns: Optional[Dict[str, str]] = None) -> Tuple[np.ndarray, np.ndarray]:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if path.lower().endswith(".npz"):
        data = np.load(path, allow_pickle=True)
        years = data["years"] if "years" in data else data[list(data.keys())[0]]
        values = data["values"] if "values" in data else data[list(data.keys())[1]]
        return np.asarray(years, dtype=float), np.asarray(values, dtype=float)
    import pandas as pd

    df = pd.read_csv(path)
    cols = columns or {"year": "year", "value": "value"}
    if cols["year"] not in df.columns or cols["value"] not in df.columns:
        raise ConfigError(f"{path} must have columns {cols['year']!r} and {cols['value']!r}; found {list(df.columns)}")
    df = df[[cols["year"], cols["value"]]].dropna().sort_values(cols["year"])
    return df[cols["year"]].to_numpy(dtype=float), df[cols["value"]].to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------
def save_field(da: xr.DataArray, path: str, name: str = "field") -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    out = da.rename(name) if da.name != name else da
    encoding = {name: {"zlib": True, "complevel": 4}}
    try:
        out.to_netcdf(path, encoding=encoding, engine="netcdf4")
    except Exception:
        out.to_netcdf(path)
    return path


def save_npz(path: str, **arrays) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    np.savez_compressed(path, **arrays)
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
    yl = _year_list(years)
    if experiment is None:
        files: List[str] = []
        for template in spec.get("paths", []) or []:
            files.extend(expand_paths(template, yl))
    else:
        files = _files_for_grid_source(spec, experiment, model, years)
    if not files:
        raise ConfigError("no daily input files matched the configured pattern")

    blocks: List[xr.DataArray] = []
    for path in files:
        try:
            da = open_dataarray(path, spec.get("variable"), spec.get("plev"))
        except ConfigError as exc:
            log.warning("skipping %s (%s)", path, exc)
            continue
        if "time" not in da.dims:
            continue
        da = select_bbox(da, bbox or spec.get("bbox"))
        da = convert_units(da, spec.get("units_transform"))
        months_arr = da.time.dt.month.values
        sel = np.isin(months_arr, list(months))
        if not sel.any():
            continue
        blocks.append(da.isel(time=sel))
    if not blocks:
        raise ConfigError("no daily timesteps inside the requested months")
    field = xr.concat(blocks, dim="time")
    field = field.sortby("time")
    times = np.asarray(field.time.values)
    years_arr = np.asarray([getattr(t, "year", int(str(t)[:4])) for t in times], dtype=int)
    keep = (years_arr >= int(years[0])) & (years_arr <= int(years[1]))
    if not keep.any():
        raise ConfigError(f"no daily timesteps inside {years[0]}-{years[1]}")
    field = field.isel(time=keep)
    # drop duplicated timestamps that overlapping files can produce
    times = np.asarray(field.time.values)
    if times.size and len(np.unique(times)) != times.size:
        _, first = np.unique(times, return_index=True)
        field = field.isel(time=np.sort(first))
    n_values = field.sizes["time"] * field.sizes.get("latitude", 1) * field.sizes.get("longitude", 1)
    if n_values > max_values:
        raise ConfigError(
            f"daily field would hold {n_values/1e6:.1f}M values; "
            f"set data.*.bbox to the study region or shorten the period"
        )
    field.name = "field"
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
    if itype in ("region_mean", "region_mean_anomaly"):
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
    raise ConfigError(f"daily extremes do not support index.type={itype!r}; use region_mean or region_diff")
