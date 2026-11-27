"""Real-spawn regressions for parallel execution and scientific equivalence.

Run with the same interpreter as the workflow. Tests compare complete model
application windows; they never split QDM by year or assert a speed-up.
"""
from __future__ import annotations

import copy
import hashlib
import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest

import numpy as np
import xarray as xr
import yaml

SKILL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL / "scripts"))
from cmipkit.config import ConfigError, load_config
from cmipkit.pipeline import Pipeline

CLI = SKILL / "scripts" / "run.py"
MODELS = ["MODEL-C", "MODEL-A", "MODEL-B"]
SCENARIOS = ["ssp245", "ssp585"]
MEMBERS = ["r1i1p1f1", "r2i1p1f1"]


# These functions must be importable under spawn, rather than local closures or
# fork-only mocks. The barrier makes concurrent workers observable without
# depending on CPU speed or on how quickly a trivial job happens to complete.
def barrier_worker(item):
    from threadpoolctl import threadpool_info
    root = Path(item["root"])
    (root / f"ready-{os.getpid()}").write_text("ready")
    deadline = time.monotonic() + 30
    while len(list(root.glob("ready-*"))) < 2:
        if time.monotonic() > deadline:
            raise RuntimeError("second worker never reached the barrier")
        time.sleep(0.01)
    return {"value": item["value"], "pid": os.getpid(),
            "threads": [entry["num_threads"] for entry in threadpool_info()]}


def crash_worker(item):
    os._exit(23)


def echo_worker(item):
    return item


def fail_initializer():
    raise ValueError("intentional initializer failure")


def sleeping_worker(directory):
    (Path(directory) / f"sleeping-{os.getpid()}").write_text("ready")
    time.sleep(120)
    return os.getpid()


def large_reference_initializer(reference):
    pass


def owned_children():
    """Use OS ownership, including zombies not registered by multiprocessing."""
    result = set()
    for path in Path("/proc").glob("[0-9]*/status"):
        try:
            contents = path.read_text()
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        for line in contents.splitlines():
            if line.startswith("PPid:") and int(line.split()[1]) == os.getpid():
                result.add(int(path.parent.name))
                break
    return result


def startup_interrupt_probe():
    """Interrupt the real spawn handoff while its reference payload is sent.

    A spawned interpreter imports this module's scientific dependencies before
    reading the large payload. No startup sleep or patched Process is used.
    """
    from multiprocessing import resource_tracker
    from cmipkit.parallel import ordered_process_map
    resource_tracker.ensure_running()
    before = owned_children()
    started, errors = [], []
    during_start = False

    def interrupt():
        nonlocal during_start
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            new = owned_children() - before
            frame = sys._current_frames().get(threading.main_thread().ident)
            in_start = False
            while frame is not None:
                if (frame.f_code.co_name == "start" and
                        Path(frame.f_code.co_filename).as_posix().endswith("/multiprocessing/process.py")):
                    in_start = True
                    break
                frame = frame.f_back
            if new and in_start:
                started.extend(sorted(new))
                during_start = True
                os.kill(os.getpid(), signal.SIGINT)
                return
            time.sleep(0.001)
        errors.append("did not observe an OS child while Process.start was active")

    timer = threading.Thread(target=interrupt, daemon=True)
    timer.start()
    interrupted = False
    try:
        list(ordered_process_map(echo_worker, [1, 2], workers=2,
                                 initializer=large_reference_initializer,
                                 initargs=(b"x" * (16 * 1024 * 1024),)))
    except KeyboardInterrupt:
        interrupted = True
    timer.join(timeout=25)
    leftover = owned_children() - before
    report = {"interrupted": interrupted, "during_start": during_start,
              "started": started, "leftover": sorted(leftover), "errors": errors}
    print(json.dumps(report), flush=True)
    # Keep a failed regression from leaving live processes or zombies behind.
    for pid in leftover:
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        try:
            os.waitpid(pid, 0)
        except ChildProcessError:
            pass
    return 0 if interrupted and during_start and started and not leftover and not errors else 1


class ParallelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cmip_parallel_")
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        self.data = self.root / "data"
        self.data.mkdir()

    def write_config(self, config, name):
        path = self.root / f"{name}.yaml"
        path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        return path

    def run_pipeline(self, config, name, workers, stages=("report",)):
        config = copy.deepcopy(config)
        config["project"]["outdir"] = str(self.root / name)
        config.setdefault("runtime", {}).update(workers=workers, worker_threads=1)
        pipe = Pipeline(load_config(str(self.write_config(config, name))))
        self.run_and_check_log_release(pipe, stages)
        return pipe

    def run_and_check_log_release(self, pipe, stages):
        # Observe the real streams opened during run(), retaining references so
        # garbage collection cannot hide a detached-but-unclosed file handler.
        streams = {}

        def observe_streams(record):
            for handler in pipe.log.handlers:
                if isinstance(handler, logging.FileHandler) and handler.stream is not None:
                    streams[handler] = handler.stream
            return True

        pipe.log.addFilter(observe_streams)
        try:
            return pipe.run(stages)
        finally:
            pipe.log.removeFilter(observe_streams)
            self.assertTrue(streams, "the run did not write a file log")
            for handler, stream in streams.items():
                self.assertTrue(stream.closed, "run() left a log stream open")
                self.assertNotIn(handler, pipe.log.handlers,
                                 "a retained append-mode handler can reopen a completed run's log")
            self.assertFalse(any(isinstance(handler, logging.FileHandler)
                                 and handler.stream is not None and not handler.stream.closed
                                 for handler in pipe.log.handlers))
            # Windows rejects renaming an open file; keep this actual filesystem
            # check as well as the stream assertions used on every platform.
            log_path = Path(pipe.outdir) / "logs" / "run.log"
            moved = log_path.with_name("run.released.log")
            log_path.rename(moved)
            moved.rename(log_path)

    def run_cli(self, *arguments, timeout=90):
        return subprocess.run([sys.executable, str(CLI), *map(str, arguments)],
                              cwd=self.root, capture_output=True, text=True, timeout=timeout)

    def assert_nested_close(self, first, second):
        if isinstance(first, (xr.Dataset, xr.DataArray)):
            xr.testing.assert_allclose(first, second, rtol=1e-12, atol=1e-12)
        elif isinstance(first, dict):
            self.assertEqual(set(first), set(second))
            for key in first:
                with self.subTest(field=key):
                    self.assert_nested_close(first[key], second[key])
        elif isinstance(first, (tuple, list)):
            self.assertEqual(len(first), len(second))
            for one, two in zip(first, second):
                self.assert_nested_close(one, two)
        elif isinstance(first, np.ndarray):
            if first.dtype.kind in "biufc":
                np.testing.assert_allclose(first, second, rtol=1e-12, atol=1e-12, equal_nan=True)
            else:
                np.testing.assert_array_equal(first, second)
        elif isinstance(first, (float, np.floating)):
            np.testing.assert_allclose(first, second, rtol=1e-12, atol=1e-12, equal_nan=True)
        else:
            self.assertEqual(first, second)

    def assert_products_equal(self, first, second, directories):
        for directory in directories:
            left = Path(first.outdir) / directory
            right = Path(second.outdir) / directory
            files = sorted(p.relative_to(left) for p in left.rglob("*") if p.is_file())
            self.assertEqual(files, sorted(p.relative_to(right) for p in right.rglob("*") if p.is_file()))
            for relative in files:
                with self.subTest(product=f"{directory}/{relative}"):
                    a, b = left / relative, right / relative
                    if a.suffix == ".nc":
                        with xr.open_dataset(a) as da, xr.open_dataset(b) as db:
                            xr.testing.assert_allclose(da, db, rtol=1e-12, atol=1e-12)
                    elif a.suffix == ".npz":
                        # Only test-generated, trusted products are opened here.
                        with np.load(a, allow_pickle=True) as aa, np.load(b, allow_pickle=True) as bb:
                            self.assertEqual(set(aa.files), set(bb.files))
                            for key in aa.files:
                                self.assert_nested_close(aa[key], bb[key])
                    elif a.suffix == ".csv":
                        self.assertEqual(a.read_text(), b.read_text())

    def assert_parallel_stage(self, pipe, stage, requested, effective):
        manifest = json.loads((Path(pipe.outdir) / "manifest.json").read_text())
        meta = manifest["parallel"]
        self.assertEqual(meta["requested_workers"], requested)
        self.assertEqual(meta["worker_threads"], 1)
        detail = meta["stages"][stage]
        self.assertEqual(detail["effective_workers"], effective)
        self.assertGreaterEqual(detail["elapsed_seconds"], 0)
        if effective > 1:
            self.assertEqual(meta["start_method"], "spawn")
            self.assertTrue(detail["worker_pids"])
            self.assertNotIn(os.getpid(), detail["worker_pids"])
        self.assertEqual([task["model"] for task in detail["tasks"]], pipe.models)
        for task in detail["tasks"]:
            self.assertGreaterEqual(task["elapsed_seconds"], 0)
            self.assertGreaterEqual(task["finished_at"], task["started_at"])
        return manifest

    def series_config(self):
        years = np.arange(1981, 2001)
        future = np.arange(2031, 2051)
        observed = 10 + 0.15 * np.arange(20) + np.sin(np.arange(20) / 2)
        np.savetxt(self.data / "observed.csv", np.column_stack((years, observed)),
                   delimiter=",", header="year,value", comments="")
        for index, model in enumerate(MODELS):
            for experiment in ["historical", *SCENARIOS]:
                yy = years if experiment == "historical" else future
                values = observed + 1 + index + (0 if experiment == "historical" else 3)
                if experiment == "ssp585":
                    values = values + 0.1 * np.arange(20)
                if model == "MODEL-A" and experiment == "ssp245":
                    keep = yy != 2035
                    yy, values = yy[keep], values[keep]
                np.savetxt(self.data / f"{model}_{experiment}.csv", np.column_stack((yy, values)),
                           delimiter=",", header="year,value", comments="")
        return {
            "project": {"name": "parallel_series", "outdir": str(self.root / "out")},
            "variable": {"name": "index", "units": "K"},
            "data": {"reference": {"kind": "series", "series_file": str(self.data / "observed.csv")},
                     "models": {"kind": "series", "series_dir": str(self.data),
                                "models": MODELS, "scenarios": SCENARIOS}},
            "periods": {"historical": [1981, 2000], "reference_climate": [1981, 2000],
                        "future": [2031, 2050], "future_windows": [[2031, 2040], [2041, 2050]]},
            "index": {"type": "series"},
            "correction": {"method": "quantile_delta_mapping", "multiplicative": False},
            "ensemble": {"min_models": 2, "weights": "skill"},
            "extremes": {"enabled": False}, "plots": {"figures": []},
            "runtime": {"cache": False, "verbose": False},
        }

    def grid_config(self):
        rng = np.random.default_rng(892)
        latitude, longitude = [20., 22.], [100., 102., 104.]
        hist = np.arange("1981-01", "2001-01", dtype="datetime64[M]").astype("datetime64[ns]")
        future = np.arange("2031-01", "2051-01", dtype="datetime64[M]").astype("datetime64[ns]")

        def write(name, dates, offset, metadata=None):
            native_lat, native_lon = latitude, longitude
            if metadata and metadata["source_id"] == "MODEL-A":
                native_lat, native_lon = [19., 21., 23.], [99., 101., 103., 105.]
            elif metadata and metadata["source_id"] == "MODEL-B":
                native_lat, native_lon = [19., 20., 21., 22., 23.], [99., 100., 101., 102., 103., 104., 105.]
            shape = (len(dates), len(native_lat), len(native_lon))
            t = np.arange(len(dates))
            values = 280 + offset + 0.007 * t[:, None, None] + rng.normal(0, 0.8, shape)
            values += np.sin(t[:, None, None] * np.pi / 6)
            values += (np.asarray(native_lat)[None, :, None] - 20) * 0.1
            values += (np.asarray(native_lon)[None, None, :] - 100) * 0.1
            values[:, -1, -1] = np.nan
            ds = xr.Dataset({"tas": (("time", "lat", "lon"), values)},
                            coords={"time": dates, "lat": native_lat, "lon": native_lon},
                            attrs=metadata or {"frequency": "mon"})
            ds.tas.attrs["units"] = "K"
            ds.to_netcdf(self.data / name)

        write("observed.nc", hist, 0)
        for i, model in enumerate(MODELS):
            for k, member in enumerate(MEMBERS):
                for j, experiment in enumerate(["historical", *SCENARIOS]):
                    write(f"{model}_{experiment}_{member}.nc", hist if j == 0 else future,
                          1 + i + k / 2 + j * 1.5,
                          {"source_id": model, "experiment_id": experiment, "variant_label": member,
                           "grid_label": "gn", "frequency": "mon", "table_id": "Amon"})
        aggregate = {"frequency": "monthly", "months": [6, 7, 8], "stat": "mean"}
        return {
            "project": {"name": "parallel_grid", "outdir": str(self.root / "out")},
            "variable": {"name": "tas", "units": "K"},
            "data": {"reference": {"kind": "grid", "paths": [str(self.data / "observed.nc")],
                                   "variable": "tas", "years": [1981, 2000], "aggregate": aggregate},
                     "models": {"kind": "grid", "root": str(self.data),
                                "pattern": "{root}/{model}_{experiment}_{member}.nc",
                                "models": MODELS, "members": MEMBERS, "member_combine": "mean",
                                "scenarios": SCENARIOS, "variable": "tas", "aggregate": aggregate}},
            # Workers must inherit the target established by the reference;
            # each model deliberately has a different native grid.
            "grid": {"lat": None, "lon": None},
            "periods": {"historical": [1981, 2000], "reference_climate": [1981, 2000],
                        "future": [2031, 2050], "future_windows": [[2031, 2040], [2041, 2050]]},
            "index": {"type": "region_mean", "region": {"lat": [20, 22], "lon": [100, 104]},
                      "anomaly": True, "standardize": True},
            "correction": {"method": "quantile_delta_mapping", "multiplicative": False, "n_quantiles": 50},
            "ensemble": {"min_models": 2, "weights": "skill"},
            "quality": {"scan": "full", "coverage": "error"},
            "extremes": {"enabled": False}, "plots": {"figures": []},
            "runtime": {"cache": True, "verbose": False},
        }

    def daily_config(self):
        rng = np.random.default_rng(417)
        latitude, longitude = [30., 32.], [100., 102.]
        bbox = {"lat": [30, 32], "lon": [100, 102]}
        hist = np.concatenate([np.datetime64(f"{year}-06-01") + np.arange(92) for year in range(1981, 2001)])
        future = np.concatenate([np.datetime64(f"{year}-06-01") + np.arange(92) for year in range(2031, 2051)])

        def write(name, dates, bias, trend, metadata=None):
            years = dates.astype("datetime64[Y]").astype(int) + 1970
            noise = rng.normal(0, 1, len(dates))
            spikes = rng.gamma(2, 2, len(dates)) * (rng.random(len(dates)) < 0.12)
            values = 28 + bias + trend * (years - years[0]) + noise + spikes
            values = values[:, None, None] + np.array([[0., 0.3], [0.1, 0.4]])[None]
            ds = xr.Dataset({"tas": (("time", "lat", "lon"), values)},
                            coords={"time": dates.astype("datetime64[ns]"), "lat": latitude, "lon": longitude},
                            attrs=metadata or {"frequency": "day"})
            ds.tas.attrs["units"] = "degC"
            ds.to_netcdf(self.data / name)

        write("observed.nc", hist, 0, 0.03)
        for i, model in enumerate(MODELS):
            for j, experiment in enumerate(["historical", *SCENARIOS]):
                write(f"{model}_{experiment}_{MEMBERS[0]}.nc", hist if j == 0 else future,
                      1 + i / 2 + j, 0.025 + j * 0.02,
                      {"source_id": model, "experiment_id": experiment, "variant_label": MEMBERS[0],
                       "grid_label": "gn", "frequency": "day", "table_id": "day"})
        aggregate = {"frequency": "daily", "months": [6, 7, 8], "stat": "mean"}
        return {
            "project": {"name": "parallel_daily", "outdir": str(self.root / "out")},
            "variable": {"name": "tas", "units": "degC"},
            "data": {"reference": {"kind": "grid", "paths": [str(self.data / "observed.nc")],
                                   "variable": "tas", "years": [1981, 2000], "aggregate": aggregate,
                                   "bbox": bbox},
                     "models": {"kind": "grid", "root": str(self.data),
                                "pattern": "{root}/{model}_{experiment}_{member}.nc",
                                "models": MODELS, "members": [MEMBERS[0]], "scenarios": SCENARIOS,
                                "variable": "tas", "aggregate": aggregate, "bbox": bbox}},
            "grid": {"lat": {"start": 30, "stop": 32, "step": 2},
                     "lon": {"start": 100, "stop": 102, "step": 2}},
            "periods": {"historical": [1981, 2000], "reference_climate": [1981, 2000],
                        "future": [2031, 2050], "future_windows": [[2031, 2050]]},
            "index": {"type": "region_mean", "region": {"lat": [30, 32], "lon": [100, 102]},
                      "anomaly": False},
            "correction": {"method": "mean_std", "multiplicative": False},
            "ensemble": {"min_models": 2}, "extremes": {"enabled": False},
            "daily": {"enabled": True, "months": [6, 7, 8], "correction": "quantile_delta_mapping",
                      "plots": False,
                      "extremes": {"kind": "pot", "threshold_mode": "quantile", "threshold_quantile": 0.9,
                                   "min_peaks": 20, "min_years": 12, "decluster_days": 3,
                                   "return_periods": [2, 10, 50]}},
            "plots": {"figures": []}, "runtime": {"cache": False, "verbose": False},
        }

    def test_grid_qdm_members_products_and_diagnostics_match_serial(self):
        config = self.grid_config()
        hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in self.data.glob("*.nc")}
        serial = self.run_pipeline(config, "serial_grid", 1)
        parallel = self.run_pipeline(config, "parallel_grid", 2)
        self.assertEqual(list(parallel.series), MODELS)
        for key in ("series", "raw_series", "scores", "params", "map_acc", "pdf_data", "ens"):
            with self.subTest(state=key):
                self.assert_nested_close(getattr(serial, key), getattr(parallel, key))
        self.assert_products_equal(serial, parallel, ["corrected", "params", "series", "ensemble"])
        manifest = self.assert_parallel_stage(parallel, "correct", 2, 2)
        self.assertEqual(manifest["status"], "completed")
        self.assertEqual(manifest["models_loaded"], MODELS)
        self.assertEqual(manifest["model_members"], {model: MEMBERS for model in MODELS})
        expected_metadata = {f"{model}/{member}" for model in MODELS for member in MEMBERS}
        self.assertTrue(expected_metadata.issubset(serial.source.model_metadata))
        self.assertEqual(serial.source.model_metadata, parallel.source.model_metadata)
        source = {item["label"]: item for item in serial.manifest["actual_inputs"]}
        self.assertEqual(source, {item["label"]: item for item in manifest["actual_inputs"]})
        output = Path(parallel.outdir)
        self.assertEqual(len(list(output.rglob("manifest.json"))), 1)
        self.assertEqual(len(list(output.rglob("index_series.npz"))), 1)
        self.assertEqual(len(list(output.glob("corrected/*.nc"))), len(MODELS) * 3)
        self.assertEqual(hashes, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in hashes})

    def test_series_qdm_alignment_and_skill_weights_match_serial(self):
        config = self.series_config()
        serial = self.run_pipeline(config, "serial_series", 1)
        parallel = self.run_pipeline(config, "parallel_series", 2)
        self.assertEqual(list(parallel.series), MODELS)
        for key in ("series", "raw_series", "scores", "params", "ens"):
            self.assert_nested_close(getattr(serial, key), getattr(parallel, key))
        self.assert_products_equal(serial, parallel, ["params", "series", "ensemble"])
        self.assert_parallel_stage(parallel, "correct", 2, 2)
        ensemble = parallel.ens["ssp245"]
        self.assertEqual(ensemble["n_models_by_year"][ensemble["years"] == 2035].item(), 2)

    def test_daily_pot_and_block_match_serial(self):
        config = self.daily_config()
        for kind in ("pot", "block"):
            with self.subTest(kind=kind):
                config["daily"]["extremes"]["kind"] = kind
                serial = self.run_pipeline(config, f"serial_{kind}", 1, ["daily"])
                parallel = self.run_pipeline(config, f"parallel_{kind}", 2, ["daily"])
                self.assert_parallel_stage(parallel, "daily", 2, 2)
                self.assertEqual(parallel.daily["kind"], kind)
                self.assert_nested_close(serial.daily["historical"].as_dict(), parallel.daily["historical"].as_dict())
                for scenario in SCENARIOS:
                    a, b = serial.daily["scenarios"][scenario], parallel.daily["scenarios"][scenario]
                    self.assert_nested_close(a["curves"], b["curves"])
                    self.assert_nested_close(a["stats"], b["stats"])
                    self.assertEqual(b["n_models"], len(MODELS))
                    for model in MODELS:
                        self.assert_nested_close(a["fits"][model].as_dict(), b["fits"][model].as_dict())
                self.assert_products_equal(serial, parallel, ["extremes", "tables"])
                self.assertEqual({x["label"]: x for x in serial.manifest["actual_inputs"]},
                                 {x["label"]: x for x in parallel.manifest["actual_inputs"]})
                self.assertEqual(parallel.manifest["model_members"], {model: [MEMBERS[0]] for model in MODELS})

    def test_worker_scientific_failure_stops_strict_run_and_records_model(self):
        config = self.series_config()
        bad_path = self.data / "MODEL-A_historical.csv"
        values = np.loadtxt(bad_path, delimiter=",", skiprows=1)
        values[:, 1] = -1
        np.savetxt(bad_path, values, delimiter=",", header="year,value", comments="")
        config["correction"].update(method="mean", multiplicative=True)
        config["runtime"].update(workers=2, worker_threads=1, strict=True)
        cfg = load_config(str(self.write_config(config, "strict_failure")))
        pipe = Pipeline(cfg)
        with self.assertRaisesRegex(Exception, "nonnegative"):
            self.run_and_check_log_release(pipe, ["correct"])
        manifest = json.loads((Path(pipe.outdir) / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertNotIn("correct", manifest["completed_stages"])
        self.assertTrue(any("MODEL-A" in item["dataset"] for item in manifest["failures"]))
        self.assertIn("Run failed:", (Path(pipe.outdir) / "logs" / "run.log").read_text(encoding="utf-8"))

    def test_worker_scientific_failure_can_skip_without_losing_other_models(self):
        config = self.series_config()
        bad_path = self.data / "MODEL-A_historical.csv"
        values = np.loadtxt(bad_path, delimiter=",", skiprows=1)
        values[:, 1] = -1
        np.savetxt(bad_path, values, delimiter=",", header="year,value", comments="")
        config["correction"].update(method="mean", multiplicative=True)
        config["runtime"]["strict"] = False
        parallel = self.run_pipeline(config, "skip_failure", 2)
        self.assertEqual(list(parallel.series), ["MODEL-C", "MODEL-B"])
        self.assertEqual(parallel.manifest["status"], "completed_with_skips")
        self.assertTrue(any("MODEL-A" in item["dataset"] for item in parallel.manifest["failures"]))
        self.assertEqual(parallel.ens["ssp245"]["n_models"], 2)

    def test_single_model_caps_worker_count(self):
        config = self.series_config()
        config["data"]["models"]["models"] = [MODELS[0]]
        config["ensemble"]["min_models"] = 1
        parallel = self.run_pipeline(config, "one_model", 8)
        self.assert_parallel_stage(parallel, "correct", 8, 1)
        self.assertEqual(list(parallel.series), [MODELS[0]])

    def test_worker_setting_can_change_in_existing_output(self):
        config = self.series_config()
        serial = self.run_pipeline(config, "same_output", 1)
        before = copy.deepcopy(serial.series)
        parallel = self.run_pipeline(config, "same_output", 2)
        self.assert_nested_close(before, parallel.series)
        self.assertIsNone(parallel.previous_run_directory)
        self.assert_parallel_stage(parallel, "correct", 2, 2)
        # The same instance must acquire a fresh log for its next invocation,
        # append to the existing history, and release that new stream as well.
        log_path = Path(parallel.outdir) / "logs" / "run.log"
        previous_log = log_path.read_text(encoding="utf-8")
        self.run_and_check_log_release(parallel, ["report"])
        current_log = log_path.read_text(encoding="utf-8")
        self.assertTrue(current_log.startswith(previous_log))
        self.assertEqual(current_log[len(previous_log):].count(f"=== {parallel.name} | stages:"), 1)
        self.assert_nested_close(before, parallel.series)
        self.assert_parallel_stage(parallel, "correct", 2, 2)

    def test_cli_overrides_workers_and_threads_before_dry_run(self):
        config = self.series_config()
        config["runtime"].update(workers=1, worker_threads=2)
        cfg_path = self.write_config(config, "cli_override")
        result = self.run_cli("run", cfg_path, "--workers", "2", "--worker-threads", "1", "--dry-run", "--json")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["parallel"]["requested_workers"], 2)
        self.assertEqual(payload["parallel"]["effective_workers"], 2)
        self.assertEqual(payload["parallel"]["worker_threads"], 1)
        self.assertFalse(Path(config["project"]["outdir"]).exists())

    def test_parallel_configuration_rejects_invalid_counts_and_duplicate_models(self):
        config = self.series_config()
        for key in ("workers", "worker_threads"):
            for bad in (0, -1, True, 1.5, "2", None):
                with self.subTest(key=key, value=bad):
                    candidate = copy.deepcopy(config)
                    candidate["runtime"][key] = bad
                    with self.assertRaises(ConfigError):
                        load_config(str(self.write_config(candidate, "invalid_count")))
        config["data"]["models"]["models"] = ["MODEL-A", "MODEL-A"]
        with self.assertRaises(ConfigError):
            load_config(str(self.write_config(config, "duplicate_models")))

    def test_cli_run_starts_workers_and_saves_effective_configuration(self):
        config = self.series_config()
        config["runtime"].update(workers=1, worker_threads=2)
        cfg_path = self.write_config(config, "cli_parallel")
        result = self.run_cli("run", cfg_path, "--workers", "2", "--worker-threads", "1", "--stages", "ensemble")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = Path(config["project"]["outdir"])
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "completed")
        self.assertEqual(manifest["parallel"]["requested_workers"], 2)
        self.assertEqual(manifest["parallel"]["stages"]["correct"]["effective_workers"], 2)
        self.assertTrue(manifest["parallel"]["stages"]["correct"]["worker_pids"])
        resolved = yaml.safe_load((output / "config.resolved.yaml").read_text())
        self.assertEqual(resolved["runtime"]["workers"], 2)
        self.assertEqual(resolved["runtime"]["worker_threads"], 1)

    def test_cli_rejects_invalid_parallel_counts(self):
        cfg_path = self.write_config(self.series_config(), "invalid_cli")
        for flag, value in (("--workers", "0"), ("--workers", "-2"),
                            ("--worker-threads", "1.5"), ("--worker-threads", "false")):
            with self.subTest(flag=flag, value=value):
                result = self.run_cli("run", cfg_path, flag, value, "--dry-run")
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)

    def test_spawn_workers_are_concurrent_ordered_and_limit_native_threads(self):
        from cmipkit.parallel import ordered_process_map
        tasks = [{"root": str(self.root), "value": value} for value in (7, 2)]
        results = list(ordered_process_map(barrier_worker, tasks, workers=2, worker_threads=1))
        self.assertEqual([result["index"] for result in results], [0, 1])
        self.assertEqual([result["value"]["value"] for result in results], [7, 2])
        self.assertEqual(len({result["pid"] for result in results}), 2)
        self.assertNotIn(os.getpid(), [result["pid"] for result in results])
        for result in results:
            self.assertEqual(result["pid"], result["value"]["pid"])
            self.assertTrue(all(threads <= 1 for threads in result["value"]["threads"]))

    def _helper_script(self, body):
        script = self.root / "helper.py"
        script.write_text("import sys\n" +
                          f"sys.path.insert(0, {str(SKILL / 'tests')!r})\n" +
                          "from test_parallel import crash_worker, echo_worker, fail_initializer, sleeping_worker, startup_interrupt_probe\n" +
                          "from cmipkit.parallel import ordered_process_map\n" +
                          "if __name__ == '__main__':\n" +
                          "\n".join("    " + line for line in body.splitlines()) + "\n")
        return script

    def _kill_group(self, process):
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        elif process.poll() is None:
            process.kill()
        process.wait(timeout=10)

    def test_worker_crash_and_initializer_failure_exit_without_hanging(self):
        for task in ("list(ordered_process_map(crash_worker, [1, 2], workers=2))",
                     "list(ordered_process_map(echo_worker, [1, 2], workers=2, initializer=fail_initializer))"):
            with self.subTest(task=task):
                script = self._helper_script("try:\n    " + task +
                                             "\nexcept Exception as exc:\n    print(type(exc).__name__, str(exc))\n"
                                             "    sys.exit(0)\nraise SystemExit('unexpected success')")
                process = subprocess.Popen([sys.executable, str(script)], stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE, text=True, start_new_session=(os.name == "posix"))
                try:
                    stdout, stderr = process.communicate(timeout=30)
                    self.assertEqual(process.returncode, 0, stdout + stderr)
                    self.assertIn("Error", stdout)
                finally:
                    self._kill_group(process)

    @unittest.skipUnless(os.name == "posix", "SIGINT subprocess lifecycle contract is POSIX-specific")
    def test_interrupt_reaps_active_workers(self):
        script = self._helper_script(
            "try:\n    list(ordered_process_map(sleeping_worker, " + repr([str(self.root)] * 2) +
            ", workers=2))\nexcept KeyboardInterrupt:\n    sys.exit(130)")
        process = subprocess.Popen([sys.executable, str(script)], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, start_new_session=True)
        worker_pids = []
        try:
            deadline = time.monotonic() + 20
            while len(list(self.root.glob("sleeping-*"))) < 2 and process.poll() is None:
                if time.monotonic() > deadline:
                    self.fail("workers did not start within the bounded setup interval")
                time.sleep(0.05)
            worker_pids = [int(path.name.split("-")[-1]) for path in self.root.glob("sleeping-*")]
            self.assertEqual(len(worker_pids), 2)
            process.send_signal(signal.SIGINT)
            stdout, stderr = process.communicate(timeout=15)
            self.assertEqual(process.returncode, 130, stdout + stderr)
            for pid in worker_pids:
                with self.assertRaises(ProcessLookupError, msg=f"worker {pid} survived parent interruption"):
                    os.kill(pid, 0)
        finally:
            self._kill_group(process)

    @unittest.skipUnless(os.name == "posix" and Path("/proc/self/status").exists(),
                         "spawn-startup OS child ownership check needs Linux /proc and SIGINT")
    def test_interrupt_during_spawn_reference_handoff_reaps_unregistered_child(self):
        script = self._helper_script("raise SystemExit(startup_interrupt_probe())")
        process = subprocess.Popen([sys.executable, str(script)], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, start_new_session=True)
        try:
            stdout, stderr = process.communicate(timeout=45)
            self.assertEqual(process.returncode, 0, stdout + stderr)
            result = json.loads(stdout.strip().splitlines()[-1])
            self.assertTrue(result["interrupted"])
            self.assertTrue(result["during_start"])
            self.assertTrue(result["started"])
            self.assertEqual(result["leftover"], [])
            self.assertEqual(result["errors"], [])
        finally:
            self._kill_group(process)


if __name__ == "__main__":
    unittest.main(verbosity=2)
