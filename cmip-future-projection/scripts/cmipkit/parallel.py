"""Bounded, ordered process execution with explicit child-process ownership.

Each spawned process receives initialization data once, and only ``workers``
results can be waiting for an earlier result.  The caller must exhaust or close
``ordered_process_map``; closing it also terminates outstanding work.
"""
from __future__ import annotations

from contextlib import contextmanager
import multiprocessing as mp
from multiprocessing.connection import wait
import os
import signal
import threading
import time
import traceback
from typing import Any, Callable, Iterable, Iterator, Optional, Sequence


class RemoteWorkerError(RuntimeError):
    """A worker failed, exited unexpectedly, or could not transport a result."""


@contextmanager
def _complete_start_before_sigint():
    """Deliver Ctrl-C only after Process.start has retained the child handle.

    The Python handler defers interrupts even if another existing parent thread
    receives SIGINT.  On POSIX the spawning thread also blocks SIGINT, so the
    child inherits protection throughout bootstrap/imports.  Neither measure
    alone covers both the parent and the not-yet-initialized child.
    """
    previous_handler = signal.getsignal(signal.SIGINT)
    defer_handler = (threading.current_thread() is threading.main_thread()
                     and callable(previous_handler))
    pending = []

    def defer(signum, frame):
        if not pending:
            pending.append((signum, frame))

    if defer_handler:
        signal.signal(signal.SIGINT, defer)
    previous_mask = None
    try:
        if hasattr(signal, "pthread_sigmask"):
            previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT})
        yield
    finally:
        # Restore the mask while the temporary Python handler is still active:
        # an OS-pending signal must be recorded before restoring the handler.
        try:
            if previous_mask is not None:
                signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        finally:
            if defer_handler:
                signal.signal(signal.SIGINT, previous_handler)
                if pending:
                    previous_handler(*pending[0])


def _worker_entry(connection, function, initializer, initargs, worker_threads):
    # The parent owns cancellation.  A Ctrl-C must not leave child processes
    # trying to write tracebacks or incomplete results into an abandoned pipe.
    if hasattr(signal, "SIGINT"):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        if hasattr(signal, "pthread_sigmask"):
            signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGINT})
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[name] = str(worker_threads)
    try:
        from threadpoolctl import threadpool_limits
        # Imports in __main__ can already have loaded numerical libraries under
        # spawn.  Apply runtime limits too, and keep the controller alive.
        with threadpool_limits(limits=worker_threads):
            if initializer is not None:
                initializer(*initargs)
            # Initializers may import a numerical library for the first time.
            with threadpool_limits(limits=worker_threads):
                connection.send({"kind": "ready", "pid": os.getpid()})
                while True:
                    task = connection.recv()
                    if task is None:
                        break
                    index, item = task
                    started_at = time.time()
                    start = time.perf_counter()
                    try:
                        value = function(item)
                        connection.send({
                            "kind": "result", "index": index, "pid": os.getpid(),
                            "started_at": started_at, "finished_at": time.time(),
                            "elapsed_seconds": time.perf_counter() - start,
                            "value": value,
                        })
                    except BaseException as exc:
                        connection.send({"kind": "error", "index": index,
                                         "pid": os.getpid(), "error_type": type(exc).__name__,
                                         "message": str(exc), "traceback": traceback.format_exc()})
                        break
    except (EOFError, BrokenPipeError):
        pass  # The parent was interrupted or rejected another worker result.
    except BaseException as exc:
        try:
            connection.send({"kind": "error", "pid": os.getpid(),
                             "error_type": type(exc).__name__, "message": str(exc),
                             "traceback": traceback.format_exc()})
        except (EOFError, BrokenPipeError, OSError):
            pass
    finally:
        connection.close()


def _stop_workers(records, successful):
    started = []
    for record in records:
        if record["process"].pid is None:
            record["connection"].close()
            record["process"].close()
        else:
            started.append(record)
    records = started
    if successful:
        for record in records:
            try:
                record["connection"].send(None)
            except (OSError, EOFError):
                pass
    else:
        for record in records:
            if record["process"].is_alive():
                record["process"].terminate()
    deadline = time.monotonic() + 3.0
    for record in records:
        record["process"].join(max(0.0, deadline - time.monotonic()))
    for record in records:
        process = record["process"]
        if process.is_alive():
            process.kill()
    deadline = time.monotonic() + 3.0
    for record in records:
        process = record["process"]
        process.join(max(0.0, deadline - time.monotonic()))
        record["connection"].close()
        if not process.is_alive():
            process.close()


def ordered_process_map(
    function: Callable[[Any], Any], items: Iterable[Any], *, workers: int,
    worker_threads: int = 1, initializer: Optional[Callable[..., None]] = None,
    initargs: Sequence[Any] = (),
) -> Iterator[dict]:
    """Yield result/timing records in input order from real spawn workers.

    ``function`` and ``initializer`` must be importable module-level callables;
    use the standard ``if __name__ == '__main__':`` guard in Python programs.
    Worker crashes are detected without waiting forever for a missing result.
    Infrastructure/task errors abort the pool.  Applications may return their
    own structured error values when selected task failures are recoverable.
    """
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive integer")
    if isinstance(worker_threads, bool) or not isinstance(worker_threads, int) or worker_threads < 1:
        raise ValueError("worker_threads must be a positive integer")
    items = list(items)
    if not items:
        return
    context = mp.get_context("spawn")
    count = min(workers, len(items))
    records = []
    successful = False
    try:
        for index in range(count):
            parent, child = context.Pipe(duplex=True)
            process = context.Process(target=_worker_entry,
                                      args=(child, function, initializer, tuple(initargs), worker_threads),
                                      name=f"cmip-worker-{index + 1}")
            # Register ownership before start(), so an interrupt immediately
            # after it launches cannot leave an untracked live child.
            records.append({"process": process, "connection": parent,
                            "ready": False, "task": None})
            # spawn can launch an OS child and block transferring initargs
            # before assigning Process._popen.  Defer SIGINT across that gap.
            with _complete_start_before_sigint():
                try:
                    process.start()
                finally:
                    child.close()
        by_connection = {r["connection"]: r for r in records}
        pending = {}
        next_submit = 0
        next_yield = 0
        while next_yield < len(items):
            connections = wait(list(by_connection), timeout=0.2)
            for connection in connections:
                record = by_connection[connection]
                try:
                    message = connection.recv()
                except (EOFError, OSError) as exc:
                    process = record["process"]
                    process.join(timeout=0.1)
                    raise RemoteWorkerError(
                        f"worker {process.pid} exited unexpectedly (exitcode={process.exitcode}, "
                        f"task={record['task']}); inspect the worker log") from exc
                if message.get("kind") == "ready":
                    record["ready"] = True
                elif message.get("kind") == "error":
                    raise RemoteWorkerError(
                        f"worker {message.get('pid')} task {message.get('index', 'initialization')} "
                        f"failed: {message.get('error_type')}: {message.get('message')}\n"
                        f"{message.get('traceback', '')}")
                elif message.get("kind") == "result":
                    if message.get("index") != record["task"]:
                        raise RemoteWorkerError("worker returned an unexpected task index")
                    pending[message["index"]] = message
                    record["task"] = None
                else:
                    raise RemoteWorkerError("worker returned an invalid protocol message")
            for record in records:
                process = record["process"]
                if process.exitcode is not None:
                    raise RemoteWorkerError(
                        f"worker {process.pid} exited unexpectedly (exitcode={process.exitcode}, "
                        f"task={record['task']}); inspect the worker log")
            # Bounded lookahead prevents a slow first model from causing every
            # later model's gridded summaries to accumulate in parent memory.
            while next_yield in pending:
                message = pending.pop(next_yield)
                message.pop("kind", None)
                yield message
                next_yield += 1
            if all(record["ready"] for record in records):
                for record in records:
                    if (record["task"] is None and next_submit < len(items)
                            and next_submit < next_yield + count):
                        record["connection"].send((next_submit, items[next_submit]))
                        record["task"] = next_submit
                        next_submit += 1
        successful = True
    finally:
        _stop_workers(records, successful)
