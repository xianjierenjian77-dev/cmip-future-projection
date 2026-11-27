# Multi-core model processing

Contents: 1. Execution scope; 2. Configure a run; 3. Size resources;
4. Verify results and timing; 5. Diagnose failures; 6. Use the Python API.

Use model parallelism when independent models provide enough work to justify
starting additional processes. Preserve the scientific configuration and first
confirm a serial pilot with the same methods and representative inputs.

## 1. Execution scope

| Work | Execution |
|---|---|
| Input discovery and quality scans | Coordinating process, before analysis. |
| Reference data and projection-template preparation | Coordinating process. |
| Standard grid or annual-series correction | One task per model; multiple models can run concurrently. |
| Daily regional correction and POT/block analysis | One task per model after the reference daily series is prepared. |
| Scenarios and time windows within a model | Process in their configured order, retaining the complete calibration/application windows. |
| Ensemble reduction, final figures and report | Coordinating process after the required model results are collected. |

Workers start with Python's `spawn` method and open their own input files. This
avoids depending on inherited open NetCDF handles. The coordinating process
merges model results in the original model order, rather than completion order.
Model/member selection and missing-model policies still apply.

Parallelism does not change the data definition. Standard gridded correction
acts on annual/seasonal fields. The daily branch first reduces each day's grid
to a regional series; it does not become full-domain daily grid-cell correction.
This implementation uses local CPU processes, not GPU or multi-node execution.

## 2. Configure a run

~~~yaml
runtime:
  workers: 4
  worker_threads: 1
~~~

Both options require positive integers; `0`, negative values, booleans and
fractional numbers are invalid. `workers` defaults to `1`, retaining serial
execution. Model worker processes are capped at the number of selected models.
One selected model runs serially even if a larger worker count was requested.
The coordinating process and any Python resource-tracker process are additional
OS processes; `workers` is not a limit on the total number shown by a system monitor.

Override either setting for a single CLI invocation:

~~~bash
python /path/to/skill/scripts/run.py run myrun/config.yaml --workers 4 --worker-threads 1
python /path/to/skill/scripts/run.py run myrun/config.yaml --workers 4 --dry-run
~~~

A dry run validates the request and prints the plan without creating a process
pool or computing projections. On an actual run, `config.resolved.yaml` records
the effective CLI overrides. Changing worker/thread counts alone does not require
`--overwrite`.

Install the project dependencies, including `threadpoolctl`, before enabling
multiple worker processes. A serial run can continue without this parallel-only
dependency. `worker_threads` limits supported BLAS/OpenMP numerical thread pools
inside workers; it is not an operating-system limit on every possible thread,
and does not change the coordinating process's thread settings.

Run one coordinating command per output directory. Launching several copies of
the same command against one directory is not the supported parallel interface.

## 3. Size resources

Start with a representative serial run. Then try 4 workers, followed by 8 or 16
only if enough models, memory and storage throughput are available. Keep
`worker_threads: 1` initially to avoid multiplying process concurrency by large
native thread pools. On a scheduled cluster, stay within the allocated CPUs and
memory; the command does not request resources from the scheduler.

Each worker needs its own reference data, current model arrays, correction
parameters and temporary arrays. Estimate memory for the coordinating process
plus every active worker. `runtime.max_input_values` is still a per-load array
guard inside each process, not a global memory budget or reservation. Smaller
source regions and suitable output grids reduce work; changing scientific time
windows or member combinations for speed changes the analysis.

Shared disks, compressed NetCDF chunks and output writes can limit throughput.
Quality scans and final reductions remain serial work, and process startup adds
overhead. Small annual-series jobs can be faster with one worker. Do not promise
linear speedup or choose all available logical CPUs without measuring.

## 4. Verify results and timing

Compare serial and parallel runs with the same model list, members, methods,
periods, source versions and quality settings. Check corrected fields, parameters,
index series, ensemble counts/maps and daily results within appropriate numerical
tolerances. Confirm actual input inventories and any exclusions, not only exit
status or figure existence.

Use the `parallel` object in `manifest.json` to inspect the execution records:

| Field below `parallel` | Interpretation |
|---|---|
| `requested_workers`, `worker_threads` | Requested process count and per-child numerical thread limit. |
| `start_method` | The configured parallel strategy, `spawn`; this alone does not prove that children ran. |
| `stages.correct` / `stages.daily` | Execution details for each model stage actually entered. |
| `stages.<stage>.effective_workers` | Model execution concurrency; one uses the coordinating process. |
| `stages.<stage>.start_method` | Actual stage mode: `serial` or `spawn`. |
| `stages.<stage>.worker_pids` | PIDs whose results were received; a serial stage records the coordinator's PID. |
| `stages.<stage>.elapsed_seconds` | Wall time for model execution and collection, including process startup/cleanup when used. |
| `stages.<stage>.tasks` | Received task records: model, PID, start/finish timestamps, elapsed seconds and `completed`/`failed` status. Spawned tasks include a relative `log` path; failed tasks include `error`. |

Task timestamps `started_at` and `finished_at` are Unix seconds. A completed task
can contain recorded scenario or fit exclusions; also inspect the top-level
`failures` and final model counts. On failure or interruption, the task list may
not include work whose result was never received by the coordinator.

Read `logs/run.log` for the coordinating process and
`logs/{stage}-worker-{pid}.log` for worker details. Worker PIDs identify processes,
not CPU affinity or a guaranteed number of continuously busy cores. `report.md`
summarizes the stage mode, worker count and execution time when a report is produced.

Measure the full command wall time on the target machine. A task's worker time
does not include all common preparation and finalization costs. Compare equivalent
cold/warm prepared-cache states and record raw-data storage and quality-scan level;
OS filesystem caches can also affect repeated measurements. Do not extrapolate
production timing from the small synthetic selftest alone.

## 5. Diagnose failures

With `runtime.strict: true`, a raised model-processing error stops the run. With explicit
`runtime.strict: false`, ordinary model/scenario failures are recorded and the
supported remaining results may be assembled; inspect the resulting model counts.
Input-quality errors still stop preflight. Worker initialization failure or an
unexpected process exit is a pool/infrastructure failure and stops the run even
in permissive mode.

On a handled exception or Ctrl-C, the coordinator stops and joins its worker processes
and records failed/interrupted status. Completed files remain available for
inspection; their presence does not certify completion of the stage or whole run.
Files written by a worker whose result was not collected may not be represented
in the coordinator's model summary. Inspect the manifest's current artifacts and
worker logs together rather than treating the summary as a complete file inventory.
Keep the manifest and logs, fix the cause and rerun using the documented stage
and cache behavior. This is not an arbitrary artifact-resume system.

## 6. Use the Python API

Use an importable script and protect process creation with the main guard on
Linux, macOS and Windows:

~~~python
from cmipkit.config import load_config
from cmipkit.pipeline import Pipeline

def main():
    cfg = load_config("myrun/config.yaml")
    cfg.set_path("runtime.workers", 4)
    Pipeline(cfg).run()

if __name__ == "__main__":
    main()
~~~

Ensure the skill's `scripts` directory is on the import path. In notebooks,
prefer launching the CLI with the notebook environment's Python interpreter;
do not start spawned workers from unguarded cell-defined application code.
Keep `spawn` rather than switching to `fork` to hide an import/entry-point error.
