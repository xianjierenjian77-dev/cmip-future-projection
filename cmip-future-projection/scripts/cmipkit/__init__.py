"""cmipkit - a general, config-driven pipeline for CMIP6 (or any model-ensemble)
bias correction, future projection, extreme-value (GEV) analysis and publication
figures.

The package is deliberately plain-python + xarray/numpy/scipy/matplotlib so that
it can run on the cluster (``/disk1/wsr/miniconda3/envs/climate/bin/python``)
without any custom services.

Typical programmatic use::

    from cmipkit import load_config, Pipeline
    cfg = load_config("config.yaml")
    pl = Pipeline(cfg)
    pl.run()

or from the command line::

    python .../scripts/run.py run config.yaml
"""

from __future__ import annotations

__version__ = "1.0.0"

from .config import Config, ConfigError, load_config, write_resolved  # noqa: F401
from .pipeline import Pipeline  # noqa: F401

__all__ = ["Config", "ConfigError", "load_config", "write_resolved", "Pipeline", "__version__"]
