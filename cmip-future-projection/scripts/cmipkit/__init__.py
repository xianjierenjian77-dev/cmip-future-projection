"""cmipkit - a general, config-driven pipeline for CMIP6 (or any model-ensemble)
bias correction, future projection, extreme-value (GEV) analysis and publication
figures.

The package uses local files and ordinary scientific Python dependencies.
Heavy modules are loaded lazily so doctor/help work before installation.

Typical programmatic use::

    from cmipkit import load_config, Pipeline
    cfg = load_config("config.yaml")
    pl = Pipeline(cfg)
    pl.run()

or from the command line::

    python .../scripts/run.py run config.yaml
"""

from __future__ import annotations

__version__ = "1.2.0"

__all__ = ["Config", "ConfigError", "load_config", "write_resolved", "Pipeline", "__version__"]


def __getattr__(name):
    if name == "Pipeline":
        from .pipeline import Pipeline
        return Pipeline
    if name in ("Config", "ConfigError", "load_config", "write_resolved"):
        from . import config
        return getattr(config, name)
    raise AttributeError(name)
