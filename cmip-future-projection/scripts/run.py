#!/usr/bin/env python
"""Entry point for the cmip-future-projection skill.

Examples
--------
    python run.py doctor
    python run.py selftest
    python run.py init myrun
    python run.py run myrun/config.yaml
    python run.py run myrun/config.yaml --stages correct,figures
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from cmipkit.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
