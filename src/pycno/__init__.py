"""PyCNO: Physiologically based radiopharmacokinetic modeling utilities."""

import ctypes
import os
import sys
from pathlib import Path

# --- Preload libpython for libroadrunner in Conda environments ---
conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix:
    libpython_path = (
        Path(conda_prefix)
        / "lib"
        / f"libpython{sys.version_info.major}.{sys.version_info.minor}.so.1.0"
    )
    if libpython_path.exists():
        try:
            ctypes.CDLL(str(libpython_path))
        except OSError as e:
            print(f"Warning: failed to preload libpython: {e}")

from .modeling.dosing import Dose  # noqa: E402
from .modeling.modeling import Model  # noqa: E402

__all__ = ["Dose", "Model"]
