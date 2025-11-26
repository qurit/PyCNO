"""PyCNO: Physiologically based radiopharmacokinetic modeling utilities."""
 
import os
from pathlib import Path
import ctypes
import sys

# --- Preload libpython for libroadrunner in Conda environments ---
conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix:
    libpython_path = Path(conda_prefix) / "lib" / f"libpython{sys.version_info.major}.{sys.version_info.minor}.so.1.0"
    if libpython_path.exists():
        try:
            ctypes.CDLL(str(libpython_path))
        except OSError as e:
            print(f"Warning: failed to preload libpython: {e}")

from .functions import Model

__all__ = ["Model"]