"""PyCNO: Physiologically based radiopharmacokinetic modeling utilities."""

import os
from pathlib import Path

# --- Compatibility fix for libroadrunner in Conda environments ---
def _fix_conda_ld_library_path():
    """Ensure Conda's libpython is discoverable."""
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        return

    lib_dir = Path(conda_prefix) / "lib"
    if not lib_dir.exists():
        return

    ld_path = os.environ.get("LD_LIBRARY_PATH", "")
    if str(lib_dir) not in ld_path.split(":"):
        os.environ["LD_LIBRARY_PATH"] = f"{lib_dir}:{ld_path}"
        try:
            import ctypes
            ctypes.CDLL(None)
        except Exception:
            pass

_fix_conda_ld_library_path()

from .functions import run_model

__all__ = ["run_model"]
