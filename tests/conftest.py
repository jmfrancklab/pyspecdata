import importlib.util
import pathlib
import sys
import types
import numpy as np
import os

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[1] / "pyspecdata"


def load_module(name: str):
    """Load a pyspecdata submodule without requiring compiled extensions."""
    # provide dummy replacements for optional compiled dependencies
    sys.modules.setdefault("_nnls", types.ModuleType("_nnls"))
    sys.modules.setdefault("tables", types.ModuleType("tables"))
    sys.modules.setdefault("h5py", types.ModuleType("h5py"))
    if "numpy.core.rec" not in sys.modules:
        rec = types.ModuleType("rec")
        rec.fromarrays = np.core.records.fromarrays
        sys.modules["numpy.core.rec"] = rec
    sys.modules.setdefault("pyspecdata.fornotebook", types.ModuleType("fornotebook"))
    fig_stub = sys.modules.setdefault("pyspecdata.figlist", types.ModuleType("figlist"))
    if not hasattr(fig_stub, "figlist"):
        fig_stub.figlist = type("figlist", (), {})
    os.environ.setdefault("pyspecdata_figures", "standard")

    pkg = sys.modules.setdefault("pyspecdata", types.ModuleType("pyspecdata"))
    if not hasattr(pkg, "__path__"):
        pkg.__path__ = [str(PACKAGE_ROOT)]

    spec = importlib.util.spec_from_file_location(
        f"pyspecdata.{name}", PACKAGE_ROOT / f"{name.replace('.', '/')}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[f"pyspecdata.{name}"] = module
    return module
