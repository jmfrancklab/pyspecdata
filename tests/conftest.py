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
    try:
        import matplotlib  # noqa: F401
    except Exception:
        mpl_stub = types.ModuleType("matplotlib")
        pylab_stub = types.ModuleType("matplotlib.pylab")

        class DummyImg:
            def get_clim(self):
                return (-1, 1)

            def set_clim(self, *_args):
                pass

        def gci():
            return DummyImg()

        pylab_stub.gci = gci
        pylab_stub.pi = np.pi
        pylab_stub.r_ = np.r_
        mpl_stub.pylab = pylab_stub
        sys.modules["matplotlib"] = mpl_stub
        sys.modules["matplotlib.pylab"] = pylab_stub
        sys.modules.setdefault("pylab", pylab_stub)
    try:
        import mpl_toolkits.mplot3d  # noqa: F401
    except Exception:
        mpl_toolkits_stub = types.ModuleType("mpl_toolkits")
        mplot3d_stub = types.ModuleType("mpl_toolkits.mplot3d")
        mplot3d_stub.axes3d = types.SimpleNamespace()
        mpl_toolkits_stub.mplot3d = mplot3d_stub
        sys.modules["mpl_toolkits"] = mpl_toolkits_stub
        sys.modules["mpl_toolkits.mplot3d"] = mplot3d_stub
    try:
        import pint  # noqa: F401
    except Exception:
        pint_stub = types.ModuleType("pint")

        class DummyQuantity:
            def __init__(self, magnitude, units=""):
                self.magnitude = magnitude
                self.units = units

            def to_base_units(self):
                return self

        class DummyUnitRegistry:
            Quantity = DummyQuantity

            def __call__(self, *args):
                if len(args) == 1:
                    mag = 1
                    units = args[0]
                else:
                    mag, units = args
                return DummyQuantity(mag, units)

            def define(self, *_args, **_kwargs):
                pass

        pint_stub.UnitRegistry = DummyUnitRegistry
        sys.modules["pint"] = pint_stub
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
