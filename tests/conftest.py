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
        ticker_stub = types.ModuleType("matplotlib.ticker")
        transforms_stub = types.ModuleType("matplotlib.transforms")

        class DummyImg:
            def get_clim(self):
                return (-1, 1)

            def set_clim(self, *_args):
                pass

        def gci():
            return DummyImg()

        class DummyAx:
            def __init__(self):
                pass

        def gca():
            return DummyAx()

        def sca(_ax):
            pass

        def imshow(*_args, **_kwargs):
            return DummyImg()

        def plot(*_args, **_kwargs):
            return DummyImg()

        def xlabel(*_args, **_kwargs):
            pass

        def ylabel(*_args, **_kwargs):
            pass

        def title(*_args, **_kwargs):
            pass

        def colorbar(*_args, **_kwargs):
            return DummyImg()

        def setp(*_args, **_kwargs):
            pass

        def rc(*_args, **_kwargs):
            pass

        pylab_stub.gci = gci
        pylab_stub.gca = gca
        pylab_stub.sca = sca
        pylab_stub.imshow = imshow
        pylab_stub.plot = plot
        pylab_stub.xlabel = xlabel
        pylab_stub.ylabel = ylabel
        pylab_stub.title = title
        pylab_stub.colorbar = colorbar
        pylab_stub.setp = setp
        pylab_stub.cm = types.SimpleNamespace(gray=None)
        pylab_stub.rc = rc
        pylab_stub.rcParams = {
            "axes.prop_cycle": types.SimpleNamespace(by_key=lambda: {"color": ["k"]})
        }
        pylab_stub.pi = np.pi
        pylab_stub.r_ = np.r_
        ticker_stub.AutoMinorLocator = type("AutoMinorLocator", (), {})
        ticker_stub.MaxNLocator = type("MaxNLocator", (), {"__init__": lambda self, *a, **k: None})
        class DummyBbox:
            @staticmethod
            def union(_bboxes):
                return DummyBbox()

            def transformed(self, *_args):
                return self

            @property
            def width(self):
                return 0

            @property
            def height(self):
                return 0

        transforms_stub.Bbox = DummyBbox
        mpl_stub.pylab = pylab_stub
        mpl_stub.pyplot = pylab_stub
        mpl_stub.ticker = ticker_stub
        mpl_stub.transforms = transforms_stub
        sys.modules["matplotlib"] = mpl_stub
        sys.modules["matplotlib.pylab"] = pylab_stub
        sys.modules["matplotlib.pyplot"] = pylab_stub
        sys.modules["matplotlib.ticker"] = ticker_stub
        sys.modules["matplotlib.transforms"] = transforms_stub
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

        class DummyUnit(str):
            def __format__(self, _spec):
                return str(self)

        class DummyQuantity:
            def __init__(self, magnitude, units=""):
                self.magnitude = magnitude
                self.units = DummyUnit(units)

            def _combine_units(self, other, op):
                if isinstance(other, DummyQuantity):
                    other_units = other.units
                    other = other.magnitude
                else:
                    other_units = ""
                if op == "*":
                    if self.units and other_units:
                        units = f"{self.units}*{other_units}"
                    else:
                        units = self.units or other_units
                    magnitude = self.magnitude * other
                elif op == "/":
                    if self.units and other_units:
                        units = f"{self.units}/{other_units}"
                    else:
                        units = self.units or other_units
                    magnitude = self.magnitude / other
                elif op == "+":
                    units = self.units
                    magnitude = self.magnitude + other
                elif op == "-":
                    units = self.units
                    magnitude = self.magnitude - other
                else:
                    raise NotImplementedError
                return DummyQuantity(magnitude, units)

            __add__ = lambda self, other: self._combine_units(other, "+")
            __sub__ = lambda self, other: self._combine_units(other, "-")
            __mul__ = lambda self, other: self._combine_units(other, "*")
            __truediv__ = lambda self, other: self._combine_units(other, "/")
            __radd__ = __add__
            __rsub__ = lambda self, other: DummyQuantity(other, self.units)._combine_units(self, "-")
            __rmul__ = __mul__
            __rtruediv__ = lambda self, other: DummyQuantity(other, self.units)._combine_units(self, "/")

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
    try:
        import sympy  # noqa: F401
    except Exception:
        sympy_stub = types.ModuleType("sympy")

        class Expr:
            pass

        class Symbol(Expr):
            def __init__(self, name):
                self.name = name

        def sqrt(x):
            return np.sqrt(x)

        def nsimplify(x, *args, **kwargs):
            return x

        sympy_stub.Expr = Expr
        sympy_stub.Symbol = Symbol
        sympy_stub.sqrt = sqrt
        sympy_stub.nsimplify = nsimplify
        sympy_stub.I = 1j
        sympy_stub.pi = np.pi
        sympy_stub.utilities = types.SimpleNamespace(lambdify=lambda *a, **k: None)
        fem = types.ModuleType("sympy.functions.elementary.miscellaneous")
        fem.sqrt = sqrt
        sys.modules["sympy"] = sympy_stub
        sys.modules["sympy.functions"] = types.ModuleType("sympy.functions")
        sys.modules["sympy.functions.elementary"] = types.ModuleType("sympy.functions.elementary")
        sys.modules["sympy.functions.elementary.miscellaneous"] = fem
    try:
        import scipy  # noqa: F401
    except Exception:
        scipy_stub = types.ModuleType("scipy")
        sparse_stub = types.ModuleType("scipy.sparse")
        interpolate_stub = types.ModuleType("scipy.interpolate")
        interpolate_stub.interp1d = lambda *a, **k: None
        scipy_stub.sparse = sparse_stub
        scipy_stub.interpolate = interpolate_stub
        sys.modules["scipy"] = scipy_stub
        sys.modules["scipy.sparse"] = sparse_stub
        sys.modules["scipy.interpolate"] = interpolate_stub
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
