import importlib.util
import importlib.machinery
import pathlib
import sys
import types
import numpy as np
import os

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[1] / "pyspecdata"


def load_module(name: str, *, use_real_pint: bool = False):
    """Load a pyspecdata submodule without requiring external packages.

    Parameters
    ----------
    name
        The submodule to import from ``pyspecdata``.
    use_real_pint
        If ``True`` the real :mod:`pint` package is required and will be
        imported.  When ``False`` (the default), a very small dummy
        implementation of :mod:`pint` is supplied so that tests can run in an
        environment without the real dependency installed.
    """
    # provide dummy replacements for optional compiled dependencies
    sys.modules.setdefault("tables", types.ModuleType("tables"))
    sys.modules.setdefault("h5py", types.ModuleType("h5py"))

    if not use_real_pint:
        try:
            import pint  # noqa: F401
        except Exception:
            pint_stub = types.ModuleType("pint")
            pint_stub.__spec__ = importlib.machinery.ModuleSpec(
                "pint", loader=None
            )
            pint_stub.__pyspec_stub__ = True

            class DummyUnit:
                def __init__(self, name=""):
                    self.name = name

                def __str__(self):
                    return self.name

                def __format__(self, _spec):
                    return self.name

                __repr__ = __str__

            class DummyQuantity:
                def __init__(self, magnitude=1.0, units=""):
                    self.magnitude = magnitude
                    self._units = units

                def __mul__(self, other):
                    if isinstance(other, DummyQuantity):
                        return DummyQuantity(
                            self.magnitude * other.magnitude,
                            f"{self._units}*{other._units}".strip("*"),
                        )
                    return DummyQuantity(self.magnitude * other, self._units)

                def __truediv__(self, other):
                    if isinstance(other, DummyQuantity):
                        units = (
                            f"{self._units}/{other._units}"
                            if other._units
                            else self._units
                        )
                        return DummyQuantity(
                            self.magnitude / other.magnitude, units
                        )
                    return DummyQuantity(self.magnitude / other, self._units)

                def __add__(self, other):
                    other_val = (
                        other.magnitude
                        if isinstance(other, DummyQuantity)
                        else other
                    )
                    return DummyQuantity(
                        self.magnitude + other_val, self._units
                    )

                def __sub__(self, other):
                    other_val = (
                        other.magnitude
                        if isinstance(other, DummyQuantity)
                        else other
                    )
                    return DummyQuantity(
                        self.magnitude - other_val, self._units
                    )

                def to_compact(self):
                    return self

                @property
                def units(self):
                    return DummyUnit(self._units)

                def __repr__(self):
                    if self._units:
                        return f"{self.magnitude} {self._units}"
                    return str(self.magnitude)

            class DummyUnitRegistry:
                def __init__(self):
                    pass

                def define(self, *_a, **_k):
                    pass

                def Quantity(self, *args):
                    if len(args) == 1:
                        magnitude = 1.0
                        units = args[0]
                    elif len(args) == 2:
                        magnitude, units = args
                    else:
                        raise TypeError(
                            "Quantity expects one or two arguments"
                        )
                    return DummyQuantity(magnitude, units)

            pint_stub.UnitRegistry = DummyUnitRegistry
            pint_stub.Quantity = DummyQuantity
            sys.modules.setdefault("pint", pint_stub)
    else:
        if "pint" in sys.modules:
            del sys.modules["pint"]
        import pint  # noqa: F401
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
            "axes.prop_cycle": types.SimpleNamespace(
                by_key=lambda: {"color": ["k"]}
            )
        }
        pylab_stub.pi = np.pi
        pylab_stub.r_ = np.r_
        ticker_stub.AutoMinorLocator = type("AutoMinorLocator", (), {})
        ticker_stub.MaxNLocator = type(
            "MaxNLocator", (), {"__init__": lambda self, *a, **k: None}
        )

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
        sympy_stub.utilities = types.SimpleNamespace(
            lambdify=lambda *a, **k: None
        )
        fem = types.ModuleType("sympy.functions.elementary.miscellaneous")
        fem.sqrt = sqrt
        sys.modules["sympy"] = sympy_stub
        sys.modules["sympy.functions"] = types.ModuleType("sympy.functions")
        sys.modules["sympy.functions.elementary"] = types.ModuleType(
            "sympy.functions.elementary"
        )
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
        try:
            rec.fromarrays = np._core.records.fromarrays
        except Exception:  # fall back for older NumPy versions
            rec.fromarrays = np.core.records.fromarrays
        sys.modules["numpy.core.rec"] = rec
    sys.modules.setdefault(
        "pyspecdata.fornotebook", types.ModuleType("fornotebook")
    )
    fig_stub = sys.modules.setdefault(
        "pyspecdata.figlist", types.ModuleType("figlist")
    )
    if not hasattr(fig_stub, "figlist"):
        fig_stub.figlist = type("figlist", (), {})
    os.environ.setdefault("pyspecdata_figures", "standard")

    pkg = sys.modules.setdefault("pyspecdata", types.ModuleType("pyspecdata"))
    existing_paths = list(getattr(pkg, "__path__", []))
    if str(PACKAGE_ROOT) not in existing_paths:
        pkg.__path__ = existing_paths + [str(PACKAGE_ROOT)]

    spec = importlib.util.spec_from_file_location(
        f"pyspecdata.{name}", PACKAGE_ROOT / f"{name.replace('.', '/')}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[f"pyspecdata.{name}"] = module
    return module
