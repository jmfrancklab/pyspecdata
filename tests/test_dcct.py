from collections import OrderedDict
import sys
import types

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as s

from conftest import load_module

sys.modules.setdefault("_nnls", types.ModuleType("_nnls"))

core = load_module("core", use_real_pint=True)
pkg = sys.modules["pyspecdata"]
pkg.nddata = core.nddata
pkg.ndshape = core.ndshape
dcct_mod = load_module("plot_funcs.DCCT_function", use_real_pint=True)
fake_data_mod = load_module("generate_fake_data", use_real_pint=True)
nddata = core.nddata


def test_getitem_preserves_none_axis_coords_after_chunk():
    idx = nddata(np.r_[0:8], [-1], ["smooshed"])
    idx.chunk("smooshed", ["ph1", "ph2"], [2, 4])
    sliced = idx["ph1", 0]
    assert sliced.dimlabels == ["ph2"]
    assert sliced.getaxis("ph2") is None
    np.testing.assert_array_equal(sliced.data, np.r_[0:4])


def test_dcct_smoke_with_phase_dimensions():
    data = np.arange(2 * 4 * 3 * 16, dtype=float).reshape(2, 4, 3, 16)
    data = data + 1j * data[::-1]
    d = nddata(data, [2, 4, 3, 16], ["ph1", "ph2", "vd", "t2"])
    d.set_axis("ph1", np.r_[0:2] / 4)
    d.set_axis("ph2", np.r_[0:4] / 4)
    d.set_axis("vd", np.linspace(0.0, 1.0, 3))
    d.set_axis("t2", np.linspace(0.0, 15e-3, 16))
    d.set_units("t2", "s")
    d.ft("t2", shift=True).ft(["ph1", "ph2"])
    ax_list, _ = dcct_mod.DCCT(d)
    assert len(ax_list) == 8
    plt.close("all")


def test_fake_data_uses_scan_based_temp_units():
    t2, vd = s.symbols("t2 vd")
    data = fake_data_mod.fake_data(
        s.exp(-vd / 0.2) * s.exp(1j * 2 * s.pi * 10 * t2),
        OrderedDict(
            [
                ("vd", nddata(np.linspace(0.0, 1.0, 4), "vd")),
                ("t2", nddata(np.linspace(0.0, 10e-3, 16), "t2")),
            ]
        ),
        {},
    )
    assert data.get_units("t2") == "s"
    assert data.dimlabels == ["vd", "t2"]
