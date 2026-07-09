import numpy as np
import pytest
from conftest import load_module

# Skip if the real Pint library isn't available. ``load_module`` supplies a
# minimal stub when Pint can't be imported, which lacks the functionality this
# test requires.
try:
    import pint
except ModuleNotFoundError:
    pytest.skip("pint not installed", allow_module_level=True)
if getattr(pint, "__pyspec_stub__", False):
    pytest.skip("pint stub in use", allow_module_level=True)

# Load dependencies via the test loader (ensures real Pint is imported)
gf = load_module("general_functions", use_real_pint=True)
core = load_module("core", use_real_pint=True)
nddata = core.nddata


def test_voltage_times_current_yields_power_units():
    v = nddata(np.ones(3), "t")
    v.set_units("V")
    i = nddata(np.ones(3), "t")
    i.set_units("A")
    p = v * i
    assert gf.Q_(1, p.get_units()).to("W").units == gf.Q_(1, "W").units


def test_fourier_transform_updates_time_units_to_frequency():
    d = nddata(np.ones(8), "t")
    d.set_axis("t", np.linspace(0.0, 7.0, 8))
    d.set_units("t", "µs")
    d.ft("t")
    assert d.get_units("t") == "MHz"


def test_scan_unit_is_registered():
    quantity = gf.Q_(1, "cyc / scan")
    assert quantity.check(gf.Q_(1, "cyc") / gf.Q_(1, "scan"))


def test_compact_square_root_unit_is_parseable():
    assert gf.Q_(1, "s√W").check(gf.Q_(1, "s * W**0.5"))


def test_pretty_unit_string_is_parseable():
    pretty = gf.Q_(1, "g⁰⋅⁵·µm/s⁰⋅⁵")
    parseable = gf.Q_(1, "g**0.5 * um / s**0.5")
    assert pretty.check(parseable)
    assert gf.det_unit_prefactor("g⁰⋅⁵·µm/s⁰⋅⁵") == -9


def test_div_units_accepts_pretty_unit_labels():
    d = nddata(np.ones(2), "beta")
    d.setaxis("beta", np.r_[0:2])
    d.set_units("beta", "g⁰⋅⁵·µm/s⁰⋅⁵")
    assert np.isclose(
        d.div_units("beta", "s * W**0.5"), np.sqrt(1e-3) * 1e-6
    )


def test_inverse_fourier_transform_accepts_cycles_per_scan_units():
    d = nddata(np.ones(8), "temp")
    d.set_axis("temp", np.linspace(-0.5, 0.5, 8, endpoint=False))
    d.set_units("temp", "cyc / scan")
    d.ift("temp")
    assert d.get_units("temp") == "scan"
