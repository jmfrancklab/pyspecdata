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
