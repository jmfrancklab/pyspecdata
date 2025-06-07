import importlib.util
import numpy as np
import pytest
from conftest import load_module

# load dependencies via test loader (ensures real pint is imported)
gf = load_module("general_functions", use_real_pint=True)

# Skip this module entirely if pint is not available
if importlib.util.find_spec("pint") is None:
    pytest.skip("pint not installed", allow_module_level=True)

core = load_module("core", use_real_pint=True)
nddata = core.nddata


def test_voltage_times_current_yields_power_units():
    v = nddata(np.ones(3), "t")
    v.set_units("V")
    i = nddata(np.ones(3), "t")
    i.set_units("A")
    p = v * i
    assert gf.Q_(1, p.get_units()).to("W").units == gf.Q_(1, "W").units
