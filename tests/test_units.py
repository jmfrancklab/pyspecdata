import numpy as np
from conftest import load_module

# load dependencies via test loader
load_module("general_functions")
core = load_module("core")
nddata = core.nddata


def test_voltage_times_current_yields_power_units():
    v = nddata(np.ones(3), "t")
    v.set_units("V")
    i = nddata(np.ones(3), "t")
    i.set_units("A")
    p = v * i
    gf = load_module("general_functions")
    Q_ = gf.Q_
    assert Q_(1, p.get_units()).to("W").units == gf.ureg.watt
