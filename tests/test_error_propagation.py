import numpy as np
from conftest import load_module

load_module("general_functions")
load_module("ndshape")
core = load_module("core")
nddata = core.nddata


def test_error_propagation_basic():
    a = nddata(np.array([2.0, 3.0]), "x")
    a.set_error(np.array([0.1, 0.1]))
    b = nddata(np.array([1.0, 2.0]), "x")
    b.set_error(np.array([0.2, 0.1]))

    c = a + b
    expected_add = np.sqrt(a.get_error() ** 2 + b.get_error() ** 2)
    assert np.allclose(c.get_error(), expected_add)

    d = a * b
    expected_mul = np.sqrt(
        (a.get_error() * b.data) ** 2 + (b.get_error() * a.data) ** 2
    )
    assert np.allclose(d.get_error(), expected_mul)

    e = a / b
    expected_div = np.sqrt(
        (a.get_error() / b.data) ** 2
        + (a.data * b.get_error() / (b.data**2)) ** 2
    )
    assert np.allclose(e.get_error(), expected_div)
