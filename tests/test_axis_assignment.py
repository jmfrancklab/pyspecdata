import numpy as np
import pytest
from conftest import load_module

load_module("general_functions")
load_module("ndshape")
core = load_module("core")
nddata = core.nddata


def build_axis_test_data():
    d = nddata(np.array([10.0, 20.0, 30.0]), "t")
    d.setaxis("t", np.array([0.0, 1.0, 2.0]))
    d.set_error("t", np.array([0.1, 0.2, 0.3]))
    d.set_units("t", "s")
    return d


def test_direct_axis_assignment_clears_axis_error():
    d = build_axis_test_data()

    d["t"] = np.array([5.0, 6.0, 7.0])

    np.testing.assert_allclose(d.getaxis("t"), np.array([5.0, 6.0, 7.0]))
    assert d.get_error("t") is None
    assert d.get_units("t") == "s"


def test_inplace_axis_assignment_preserves_axis_error():
    d = build_axis_test_data()

    d["t"] += 10.0

    np.testing.assert_allclose(d.getaxis("t"), np.array([10.0, 11.0, 12.0]))
    np.testing.assert_allclose(d.get_error("t"), np.array([0.1, 0.2, 0.3]))
    assert d.get_units("t") == "s"


def test_nddata_axis_assignment_copies_error_and_units():
    d = build_axis_test_data()
    rhs = nddata(np.array([100.0, 200.0, 300.0]), "t")
    rhs.set_error(np.array([1.0, 2.0, 3.0]))
    rhs.set_units("ms")

    d["t"] = rhs

    np.testing.assert_allclose(d.getaxis("t"), rhs.data)
    np.testing.assert_allclose(d.get_error("t"), np.array([1.0, 2.0, 3.0]))
    assert d.get_units("t") == "ms"


@pytest.mark.parametrize(
    "rhs,match",
    [
        pytest.param(
            lambda: nddata(np.ones((3, 1)), [3, 1], ["t", "u"]),
            "exactly one dimension",
            id="too-many-dimensions",
        ),
        pytest.param(
            lambda: nddata(np.array([1.0, 2.0, 3.0]), "u"),
            "rhs dimension must also be",
            id="wrong-dimension-name",
        ),
    ],
)
def test_nddata_axis_assignment_validates_dimension(rhs, match):
    d = build_axis_test_data()

    with pytest.raises(ValueError, match=match):
        d["t"] = rhs()
