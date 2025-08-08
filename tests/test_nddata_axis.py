import numpy as np
from conftest import load_module

axis_mod = load_module("axis_class")
ns = load_module("core")


def build_axis():
    """Helper constructing a time axis."""
    return axis_mod.ax_["t", 0:10:11j]


def test_axis_builder_sets_name():
    axis = axis_mod.ax_["t", 0:5:6j]
    assert getattr(axis, "name", None) == "t"


def test_indexing_returns_value():
    a = build_axis()
    assert np.isclose(a[1], 1.0)


def test_slicing_returns_axis():
    a = build_axis()
    sliced = a[0:5:2]
    assert np.allclose(sliced.to_array(), np.r_[0:5:2])


def test_nddata_like_slicing():
    a = build_axis()
    between = a[(0, 5.5)]
    assert np.allclose(between.to_array(), np.r_[0:6])
    above = a[(3, None)]
    assert np.allclose(above.to_array(), np.r_[3:11])


def test_multiplication_with_numpy_array():
    a = build_axis()
    array = np.ones(a.size)
    result = a * array
    assert isinstance(result, np.ndarray)
    assert np.allclose(result, a.to_array())


def test_multiplication_with_nddata():
    a = build_axis()
    data = ns.nddata(np.ones(a.size), ["t"])
    result = a * data
    assert isinstance(result, ns.nddata)
    assert np.allclose(result.data, a.to_array())


def test_multiplication_with_nddata_diffdim():
    a = build_axis()
    data = ns.nddata(np.ones(3), ["s"])
    result = a * data
    assert isinstance(result, ns.nddata)
    assert result.dimlabels == ["t", "s"]
    expected = np.outer(a.to_array(), data.data)
    assert np.allclose(result.data, expected)


def test_multiplication_between_axes_diffdim():
    a = build_axis()
    b = axis_mod.ax_["s", 0:5:6j]
    result = a * b
    assert isinstance(result, ns.nddata)
    assert result.dimlabels == ["t", "s"]
    expected = np.outer(a.to_array(), b.to_array())
    assert np.allclose(result.data, expected)


def test_function_application_returns_axis():
    a = build_axis()
    new_axis = np.exp(a)
    assert isinstance(new_axis, axis_mod.nddata_axis)
    assert np.allclose(new_axis.to_array(), np.exp(a.to_array()))


def test_interpolation_operator():
    a = build_axis()
    data = ns.nddata(np.sin(a.to_array()), ["t"])
    new_axis = axis_mod.ax_["t", 0:10:21j]
    interp = data @ new_axis
    assert isinstance(interp, ns.nddata)
    assert np.allclose(interp.axis("t").to_array(), new_axis.to_array())
    expected = np.interp(new_axis.to_array(), a.to_array(), data.data)
    assert np.allclose(interp.data, expected)
