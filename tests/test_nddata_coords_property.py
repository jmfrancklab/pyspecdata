import numpy as np
import pytest
from conftest import load_module

axis_mod = load_module("axis_class")
core_mod = load_module("core")


def test_setaxis_with_numpy_array_creates_axis_object():
    # setting an axis from raw coordinates should wrap them in an axis object
    data = core_mod.nddata(np.zeros(5), ["t"])
    coord_values = np.linspace(0.0, 1.0, 5)
    data.setaxis("t", coord_values)
    assert isinstance(data.coords["t"], axis_mod.nddata_axis)
    assert np.allclose(data.coords["t"].to_array(), coord_values)


def test_setaxis_accepts_axis_object_directly():
    # passing an axis object to setaxis should install it without conversion
    axis = axis_mod.ax_["t", 0:2:5j]
    data = core_mod.nddata(np.zeros(5), ["t"])
    data.setaxis("t", axis)
    assert data.coords["t"] is axis


def test_coords_property_accepts_dict_of_arrays_and_builds_collection():
    # assigning a dict of arrays to coords should create an axis collection
    data = core_mod.nddata(np.zeros((2, 3)), ["x", "y"])
    axes_dict = {"x": np.array([1.0, 2.0]), "y": np.array([3.0, 4.0, 5.0])}
    data.coords = axes_dict
    assert isinstance(data.coords, axis_mod.axis_collection)
    assert np.allclose(data.coords["x"].to_array(), axes_dict["x"])
    assert np.allclose(data.coords["y"].to_array(), axes_dict["y"])


def test_labels_uses_coords_property_for_backwards_support():
    # legacy labels API should still populate the coords collection
    data = core_mod.nddata(np.zeros((2, 2)), ["a", "b"])
    data.labels({"a": np.array([0.0, 1.0]), "b": np.array([10.0, 20.0])})
    assert isinstance(data.coords, axis_mod.axis_collection)
    assert np.allclose(data.coords["a"].to_array(), np.array([0.0, 1.0]))
    assert np.allclose(data.coords["b"].to_array(), np.array([10.0, 20.0]))


def test_getaxis_reads_from_coords_collection():
    # getaxis should reflect the coordinate object stored in coords
    axis = axis_mod.ax_["x", 0:3:3j]
    data = core_mod.nddata(np.zeros(3), ["x"])
    data.coords = {"x": axis}
    assert np.allclose(data.getaxis("x"), axis.to_array())


def test_coords_units_reference_axis_units_object():
    # units reported by the collection should point back to the axis units
    axis = axis_mod.ax_["x", 0:1:2j]
    axis.units = object()
    axes = axis_mod.axis_collection()
    axes += axis
    assert axes.units["x"] is axis.units


def test_axis_coords_error_replaced_by_coords_property():
    # explicit axis_coords_error should be unnecessary once coords stores axis objects
    data = core_mod.nddata(np.zeros(3), ["x"])
    data.coords = {"x": np.arange(3)}
    with pytest.raises(AttributeError):
        _ = data.axis_coords_error
