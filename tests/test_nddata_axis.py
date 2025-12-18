import numpy as np
import pytest
from conftest import load_module

axis_mod = load_module("axis_class")
core_mod = load_module("core")


def test_ax_slice_notation_with_complex_step_matches_linspace():
    # complex step should be interpreted as a number of points like numpy.r_
    axis = axis_mod.ax_["t", 0:1.2:100j]
    assert axis.size == 100
    assert axis.start == 0
    assert np.isclose(axis.dx, 1.2 / 99)
    assert np.allclose(axis.to_array(), np.linspace(0, 1.2, 100))


def test_ax_slice_notation_with_float_step_matches_arange():
    # real-valued step should mirror numpy.arange spacing
    axis = axis_mod.ax_["t", 0:2:0.5]
    assert axis.size == 4
    assert axis.start == 0
    assert np.isclose(axis.dx, 0.5)
    assert np.allclose(axis.to_array(), np.arange(0, 2, 0.5))


def test_nddata_axis_supports_position_indexing():
    # indexing by integer should return the coordinate at that position
    axis = axis_mod.ax_["pos", 0:2:5j]
    base = axis.to_array()
    assert axis[2] == base[2]


def test_nddata_axis_slice_by_coordinate_range():
    # slicing with a tuple should keep coordinates within the inclusive bounds
    axis = axis_mod.ax_["pos", 0:10:6j]
    sliced = axis[(0, 5.5)]
    assert np.allclose(sliced.to_array(), np.array([0, 2, 4]))


def test_nddata_axis_multiplies_numpy_arrays():
    # arithmetic with numpy arrays should behave like elementwise coordinate
    # scaling
    axis = axis_mod.ax_["x", 1:5:1]
    doubled = axis * 2
    assert np.allclose(doubled, np.array([2, 4, 6, 8]))
    ones = np.ones(4)
    product_left = axis * ones
    product_right = ones * axis
    assert np.allclose(product_left, np.array([1, 2, 3, 4]))
    assert np.allclose(product_right, np.array([1, 2, 3, 4]))


def test_nddata_axis_multiplies_nddata_by_axis_coordinates():
    # multiplying an nddata by an axis should scale along the matching
    # dimension
    data = core_mod.nddata(np.ones(4), "time")
    axis = axis_mod.ax_["time", 10:14]
    axis.names = ["time"]
    scaled = axis * data
    assert isinstance(scaled, core_mod.nddata)
    assert np.allclose(scaled.data, np.arange(10, 14))


def test_nddata_axis_interpolates_nddata_with_matmul():
    # using @ should interpolate the nddata onto the new axis coordinates
    data = core_mod.nddata(np.array([0.0, 1.0, 2.0, 3.0]), ["time"])
    data.setaxis("time", np.array([0.0, 1.0, 2.0, 3.0]))
    target_axis = axis_mod.ax_["time", 0:3:6j]
    target_axis.names = ["time"]
    interpolated = data @ target_axis
    assert isinstance(interpolated, core_mod.nddata)
    assert np.allclose(interpolated.getaxis("time"), target_axis.to_array())
    assert interpolated.data.shape[0] == target_axis.size
    direct_interp = data.interp("time", target_axis.to_array())
    assert np.allclose(interpolated.data, direct_interp.data)


def test_function_application_stores_transformation_on_axis():
    # passing the axis to a numpy ufunc should record the transformation used
    axis = axis_mod.ax_["t", 0:2:5j]
    transformed = np.exp(axis)
    assert isinstance(transformed, axis_mod.nddata_axis)
    assert transformed.transf_func is np.exp
    assert np.allclose(transformed.to_array(), np.exp(axis.to_array()))


def test_axis_arithmetic_two_names_yields_grid_ndarray():
    # combining two differently named axes should form a 2D grid (x first, y
    # second) matching fromaxis behavior
    x_axis = axis_mod.ax_["x", 0:2:3j]
    y_axis = axis_mod.ax_["y", 1:3:3j]
    product = x_axis * y_axis
    expected = x_axis.data[:, None] * y_axis.data[None, :]
    assert isinstance(product, np.ndarray)
    assert product.shape == expected.shape
    assert np.allclose(product, expected)
    quad_expr = x_axis**2 + 2 * y_axis + 1
    quad_expected = (x_axis.data[:, None] ** 2) + 2 * y_axis.data[None, :] + 1
    assert np.allclose(quad_expr, quad_expected)
    data = core_mod.nddata(np.zeros(expected.shape), ["x", "y"])
    data.setaxis("x", x_axis.data)
    data.setaxis("y", y_axis.data)
    assert np.allclose(
        product, data.fromaxis(["x", "y"], lambda x, y: x * y, as_array=True)
    )
    assert np.allclose(
        quad_expr,
        data.fromaxis(
            ["x", "y"], lambda x, y: x**2 + 2 * y + 1, as_array=True
        ),
    )


def test_axis_arithmetic_same_name_returns_axis_object():
    # combining axes with the same name should return a new axis coordinate
    # object
    first = axis_mod.ax_["t", 0:1:4j]
    second = axis_mod.ax_["t", 1:2:4j]
    summed = first + second
    assert isinstance(summed, axis_mod.nddata_axis)
    assert summed.names[0] == "t"
    assert np.allclose(summed.to_array(), first.to_array() + second.to_array())


def test_ax_helper_without_name_raises_error():
    # ax_ must always be called with a name before the slice
    with pytest.raises(ValueError):
        _ = axis_mod.ax_[0:1:3j]
