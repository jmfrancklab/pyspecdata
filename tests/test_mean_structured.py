import numpy as np
from conftest import load_module

load_module("general_functions")
load_module("ndshape")
core = load_module("core")
nddata = core.nddata


def _structured_test_data():
    dtype = np.dtype([("amp", "f8"), ("count", "i4")])
    data = np.array(
        [
            [(1.0, 1), (2.0, 3)],
            [(3.0, 5), (4.0, 7)],
            [(5.0, 9), (6.0, 11)],
        ],
        dtype=dtype,
    )
    return data


def test_mean_structured_sets_std_per_field():
    data = _structured_test_data()
    x = nddata(data.copy(), ["rep", "obs"])

    x.mean("rep", std=True)

    expected_mean = np.empty((2,), dtype=data.dtype)
    expected_std = np.empty((2,), dtype=data.dtype)
    for field in data.dtype.names:
        expected_mean[field] = np.mean(data[field], axis=0)
        expected_std[field] = np.std(data[field], axis=0)

    assert x.data.dtype == data.dtype
    assert x.get_error().dtype == data.dtype
    np.testing.assert_array_equal(x.data, expected_mean)
    np.testing.assert_array_equal(x.get_error(), expected_std)


def test_mean_structured_sets_stderr_per_field():
    data = _structured_test_data()
    x = nddata(data.copy(), ["rep", "obs"])

    x.mean("rep", stderr=True)

    expected_stderr = np.empty((2,), dtype=data.dtype)
    for field in data.dtype.names:
        expected_stderr[field] = np.std(data[field], axis=0) / np.sqrt(
            data.shape[0]
        )

    assert x.get_error().dtype == data.dtype
    np.testing.assert_array_equal(x.get_error(), expected_stderr)


def test_mean_structured_scalar_str():
    dtype = np.dtype([("field", "f8"), ("power", "f8")])
    data = np.zeros(10, dtype=dtype)
    data["field"] = np.arange(10, dtype=float)
    data["power"] = np.arange(10, 20, dtype=float)
    x = nddata(data, ["t"]).setaxis("t", np.r_[0:1:0.1])

    x.mean("t")

    assert x.data.dtype == dtype
    assert x.data.shape == ()
    assert str(x) == "(field=4.5000, power=14.500)"
