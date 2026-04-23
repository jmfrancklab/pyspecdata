import numpy as np
from conftest import load_module

load_module("general_functions")
load_module("ndshape")
core = load_module("core")
nddata = core.nddata


def test_mean_structured_sets_std_per_field():
    data = np.array(
        [
            [(1.0, 1), (2.0, 3)],
            [(3.0, 5), (4.0, 7)],
            [(5.0, 9), (6.0, 11)],
        ],
        dtype=[("amp", "f8"), ("count", "i4")],
    )
    x = nddata(data.copy(), ["rep", "obs"])

    x.mean("rep", std=True)

    expected_dtype = np.dtype([("amp", "f8"), ("count", "f8")])
    expected_mean = np.empty((2,), dtype=expected_dtype)
    expected_std = np.empty((2,), dtype=expected_dtype)
    for field in data.dtype.names:
        expected_mean[field] = np.mean(data[field], axis=0)
        expected_std[field] = np.std(data[field], axis=0)

    assert x.data.dtype == expected_dtype
    assert x.get_error().dtype == expected_dtype
    np.testing.assert_array_equal(x.data, expected_mean)
    np.testing.assert_array_equal(x.get_error(), expected_std)


def test_mean_structured_sets_stderr_per_field():
    data = np.array(
        [
            [(1.0, 1), (2.0, 3)],
            [(3.0, 5), (4.0, 7)],
            [(5.0, 9), (6.0, 11)],
        ],
        dtype=[("amp", "f8"), ("count", "i4")],
    )
    x = nddata(data.copy(), ["rep", "obs"])

    x.mean("rep", stderr=True)

    expected_dtype = np.dtype([("amp", "f8"), ("count", "f8")])
    expected_stderr = np.empty((2,), dtype=expected_dtype)
    for field in data.dtype.names:
        expected_stderr[field] = np.std(data[field], axis=0) / np.sqrt(
            data.shape[0]
        )

    assert x.get_error().dtype == expected_dtype
    np.testing.assert_array_equal(x.get_error(), expected_stderr)


def test_mean_structured_preserves_complex_fields():
    data = np.array(
        [
            [(1.0 + 1.0j, 1), (2.0 + 2.0j, 3)],
            [(3.0 + 3.0j, 5), (4.0 + 4.0j, 7)],
            [(5.0 + 5.0j, 9), (6.0 + 6.0j, 11)],
        ],
        dtype=[("signal", "c16"), ("count", "i4")],
    )
    x = nddata(data.copy(), ["rep", "obs"])

    x.mean("rep", std=True)

    expected_dtype = np.dtype([("signal", "c16"), ("count", "f8")])
    expected_mean = np.empty((2,), dtype=expected_dtype)
    expected_std = np.empty((2,), dtype=expected_dtype)
    for field in data.dtype.names:
        expected_mean[field] = np.mean(data[field], axis=0)
        expected_std[field] = np.std(data[field], axis=0)

    assert x.data.dtype == expected_dtype
    assert x.get_error().dtype == expected_dtype
    np.testing.assert_array_equal(x.data, expected_mean)
    np.testing.assert_array_equal(x.get_error(), expected_std)
