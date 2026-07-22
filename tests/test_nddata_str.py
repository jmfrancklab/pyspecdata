import numpy as np
import pytest
from conftest import load_module

load_module("general_functions")
load_module("ndshape")
nddata = load_module("core").nddata


def build_0d_plain():
    return nddata(np.array(1.25))


def build_0d_plain_with_error():
    x = nddata(np.array(1.25))
    x.set_error(np.array(0.05))
    return x


def build_0d_structured():
    data = np.zeros((), dtype=[("a", "f8"), ("b", "f8")])
    data["a"] = 1.25
    data["b"] = 2.5
    return nddata(data)


def build_0d_structured_with_error():
    data = np.zeros((), dtype=[("a", "f8"), ("b", "f8")])
    data["a"] = 1.25
    data["b"] = 2.5
    x = nddata(data)
    err = np.zeros((), dtype=[("a", "f8"), ("b", "f8")])
    err["a"] = 0.05
    err["b"] = 0.2
    x.set_error(err)
    return x


def build_1d_plain():
    x = nddata(np.array([1.0, 2.0, 3.0]), "t")
    x.setaxis("t", np.array([0.0, 1.0, 2.0]))
    return x


def build_1d_plain_with_errors():
    x = build_1d_plain()
    x.set_error(np.array([0.1, 0.2, 0.3]))
    x.set_error("t", np.array([0.01, 0.02, 0.03]))
    return x


def build_1d_structured_data():
    data = np.zeros(3, dtype=[("a", "f8"), ("b", "f8")])
    data["a"] = [1.0, 2.0, 3.0]
    data["b"] = [4.0, 5.0, 6.0]
    x = nddata(data, "t")
    x.setaxis("t", np.array([0.0, 1.0, 2.0]))
    return x


def build_1d_structured_data_with_errors():
    x = build_1d_structured_data()
    err = np.zeros(3, dtype=[("a", "f8"), ("b", "f8")])
    err["a"] = [0.1, 0.2, 0.3]
    err["b"] = [0.4, 0.5, 0.6]
    x.set_error(err)
    x.set_error("t", np.array([0.01, 0.02, 0.03]))
    return x


def build_1d_structured_axis():
    x = nddata(np.array([1.0, 2.0, 3.0]), "t")
    axis = np.zeros(3, dtype=[("f", "f8"), ("p", "f8")])
    axis["f"] = [10.0, 20.0, 30.0]
    axis["p"] = [1.0, 2.0, 3.0]
    x.setaxis("t", axis)
    return x


def build_1d_structured_axis_with_error():
    x = build_1d_structured_axis()
    axiserr = np.zeros(3, dtype=[("f", "f8"), ("p", "f8")])
    axiserr["f"] = [0.1, 0.2, 0.3]
    axiserr["p"] = [0.01, 0.02, 0.03]
    x.set_error("t", axiserr)
    return x


@pytest.mark.parametrize(
    "builder,expected",
    [
        pytest.param(build_0d_plain, "1.2500", id="0d-plain"),
        pytest.param(
            build_0d_plain_with_error,
            "1.25 ± 0.050",
            id="0d-plain-data-error",
        ),
        pytest.param(
            build_0d_structured,
            "(a=1.2500, b=2.5000)",
            id="0d-structured-data",
        ),
        pytest.param(
            build_0d_structured_with_error,
            "(a=1.25 ± 0.050, b=2.5 ± 0.20)",
            id="0d-structured-data-error",
        ),
        pytest.param(
            build_1d_plain,
            "array([1., 2., 3.])\n"
            "\tdimlabels=['t']\n"
            "\taxes={`t':array([0., 1., 2.])\n"
            "\t\t\t+/-None}\n",
            id="1d-plain",
        ),
        pytest.param(
            build_1d_plain_with_errors,
            "array([1., 2., 3.])\n"
            "\t\t±array([0.1, 0.2, 0.3])\n"
            "\tdimlabels=['t']\n"
            "\taxes={`t':array([0., 1., 2.])\n"
            "\t\t\t+/-array([0.01, 0.02, 0.03])}\n",
            id="1d-plain-data-and-coord-error",
        ),
        pytest.param(
            build_1d_structured_data,
            "array([(1., 4.), (2., 5.), (3., 6.)], "
            "dtype=[('a', '<f8'), ('b', '<f8')])\n"
            "\tdimlabels=['t']\n"
            "\taxes={`t':array([0., 1., 2.])\n"
            "\t\t\t+/-None}\n",
            id="1d-structured-data",
        ),
        pytest.param(
            build_1d_structured_data_with_errors,
            "array([(1., 4.), (2., 5.), (3., 6.)], "
            "dtype=[('a', '<f8'), ('b', '<f8')])\n"
            "\t\t±array([(0.1, 0.4), (0.2, 0.5), (0.3, 0.6)],\n"
            "      dtype=[('a', '<f8'), ('b', '<f8')])\n"
            "\tdimlabels=['t']\n"
            "\taxes={`t':array([0., 1., 2.])\n"
            "\t\t\t+/-array([0.01, 0.02, 0.03])}\n",
            id="1d-structured-data-and-errors",
        ),
        pytest.param(
            build_1d_structured_axis,
            "array([1., 2., 3.])\n"
            "\tdimlabels=['t']\n"
            "\taxes={`t':array([(10., 1.), (20., 2.), (30., 3.)],\n"
            "      dtype=[('f', '<f8'), ('p', '<f8')])\n"
            "\t\t\t+/-None}\n",
            id="1d-structured-axis",
        ),
        pytest.param(
            build_1d_structured_axis_with_error,
            "array([1., 2., 3.])\n"
            "\tdimlabels=['t']\n"
            "\taxes={`t':array([(10., 1.), (20., 2.), (30., 3.)],\n"
            "      dtype=[('f', '<f8'), ('p', '<f8')])\n"
            "\t\t\t+/-array([(0.1, 0.01), (0.2, 0.02), (0.3, 0.03)],\n"
            "      dtype=[('f', '<f8'), ('p', '<f8')])}\n",
            id="1d-structured-axis-error",
        ),
    ],
)
def test_nddata_str(builder, expected):
    assert str(builder()) == expected
