import numpy as np
import pytest
from conftest import load_module

axis_mod = load_module("axis_class")


def test_add_axis_via_builder():
    axes = axis_mod.axis_collection()
    t_axis = axis_mod.ax_["t2", 0:1:10j]
    axes += t_axis
    arr = axes["t2"].to_array()
    assert arr.size == 10
    assert np.isclose(arr[0], 0)
    assert np.isclose(arr[-1], 1)


def test_getitem_returns_axis():
    axes = axis_mod.axis_collection()
    t_axis = axis_mod.ax_["t", 0:5:6j]
    axes += t_axis
    assert axes["t"] is t_axis


def test_no_duplicate_axis_names():
    axes = axis_mod.axis_collection()
    axes += axis_mod.ax_["t", 0:1:10j]
    with pytest.raises(ValueError):
        axes += axis_mod.ax_["t", 0:1:10j]


def test_alias_conflicts_with_name():
    axes = axis_mod.axis_collection()
    a = axis_mod.ax_["t", 0:1:10j]
    axes += a
    b = axis_mod.ax_["f", 0:1:10j]
    b.aliases = ["t"]
    with pytest.raises(ValueError):
        axes += b


def test_no_duplicate_aliases():
    axes = axis_mod.axis_collection()
    a = axis_mod.ax_["t", 0:1:10j]
    a.aliases = ["time"]
    axes += a
    b = axis_mod.ax_["f", 0:1:10j]
    b.aliases = ["time"]
    with pytest.raises(ValueError):
        axes += b
