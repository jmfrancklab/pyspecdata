import numpy as np
import pytest
from conftest import load_module

axis_mod = load_module("axis_class")


def test_add_axis_via_builder():
    axes = axis_mod.axis_collection(["t2"])
    axes["t2"] = axis_mod.ax_[0:1:10j]
    arr = axes["t2"].to_array()
    assert arr.size == 10
    assert np.isclose(arr[0], 0)
    assert np.isclose(arr[-1], 1)


def test_rename_axis_updates_dimlabels():
    axes = axis_mod.axis_collection(["t2"])
    axes["t2"] = axis_mod.ax_[0:1:10j]
    axes.rename("t2", "time")
    assert "time" in axes.dimlabels
    assert "t2" not in axes.dimlabels


def test_add_existing_axis_requires_matching_dimlabel():
    axes = axis_mod.axis_collection(["t2"])
    existing = axis_mod.ax_[0:1:10j]
    existing.names = ["t2"]
    axes += existing
    assert axes["t2"] is existing
    other = axis_mod.ax_[0:1:10j]
    other.names = ["other"]
    with pytest.raises(ValueError):
        axes += other


def test_no_duplicate_axis_names():
    axes = axis_mod.axis_collection(["t2", "f2"])
    axes["t2"] = axis_mod.ax_[0:1:10j]
    duplicate = axis_mod.ax_[0:1:10j]
    duplicate.names = ["t2"]
    with pytest.raises(ValueError):
        axes["f2"] = duplicate
