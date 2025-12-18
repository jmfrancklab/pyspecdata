import numpy as np
import pytest
from conftest import load_module

axis_mod = load_module("axis_class")


def test_axis_collection_allows_slice_notation_assignment():
    # adding an axis with numpy r_ style notation should populate the collection
    axes = axis_mod.axis_collection()
    axes += axis_mod.ax_["t2", 0:1.2:5j]
    assert axes.dimlabels == ["t2"]
    assert axes["t2"].names[0] == "t2"
    assert np.allclose(axes["t2"].to_array(), np.linspace(0, 1.2, 5))


def test_axis_collection_rejects_duplicate_aliases():
    # collection should stop two axes from sharing any alias
    axes = axis_mod.axis_collection()
    axes += axis_mod.ax_["f", 0:1:3j]
    axes["f"].names = ["f", "freq"]
    new_axis = axis_mod.ax_["t", 0:1:3j]
    new_axis.names = ["t", "freq"]
    with pytest.raises(ValueError):
        axes["t"] = new_axis


def test_axis_collection_requires_matching_dimlabel_on_add():
    # adding an existing axis should require one of its names to match a dimlabel
    axes = axis_mod.axis_collection()
    axes.dimlabels = ["t", "f"]
    existing_axis = axis_mod.ax_["t", 0:2:5j]
    existing_axis.names = ["t", "time"]
    axes += existing_axis
    mismatched_axis = axis_mod.ax_["z", 0:1:3j]
    mismatched_axis.names = ["z"]
    with pytest.raises(ValueError):
        axes += mismatched_axis
    matched_axis = axis_mod.ax_["f", 0:1:3j]
    matched_axis.names = ["f"]
    axes += matched_axis
    assert axes["f"] is matched_axis


def test_axis_collection_rename_updates_dimlabels_and_axis():
    # rename should relabel both dimlabels and the preferred axis name
    axes = axis_mod.axis_collection()
    axes += axis_mod.ax_["t2", 0:1:4j]
    axes["t2"].names = ["t2", "direct"]
    renamed = axes.rename("t2", "time")
    assert renamed is axes
    assert axes.dimlabels == ["time"]
    assert axes["time"].names[0] == "time"


def test_axis_collection_supports_append_api():
    # append should be equivalent to using the in-place add syntax
    axes = axis_mod.axis_collection()
    axes.append(axis_mod.ax_["f", 0:1:3j])
    assert "f" in axes
    axes.append(axis_mod.ax_["t", 0:2:2j])
    assert set(axes.dimlabels) == {"f", "t"}


def test_ax_helper_requires_name_when_slicing():
    # constructing an axis without a name should raise a clear error
    with pytest.raises(ValueError):
        _ = axis_mod.ax_[0:1:3j]
