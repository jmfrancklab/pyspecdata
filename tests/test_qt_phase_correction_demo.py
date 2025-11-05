"""Tests for the Qt phase correction example using generated data."""
import os
import sys
from pathlib import Path

import importlib.util

import numpy as np
import pytest
from PyQt5.QtWidgets import QApplication

# Ensure the repository's pyspecdata package is importable without installation.
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

# Provide a minimal stub for the optional compiled nnls extension.
sys.modules.setdefault("_nnls", type(sys)("_nnls"))

from pyspecdata import nddata

# Load the Qt demo module directly from its file path.
demo_path = repo_root / "examples" / "qt_phase_correction_demo.py"
spec = importlib.util.spec_from_file_location(
    "qt_phase_correction_demo", demo_path
)
qt_demo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qt_demo)
PhaseCorrectionWidget = qt_demo.PhaseCorrectionWidget
apply_phase_corrections = qt_demo.apply_phase_corrections


def generate_fake_dataset():
    """Create a simple hypercomplex dataset with μs axes for testing."""
    num_points = 1000
    t_axis = np.linspace(0.0, 4.0, num_points)
    # Construct a smooth complex surface so the plots have meaningful content.
    t1_grid, t2_grid = np.meshgrid(t_axis, t_axis, indexing="ij")
    base_surface = np.exp(-((t1_grid - 2.0) ** 2 + (t2_grid - 2.0) ** 2))
    data = np.empty((2, num_points, num_points), dtype=np.complex128)
    data[0] = base_surface * np.exp(1j * (t1_grid + t2_grid))
    data[1] = base_surface * np.exp(-1j * (t1_grid - t2_grid))
    dataset = nddata(data, ["ph1", "t1", "t2"])
    dataset.setaxis("ph1", np.r_[0.0, 0.5])
    dataset.setaxis("t1", t_axis).set_units("t1", "μs")
    dataset.setaxis("t2", t_axis).set_units("t2", "μs")
    return dataset


@pytest.fixture(scope="module")
def qapp():
    """Provide a QApplication in offscreen mode for the Qt tests."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    yield app


def test_apply_phase_corrections_returns_expected_shapes():
    """The helper should keep axis dimensions intact after corrections."""
    dataset = generate_fake_dataset()
    time_domain, corrected = apply_phase_corrections(dataset, 0.1, -0.15)
    assert time_domain.shape["t1"] >= dataset.shape["t1"]
    assert corrected.dimlabels == dataset.dimlabels
    for axis in dataset.dimlabels:
        assert corrected.shape[axis] == dataset.shape[axis]


def test_widget_updates_with_fake_data(qapp):
    """Changing slider positions should update the labels without errors."""
    dataset = generate_fake_dataset()
    widget = PhaseCorrectionWidget(base_dataset=dataset)
    try:
        initial_time_image = np.array(widget.time_image.get_array())
        start_axes_count = len(widget.figure.axes)
        widget.diag_slider.setValue(25)
        widget.anti_slider.setValue(-30)
        qapp.processEvents()
        assert np.array_equal(initial_time_image, np.array(widget.time_image.get_array()))
        widget.diag_slider.sliderReleased.emit()
        qapp.processEvents()
        assert widget.diag_label.text() == "0.25 μs"
        assert widget.anti_label.text() == "-0.30 μs"
        updated_time_image = np.array(widget.time_image.get_array())
        assert not np.array_equal(initial_time_image, updated_time_image)
        assert len(widget.figure.axes) == start_axes_count
        assert widget.ax_time.images
        assert widget.ax_scp.images
        assert widget.ax_scm.images
    finally:
        widget.close()


def test_sensitivity_slider_adjusts_phase_ranges(qapp):
    """The sensitivity slider should expand and contract the phase ranges."""
    dataset = generate_fake_dataset()
    widget = PhaseCorrectionWidget(base_dataset=dataset)
    try:
        center_before = widget.diag_slider.value()
        range_before = widget.diag_slider.maximum() - widget.diag_slider.minimum()
        widget.sensitivity_slider.setValue(widget.sensitivity_slider.value() + 50)
        widget.sensitivity_slider.sliderReleased.emit()
        qapp.processEvents()
        center_after = widget.diag_slider.value()
        range_after = widget.diag_slider.maximum() - widget.diag_slider.minimum()
        assert center_after == center_before
        assert range_after > range_before
        assert (
            widget.diag_slider.maximum() - center_after
            == center_after - widget.diag_slider.minimum()
        )
    finally:
        widget.close()
