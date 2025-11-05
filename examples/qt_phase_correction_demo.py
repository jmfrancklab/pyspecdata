"""Qt example showing diagonal and anti-diagonal corrections with sliders."""
import configparser
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import (
    NavigationToolbar2QT,
    new_figure_manager_given_figure,
)
from matplotlib.figure import Figure
from numpy import pi
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QApplication,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from pyspecdata import find_file, image
from pyspecdata.datadir import pyspec_config
from pyspecdata.plot_funcs.image import imagehsv


mpl.rcParams.update(
    {
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
    }
)


mpl.colormaps.register(
    name="custom_diverging",
    cmap=mpl.colors.LinearSegmentedColormap.from_list(
        "custom_diverging",
        [
            (0.0, "black"),
            (0.25, "blue"),
            (0.5, "white"),
            (0.75, "red"),
            (1.0, "black"),
        ],
    ),
    force=True,
)
mpl.colormaps.register(
    name="single_sided",
    cmap=mpl.colors.LinearSegmentedColormap.from_list(
        "single_sided",
        [
            (0.0, "white"),
            (0.2, "black"),
            (0.4, "blue"),
            (0.6, "white"),
            (0.8, "red"),
            (1.0, "black"),
        ],
    ),
    force=True,
)
mpl.rcParams["image.cmap"] = "custom_diverging"


T1_RANGE = (-60, 60)
T2_RANGE = (-60, 60)
DEFAULT_EXP_TYPE = "2DELDOR/251027_DiPhNO_500uM/Hypercomplex"
DEFAULT_SEARCH = "1500ns\\.DSC"


# {{{ data preparation helpers
# Load the dataset once so that slider updates only adjust phase factors.
def load_base_dataset(exp_type, search_string):
    print(f"load_base_dataset called with {exp_type} {search_string}")
    d = find_file(search_string, exp_type=exp_type)
    d.chunk("t1", ["t1", "ph1"], [-1, 2]).reorder(["ph1", "t1", "t2"])
    d.set_units("t1", "ns")
    for axis in ["t1", "t2"]:
        d[axis] *= d.div_units(axis, "μs")
        d.set_units(axis, "μs")
    d = d["t2":T2_RANGE]
    d["ph1", 1] *= 1j
    d.setaxis("ph1", np.r_[0, 0.5]).ft("ph1")
    d.ft(["t1", "t2"], shift=True)
    d *= -1
    return d


def apply_phase_corrections(base_dataset, diag_corr, anti_corr):
    corrected = base_dataset.copy()
    # Apply diagonal phase adjustment shared by both frequency axes.
    corrected *= np.exp(
        -1j
        * 2
        * pi
        * diag_corr
        * (corrected.fromaxis("t1") + corrected.fromaxis("t2"))
        / 2
    )
    # Apply anti-diagonal phase adjustment for the relative phase between axes.
    corrected *= np.exp(
        -1j
        * 2
        * pi
        * anti_corr
        * (corrected.fromaxis("t1") - corrected.fromaxis("t2"))
        / 2
    )
    # Reconstruct the time-domain Hermitian combination.
    together = corrected.C.ift(["t1", "t2"])
    scm = together["ph1", 1]
    t1_len = scm.shape["t1"]
    time_domain = together["ph1", 0]
    time_domain.extend("t1", -time_domain["t1"][-1])
    negative_startpoint = time_domain["t1"][0]
    offset = negative_startpoint + scm["t1"][-1]
    scm.ft("t1")
    scm *= np.exp(-1j * 2 * pi * offset * scm.fromaxis("t1"))
    scm.ift("t1")
    negative_region = time_domain["t1":(None, 0)]
    mirror_len = min(negative_region.shape["t1"], t1_len)
    # Copy the mirrored Sc- component into the negative time region.
    negative_region["t1", :mirror_len] = scm["t1", ::-1]["t1", :mirror_len]
    # Prepare the frequency-domain views used for the two right-hand panels.
    corrected["ph1", 1]["t1", :] = corrected["ph1", 1]["t1", ::-1]
    corrected = corrected["t1":T1_RANGE]["t2":T2_RANGE]
    return time_domain, corrected


# }}}


class PhaseCorrectionWidget(QWidget):
    def __init__(self, base_dataset=None, exp_types=None, search_string=None):
        QWidget.__init__(self)
        self.setWindowTitle("Diagonal and Anti-diagonal Corrections")
        self.figure = Figure(figsize=(10, 8))
        self.figure.patch.set_alpha(0)
        manager = new_figure_manager_given_figure(id(self), self.figure)
        self.manager = manager
        self.canvas = manager.canvas
        self.canvas.setParent(self)
        self.canvas.setStyleSheet("background: transparent")
        if not hasattr(self.manager, "_cidgcf"):
            self.manager._cidgcf = None
        plt.figure(manager.num)
        grid = self.figure.add_gridspec(
            2, 3, width_ratios=[1, 1.2, 1], height_ratios=[1, 1]
        )
        self.ax_hermitian_ft = self.figure.add_subplot(grid[:, 0])
        self.ax_time = self.figure.add_subplot(grid[:, 1])
        self.ax_scp = self.figure.add_subplot(grid[0, 2])
        self.ax_scm = self.figure.add_subplot(grid[1, 2])
        for axis in [self.ax_time, self.ax_scp, self.ax_scm]:
            axis.set_facecolor((1.0, 1.0, 1.0, 0.0))
        self.ax_hermitian_ft.set_facecolor((1.0, 1.0, 1.0, 0.0))
        self.ax_hermitian_ft.set_title("FT of Hermitian Time")
        self.ax_time.set_title("Time Domain Hermitian Time, no apo")
        self.ax_scp.set_title("$S_{c+}$")
        self.ax_scm.set_title("$S_{c-}$")
        self.slider_scale = 100
        self.sensitivity_scale = 100
        self.diag_value = 0.0
        self.anti_value = 0.0
        self.sensitivity_value = 0.5
        self.diag_label = QLabel()
        self.diag_slider = QSlider(Qt.Horizontal)
        self.diag_slider.setTickPosition(QSlider.TicksBelow)
        self.anti_label = QLabel()
        self.anti_slider = QSlider(Qt.Horizontal)
        self.anti_slider.setTickPosition(QSlider.TicksBelow)
        self.sensitivity_label = QLabel()
        self.sensitivity_slider = QSlider(Qt.Horizontal)
        self.sensitivity_slider.setTickPosition(QSlider.TicksBelow)
        self.sensitivity_slider.setMinimum(int(0.05 * self.sensitivity_scale))
        self.sensitivity_slider.setMaximum(int(3.0 * self.sensitivity_scale))
        self.sensitivity_slider.setSingleStep(1)
        self.sensitivity_slider.setPageStep(int(0.1 * self.sensitivity_scale))
        self.sensitivity_slider.setTickInterval(max(int(0.1 * self.sensitivity_scale), 1))
        self.diag_slider.setSingleStep(1)
        self.anti_slider.setSingleStep(1)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.hermitian_ft_image = None
        slider_layout = QVBoxLayout()
        diag_layout = QHBoxLayout()
        diag_layout.addWidget(QLabel("Diagonal correction"))
        diag_layout.addWidget(self.diag_label)
        slider_layout.addLayout(diag_layout)
        slider_layout.addWidget(self.diag_slider)
        anti_layout = QHBoxLayout()
        anti_layout.addWidget(QLabel("Anti-diagonal correction"))
        anti_layout.addWidget(self.anti_label)
        slider_layout.addLayout(anti_layout)
        slider_layout.addWidget(self.anti_slider)
        sensitivity_layout = QHBoxLayout()
        sensitivity_layout.addWidget(QLabel("Sensitivity"))
        sensitivity_layout.addWidget(self.sensitivity_label)
        slider_layout.addLayout(sensitivity_layout)
        slider_layout.addWidget(self.sensitivity_slider)
        # Build the data-loading controls so the user can switch datasets.
        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel("Experiment type"))
        self.exp_type_dropdown = QComboBox()
        if exp_types is None:
            if pyspec_config._config_parser is None:
                parser = configparser.ConfigParser()
                parser.read(pyspec_config.config_location)
            else:
                parser = pyspec_config._config_parser
            if parser.has_section("ExpTypes"):
                exp_types = sorted(name for name, _ in parser.items("ExpTypes"))
            else:
                exp_types = []
        if not exp_types:
            exp_types = [DEFAULT_EXP_TYPE]
        for name in exp_types:
            self.exp_type_dropdown.addItem(name)
        top_layout.addWidget(self.exp_type_dropdown)
        top_layout.addWidget(QLabel("Search"))
        self.search_entry = QLineEdit()
        if search_string is None:
            search_string = DEFAULT_SEARCH
        self.search_entry.setText(search_string)
        top_layout.addWidget(self.search_entry)
        self.go_button = QPushButton("Go")
        top_layout.addWidget(self.go_button)
        # Provide text boxes for truncating the time-domain points.
        truncation_layout = QHBoxLayout()
        truncation_layout.addWidget(QLabel("t1 points"))
        self.t1_entry = QLineEdit()
        truncation_layout.addWidget(self.t1_entry)
        truncation_layout.addWidget(QLabel("t2 points"))
        self.t2_entry = QLineEdit()
        truncation_layout.addWidget(self.t2_entry)
        layout = QVBoxLayout()
        layout.addLayout(top_layout)
        layout.addLayout(truncation_layout)
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        layout.addLayout(slider_layout)
        self.setLayout(layout)
        if base_dataset is None:
            base_dataset = load_base_dataset(
                    DEFAULT_EXP_TYPE, DEFAULT_SEARCH
            )
        self.full_dataset = base_dataset
        self.base_dataset = base_dataset
        self.current_t1_points = self.base_dataset.shape["t1"]
        self.current_t2_points = self.base_dataset.shape["t2"]
        self.t1_entry.setText(str(self.current_t1_points))
        self.t2_entry.setText(str(self.current_t2_points))
        self.time_image = None
        self.time_colorbar = None
        self.scp_image = None
        self.scm_image = None
        self.hermitian_ft_colorbar = None
        self.current_time_shape = None
        self.current_freq_shape = None
        self.last_time_domain = None
        self.go_button.clicked.connect(self.load_new_dataset)
        self.t1_entry.editingFinished.connect(self.update_truncation)
        self.t2_entry.editingFinished.connect(self.update_truncation)
        for slider in [self.diag_slider, self.anti_slider, self.sensitivity_slider]:
            slider.blockSignals(True)
        self.diag_slider.setPageStep(10)
        self.anti_slider.setPageStep(10)
        self.sensitivity_slider.setValue(int(self.sensitivity_value * self.sensitivity_scale))
        self.diag_slider.setValue(0)
        self.anti_slider.setValue(0)
        for slider in [self.diag_slider, self.anti_slider, self.sensitivity_slider]:
            slider.blockSignals(False)
        self.diag_slider_moved(self.diag_slider.value())
        self.anti_slider_moved(self.anti_slider.value())
        self.sensitivity_slider_moved(self.sensitivity_slider.value())
        self.diag_slider.valueChanged.connect(self.diag_slider_moved)
        self.diag_slider.sliderReleased.connect(self.handle_slider_release)
        self.anti_slider.valueChanged.connect(self.anti_slider_moved)
        self.anti_slider.sliderReleased.connect(self.handle_slider_release)
        self.sensitivity_slider.valueChanged.connect(self.sensitivity_slider_moved)
        self.sensitivity_slider.sliderReleased.connect(self.update_phase_slider_limits)
        self.update_phase_slider_limits()
        self.update_plots()

    def load_new_dataset(self):
        """Reload the raw data using the search controls."""
        exp_type = self.exp_type_dropdown.currentText()
        search_string = self.search_entry.text()
        try:
            dataset = load_base_dataset(exp_type, search_string)
        except Exception as error:  # noqa: BLE001
            print("Failed to load dataset:", error)
            return
        self.full_dataset = dataset
        self.base_dataset = dataset
        self.current_t1_points = self.base_dataset.shape["t1"]
        self.current_t2_points = self.base_dataset.shape["t2"]
        self.t1_entry.setText(str(self.current_t1_points))
        self.t2_entry.setText(str(self.current_t2_points))
        self.reset_images()
        self.update_plots()

    def update_truncation(self):
        """Slice the dataset before corrections when the limits change."""
        max_t1 = self.full_dataset.shape["t1"]
        max_t2 = self.full_dataset.shape["t2"]
        try:
            requested_t1 = int(self.t1_entry.text())
        except ValueError:
            requested_t1 = self.current_t1_points
        try:
            requested_t2 = int(self.t2_entry.text())
        except ValueError:
            requested_t2 = self.current_t2_points
        requested_t1 = max(1, min(requested_t1, max_t1))
        requested_t2 = max(1, min(requested_t2, max_t2))
        self.current_t1_points = requested_t1
        self.current_t2_points = requested_t2
        self.t1_entry.setText(str(self.current_t1_points))
        self.t2_entry.setText(str(self.current_t2_points))
        self.base_dataset = self.full_dataset["t1", :self.current_t1_points][
            "t2", :self.current_t2_points
        ]
        self.reset_images()
        self.update_plots()

    def reset_images(self):
        """Clear cached images when the dataset shape changes."""
        if self.time_colorbar is not None:
            self.time_colorbar.remove()
            self.time_colorbar = None
        if self.hermitian_ft_colorbar is not None:
            self.hermitian_ft_colorbar.remove()
            self.hermitian_ft_colorbar = None
        self.ax_time.cla()
        self.ax_time.set_title("Time Domain Hermitian Time, no apo")
        self.ax_hermitian_ft.cla()
        self.ax_hermitian_ft.set_title("FT of Hermitian Time")
        self.ax_scp.cla()
        self.ax_scp.set_title("$S_{c+}$")
        self.ax_scm.cla()
        self.ax_scm.set_title("$S_{c-}$")
        for axis in [self.ax_hermitian_ft, self.ax_time, self.ax_scp, self.ax_scm]:
            axis.set_facecolor((1.0, 1.0, 1.0, 0.0))
        self.time_image = None
        self.scp_image = None
        self.scm_image = None
        self.hermitian_ft_image = None
        self.current_time_shape = None
        self.current_freq_shape = None

    def diag_slider_moved(self, value):
        self.diag_value = value / self.slider_scale
        self.diag_label.setText(f"{self.diag_value:.2f} μs")

    def anti_slider_moved(self, value):
        self.anti_value = value / self.slider_scale
        self.anti_label.setText(f"{self.anti_value:.2f} μs")

    def sensitivity_slider_moved(self, value):
        self.sensitivity_value = value / self.sensitivity_scale
        self.sensitivity_label.setText(f"±{self.sensitivity_value:.2f} μs")

    def update_phase_slider_limits(self):
        range_steps = max(int(round(self.sensitivity_value * self.slider_scale)), 1)
        for slider, current in [
            (self.diag_slider, self.diag_value),
            (self.anti_slider, self.anti_value),
        ]:
            center = int(round(current * self.slider_scale))
            slider.blockSignals(True)
            slider.setMinimum(center - range_steps)
            slider.setMaximum(center + range_steps)
            slider.setPageStep(max(range_steps // 5, 1))
            slider.setTickInterval(max(range_steps // 5, self.slider_scale // 10, 1))
            slider.setValue(center)
            slider.blockSignals(False)
        self.diag_slider_moved(self.diag_slider.value())
        self.anti_slider_moved(self.anti_slider.value())

    def handle_slider_release(self):
        self.update_plots()

    def update_plots(self):
        diag_corr = self.diag_value
        anti_corr = self.anti_value
        plt.figure(self.manager.num)
        truncated = self.base_dataset
        time_domain, corrected = apply_phase_corrections(
            truncated, diag_corr, anti_corr
        )
        self.last_time_domain = time_domain
        time_magnitude_nd = abs(time_domain)
        time_magnitude = time_magnitude_nd.data
        scp = corrected["ph1", 0]
        scm = corrected["ph1", 1]
        scale = abs(corrected).max()
        hermitian_source = time_domain.copy()
        # Reset Fourier bookkeeping before transforming the Hermitian data again.
        for axis in ["t1", "t2"]:
            hermitian_source.set_ft_prop(axis, None, None)
            hermitian_source.set_ft_prop(axis, ["start", "freq"], None)
        hermitian_ft = hermitian_source.ft(["t1", "t2"], shift=True)
        hermitian_magnitude_nd = abs(hermitian_ft)
        hermitian_magnitude = hermitian_magnitude_nd.data
        time_shape = time_magnitude.shape
        freq_shape = scp.data.shape
        if self.current_time_shape != time_shape or self.current_freq_shape != freq_shape:
            self.reset_images()
        if self.time_image is None:
            self.time_image = image(time_magnitude_nd, cmap="single_sided", ax=self.ax_time)
            self.time_colorbar = self.time_image.colorbar
            self.ax_time.set_aspect("equal")
            self.hermitian_ft_image = image(
                hermitian_magnitude_nd, cmap="single_sided", ax=self.ax_hermitian_ft
            )
            self.hermitian_ft_colorbar = self.hermitian_ft_image.colorbar
            self.ax_hermitian_ft.set_aspect("equal")
            self.scp_image = image(scp, scaling=scale, ax=self.ax_scp)
            self.scm_image = image(scm, scaling=scale, ax=self.ax_scm)
            self.current_time_shape = time_shape
            self.current_freq_shape = freq_shape
            self.canvas.draw_idle()
            return
        self.time_image.set_data(time_magnitude)
        if self.time_colorbar is not None:
            vmax = float(np.nanmax(time_magnitude))
            if vmax <= 0:
                vmax = 1.0
            self.time_image.set_clim(0, vmax)
            self.time_colorbar.update_normal(self.time_image)
        if self.hermitian_ft_image is not None:
            hermitian_vmax = float(np.nanmax(hermitian_magnitude))
            if hermitian_vmax <= 0:
                hermitian_vmax = 1.0
            self.hermitian_ft_image.set_data(hermitian_magnitude)
            self.hermitian_ft_image.set_clim(0, hermitian_vmax)
            if self.hermitian_ft_colorbar is not None:
                self.hermitian_ft_colorbar.update_normal(self.hermitian_ft_image)
        scp_colors = imagehsv(scp.data, scaling=scale)
        scm_colors = imagehsv(scm.data, scaling=scale)
        self.scp_image.set_data(scp_colors)
        self.scm_image.set_data(scm_colors)
        self.current_time_shape = time_shape
        self.current_freq_shape = freq_shape
        self.canvas.draw_idle()


def main():
    app = QApplication(sys.argv)
    widget = PhaseCorrectionWidget()
    widget.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
