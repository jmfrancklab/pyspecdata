"""Qt example showing diagonal and anti-diagonal corrections with sliders."""
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import new_figure_manager_given_figure
from matplotlib.figure import Figure
from numpy import pi
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QLabel, QHBoxLayout, QSlider, QVBoxLayout, QWidget

from pyspecdata import find_file, image


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


# {{{ data preparation helpers
# Load the dataset once so that slider updates only adjust phase factors.
def load_base_dataset():
    d = find_file(
        "1500ns\\.DSC", exp_type="2DELDOR/251027_DiPhNO_500uM/Hypercomplex"
    )
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
    def __init__(self, base_dataset=None):
        QWidget.__init__(self)
        self.setWindowTitle("Diagonal and Anti-diagonal Corrections")
        # Use a pyplot-managed figure so pyspecdata.image can locate it.
        self.figure = Figure(figsize=(10, 8))
        manager = new_figure_manager_given_figure(id(self), self.figure)
        self.manager = manager
        self.canvas = manager.canvas
        self.canvas.setParent(self)
        if not hasattr(self.manager, "_cidgcf"):
            self.manager._cidgcf = None
        plt.figure(manager.num)
        grid = self.figure.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
        self.ax_time = self.figure.add_subplot(grid[:, 0])
        self.ax_scp = self.figure.add_subplot(grid[0, 1])
        self.ax_scm = self.figure.add_subplot(grid[1, 1])
        for axis in [self.ax_time, self.ax_scp, self.ax_scm]:
            axis.set_facecolor((1.0, 1.0, 1.0, 0.9))
        self.ax_time.set_title("Time Domain Hermitian Time, no apo")
        self.ax_scp.set_title("$S_{c+}$")
        self.ax_scm.set_title("$S_{c-}$")
        self.diag_label = QLabel()
        self.diag_slider = QSlider(Qt.Horizontal)
        self.diag_slider.setMinimum(-200)
        self.diag_slider.setMaximum(200)
        self.diag_slider.setSingleStep(1)
        self.diag_slider.setValue(0)
        self.anti_label = QLabel()
        self.anti_slider = QSlider(Qt.Horizontal)
        self.anti_slider.setMinimum(-200)
        self.anti_slider.setMaximum(200)
        self.anti_slider.setSingleStep(1)
        self.anti_slider.setValue(0)
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
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addLayout(slider_layout)
        self.setLayout(layout)
        if base_dataset is None:
            # Store the loaded dataset so repeated slider updates avoid disk I/O.
            base_dataset = load_base_dataset()
        self.base_dataset = base_dataset
        self.diag_slider.valueChanged.connect(self.slider_changed)
        self.anti_slider.valueChanged.connect(self.slider_changed)
        self.update_plots()

    def slider_changed(self, value):
        # Trigger a refresh whenever either slider reports a new value.
        self.update_plots()

    def update_plots(self):
        diag_corr = self.diag_slider.value() / 100.0
        anti_corr = self.anti_slider.value() / 100.0
        self.diag_label.setText(f"{diag_corr:.2f}")
        self.anti_label.setText(f"{anti_corr:.2f}")
        # Re-select the figure so pyspecdata.image can find the active canvas.
        plt.figure(self.manager.num)
        time_domain, corrected = apply_phase_corrections(
            self.base_dataset, diag_corr, anti_corr
        )
        self.ax_time.clear()
        image(abs(time_domain), cmap="single_sided", ax=self.ax_time)
        self.ax_time.set_aspect("equal")
        self.ax_time.set_title("Time Domain Hermitian Time, no apo")
        self.ax_scp.clear()
        self.ax_scm.clear()
        scale = abs(corrected).max()
        image(corrected["ph1", 0], scaling=scale, ax=self.ax_scp)
        image(corrected["ph1", 1], scaling=scale, ax=self.ax_scm)
        self.ax_scp.set_title("$S_{c+}$")
        self.ax_scm.set_title("$S_{c-}$")
        self.figure.tight_layout()
        self.canvas.draw_idle()


def main():
    app = QApplication(sys.argv)
    widget = PhaseCorrectionWidget()
    widget.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
