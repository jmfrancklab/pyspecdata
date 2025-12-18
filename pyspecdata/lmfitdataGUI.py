from .core import plot as psd_plot

try:
    from PySide6 import QtWidgets, QtCore
    from matplotlib.backends.backend_qtagg import (
        FigureCanvasQTAgg as FigureCanvas,
    )
    from matplotlib.backends.backend_qtagg import (
        NavigationToolbar2QT as NavigationToolbar,
    )
except Exception:
    from PyQt5 import QtWidgets, QtCore  # type: ignore
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas  # type: ignore
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar  # type: ignore
import numpy as np
import matplotlib.pyplot as plt


class lmfitdataGUI(QtWidgets.QWidget):
    """Qt UI for exploring lmfit parameters on a pyspecdata.lmfitdata.

    Pass the lmfitdata object ``d`` to ``__init__``. The top area embeds a
    Matplotlib Axes plotting ``d.settoguess().eval()`` via ``pyspecdata.plot``.

    Below that:
      • A *parameter checkbox grid* (alpha-sorted). Checking a box sets
        ``param.vary=True`` and shows its slider; unchecking hides the slider
        and sets ``param.vary=False``.
      • Optionally, a *control slider* provided via ``control=(func, (lo, hi, init))``.
        ``func`` is called as ``func(d, value)`` on change; the row is styled to
        stand out and is placed above the parameter sliders.
      • A scrollable list of sliders for varying parameters. Slider ranges use
        (min,max) when finite; otherwise ±5% around the current value (expanded
        to include the current value).

    This class always creates its own ``QApplication`` in ``__init__``. The
    event loop is not started automatically; call ``.exec()`` to show the
    window and start the loop.
    """

    def __init__(self, d, control=None):
        # Always create our own QApplication (script-owned GUI)
        self._app = QtWidgets.QApplication([])

        super().__init__()
        self.d = d
        self.params = d.guess_parameters  # lmfit.Parameters
        self.setWindowTitle("lmfit parameter explorer")

        # --- main vertical layout ---
        vbox = QtWidgets.QVBoxLayout(self)

        # --- Matplotlib canvas + toolbar ---
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.fig)
        # Transparent figure/bg so Qt widget shows through
        try:
            self.fig.patch.set_alpha(0)
            self.fig.patch.set_facecolor("none")
            self.ax.set_facecolor("none")
            self.canvas.setStyleSheet("background: transparent")
        except Exception:
            pass
        # Toolbar row (top)
        toolbar_row = QtWidgets.QWidget()
        hbar = QtWidgets.QHBoxLayout(toolbar_row)
        hbar.setContentsMargins(0, 0, 0, 0)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.zoomfit_btn = QtWidgets.QPushButton("Zoom to fit")
        self.zoomfit_btn.clicked.connect(self._zoom_to_fit)
        hbar.addWidget(self.toolbar)
        hbar.addStretch(1)
        hbar.addWidget(self.zoomfit_btn)
        vbox.addWidget(toolbar_row)
        vbox.addWidget(self.canvas)

        # --- Scrollable center area ---
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        vbox.addWidget(scroll)
        inner = QtWidgets.QWidget()
        scroll.setWidget(inner)
        inner_v = QtWidgets.QVBoxLayout(inner)
        inner_v.setContentsMargins(0, 0, 0, 0)

        # ===== Parameter checkbox grid (alpha-sorted) =====
        chk_container = QtWidgets.QGroupBox("Parameters (toggle vary on/off)")
        grid = QtWidgets.QGridLayout(chk_container)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(6)
        self._check_boxes = {}
        p_names = sorted(list(self.params))
        ncol = 4  # compact grid; adjust if you like
        for idx, name in enumerate(p_names):
            row, col = divmod(idx, ncol)
            cb = QtWidgets.QCheckBox(name)
            vary = bool(getattr(self.params[name], "vary", True))
            cb.setChecked(vary)
            cb.stateChanged.connect(self._make_toggle_vary_handler(name))
            self._check_boxes[name] = cb
            grid.addWidget(cb, row, col)
        inner_v.addWidget(chk_container)

        # ===== Optional control slider (distinct style) =====
        self._control = None
        if control is not None:
            func, bounds = control
            if not callable(func):
                raise TypeError("control[0] must be callable: func(d, value)")
            if not (isinstance(bounds, (tuple, list)) and len(bounds) == 3):
                raise TypeError("control[1] must be (lo, hi, init)")
            clo, chi, cinit = map(float, bounds)
            cinit = max(min(cinit, chi), clo)
            self._control = (func, (clo, chi))

            ctrl_row = QtWidgets.QWidget()
            ch = QtWidgets.QHBoxLayout(ctrl_row)
            ch.setContentsMargins(0, 0, 0, 0)
            ctrl_label = QtWidgets.QLabel("control")
            ctrl_label.setMinimumWidth(120)
            self._ctrl_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            self._ctrl_slider.setRange(0, 1000)
            self._ctrl_box = QtWidgets.QDoubleSpinBox()
            self._ctrl_box.setDecimals(10)
            self._ctrl_box.setRange(min(clo, chi), max(clo, chi))
            self._ctrl_box.setSingleStep(
                (chi - clo) / 1000.0 if chi > clo else 1.0
            )
            # distinctive color via stylesheet
            self._ctrl_slider.setStyleSheet(
                "QSlider::groove:horizontal{height:6px;background:#2255cc;}QSlider::handle:horizontal{background:#66a3ff;border:1px"
                " solid #1a4ca3;width:12px;margin:-4px 0;}"
            )
            # init values
            self._ctrl_slider.setValue(self._to_slider(cinit, clo, chi))
            self._ctrl_box.setValue(cinit)

            updating = {"flag": False}

            def on_ctrl_slider(i):
                if updating["flag"]:
                    return
                updating["flag"] = True
                val = self._from_slider(i, clo, chi)
                self._ctrl_box.setValue(val)
                try:
                    func(self.d, val)
                finally:
                    self._replot()
                    updating["flag"] = False

            def on_ctrl_box(val):
                if updating["flag"]:
                    return
                updating["flag"] = True
                v = float(np.clip(val, clo, chi))
                self._ctrl_slider.setValue(self._to_slider(v, clo, chi))
                try:
                    func(self.d, v)
                finally:
                    self._replot()
                    updating["flag"] = False

            self._ctrl_slider.valueChanged.connect(on_ctrl_slider)
            self._ctrl_box.valueChanged.connect(on_ctrl_box)

            ch.addWidget(ctrl_label)
            ch.addWidget(self._ctrl_slider, 1)
            ch.addWidget(self._ctrl_box)
            # put control row ABOVE parameter sliders
            inner_v.addWidget(ctrl_row)

        # ===== Sliders list (form) =====
        form_container = QtWidgets.QWidget()
        self.form = QtWidgets.QFormLayout(form_container)
        self.form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.ExpandingFieldsGrow
        )
        inner_v.addWidget(form_container)

        # bookkeeping: name -> {row, slider, box, bounds}
        self._widgets = {}
        for name in p_names:
            p = self.params[name]
            lo, hi = self._bounds_for_param(p)
            row = self._make_param_row(name, p, lo, hi)
            self.form.addRow(row)
            # show only if vary==True
            row.setVisible(bool(getattr(p, "vary", True)))

        self.resize(1000, 750)
        self._replot()

    # ---------------- helpers -----------------
    def _make_toggle_vary_handler(self, name):
        def _handler(state):
            checked = state == QtCore.Qt.Checked
            # update lmfit parameter
            if name in self.params:
                self.params[name].vary = bool(checked)
            # ensure slider row exists
            if name not in self._widgets:
                p = self.params[name]
                lo, hi = self._bounds_for_param(p)
                row = self._make_param_row(name, p, lo, hi)
                self.form.addRow(row)
                row.setVisible(checked)
            else:
                row = self._widgets[name]["row"]
                row.setVisible(checked)
            self._replot()

        return _handler

    def _bounds_for_param(self, p):
        val = float(p.value)
        pmin = getattr(p, "min", -np.inf)
        pmax = getattr(p, "max", np.inf)
        if np.isfinite(pmin) and np.isfinite(pmax) and pmin < pmax:
            lo, hi = float(pmin), float(pmax)
        else:
            delta = 0.05 * (abs(val) if abs(val) > 0 else 1.0)
            lo, hi = val - delta, val + delta
        lo = min(lo, val)
        hi = max(hi, val)
        if lo == hi:
            lo, hi = lo - 1.0, hi + 1.0
        return lo, hi

    @staticmethod
    def _to_slider(v, lo, hi):
        return (
            int(round(1000 * (float(v) - lo) / (hi - lo))) if hi != lo else 0
        )

    @staticmethod
    def _from_slider(i, lo, hi):
        return lo + (hi - lo) * (i / 1000.0) if hi != lo else lo

    def _make_param_row(self, name, p, lo, hi):
        row = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(row)
        h.setContentsMargins(0, 0, 0, 0)

        name_label = QtWidgets.QLabel(name)
        name_label.setMinimumWidth(120)
        slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        slider.setRange(0, 1000)
        box = QtWidgets.QDoubleSpinBox()
        box.setDecimals(10)
        box.setRange(min(lo, hi), max(lo, hi))
        box.setSingleStep((hi - lo) / 1000.0 if hi > lo else 1.0)

        slider.setValue(self._to_slider(p.value, lo, hi))
        box.setValue(float(p.value))

        updating = {"flag": False}

        def on_slider(i):
            if updating["flag"]:
                return
            updating["flag"] = True
            val = self._from_slider(i, lo, hi)
            self.params[name].set(value=val)
            box.setValue(val)
            self._replot()
            updating["flag"] = False

        def on_box(val):
            if updating["flag"]:
                return
            updating["flag"] = True
            v = float(np.clip(val, lo, hi))
            self.params[name].set(value=v)
            slider.setValue(self._to_slider(v, lo, hi))
            self._replot()
            updating["flag"] = False

        slider.valueChanged.connect(on_slider)
        box.valueChanged.connect(on_box)

        h.addWidget(name_label)
        h.addWidget(slider, 1)
        h.addWidget(box)
        self._widgets[name] = {
            "row": row,
            "slider": slider,
            "box": box,
            "bounds": (lo, hi),
        }
        return row

    def _zoom_to_fit(self):
        try:
            self.ax.relim()
            self.ax.autoscale_view()
            self.canvas.draw_idle()
        except Exception:
            pass

    def exec(self):
        """Show the window and start the Qt event loop we created in `__init__`.
        Returns the event-loop exit code.
        """
        self.show()
        try:
            self.raise_()
            self.activateWindow()
        except Exception:
            pass
        # Start the Qt event loop (we always own the QApplication)
        if hasattr(self._app, "exec"):
            return self._app.exec()
        else:  # PyQt5
            return self._app.exec_()

    def _replot(self):
        # Preserve current view limits if we've plotted before
        preserve = hasattr(self, "_had_first_plot") and self._had_first_plot
        if preserve:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
        self.ax.clear()
        # keep background transparent after clear
        try:
            self.ax.set_facecolor("none")
        except Exception:
            pass
        y = self.d.settoguess().eval()
        if not hasattr(self, "data_forcomp"):
            # NOTE: this should be converted to a function that we feed
            #       NO -- solved by assymetric transform
            self.data_forcomp = self.d.copy(data=False)
            self.data_forcomp.data = self.transformed_data
        try:
            psd_plot(self.data_forcomp, ax=self.ax, alpha=0.5)
            psd_plot(y, ax=self.ax, alpha=0.5)
        except TypeError:
            plt.sca(self.ax)
            psd_plot(self.data_forcomp, ax=self.ax, alpha=0.5)
            psd_plot(y, ax=self.ax, alpha=0.5)
        self.ax.relim()
        if preserve:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            self.ax.autoscale_view()
        self.canvas.draw_idle()
        self._had_first_plot = True
