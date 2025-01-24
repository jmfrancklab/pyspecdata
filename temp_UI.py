import sys
from collections import OrderedDict
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QWidget,
    QComboBox,
)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
)
from matplotlib.figure import Figure
import numpy as np
import pyspecdata as psd
import threading


class PlotApp(object):
    def __init__(self):
        self.app = QApplication(sys.argv)
        self.main_window = QMainWindow()
        self.main_window.setWindowTitle("Polynomial Plotter")
        self.main_window.setGeometry(100, 100, 800, 600)
        self.figures = OrderedDict()  # Ordered storage for figures
        self.canvas = None
        # Canvas placeholder (set later in `first_screen`)
        self.preloaded_canvases = {}

    def __getitem__(self, key):
        if key not in self.figures:
            fig = Figure()
            fig.patch.set_facecolor("none")  # Make the figure box transparent
            self.figures[key] = fig  # Add figure to the ordered dictionary
        return self.figures[key]

    def __enter__(self):
        self.initUI()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.first_screen()  # Populate combo box after all figures are added
        self.preload_figures()  # Preload other figures in a separate thread
        self.main_window.show()  # Show the main window
        self.thread.join()  # Ensure the preload thread completes before exiting
        sys.exit(self.app.exec_())  # Start the event loop

    def initUI(self):
        # Instantiate QApplication at the end of the UI setup
        # Create central widget and layout
        self.central_widget = QWidget(self.main_window)
        self.main_window.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)

        # Create combo box for figure selection
        self.combo_box = QComboBox(self.central_widget)
        self.layout.addWidget(self.combo_box)

        # Connect combo box to figure update logic
        self.combo_box.currentIndexChanged.connect(self.update_figure)

    def generate_canvas(self, key):
        print("generating", key)
        canvas = FigureCanvas(self.figures[key])
        canvas.setStyleSheet(
            "background: transparent;"
        )  # Make the canvas background transparent
        self.preloaded_canvases[key] = canvas
        return

    def preload_figures(self):
        for key in list(self.figures.keys())[1:]:  # Skip the first figure
            self.generate_canvas(key)

        def render_figures():
            for key in list(self.figures.keys())[1:]:  # Skip the first figure
                self.preloaded_canvases[key].draw_idle()

        self.thread = threading.Thread(target=render_figures)
        self.thread.start()

    def first_screen(self):
        self.combo_box.addItems(
            self.figures.keys()
        )  # Add figure names to the combo box
        # note that this seems to call update_figure

    def update_figure(self, index):
        figure_name = list(self.figures.keys())[
            index
        ]  # Access figure key by index
        print("update figure", index, figure_name, self.canvas)
        if self.canvas is None:
            # no figure set yet, so use the first key
            figure_name = next(iter(self.figures.keys()))
            # and preload the first figure
            self.generate_canvas(figure_name)
            self.preloaded_canvases[figure_name].draw_idle()
        else:
            # Remove the old canvas from the layout
            self.canvas.hide()
            self.layout.removeWidget(self.canvas)
            self.thread.join()  # Ensure the preload thread completes before exiting
        self.canvas = self.preloaded_canvases[figure_name]
        self.canvas.show()
        # Use the preloaded canvas for the selected figure
        self.layout.addWidget(self.canvas)


def main():
    with PlotApp() as main_window:
        # Create and plot figures
        for j in range(1, 6):  # Polynomial orders 1 to 5
            fig = main_window[f"Polynomial {j}"]
            ax = fig.add_subplot(111)
            data = psd.nddata(np.linspace(-10, 10, 500), "x")  # Smooth range of x values
            data.run(lambda x: x**j)  # NumPy for polynomial calculation
            psd.plot(data, ax=ax)
            ax.set_title(f"$y = x^{j}$")


if __name__ == "__main__":
    main()
