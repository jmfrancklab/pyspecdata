import sys
from collections import OrderedDict
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QComboBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np


class PlotApp(object):
    def __init__(self):
        self.app = QApplication(sys.argv)
        self.main_window = QMainWindow()
        self.main_window.setWindowTitle("Polynomial Plotter")
        self.main_window.setGeometry(100, 100, 800, 600)
        self.figures = OrderedDict()  # Ordered storage for figures

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
        self.canvas = FigureCanvas(self[list(self.figures.keys())[0]])
        self.layout.addWidget(self.canvas)
        self.populate_combo_box()  # Populate combo box after all figures are added
        self.main_window.show()  # Show the main window
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

        # Canvas placeholder (set later in `populate_combo_box`)
        self.canvas = None

        # Connect combo box to figure update logic
        self.combo_box.currentIndexChanged.connect(self.update_figure)

    def populate_combo_box(self):
        self.combo_box.addItems(self.figures.keys())  # Add figure names to the combo box

        # Initialize the canvas with the first figure
        if self.figures:
                        self.update_figure(0)  # Display the first figure

    def update_figure(self, index):
        figure_name = list(self.figures.keys())[index]  # Access figure key by index

        # Remove the old canvas from the layout
        if self.canvas is not None:
            self.layout.removeWidget(self.canvas)
            self.canvas.deleteLater()

        # Create a new canvas for the selected figure
        self.canvas = FigureCanvas(self[figure_name])
        self.canvas.setStyleSheet("background: transparent;")  # Make the canvas background transparent
        self.layout.addWidget(self.canvas)
        self.canvas.draw_idle()


def main():
    with PlotApp() as main_window:
        # Create and plot figures
        for i in range(1, 6):  # Polynomial orders 1 to 5
            fig = main_window[f"Polynomial {i}"]
            ax = fig.add_subplot(111)
            x = np.linspace(-10, 10, 500)  # Smooth range of x values
            y = np.power(x, i)  # NumPy for polynomial calculation
            ax.plot(x, y)
            ax.set_title(f"$y = x^{i}$")


if __name__ == "__main__":
    main()
