import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyspecdata as psd  # noqa: E402
from numpy import r_  # noqa: E402


def test_plot_spreading_box_spreads_fixed_pigment():
    """Boxes are 3-pixel squares unless the error spreads them out."""
    d = (
        psd.nddata(r_[1.0, 2.0, 1.5, 3.0, 2.5], "t")
        .setaxis("t", r_[0:4:5j])
        .set_units("t", "s")
        .set_error(r_[0.01, 2.0, 0.01, 0.01, 0.01])
        .set_error("t", r_[0.0, 0.0, 0.0, 3.0, 0.0])
    )
    fig, ax = plt.subplots()
    boxes = d.plot_spreading_box(ax=ax)
    fig.canvas.draw()
    # width and height (data units) and opacity of each drawn box
    corners = np.array([path.vertices[:4] for path in boxes.get_paths()])
    width = np.ptp(corners[:, :, 0], axis=1)
    height = np.ptp(corners[:, :, 1], axis=1)
    alpha = boxes.get_facecolor()[:, 3]
    # the baseline is a 3 × 3 pixel square
    pixels = np.abs(
        np.diff(
            ax.transData.transform([[0.0, 0.0], [width[0], height[0]]]),
            axis=0,
        )[0]
    )
    np.testing.assert_allclose(pixels, 3.0)
    np.testing.assert_allclose(width[[0, 1, 2, 4]], width[0])
    np.testing.assert_allclose(height[[0, 2, 3, 4]], height[0])
    np.testing.assert_allclose(alpha[[0, 2, 4]], 0.5)
    # a large error in the data or the coordinate spreads the box to 2σ
    np.testing.assert_allclose(height[1], 4.0)
    np.testing.assert_allclose(width[3], 6.0)
    # ...with the same amount of pigment (area × opacity)
    np.testing.assert_allclose(
        alpha * width * height, 0.5 * width[0] * height[0]
    )
    assert np.all(boxes.get_linewidths() == 0)
    plt.close(fig)
