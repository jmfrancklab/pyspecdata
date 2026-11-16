import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from numpy import r_  # noqa: E402
from conftest import load_module  # noqa: E402


def test_plot_spreading_box_spreads_fixed_pigment():
    """Boxes are opaque 3-pixel squares unless the error spreads them
    out."""
    # like the other tests, load core through conftest (and do it here,
    # rather than when the tests are collected), so that we use the copy of
    # the nddata class that is current, regardless of which tests have run
    # before this one
    nddata = load_module("core", use_real_pint=True).nddata
    d = (
        nddata(r_[1.0, 2.0, 1.5, 3.0, 2.5], "t")
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
    # the small boxes are loaded with more pigment than they can show, so
    # they are thresholded to fully opaque
    np.testing.assert_allclose(alpha[[0, 2, 4]], 1.0)
    # a large error in the data or the coordinate spreads the box to 2σ
    np.testing.assert_allclose(height[1], 4.0)
    np.testing.assert_allclose(width[3], 6.0)
    # ...and spreads out all of the pigment (area × opacity), which is the
    # default base_alpha of 5 times what covers the baseline square
    np.testing.assert_allclose(
        (alpha * width * height)[[1, 3]], 5 * width[0] * height[0]
    )
    assert np.all(boxes.get_linewidths() == 0)
    # the amount of pigment belongs to each plot, so a later plot with a
    # different base_alpha doesn't change this one
    fig_other, ax_other = plt.subplots()
    boxes_other = d.plot_spreading_box(ax=ax_other, base_alpha=0.5)
    fig_other.canvas.draw()
    fig.canvas.draw()
    np.testing.assert_allclose(boxes_other.get_facecolor()[[0, 2, 4], 3], 0.5)
    np.testing.assert_allclose(boxes.get_facecolor()[:, 3], alpha)
    plt.close("all")
