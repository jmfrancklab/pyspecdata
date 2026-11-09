"""Plot points as boxes that spread a fixed amount of pigment over their
error."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgba


class SpreadingBoxes(PolyCollection):
    """Boxes that keep the amount of pigment per point fixed.

    Each point starts from a baseline square whose side is *min_pixels*
    pixels on screen (found from the display transform).  Along any
    direction with an error σ, the box is widened to 2σ when that is larger
    than the baseline.  The opacity is
    then lowered so that (area × opacity) equals the baseline square at
    50 % opacity: the same amount of pigment, spread over an area governed
    by the error.  The boxes are recomputed at every draw, so they stay
    correct as the axes rescale, zoom or resize.
    """

    def _update_boxes(self):
        ax = self.axes
        assert (
            ax.get_xscale() == "linear" and ax.get_yscale() == "linear"
        ), "spreading boxes need linear axes"
        # the baseline square, converted from pixels to data units along x
        # and y at the same time
        baseline = np.abs(
            np.diff(
                ax.transData.inverted().transform(
                    [[0.0, 0.0], [self._min_pixels, self._min_pixels]]
                ),
                axis=0,
            )[0]
        )
        width = np.maximum(baseline[0], 2 * self._x_error)
        height = np.maximum(baseline[1], 2 * self._y_error)
        corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) / 2
        self.set_verts(
            self._centers[:, None, :]
            + corners[None, :, :] * np.stack([width, height], axis=-1)[:, None]
        )
        rgba = np.tile(self._base_rgba, (len(self._centers), 1))
        rgba[:, 3] =  self.base_alpha * baseline[0] * baseline[1] / (width * height)
        mask = rgba[:,3] > 1
        rgba[mask,3] = 1
        self.set_facecolor(rgba)

    def draw(self, renderer):
        # legend handles are plain copies that carry no box data
        if hasattr(self, "_centers"):
            self._update_boxes()
        super().draw(renderer)


def plot_spreading_box(self, *args, min_pixels=3, ax=None, base_alpha=5, **kwargs):
    r"""Plot each point as a box whose size shows its error.

    Think of using the same amount of pigment for every point, and spreading
    it out over an area that is governed by the error.  A point with no (or
    small) error is drawn as a small square, 50 % opaque, *min_pixels*
    pixels on a side.  If the data have an error σ, the box grows to a
    height of 2σ when that exceeds the square; if the coordinate has an
    error, the box grows to a width of 2σ in the same way.  A bigger box is
    proportionally more transparent, so (area × opacity) stays the same, and
    points with large errors fade out instead of dominating the plot.  The
    boxes have no edges, since an edge would add pigment and spoil the
    effect.

    The data are first drawn with :func:`~pyspecdata.plot`, so the axes are
    autoscaled (including the error bars) and labelled as usual, from the
    axis names and units; that drawing is then replaced by the boxes, and
    the limits are widened to show whole boxes.

    Parameters
    ==========
    *args :
        Passed to :func:`~pyspecdata.plot` (*e.g.* a format string).
    min_pixels : float
        Side of the baseline square, in pixels.
    ax : matplotlib Axes, optional
        Axes to draw on (default: current axes).
    **kwargs :
        Passed to :func:`~pyspecdata.plot` (*e.g.* ``color``, ``label`` or
        ``human_units``).

    Returns
    =======
    SpreadingBoxes
        The collection of boxes that was added to *ax*.
    """
    from ..core import plot

    assert len(self.dimlabels) == 1, (
        "plot_spreading_box needs one-dimensional data -- loop over (or"
        " select along) the other dimensions"
    )
    if ax is None:
        ax = plt.gca()
    (dimname,) = self.dimlabels
    drawn = plot(self, *args, ax=ax, **kwargs)
    # the drawn line tells us the color, the label, and the coordinates as
    # actually plotted (human units may have rescaled them)
    if hasattr(drawn, "lines"):  # an errorbar container
        line = drawn.lines[0]
        drawn = [drawn]
    else:  # a list of lines
        line = drawn[0]
    color = line.get_color()
    centers = np.column_stack([line.get_xdata(), line.get_ydata()])
    x_scale = _scale_between(self.getaxis(dimname), centers[:, 0])
    y_scale = _scale_between(self.data.real, centers[:, 1])
    for artist in drawn:
        artist.remove()
    ax.autoscale_view()
    boxes = SpreadingBoxes(
        [], edgecolors="none", linewidths=0, label=kwargs.get("label")
    )
    SpreadingBoxes.base_alpha = base_alpha
    boxes._centers = centers
    boxes._base_rgba = np.array(to_rgba(color))
    boxes._min_pixels = min_pixels
    boxes._x_error = _error_or_zero(self.get_error(dimname), len(centers))
    boxes._x_error = boxes._x_error * x_scale
    boxes._y_error = _error_or_zero(self.get_error(), len(centers))
    boxes._y_error = boxes._y_error * y_scale
    # the removed drawing set the data limits and so the pixel scale; now
    # widen the limits once so that whole boxes are visible (the boxes are
    # recomputed for the final scale when drawn)
    ax.add_collection(boxes, autolim=False)
    boxes._update_boxes()
    corners = np.concatenate([path.vertices for path in boxes.get_paths()])
    ax.update_datalim(corners)
    ax.autoscale_view()
    # plot may have switched autoscaling off along an axis -- then extend
    # its current limits instead
    for autoscaled, get_limits, set_limits, column in [
        (ax.get_autoscalex_on(), ax.get_xlim, ax.set_xlim, 0),
        (ax.get_autoscaley_on(), ax.get_ylim, ax.set_ylim, 1),
    ]:
        if autoscaled:
            continue
        first, last = get_limits()
        low = min(first, last, corners[:, column].min())
        high = max(first, last, corners[:, column].max())
        set_limits((low, high) if first <= last else (high, low))
    return boxes


def _scale_between(original, plotted):
    """Factor that human units applied to go from *original* to *plotted*."""
    reference = np.abs(original).max()
    if reference == 0:
        return 1.0
    return np.abs(plotted).max() / reference


def _error_or_zero(error, length):
    if error is None:
        return np.zeros(length)
    return np.broadcast_to(np.abs(error), (length,)).astype(float)
