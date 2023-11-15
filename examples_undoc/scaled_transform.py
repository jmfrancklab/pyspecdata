"""
Scaled Transform
================

Transforms added figures or text from the display coordinates of the 
original plot to a new sizing using the IdentityTransform.
"""
from pylab import *
from matplotlib.transforms import ScaledTranslation, IdentityTransform
from matplotlib.patches import FancyArrow, FancyArrowPatch, Circle
import matplotlib.ticker as mticker

fig = figure()
plot(r_[0:10])
ax = gca()
axis_to_figure = ax.transAxes + fig.transFigure.inverted()
ax_x, ax_y = axis_to_figure.transform(
    r_[0, 0]
)  # get the location of the axis in figure coordinates
a = Circle(
    (-10, -10),
    3,  # these are display coords relative to the axis corner
    clip_on=False,
    transform=(
        IdentityTransform()  # use display coords
        + ScaledTranslation(
            ax_x, ax_y, fig.transFigure
        )  # but apply a translation that would put 0,0 at the corner of the axis
    ),
    color="b",
)
fig.add_artist(a)
majorLocator = lambda: mticker.MaxNLocator(nbins="auto", steps=[1, 2, 2.5, 5, 10])
minorLocator = lambda: mticker.AutoMinorLocator(n=5)
ax.xaxis.set_major_locator(majorLocator())
ax.xaxis.set_minor_locator(minorLocator())
show()
