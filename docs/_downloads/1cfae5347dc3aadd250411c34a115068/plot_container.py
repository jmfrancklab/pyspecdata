"""
Plot Container
==============

This examples shows how you would *currently* generate a plot that shows
projections of a 2D image. 

While it does work, this is for development purposes -- while the figure list
currently works well for collecting a series of individual matplotlib plots
along with text and generating a notebook presentation from them, things that
involve multiple matplotlib Axes, like this one, are unnecessarily complicated.

The comments describe how the new proposed objects will work.
"""
from numpy import *
from pyspecdata import *
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib import transforms

# {{{ from spine_placement_demo
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(("outward", 10))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color("none")  # don't draw spine

    # turn off ticks where there is no spine
    if "left" in spines:
        ax.yaxis.set_ticks_position("left")
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if "bottom" in spines:
        ax.xaxis.set_ticks_position("bottom")
    if "top" in spines:
        ax.xaxis.set_ticks_position("top")
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


# }}}

# based off an example from matplotlib gallery

# the random data
test_data = nddata(
    random.normal(size=100 * 100) + 1j * random.normal(size=100 * 100),
    [100, 100],
    ["x", "y"],
)
test_data.setaxis("x", "#").setaxis("y", "#")
test_data.set_units("x", "s").set_units("y", "m")
test_data.reorder("y")

# definitions for the axes
left, width = 0.12, 0.65
bottom, height = 0.15, 0.65
padding = 0.05
bottom_side = bottom + height + padding
left_side = left + width + padding

rect_scatter = [left, bottom, width, height]
rect_top = [left, bottom_side, width, 1.0 - bottom_side - padding]
rect_right = [left_side, bottom, 1.0 - left_side - padding, height]

with figlist_var() as fl:
    # A lot of extra junk here b/c figure list not set up for multiple matplotlib axes.
    #
    # The new idea is to create a container, which in this case would be all
    # three of these  axes together.
    #
    # When we drop data into the container, it would automatically decide what
    # to do with it, based on its shape (here, 2D data would go in the center,
    # while 1D data would go to one of the outside plots), and would directly
    # use matplotlib commands
    # -- so, there would be no need to pass the figure to "next", set the
    # various titles to 0, etc.
    #
    # also, the container would just hold a list of data until its actually
    # read to render the plots -- it would *then* do the human units thing, so
    # that there were no human units error
    #
    # the domain overview plot would be 1 matplotlib axes object for each
    # coherence pathway and would probably take advantage of this command --
    # https://matplotlib.org/examples/axes_grid/simple_axesgrid2.html

    # start with a rectangular Figure
    fig = plt.figure(1, figsize=(9, 5.56))
    axCentral = plt.axes(rect_scatter)
    axRight = plt.axes(rect_right)
    axTop = plt.axes(rect_top)
    fl.next("test figure", fig=fig, ax=axCentral)

    fl.image(test_data, ax=axCentral, human_units=False)
    axCentral.set_title("")
    fl.plot(test_data.C.sum("y"), ax=axTop, human_units=False)
    axTop.autoscale(enable=True, tight=True)  # axis tight
    axTop.set_ylabel("")
    axTop.set_xlabel("")
    adjust_spines(axTop, ["left"])
    base = axRight.transData
    rot = transforms.Affine2D().rotate_deg(90)
    fl.plot(test_data.C.sum("x"), transform=(rot + base), ax=axRight, human_units=False)
    axRight.set_ylabel("")
    axRight.set_xlabel("")
    axRight.autoscale(enable=True, tight=True)  # axis tight
    adjust_spines(axRight, ["bottom"])
