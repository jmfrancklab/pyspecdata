"""Domain Colored Coherence Transfer (DCCT) function allows us to 
visualize the complex-valued data, as well as the formalization of the 
coherence transfer dimensions using domain coloring plotting. 
"""

from numpy import r_, nan
import numpy as np
from .core import ndshape, nddata
from .general_functions import strm, process_kwargs
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import FancyArrow, Circle
from matplotlib.lines import Line2D
from matplotlib.transforms import (
    ScaledTranslation,
    IdentityTransform,
    blended_transform_factory,
)
from pyspecdata.plot_funcs.image import imagehsv
import matplotlib.ticker as mticker
import logging


def DCCT(
    this_nddata,
    this_fig_obj,
    custom_scaling=False,
    grid_bottom=0.0,
    bbox=[0.02, 0.1, 0.94],
    grid_top=1.0,
    top_pad=0.05,
    total_spacing=0.055,
    vert_label_spacer=50,
    shareaxis=False,
    diagnostic=False,
    cmap=None,
    pass_frq_slice=False,
    just_2D=False,
    scaling_factor=1,
    max_coh_jump={"ph1": 1, "ph2": 2},
    direct="t2",
    plot_title="DCCT",
    **kwargs,
):
    """DCCT plot

    Parameters
    ==========
    this_nddata : nddata
        data being plotted
    this_fig_obj : figure
        size/type of figure to be plotted on
    custom_scaling : boolean
        allows user to scale the intensity of data presented
    bbox : list
        contains the following:
        :bbox[0]: int
            Left hand side padding between the left side of the figure
            and the left side of the decorations.
        :bbox[1]: int
            Distance between the bottom of the figure and the bottom of the
            lowest axes object (in Figure coordinates)
        :bbox[2]: int
            Distance between the right most side of the axes objects and the
            left most side of the decorations (in Figure coordinates)
    grid_bottom : float
        Figure coordinates
        y-coordinate for bottom of DCCT plot.
    grid_top : float
        Figure coordinates
        y-coordinate top of grid in figure
    top_pad : float
        Figure coordinates
        Distance between top of figure and top of grid
    total_spacing : float
        Figure coordinates
        Spacing between phase cycle dimensions
    vert_label_spacer : int
        Display coordinates
        Spacing between vertical lines that label the coherence pathways
    shareaxis : boolean
        Subplots scale together, but currently, this means there
        must be tick labels on both top and bottom
    cmap : str
        Color mapping if specified
    pass_frq_slice : boolean
        If true will show the frequency sliced out with hatching
    just_2D : boolean
        If true will only return axis coordinates/shape NOT ax_list
    scaling_factor : float
        If using custom scaling this allows user to set the scaling
        factor
    max_coh_jump : dict
        Maximum allowed transitions for each phase cycle
    direct : str
        Name of direct axis
    plot_title : str
        Title for DCCT plot
    """
    x = []  # List to put direct dim into
    y = []  # List to put indirect dims into
    my_data = this_nddata.C
    ordered_labels = {}

    def ax_to_fig(thisax, thisfig):
        """Transformation to convert axes coordinates of thisax into figure
        coordinates of thisfig"""
        return thisax.transAxes + thisfig.transFigure.inverted()

    # {{{ Generate alias labels - goes to scientific fn
    for this_dim in [j for j in my_data.dimlabels if j.startswith("ph")]:
        if my_data.get_ft_prop(this_dim):
            n_ph = ndshape(my_data)[this_dim]
            this_max_coh_jump = max_coh_jump[this_dim]
            all_possibilities = np.empty(
                (int((2 * this_max_coh_jump + 1) / n_ph) + 1) * n_ph
            )  # on reviewing, I *believe* this this is designed to fit the
            #    array from -this_max_coh_jump to +this_max_coh_jump, but it
            #    needs to round up to the closest multiple of n_ph
            all_possibilities[:] = nan
            all_possibilities[: this_max_coh_jump + 1] = r_[
                0 : this_max_coh_jump + 1
            ]  # label the positive jumps in order
            all_possibilities[-this_max_coh_jump:] = r_[
                -this_max_coh_jump:0
            ]  # and alias the negative ones into the correct locations
            all_possibilities = all_possibilities.reshape(
                (-1, n_ph)
            )  # now, reshape according to the number of dimensions we actually
            #    have for discriminating
            labels_in_order = []
            for j in range(n_ph):
                temp = all_possibilities[
                    :, j
                ]  # grab the columns, which are the labels for all aliases
                #    that belong at this index
                if j == 0:
                    temp = ", ".join([
                        "%d" % j for j in np.sort(temp[np.isfinite(temp)])
                    ])
                else:
                    temp = ", ".join([
                        "%+d" % j for j in np.sort(temp[np.isfinite(temp)])
                    ])
                if len(temp) == 0:
                    temp = "X"
                labels_in_order.append(temp)
            # }}}
            ordered_labels[this_dim] = labels_in_order
        else:
            ordered_labels[this_dim] = [
                "0" if j == 0.0 else f"{j}" for j in my_data.getaxis(this_dim)
            ]
        # ordered_labels now contains a list of the labels for each index
        # of the dimension, in order
    # }}}
    real_data = False
    if cmap is not None:
        assert all(
            np.isclose(my_data.data.imag, 0)
        ), "In order to use a color map, you must pass real data"
        if type(cmap) is str:
            cmap = plt.get_cmap(cmap)
            my_data.data = my_data.data.real
            real_data = True
    my_data.human_units()
    grid_bottom += bbox[1]  # y coord for bottom of bottom axein figure coords
    grid_top -= top_pad  # y coord for top of top ax in figure coords
    a_shape = ndshape(this_nddata)
    num_dims = len(a_shape.dimlabels[:-2])
    divisions = []

    # should be looping in backward order from printed shape
    for j, thisdim in enumerate(a_shape.dimlabels[::-1][2:]):
        old = [j / 2.0 for j in divisions]
        divisions = (old + [1]) * (a_shape[thisdim] - 1) + old
        logging.debug(strm("for", thisdim, "I get", divisions))
    divisions = [j * total_spacing / sum(divisions) for j in divisions]
    axes_height = (grid_top - grid_bottom - total_spacing) / np.prod(
        a_shape.shape[:-2]
    )
    axes_bottom = np.cumsum([
        axes_height + j for j in divisions
    ])  # becomes ndarray
    axes_bottom = r_[0, axes_bottom]
    axes_bottom += grid_bottom
    fig = this_fig_obj
    ax_list = []
    majorLocator = lambda: mticker.MaxNLocator(
        nbins="auto", steps=[1, 2, 2.5, 5, 10]
    )
    minorLocator = lambda: mticker.AutoMinorLocator(n=5)
    LHS_labels, _ = fig.transFigure.inverted().transform(
        (vert_label_spacer * num_dims, 0)
    )
    width = bbox[2] - LHS_labels
    for j, b in enumerate(axes_bottom):
        if j != 0 and shareaxis:
            ax_list.append(
                plt.axes(
                    [LHS_labels + bbox[0], b, width, axes_height],
                    sharex=ax_list[0],
                    sharey=ax_list[0],
                )
            )  # in figure coords: x, y, width, height
        else:
            ax_list.append(
                plt.axes([LHS_labels + bbox[0], b, width, axes_height])
            )  # lbwh
    # {{{ make blended transform for plotting coherence transfer labels
    ax_x0, ax_y0 = ax_to_fig(ax_list[0], fig).transform(
        r_[0, 0]
    )  # bottom left corner of bottom ax in fig coord
    dx = LHS_labels + bbox[0]  # x coord in figure coord
    total_scale_transform = IdentityTransform() + ScaledTranslation(
        dx, ax_y0, fig.transFigure
    )
    total_trans = blended_transform_factory(
        total_scale_transform, fig.transFigure
    )
    # }}}
    # {{{ transform for arrows
    arrow_trans = IdentityTransform() + ScaledTranslation(
        ax_x0, ax_y0, fig.transFigure
    )
    # }}}
    # {{{ adjust tick settings -- AFTER extents are set
    # {{{ bottom subplot
    ax_list[0].xaxis.set_major_locator(majorLocator())
    ax_list[0].xaxis.set_minor_locator(minorLocator())
    ax_list[0].set_ylabel(None)
    ax_list[0].set_xlabel(
        my_data.unitify_axis(my_data.dimlabels[-1]), labelpad=20
    )
    # }}}
    # {{{ intermediate subplots
    for j in range(1, len(axes_bottom) - 1):
        ax_list[j].xaxis.set_ticks([])
        ax_list[j].get_xaxis().set_visible(False)
        ax_list[j].set_xlabel(None)
    # }}}
    # {{{ top subplot
    ax_list[-1].xaxis.set_major_locator(majorLocator())
    # for the minor ticks, use no labels; default NullFormatter
    ax_list[-1].xaxis.set_minor_locator(minorLocator())
    ax_list[-1].xaxis.tick_top()
    if not shareaxis:
        ax_list[-1].set_xticklabels([])
    # }}}
    # {{{ all subplots
    for j in range(0, len(axes_bottom)):
        ax_list[j].set_ylabel(a_shape.dimlabels[-2])
        inner_dim = a_shape.dimlabels[-2]
        inner_dim = str(inner_dim)
        if inner_dim == "ph2":
            logging.debug("Inner dimension is phase cycling dimension")
            ax_list[j].yaxis.set_major_formatter("ph2")
            ax_list[j].yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        else:
            ax_list[j].yaxis.set_minor_locator(minorLocator())
            ax_list[j].yaxis.set_ticks_position("both")
        for tick in ax_list[j].get_yticklabels():
            tick.set_rotation(0)
            tick.set_va("center")
    # }}}
    # }}}
    # {{{ smoosh dataset if needed
    if len(a_shape.dimlabels) > 3:
        A = this_nddata.C.smoosh(
            a_shape.dimlabels[:-2], "smooshed", noaxis=True
        )
        A.reorder("smooshed", first=True)
    else:
        A = this_nddata.C
        A.rename(a_shape.dimlabels[:-2][0], "smooshed")

    # }}}
    # {{{ vertical lines for coh paths
    def draw_span(
        ax1,
        ax2,
        label,
        this_label_num,
    ):
        """Place the vertical lines and labels of the coherence transfer
        pathways"""
        label_spacing = (
            this_label_num
            + 1  # plus one so the first horizontal isn't placed at 0
            #      (overlapping with the spine of the indirect axis)
        ) * vert_label_spacer  # will space the vertical lines along x approp.
        x1_disp = LHS_labels + bbox[0] - label_spacing  # x coord is the left
        #                                                 side of the axis
        #                                                 minus the spacing for
        #                                                 text/ticks
        #                                                 (label_spacing)
        # {{{ Take y coordinate of top and bottom of axes objects to get the 2
        # points for drawing the lines. To be exact I pull this from the axes
        # objects themselves.
        _, y1_fig = ax_to_fig(ax1, fig).transform(r_[0, 0.95])
        _, y2_fig = ax_to_fig(ax2, fig).transform(r_[0, 0.05])
        # }}}
        # for scaled translation, x coords should be display and y coords
        # should be figure
        lineA = lines.Line2D(
            [x1_disp, x1_disp],
            [y1_fig, y2_fig],
            linewidth=1,
            color="k",
            transform=total_trans,
            clip_on=False,
        )
        plt.text(
            x1_disp,  # same x coord as the line to keep simple
            (y2_fig + y1_fig) / 2,  # put in middle of line
            label,
            va="center",
            ha="right",
            rotation=90,
            transform=total_trans,
            color="k",
        )
        fig.add_artist(lineA)

    # }}}
    # {{{ place arrows and dimname labels
    def place_labels(
        ax1,
        label,
        label_placed,
        this_label_num,
        check_for_label_num=True,
        arrow_width_px=4.0,
    ):
        x_axorigindisp, y_axorigindisp = ax1.transAxes.transform(r_[0, 0])
        if not check_for_label_num or not label_placed[this_label_num]:
            # {{{ determine the x and y position of the label in display coords
            if check_for_label_num:
                # the labels of the outer dimensions
                label_spacing = (this_label_num + 1) * vert_label_spacer
                # Calculate coord for base of arrow
                x_textdisp = LHS_labels + bbox[0] - label_spacing
            else:
                # same as above, but determine text
                # position based on tick labels
                label = my_data.unitify_axis(my_data.dimlabels[-2])
                # from here https://stackoverflow.com/questions/44012436/pytho\
                # n-matplotlib-get-position-of-xtick-labels
                # then searching for BBox docs
                logging.debug(
                    strm(
                        "tick locations",
                        [
                            j.get_window_extent().bounds
                            for j in ax1.get_yticklabels()
                        ],
                    )
                )
                x_textdisp = (
                    LHS_labels + bbox[0]
                )  # bottom left corner of bottom axes in fig
            # tick length is a nice consistent distance to push the arrows out
            # slightly to avoid confusion
            tick_length = [
                j.get_window_extent().bounds for j in ax1.get_yticklabels()
            ][-1][-1]
            x_textdisp -= tick_length
            x_textdisp -= arrow_width_px
            y_textdisp = -25.0
            if diagnostic:
                a = Circle(
                    (x_textdisp, y_axorigindisp - y_textdisp),
                    3,
                    clip_on=False,
                    transform=None,
                    color="b",
                )
                fig.add_artist(a)
            # }}}
            AnArrow = FancyArrow(
                x_textdisp,
                y_textdisp,
                2,
                5,
                width=arrow_width_px,
                clip_on=False,
                transform=arrow_trans,
                alpha=0.1,
                color="k",
            )
            # could do fancier w/ the following, but need to mess w/ width
            # parameters
            # arrow_base = r_[x_arrowbase_fig-arrow_width/2, y_arrowbase_fig]
            # a = FancyArrowPatch(arrow_base, arrow_base+r_[dx, dy],
            #        arrowstyle='|-|',
            #        alpha=0.1,  color='k')
            fig.add_artist(AnArrow)
            if diagnostic:
                a = Line2D(
                    [x_textdisp, x_textdisp + dx],
                    [y_textdisp, y_textdisp + dy],
                    transform=fig.transFigure,
                    color="r",
                )
                fig.add_artist(a)
                a = Circle(
                    (x_textdisp, y_textdisp),
                    0.01,
                    clip_on=False,
                    transform=fig.transFigure,
                    color="r",
                )
                fig.add_artist(a)
            x_textfig = x_textdisp + arrow_width_px
            y_textfig = y_textdisp - 5.0
            plt.text(
                x_textfig,
                y_textfig,
                label,
                va="top",
                ha="right",
                rotation=45,
                clip_on=False,
                transform=arrow_trans,
                color="k",
            )
            if check_for_label_num:
                label_placed[this_label_num] = 1

    # }}}
    label_placed = np.zeros(num_dims)
    imagehsvkwargs = {}
    for k, v in list(kwargs.items()):
        if k in ["black", "logscale"]:
            imagehsvkwargs[k] = kwargs.pop(k)

    # {{{ Place actual data into appropriate axes
    for j in range(len(ax_list)):
        spacing, ax, x_first, origin, renumber = process_kwargs(
            [
                ("spacing", 1),
                ("ax", ax_list[j]),
                ("x_first", False),
                ("origin", "lower"),
                ("renumber", None),
            ],
            kwargs,
            pass_through=True,
        )
        if isinstance(x, list):
            x = np.array(my_data.getaxis(direct))
        if isinstance(y, list):
            y = np.array(my_data.getaxis(my_data.dimlabels[-2]))
        if len(x) == 0:
            x = [1, A.data.shape[1]]
        else:
            x = x.flatten()
        if len(y) == 0:
            y = [1, A.data.shape[0]]
        else:
            y = y.flatten()
        dx = (x[-1] - x[0]) / len(x)
        dy = (y[-1] - y[0]) / len(y)
        if origin == "lower":
            myext = (
                x[0] - dx / 2.0,
                x[-1] + dx / 2.0,
                y[0] - dy / 2.0,
                y[-1] + dy / 2.0,
            )
        elif origin == "upper":
            myext = (
                x[0] - dx / 2.0,
                x[-1] + dx / 2.0,
                y[-1] + dy / 2.0,
                y[0] - dy / 2.0,
            )
        elif origin == "flip":
            # {{{ need to flip
            myext = (
                x[-1] + dx / 2.0,
                x[0] - dx / 2.0,
                y[-1] + dy / 2.0,
                y[0] - dy / 2.0,
            )
            # }}}
        else:
            raise ValueError(
                "I don't understand the value you've set for the origin"
                " keyword argument"
            )
        kwargs["origin"] = (
            origin  # required so that imshow now displays the image correctly
        )

        if real_data:
            kwargs["cmap"] = cmap
            K = A["smooshed", j].data / abs(A).data.max()
        else:
            if custom_scaling:
                K = imagehsv(
                    A["smooshed", j].data,
                    **imagehsvkwargs,
                    scaling=scaling_factor,
                )
            if not custom_scaling:
                K = imagehsv(
                    A["smooshed", j].data,
                    **imagehsvkwargs,
                    scaling=abs(A).data.max(),
                )
        plt.sca(ax_list[j])
        plt.imshow(K, extent=myext, **kwargs)
        ax_list[j].set_ylabel(None)
        if pass_frq_slice:
            frq_slice = []
            start_y = A.getaxis(A.dimlabels[1])[0]
            stop_y = A.getaxis(A.dimlabels[1])[-1]
            ax_list[j].fill(
                [x[0], frq_slice[0], frq_slice[0], x[0]],
                [start_y, start_y, stop_y, stop_y],
                fill=None,
                alpha=0.7,
                hatch="//",
            )
            ax_list[j].fill(
                [frq_slice[-1], x[-1], x[-1], frq_slice[-1]],
                [start_y, start_y, stop_y, stop_y],
                fill=None,
                alpha=0.7,
                hatch="//",
            )
    # }}}
    # to drop into ax_list, just do
    # A.smoosh(a_shape.dimlabels, 'smooshed', noaxis=True)
    # in ax_list[0] put A['smooshed',0], etc
    idx = nddata(r_[0 : np.prod(a_shape.shape[:-2])], [-1], ["smooshed"])
    idx.chunk("smooshed", a_shape.dimlabels[:-2], a_shape.shape[:-2])
    remaining_dim = a_shape.dimlabels[:-2]
    depth = num_dims

    # {{{ Place decorations on figure
    def decorate_axes(idx, remaining_dim, depth):
        thisdim = remaining_dim[0]
        logging.debug(strm("This is remaining dim", remaining_dim))
        logging.debug(strm("This dim is", thisdim))
        depth -= 1
        for j in range(a_shape[thisdim]):
            idx_slice = idx[thisdim, j]
            logging.debug(
                strm("For", thisdim, "element", j, idx_slice.data.ravel())
            )
            first_axes = ax_list[idx_slice.data.ravel()[0]]
            last_axes = ax_list[idx_slice.data.ravel()[-1]]
            if j == 0:
                draw_span(
                    last_axes,
                    first_axes,
                    "%s" % ordered_labels[thisdim][0],
                    this_label_num=depth,
                )
            else:
                draw_span(
                    last_axes,
                    first_axes,
                    "%s" % ordered_labels[thisdim][j],
                    this_label_num=depth,
                )
            place_labels(
                ax_list[0],
                "%s" % my_data.unitify_axis("%s" % thisdim),
                label_placed,
                this_label_num=depth,
            )
            new_remaining_dim = remaining_dim[1:]
            if len(remaining_dim) > 1:
                decorate_axes(idx_slice, new_remaining_dim, depth)

    decorate_axes(idx, remaining_dim, depth)
    place_labels(
        ax_list[0],
        "%s" % (a_shape.dimlabels[-2]),
        label_placed,
        this_label_num=depth - 1,
        check_for_label_num=False,
    )
    plt.title(plot_title)
    if just_2D:
        return (
            bbox[0] + LHS_labels,
            axes_bottom[0],
            width,
            axes_bottom[-1] - top_pad,
        )
    else:
        return (
            ax_list,
            bbox[0] + LHS_labels,
            axes_bottom[-1] + axes_height,
            width,
            top_pad - (1 - bbox[2] - bbox[0]),
        )
