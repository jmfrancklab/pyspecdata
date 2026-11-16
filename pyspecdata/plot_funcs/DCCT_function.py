"""Domain Colored Coherence Transfer (DCCT) plots coherence pathways.

See :func:`pyspecdata.plot_funcs.image.imagehsv` for examples of domain
colored spin-echo data.  This routine displays coherence pathways using
domain coloring.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import SubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.transforms import (
    ScaledTranslation,
    IdentityTransform,
    blended_transform_factory,
)
from pyspecdata.plot_funcs.image import imagehsv, image
from ..general_functions import strm


def _label_arrow(fig, label, x_tip, arrow_dx, arrow_dy, width, transform):
    """Draw a grey arrow whose tip is `x_tip` display units to the right of
    the lower left corner of the bottom axes (and just below it), with
    `label` at its tail.  `width` is the shaft width in points."""
    y_tip = -2
    tail = (x_tip - arrow_dx, y_tip - arrow_dy)
    fig.add_artist(
        FancyArrowPatch(
            tail,
            (x_tip, y_tip),
            clip_on=False,
            transform=transform,
            alpha=0.1,
            color="k",
            # the shaft of the default ("simple") arrow style is 0.2
            # times the mutation scale
            mutation_scale=width / 0.2,
        )
    )
    fig.text(
        *tail,
        label,
        va="top",
        ha="right",
        rotation=np.degrees(np.arctan2(arrow_dy, arrow_dx)),
        clip_on=False,
        transform=transform,
        color="k",
    )


# SINGLE_USE_EXCEPTION -- exported API
def DCCT(
    this_nddata,
    fig=None,
    custom_scaling=False,
    bbox=[0.05, 0.1, 0.90, 0.8],
    gap=0.1,
    horiz_label_spacer=35,
    cmap=None,
    shareaxis=True,
    coherence_unmask_fn=None,
    scaling_factor=1,
    max_coh_jump=None,
    direct="t2",
    title=None,
    arrow_dx=30,
    arrow_dy=30,
    arrow_width_px=None,
    no_aliasing_calc=False,
    interpolation=None,
    **kwargs,
):
    """DCCT plot.

    The last dimension is the direct dimension (the *x* axis of each image),
    the next-to-last ("inner") dimension is the *y* axis of each image, and
    one image is generated for each element of the remaining ("outer")
    dimensions.
    The images are stacked vertically, and the outer dimensions are labeled
    by nested brackets to the left of the images.

    Parameters
    ==========
    this_nddata : nddata
        data being plotted
    fig : figure
        figure to be plotted on
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
        :bbox[3]: int
            Distance from bottom of the bottom axes object to top of top
            axes object
    bbox : matplotlib.girdspec.SubplotSpec
        If you specify a gridspec element, it will use this to generate a
        bbox, as above.
    gap : float
        Figure coordinates
        Spacing between coherence transfer pathways
    horiz_label_spacer : int
        Display coordinates
        Spacing between vertical lines that label the coherence pathways
    shareaxis : boolean
        If True (default), the images share their axes, so that zooming in
        on one image zooms in on all of them.
        If False, each image can be zoomed separately.
    coherence_unmask_fn : function
        Used to show a coherence mask:
        the coherence pathways that are masked out are covered with a hatch.
        This is the same type of function that is used by
        `pyspecProcScripts.correl_align`:
        it takes an nddata that contains only the phase cycling (coherence
        transfer) dimensions and is filled with zeros,
        and it returns that nddata with the pathways that are *not* masked
        set to 1.
        (The nddata that is passed has the same `coherence_pathway`
        property as the data.)
    cmap : str
        Color mapping if specified -- only allowed for real data.
    scaling_factor : float
        If using custom scaling this allows user to set the scaling
        factor
    max_coh_jump : dict
        Maximum allowed transitions for each phase cycle
    direct : str
        Name of direct axis
    interpolation: str (default None)
        Override the auto-chosen interpolation.
    title : str
        Title for DCCT plot
    arrow_dx, arrow_dy : float
        Display coordinates
        Horizontal and vertical extent of the grey arrows that label the
        dimensions.
    arrow_width_px : float
        Display coordinates
        Width of the shaft of the grey arrows that label the dimensions
        (the arrowhead scales with it).
        Defaults to 4 points.
    no_aliasing_calc : boolean
        Label the phase cycling dimensions with their axis coordinates,
        rather than with the (aliased) coherence changes.

    Returns
    =======
    ax_list: list
        list of axes objects generated in DCCT plot
    transform_dict: dict
        dictionary containing the unique transforms generated in
        this function for scaled translations
    """
    # {{{ explain what to do about options that have been removed
    for thisoption, explanation in [
        (
            "pass_frq_slice",
            "It was intended to hatch out part of the plot, but it never"
            " worked.  To hatch out the coherence pathways that are masked,"
            " pass coherence_unmask_fn.",
        ),
        (
            "diagnostic",
            "It was intended to show extra diagnostic plots, but it never"
            " did anything, so just leave it out.",
        ),
        (
            "ax",
            "DCCT always generates its own axes, and this was ignored.  To"
            " control where the plot goes, pass fig and/or bbox (which can"
            " be a gridspec element).",
        ),
        (
            "spacing",
            "It was ignored.  To control the spacing between the images,"
            " use gap, and to control the spacing between the labels, use"
            " horiz_label_spacer.",
        ),
        (
            "x_first",
            "It was ignored.  The last dimension is always the x axis, so"
            " reorder your data, or use direct to say which dimension"
            " should be the x axis.",
        ),
        (
            "renumber",
            "It was ignored.  If you want a dimension to be labeled by"
            " index, set the coordinates of that axis before calling DCCT.",
        ),
    ]:
        if thisoption in kwargs:
            raise ValueError(
                f"The {thisoption} option of DCCT is deprecated, and has"
                " been removed.  " + explanation
            )
    # }}}
    # {{{ determine the figure and the bbox
    if fig is None:
        if hasattr(bbox, "get_gridspec") and hasattr(
            bbox.get_gridspec(), "figure"
        ):
            fig = bbox.get_gridspec().figure
        else:
            raise RuntimeError(
                "Behavior of DCCT has been changed!  You either need to"
                " (1) use the figlist_var DCCT method (2) pass a bbox or"
                " (3) pass a fig kwarg"
            )
    if isinstance(bbox, SubplotSpec):
        bbox = bbox.get_position(fig).bounds
    bbox_left, bbox_bottom, bbox_width, bbox_height = bbox
    # }}}
    if len(this_nddata.dimlabels) < 3:
        raise ValueError(
            "DCCT needs at least 3 dimensions -- for 2D data, just use image"
        )
    # we need a copy because human_units changes the axis coordinates
    my_data = this_nddata.C
    my_data.reorder(direct, first=False)
    my_data.human_units()
    outer_dims = my_data.dimlabels[:-2]
    inner_dim, direct = my_data.dimlabels[-2:]
    outer_shape = [my_data.shape[j] for j in outer_dims]
    n_axes = int(np.prod(outer_shape))
    # {{{ all phase cycles need to go exactly once around the circle, in
    #     evenly spaced steps
    for thisdim in my_data.dimlabels:
        if thisdim.startswith("ph") and my_data.shape[thisdim] > 1:
            if my_data.getaxis(thisdim) is None:
                raise ValueError(
                    f"The phase cycling dimension {thisdim} needs an axis"
                    " that gives the phases, in cycles"
                )
            # in the phase cycling domain, the steps need to add up to one
            # cycle, which is equivalent to a spacing of 1 in the coherence
            # domain
            if my_data.get_ft_prop(thisdim):
                required_step = 1
            else:
                required_step = 1 / my_data.shape[thisdim]
            if not np.allclose(
                np.diff(my_data.getaxis(thisdim)), required_step
            ):
                raise ValueError(
                    f"The phase cycle along {thisdim} needs to go exactly"
                    " once around the circle in evenly spaced steps (e.g. 0,"
                    " 0.25, 0.5, 0.75 cycles), which corresponds to a"
                    " spacing of 1 in the coherence domain, but along"
                    f" {thisdim} your data has "
                    + (
                        "coherence changes"
                        if my_data.get_ft_prop(thisdim)
                        else "phases"
                    )
                    + f" of {my_data.getaxis(thisdim)}"
                )
    # }}}
    # {{{ Generate the labels for each element along each indirect dimension
    #     -- for phase cycling dimensions, these are the alias labels
    dim_labels = {}
    for thisdim in my_data.dimlabels[:-1]:
        if not thisdim.startswith("ph"):
            if my_data.getaxis(thisdim) is None:
                my_data.set_axis(thisdim, "#")
            dim_labels[thisdim] = [
                f"{j:.4g}" if isinstance(j, (float, np.floating)) else f"{j}"
                for j in my_data.getaxis(thisdim)
            ]
        elif no_aliasing_calc:
            dim_labels[thisdim] = [f"{j}" for j in my_data.getaxis(thisdim)]
        elif my_data.get_ft_prop(thisdim):
            if max_coh_jump is None:
                all_ph = sorted(
                    j for j in my_data.dimlabels if j.startswith("ph")
                )
                if all_ph != [f"ph{j+1}" for j in range(len(all_ph))]:
                    raise ValueError(
                        "If your phases don't label by pulse number, I can't"
                        " automatically determine max_coh_jump"
                    )
                max_coh_jump = {f"ph{j+1}": j + 1 for j in range(len(all_ph))}
            # {{{ Because the phase cycle goes exactly once around the
            #     circle, the spectral width of the coherence domain is
            #     equal to the number of phase cycle steps.  All coherence
            #     changes that differ by a multiple of the spectral width
            #     alias onto the same element, so label each element with
            #     all the allowed changes (i.e. those from -max_coh_jump to
            #     +max_coh_jump) that alias onto it.  An element that no
            #     allowed change aliases onto is labeled "X".
            dp_axis = my_data.getaxis(thisdim)
            dp_sw = len(dp_axis)
            dim_labels[thisdim] = [
                ", ".join(
                    ("%d" if dp == 0 else "%+d") % dp
                    for dp in range(
                        -max_coh_jump[thisdim], max_coh_jump[thisdim] + 1
                    )
                    if np.isclose(
                        (dp - this_dp + dp_sw / 2) % dp_sw, dp_sw / 2
                    )
                )
                or "X"
                for this_dp in dp_axis
            ]
            # }}}
        else:
            dim_labels[thisdim] = [
                "0" if j == 0.0 else f"{j}" for j in my_data.getaxis(thisdim)
            ]
    # }}}
    real_data = not (np.iscomplexobj(my_data.data) and my_data.data.imag.any())
    if cmap is not None:
        if not real_data:
            raise ValueError(
                "In order to use a color map, you must pass real data"
            )
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
    # {{{ Determine the relative sizes of the gaps between the axes.
    #     Move from the innermost outer dimension to the outermost, zooming
    #     out with each step: the existing set of gaps is halved, then
    #     replicated once for each element of the new dimension, with a gap
    #     of size 1 between the replicas.
    divisions = []
    for thisdim in outer_dims[::-1]:
        halved = [j / 2.0 for j in divisions]
        divisions = (halved + [1]) * (my_data.shape[thisdim] - 1) + halved
        logging.debug(strm("for", thisdim, "I get", divisions))
    # the relative sizes are now correct, but they need to add up to "gap"
    divisions = [j * gap / sum(divisions) for j in divisions]
    # }}}
    # {{{ Generate the axes, from bottom to top.  The labels of the outer
    #     dimensions occupy a width of label_width on the left of the bbox.
    label_width = horiz_label_spacer * len(outer_dims) / fig.bbox.width
    axes_left = bbox_left + label_width
    axes_height = (bbox_height - gap) / n_axes
    ax_list = [
        fig.add_axes(
            [axes_left, thisbottom, bbox_width - label_width, axes_height]
        )
        for thisbottom in bbox_bottom + np.cumsum(
            [0] + [axes_height + j for j in divisions]
        )
    ]
    # }}}
    # {{{ Place the data into the axes
    x = np.ravel(my_data.getaxis(direct))
    y = np.ravel(my_data.getaxis(inner_dim))
    # each image is centered on its coordinate, so the image extends half a
    # pixel past the first and last coordinates
    dx = (x[-1] - x[0]) / max(len(x) - 1, 1) / 2
    dy = (y[-1] - y[0]) / max(len(y) - 1, 1) / 2
    origin = kwargs.pop("origin", "lower")
    if origin == "lower":
        extent = (x[0] - dx, x[-1] + dx, y[0] - dy, y[-1] + dy)
    elif origin == "upper":
        extent = (x[0] - dx, x[-1] + dx, y[-1] + dy, y[0] - dy)
    elif origin == "flip":
        extent = (x[-1] + dx, x[0] - dx, y[-1] + dy, y[0] - dy)
    else:
        raise ValueError(
            "I don't understand the value you've set for the origin"
            " keyword argument"
        )
    if interpolation is None:
        # if there are more data points than pixels, we need to interpolate
        # in order to avoid aliasing
        width_px, height_px = ax_list[0].get_window_extent().size
        if len(x) > width_px or len(y) > height_px:
            logging.debug(
                strm(
                    "using bilinear interpolation because data is",
                    len(x),
                    "x",
                    len(y),
                    "but each image block has only",
                    width_px,
                    "x",
                    height_px,
                    "pixels",
                )
            )
            interpolation = "bilinear"
    if interpolation is not None:
        kwargs["interpolation"] = interpolation
    imagehsvkwargs = {
        j: kwargs.pop(j) for j in ["black", "logscale"] if j in kwargs
    }
    # one image for each element of the (flattened) outer dimensions
    all_images = my_data.data.reshape((n_axes, len(y), len(x)))
    data_max = abs(all_images).max()
    for thisax, thisimage in zip(ax_list, all_images):
        if real_data:
            thisimage = thisimage.real / data_max
            kwargs["cmap"] = cmap
        else:
            thisimage = imagehsv(
                thisimage,
                **imagehsvkwargs,
                scaling=scaling_factor if custom_scaling else data_max,
            )
        thisax.imshow(thisimage, extent=extent, origin=origin, **kwargs)
    # }}}
    # {{{ hatch out the coherence pathways that are masked
    if coherence_unmask_fn is not None:
        # {{{ a slice of the data that has only the phase cycling
        #     dimensions provides the shape (and axes) for the mask
        coh_mask = my_data
        for thisdim in my_data.dimlabels:
            if not (thisdim.startswith("ph") and my_data.get_ft_prop(thisdim)):
                coh_mask = coh_mask[thisdim, 0]
        coh_mask = coh_mask.C
        coh_mask.data = np.zeros(coh_mask.data.shape)
        coh_mask.set_prop(
            "coherence_pathway", this_nddata.get_prop("coherence_pathway")
        )
        coh_mask = coherence_unmask_fn(coh_mask)
        # }}}
        # {{{ expand the mask so that we have a boolean for each element
        #     of the inner dimension inside each of the axes
        coh_mask.reorder(
            [j for j in my_data.dimlabels if j in coh_mask.dimlabels]
        )
        is_masked = np.broadcast_to(
            coh_mask.data.reshape(
                [
                    my_data.shape[j] if j in coh_mask.dimlabels else 1
                    for j in my_data.dimlabels[:-1]
                ]
            )
            == 0,
            outer_shape + [len(y)],
        ).reshape((n_axes, len(y)))
        # }}}
        for thisax, masked_rows in zip(ax_list, is_masked):
            # find where each contiguous set of masked rows starts and stops
            run_edges = np.flatnonzero(
                np.diff(np.r_[False, masked_rows, False])
            )
            for run_start, run_stop in zip(run_edges[::2], run_edges[1::2]):
                thisax.axhspan(
                    y[run_start] - dy,
                    y[run_stop - 1] + dy,
                    fill=False,
                    hatch="//",
                    linewidth=0,
                    alpha=0.5,
                )
    # }}}
    # {{{ adjust tick settings
    if shareaxis:
        for thisax in ax_list[1:]:
            thisax.sharex(ax_list[0])
            thisax.sharey(ax_list[0])
    for thisax in ax_list:
        thisax.xaxis.set_major_locator(
            mticker.MaxNLocator(nbins="auto", steps=[1, 2, 2.5, 5, 10])
        )
        thisax.xaxis.set_minor_locator(mticker.AutoMinorLocator(n=5))
        if inner_dim.startswith("ph"):
            # label the elements of the phase cycle
            thisax.set_yticks(y, dim_labels[inner_dim])
        else:
            thisax.yaxis.set_minor_locator(mticker.AutoMinorLocator(n=5))
            thisax.yaxis.set_ticks_position("both")
        for thislabel in thisax.get_yticklabels():
            thislabel.set_va("center")
    ax_list[0].set_xlabel(my_data.unitify_axis(direct), labelpad=12)
    # only the bottom axes has x tick labels, and the top axes has
    # (unlabeled) ticks on top, while the ones in between have no x ticks
    for thisax in ax_list[1:]:
        thisax.xaxis.set_visible(False)
    if n_axes > 1:
        ax_list[-1].xaxis.set_visible(True)
        ax_list[-1].tick_params(
            axis="x",
            which="both",
            top=True,
            bottom=False,
            labeltop=False,
            labelbottom=False,
        )
    # }}}
    # {{{ Find where the tick labels of the inner dimension are (in display
    #     coordinates, relative to the left side of the axes).  If they are
    #     too wide to fit between the axes and the first column of brackets,
    #     make room by moving the left side of the axes to the right.
    #     (Matplotlib also generates labels for ticks that fall outside
    #     the limits of the axes, and we need to ignore those.)
    tick_extents = [
        thislabel.get_window_extent()
        for thistick, thislabel in zip(
            ax_list[0].get_yticks(), ax_list[0].get_yticklabels()
        )
        if min(extent[2:]) <= thistick <= max(extent[2:])
    ]
    ticks_left = (
        min(j.x0 for j in tick_extents) - ax_list[0].get_window_extent().x0
    )
    ticks_right = (
        max(j.x1 for j in tick_extents) - ax_list[0].get_window_extent().x0
    )
    tick_overhang = max(0, 4 - ticks_left - horiz_label_spacer)
    if tick_overhang > 0:
        axes_left += tick_overhang / fig.bbox.width
        for thisax in ax_list:
            thisax.set_position(
                [
                    axes_left,
                    thisax.get_position().y0,
                    bbox_left + bbox_width - axes_left,
                    axes_height,
                ]
            )
    # }}}
    # {{{ transforms for the decorations:
    #     transDispTranslated is in display coordinates, with the origin at
    #     the bottom left corner of the bottom axes
    transDispTranslated = IdentityTransform() + ScaledTranslation(
        axes_left, bbox_bottom, fig.transFigure
    )
    #     transXdispYfig takes x in the same display coordinates, and y in
    #     figure coordinates
    transXdispYfig = blended_transform_factory(
        transDispTranslated, fig.transFigure
    )
    # }}}
    # {{{ Place decorations on figure
    if arrow_width_px is None:
        arrow_width = 4.0  # in points
    else:
        arrow_width = arrow_width_px * 72 / fig.dpi
    # {{{ Each outer dimension gets a column of vertical brackets, with the
    #     innermost of the outer dimensions closest to the axes.  Because
    #     the axes are numbered with the innermost outer dimension varying
    #     fastest, each element of a dimension spans a contiguous block of
    #     axes.
    for j, thisdim in enumerate(outer_dims):
        x_disp = -tick_overhang - (len(outer_dims) - j) * horiz_label_spacer
        block_size = int(np.prod(outer_shape[j + 1 :]))
        for first_ax in range(0, n_axes, block_size):
            # the bracket goes from just above the bottom of the first axes
            # to just below the top of the last axes in the block
            y_bottom = ax_list[first_ax].get_position().y0 + 0.05 * axes_height
            y_top = (
                ax_list[first_ax + block_size - 1].get_position().y0
                + 0.95 * axes_height
            )
            fig.add_artist(
                Line2D(
                    [x_disp, x_disp],
                    [y_bottom, y_top],
                    linewidth=1,
                    color="k",
                    transform=transXdispYfig,
                    clip_on=False,
                )
            )
            fig.text(
                x_disp,
                (y_bottom + y_top) / 2,
                dim_labels[thisdim][
                    first_ax // block_size % my_data.shape[thisdim]
                ],
                va="center",
                ha="right",
                rotation=90,
                transform=transXdispYfig,
                color="k",
            )
        # an arrow points to the bottom of the column of brackets
        _label_arrow(
            fig,
            my_data.unitify_axis(thisdim),
            x_disp,
            arrow_dx,
            arrow_dy,
            arrow_width,
            transDispTranslated,
        )
    # }}}
    # {{{ the arrow for the inner dimension points to its tick labels
    _label_arrow(
        fig,
        my_data.unitify_axis(inner_dim),
        (ticks_left + ticks_right) / 2,
        arrow_dx,
        arrow_dy,
        arrow_width,
        transDispTranslated,
    )
    # }}}
    # }}}
    if title is not None:
        ax_list[-1].set_title(title)
    return (
        ax_list,
        {
            "transDispTranslated": transDispTranslated,
            "transXdispYfig": transXdispYfig,
        },
    )


def fl_DCCT(self, this_nddata, **kwargs):
    this_nddata.squeeze()
    if len(this_nddata.dimlabels) < 3:
        image(this_nddata, **kwargs)
        plt.title(self.current)
    else:
        kwargs["fig"] = self.figdict[self.current]
        kwargs["title"] = self.current
        DCCT(this_nddata, **kwargs)
