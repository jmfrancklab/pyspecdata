import matplotlib.pylab as plt
import numpy as np


def pcolormesh(
    self,
    fig=None,
    cmap="viridis",
    shading="nearest",
    ax1=None,
    ax2=None,
    ax=None,
    scale_independently=False,
    vmin=None,
    vmax=None,
    human_units=True,
    force_balanced_cmap=False,
    handle_axis_sharing=True,
    mappable_list=None,
):
    """generate a pcolormesh and label it with the axis coordinate available
    from the nddata

    Parameters
    ==========
    fig : matplotlib figure object
    cmap: str (default 'viridis')
        cmap to pass to matplotlib pcolormesh
    shading: str (default 'nearest')
        the type of shading to pass to matplotlib pcolormesh
    ax1 : matplotlib axes object
        where do you want the left plot to go?
    ax2 : matplotlib axes object
        where do you want the right plot to go?
    ax : matplotlib axes object
        if passed, this is just used for ax1
    scale_independently: boolean (default False)
        Do you want each plot to be scaled independently?
        (If false, the colorbar will have the same limits for all plots)
    handle_axis_sharing: boolean (default True)
        Typically, you want the axes to scale together when you zoom
        -- *e.g.* especially when you are plotting a real and imaginary
        together.
        So, this defaults to true to do that.
        But sometimes, you want to get fancy and, *e.g.* bind the sharing of
        many plots together
        because matplotlib doesn't let you call sharex/sharey more than once,
        you need then to tell it not to handle the axis sharing, and to it
        yourself
        outside this routine.
    mappable_list : list, default []
        used to scale multiple plots along the same color
        axis. Used to make all 3x2 plots under a uniform color scale.

        List of QuadMesh objects returned by this function.

    Returns
    =======
    mappable_list : list
        list of field values for scaling color axis, used to make all 3x2 plots
        under a uniform color scale
    """
    assert len(self.dimlabels) == 2, "currently, this only supports 2D data"
    if mappable_list is None:
        # settings the default value to just [] above does NOT create a
        # new object.  This does
        mappable_list = list([])
    if human_units:
        forplot = self.C.human_units()
    else:
        forplot = self.C
    if ax is not None:
        ax1 = ax
    if forplot.data.dtype == plt.complex128:
        if ax1 is None:
            fig = plt.gcf()
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)
        ax_list = [
            (ax1, lambda x: x.real, "real"),
            (ax2, lambda x: x.imag, "imag"),
        ]
    else:
        if ax1 is None:
            fig = plt.gcf()
            ax1 = fig.add_subplot(1, 1, 1)
        ax_list = [(ax1, lambda x: x.real, "real")]
    X, Y = np.meshgrid(
        forplot.getaxis(forplot.dimlabels[1]),
        forplot.getaxis(forplot.dimlabels[0]),
    )
    Z = forplot.data
    for j, (thisax, thisfun, thislabel) in enumerate(ax_list):
        Zdata = thisfun(Z)
        mappable = thisax.pcolormesh(X, Y, Zdata, cmap=cmap, shading=shading)
        if handle_axis_sharing and thisax != ax_list[0][0]:
            thisax.sharex(ax_list[0][0])
            thisax.sharey(ax_list[0][0])
        mappable_list.append(mappable)
        thisax.set_ylabel(forplot.unitify_axis(forplot.dimlabels[0]))
        thisax.set_xlabel(forplot.unitify_axis(forplot.dimlabels[1]))
        thisax.set_title(f"({thislabel})")
        # {{{ handle the creation of colorbars here, since adjusting the
        #     mappable clim later changes the colorbars just fine
        if j == 0:
            if scale_independently:
                plt.colorbar(mappable=mappable, ax=thisax)
            else:
                pass  # b/c no use for extra colorbar if locked together
        elif j == 1:
            plt.colorbar(mappable=mappable, ax=thisax)
        # }}}
    # {{{ overall scaling
    # is manually specified:
    if vmin is not None:
        assert vmax is not None, "if vmin is specified, vmax must be too!"
        assert not scale_independently, (
            "scale_independently is True but you've manually set vmin and"
            " vmax, this doesn't make sense! If they share vmax and vmin, then"
            " they are scaled together!!"
        )
        assert not force_balanced_cmap, (
            "you're trying to force the colormap to have a balanced scale"
            " while also manually setting its limits, this doesn't make sense!"
        )
        overall_min = vmin
        overall_max = vmax
    if not scale_independently:
        if vmin is None:
            # {{{ we only need to determine the overall min and max if we
            #    haven't explicitly set them
            vmin_list = []
            vmax_list = []
            for this_mappable in mappable_list:
                thismin, thismax = this_mappable.get_clim()
                vmin_list.append(thismin)
                vmax_list.append(thismax)
            overall_min = np.min(vmin_list)
            overall_max = np.max(vmax_list)
            if force_balanced_cmap:
                if overall_max > -overall_min:
                    overall_min = -overall_max
                else:
                    overall_max = -overall_min
            # }}}
        for thismappable in mappable_list:
            thismappable.set_clim(vmin=overall_min, vmax=overall_max)
    # }}}
    return mappable_list
