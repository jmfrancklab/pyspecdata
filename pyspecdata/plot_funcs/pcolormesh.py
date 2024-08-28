import matplotlib.pylab as plt
import numpy as np


def pcolormesh(
    self,
    fig=None,
    shading="nearest",
    ax1=None,
    ax2=None,
    ax=None,
    scale_independently=False,
    human_units=True,
    force_balanced_cmap=False,
    handle_axis_sharing=True,
    mappable_list=[],
):
    """generate a pcolormesh and label it with the axis coordinate available from the nddata

    Parameters
    ==========
    fig: matplotlib figure object
    shading: str (default 'nearest')
        the type of shading to pass to matplotlib pcolormesh
    ax1: matplotlib axes object
        where do you want the left plot to go?
    ax2: matplotlib axes object
        where do you want the right plot to go?
    scale_independently: boolean (default False)
        Do you want each plot to be scaled independently?
        (If false, the colorbar will have the same limits for all plots)
    handle_axis_sharing: boolean (default True)
        Typically, you want the axes to scale together when you zoom
        -- *e.g.* especially when you are plotting a real and imaginary together.
        So, this defaults to true to do that.
        But sometimes, you want to get fancy and, *e.g.* bind the sharing of many plots together
        because matplotlib doesn't let you call sharex/sharey more than once,
        you need then to tell it not to handle the axis sharing, and to it yourself
        outside this routine.
    mappable_list: list of

    Returns
    =======
    nothing for now
    mappable_list
    """
    assert len(self.dimlabels) == 2, "currently, this only supports 2D data"
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
    # {{{ store these so that I can set the color scale for the plots together,
    #     at the end
    vmin_list = []
    vmax_list = []
    # }}}
    for thisax, thisfun, thislabel in ax_list:
        Zdata = thisfun(Z)
        mappable = thisax.pcolormesh(X, Y, Zdata, shading=shading)
        if handle_axis_sharing and thisax != ax_list[0][0]:
            thisax.sharex(ax_list[0][0])
            thisax.sharey(ax_list[0][0])
        mappable_list.append(mappable)
        thisax.set_ylabel(forplot.unitify_axis(forplot.dimlabels[0]))
        thisax.set_xlabel(forplot.unitify_axis(forplot.dimlabels[1]))
        thisax.set_title(f"({thislabel})")
    for this_mappable in mappable_list:
        thismin, thismax = this_mappable.get_clim()
        vmin_list.append(thismin)
        vmax_list.append(thismax)
    if not scale_independently:
        overall_min = np.min(vmin_list)
        overall_max = np.max(vmax_list)
        if force_balanced_cmap:
            if overall_max > -overall_min:
                overall_min = -overall_max
            else:
                overall_max = -overall_min
    for j, (thisax, thisfun, thislabel) in enumerate(ax_list):
        if not scale_independently:
            mappable_list[j].set_clim(overall_min, overall_max)
        if scale_independently or j > 0:
            plt.colorbar(mappable=mappable_list[j], ax=thisax)
    return mappable_list
