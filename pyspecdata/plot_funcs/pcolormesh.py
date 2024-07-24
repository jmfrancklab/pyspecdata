import matplotlib.pylab as plt
<<<<<<< HEAD
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
):
||||||| b2b39436
import logging
def pcolormesh(self, fig=None, shading='nearest',ax1=None,ax2=None,ax=None, scale_independently=False, human_units=True):
=======
import logging


def pcolormesh(
    self,
    fig=None,
    shading="nearest",
    ax1=None,
    ax2=None,
    ax=None,
    scale_independently=False,
    human_units=True,
):
>>>>>>> master
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

    Returns
    =======
    nothing for now
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
    X, Y = plt.meshgrid(
        forplot.getaxis(forplot.dimlabels[1]),
        forplot.getaxis(forplot.dimlabels[0]),
    )
    Z = forplot.data
    # {{{ store these so that I can set the color scale for the plots together,
    #     at the end
    vmin_list = []
    vmax_list = []
    mappable_list = []
    # }}}
    for thisax, thisfun, thislabel in ax_list:
        Zdata = thisfun(Z)
        mappable = thisax.pcolormesh(X, Y, Zdata, shading=shading)
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
    return
