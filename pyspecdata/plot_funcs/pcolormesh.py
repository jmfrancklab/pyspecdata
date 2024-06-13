from ..general_functions import *
import matplotlib.pylab as plt
import logging
def pcolormesh(self, fig=None, shading='nearest',ax1=None,ax2=None,ax=None, scale_independently=False):
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
    is_complex = False
    if ax is not None:
        ax1 = ax
    if self.data.dtype == plt.complex128:
        if ax1 is None:
            fig = plt.gcf()
            ax1 = fig.add_subplot(1,2,1)
            ax2 = fig.add_subplot(1,2,2)
        ax_list = [(ax1,lambda x: x.real,'real'),(ax2,lambda x: x.imag,'imag')]
    else:
        if ax1 is None:
            fig = plt.gcf()
            ax1 = fig.add_subplot(1,1,1)
        ax_list = [(ax1,lambda x: x.real,'real')]
    X,Y = np.meshgrid(
            self.getaxis(self.dimlabels[1]),
            self.getaxis(self.dimlabels[0]))
    Z = self.data
    # {{{ store these so that I can set the color scale for the plots together,
    #     at the end
    vmin_list = []
    vmax_list = []
    mappable_list = []
    # }}}
    for thisax, thisfun, thislabel in ax_list:
        Zdata = thisfun(Z)
        vmin_list.append(Zdata.min())
        vmax_list.append(Zdata.max())
        mappable = thisax.pcolormesh(X,Y,Zdata,shading=shading)
        mappable_list.append(mappable)
        thisax.set_ylabel(self.unitify_axis(self.dimlabels[0]))
        thisax.set_xlabel(self.unitify_axis(self.dimlabels[1]))
        thisax.set_title(f"({thislabel})")
    for j,(thisax, thisfun, thislabel) in enumerate(ax_list):
        if not scale_independently:
            mappable_list[j].set_clim(np.min(vmin_list),np.max(vmax_list))
        cbar = plt.colorbar(mappable=mappable_list[j], ax=thisax)
    return
