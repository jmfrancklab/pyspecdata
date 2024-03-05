from ..general_functions import *
import matplotlib.pylab as plt
import logging
def pcolormesh(self, fig=None, shading='nearest',ax1=None,ax2=None,ax=None):
    "generate a pcolormesh"
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
    for thisax, thisfun, thislabel in ax_list:
        mappable = thisax.pcolormesh(X,Y,thisfun(Z),shading=shading)
        thisax.set_ylabel(self.unitify_axis(self.dimlabels[0]))
        thisax.set_xlabel(self.unitify_axis(self.dimlabels[1]))
        thisax.set_title(f"({thislabel})")
        plt.colorbar(mappable=mappable, ax=thisax)
    return
