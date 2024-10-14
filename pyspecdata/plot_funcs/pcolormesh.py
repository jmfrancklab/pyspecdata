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
    vmin = None,
    vmax = None,
    human_units=True,
    force_balanced_cmap=False,
    handle_axis_sharing=True,
    mappable_list=[],
):
    """generate a pcolormesh and label it with the axis coordinate available from the nddata

    Parameters
    ==========
    fig : matplotlib figure object
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
        -- *e.g.* especially when you are plotting a real and imaginary together.
        So, this defaults to true to do that.
        But sometimes, you want to get fancy and, *e.g.* bind the sharing of many plots together
        because matplotlib doesn't let you call sharex/sharey more than once,
        you need then to tell it not to handle the axis sharing, and to it yourself
        outside this routine.
    mappable_list : list, default []
        empty list which fills with field values from color axis used for
        initial subplot, used to scale multiple plots along the same color
        axis. Used to make all 3x2 plots under a uniform color scale

    Returns
    =======
    mappable_list : list
        list of field values for scaling color axis, used to make all 3x2 plots
        under a uniform color scale
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
        # TODO ☐: I don't understand what these are for
        vmin = overall_min
        vmax = overall_max
    if scale_independently:
        vmin, vmax = mappable.get_clim()
        print("independently scaled vmin is ", vmin, "and vmax is ", vmax)
    # delete when read: I modified following to include axis
    # I changed the name of modify_colorbar_boundaries
    def adjust_colorbar(mappable, vmin, vmax, ax):
        mappable.norm.vmin = vmin
        mappable.norm.vmax = vmax
        mappable.set_norm(mappable.norm)
        # (delete when read) note that I get the following from chatGPT
        # -- we need to get the colorbar corresponding to the mappable of
        # interest
        cbar = mappable.colorbar
        if cbar is None:
            # delete when read -- I tested as follows:
            # 
            # In [53]: m = pcolormesh(r_[1,2,3],r_[1,2,3],r_[0:9].reshape(3,3))
            # 
            # In [54]: m.colorbar
            # 
            # In [55]: print(m.colorbar)
            # None
            # 
            # In [56]: colorbar()
            # Out[56]: <matplotlib.colorbar.Colorbar at 0x7f96751108d0>
            # 
            # In [57]: print(m.colorbar)
            # <matplotlib.colorbar.Colorbar object at 0x7f96751108d0>
            cbar = plt.colorbar(mappable=mappable, ax=ax)
        cbar.update_normal(mappable)
    for j, (thisax, thisfun, thislabel) in enumerate(ax_list):
        if not scale_independently:
            # if we are not scaling our plots independently,
            # then we want to set the color limits for each mappable to
            # the *overall* max z value and the *overall* min z value
            # (where z value means the value of the plot at x and y
            # (regardless of what our axis names are) that determines the
            # color
            mappable_list[j].set_clim(overall_min, overall_max)
            # TODO ☐: spinning up new color bars is exactly how we git
            #         into this trouble  -- I'm not modifying, but the
            #         following is likely going in the wrong direction
            #         so, you want to delete the following
            cbar = plt.colorbar(mappable=mappable_list[j], ax=thisax)
        elif scale_independently or j > 0:
            # if we are scaling independently, or if we are on the right
            # plot, we want to show a colorbar
            # TODO ☐: not changing, but I don't see how renaming
            #         overall_min and overall_max to vmin and vmax helps
            #         anything
            mappable_list[j].set_clim(vmin, vmax)
            # TODO ☐: I would guess that the solution is in the
            #         following.  If this is our first time touching this
            #         element of the mappable list, then you want to call
            #         plt.colorbar.  But, if the colorbar already exists,
            #         we want to adjust it.  I achieve this by modifying
            #         your existing function to grab the colorbar.  I
            #         think this should do what we want to please read +
            #         understand what I'm doing here.
            adjust_colorbar(mappable=mappable_list[j], vmin=vmin, vmax=vmax, ax=thisax)
        else: 
            # TODO ☐: this section is when you are not scaling
            #         independently, and you are in the left plot just
            #         for the sake of space/style, we don't want to
            #         generate a plot there, so this section should just not exist
            cbar = plt.colorbar(mappable)
            print("this section is being called")
    return mappable_list
