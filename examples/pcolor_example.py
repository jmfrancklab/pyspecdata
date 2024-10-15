"""
Unevenly spaced data
====================

Basic demonstration of pcolor, which can deal with unevenly spaced data

.. note::
    Right now, we just do this with real/imaginary,
    but in principal, it should be easily
    possible to extend this to use domain
    coloring (and to use it in the main DCCT
    method)
"""

import pyspecdata as psp
import matplotlib.pylab as plt
from numpy import r_

run_to_checkpoint = 1 # allows us to run to different checkpoints.  If
#                       everything is working correctly, this should go up to 4

x = psp.nddata(r_[-5, -2, -1, -0.5, 0, 0.5, 5], "x")
y = psp.nddata(3 * r_[-5, -2, -1, -0.5, 0, 0.5, 5], "y")
z1 = plt.exp(-((y - 2) ** 2) - (x - 0) ** 2 / 2) + 1j * x
z2 = 10 * z1
# {{{ plot the smaller data
# TODO ☐: note that this is currently NOT scaling independently.  two
#         colorbars are generated that are both the same.  However, when I run this
#         script on the master branch, this step runs fine!
plt.figure()
z1.C.pcolor(scale_independently=True)
# }}}
if run_to_checkpoint > 1:
    # {{{ plot the larger data
    # TODO ☐: for reasons that are not clear, the previous call messes
    #         up this call.  On the current version of the master
    #         branch, the colors look OK, but the color bars are
    #         incorrectly labeled (b/c z2 is larger than z1)
    plt.figure()
    z2.pcolor(scale_independently=True)
    # }}}
if run_to_checkpoint > 2:
    def new_figure_and_grid():
        fig = plt.figure()
        gs = plt.GridSpec(2, 2, hspace=0.5)
        ax_list = []
        for j in range(2):
            for k in range(2):
                ax_list.append(fig.add_subplot(gs[j, k]))
        return ax_list
    # {{{ independent
    ax_list = new_figure_and_grid()
    plt.suptitle("colorscales independent")
    z1.pcolor(scale_independently=True, ax1=ax_list[0], ax2=ax_list[1])
    z2.pcolor(scale_independently=True, ax1=ax_list[2], ax2=ax_list[3])
    # }}}
if run_to_checkpoint > 3:
    # {{{ small first, then large
    ax_list = new_figure_and_grid()
    plt.suptitle("colorscales dependent")
    mpb = z1.C.pcolor(ax1=ax_list[0], ax2=ax_list[1])
    z2.C.pcolor(mappable_list=mpb, ax1=ax_list[2], ax2=ax_list[3])
    # }}}
if run_to_checkpoint > 4:
    # {{{ large in first row, then small in second row
    ax_list = new_figure_and_grid()
    plt.suptitle("colorscales dependent")
    mpb = z2.C.pcolor(ax1=ax_list[0], ax2=ax_list[1])
    z1.C.pcolor(mappable_list=mpb, ax1=ax_list[2], ax2=ax_list[3])
    # }}}
plt.show()
