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

x = psp.nddata(r_[-5, -2, -1, -0.5, 0, 0.5, 5], "x")
y = psp.nddata(3 * r_[-5, -2, -1, -0.5, 0, 0.5, 5], "y")
z = plt.exp(-((y - 2) ** 2) - (x - 0) ** 2 / 2) + 1j * x
z1 = plt.exp(-((y - 2) ** 2) - (x - 0) ** 2 / 2) + 1j * x
z2 = 10 * plt.exp(-((y - 2) ** 2) - (x - 1) ** 2 / 2) + 10j * x
fig = plt.figure()
plt.suptitle("colorscales dependent")
z1.dimlabels
gs = plt.GridSpec(2, 2, hspace=0.5)
ax_list = []
for j in range(2):
    for k in range(2):
        ax_list.append(fig.add_subplot(gs[j, k]))
mpb = z2.pcolor(ax1=ax_list[0], ax2=ax_list[1])
z1.pcolor(mappable_list=mpb, ax1=ax_list[2], ax2=ax_list[3])
plt.figure()
plt.suptitle("colorscales independent")
z1.pcolor(scale_independently=True)
plt.show()
