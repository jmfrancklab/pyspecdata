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
plt.figure()
plt.suptitle("colorscales dependent")
z.dimlabels
z.pcolor()
plt.figure()
plt.suptitle("colorscales independent")
z.pcolor(scale_independently=True)
plt.show()
