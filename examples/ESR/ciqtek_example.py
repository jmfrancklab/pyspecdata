"""
CIQTEK cw ESR Data
==================

Load CIQTEK JSON EPR data through :func:`pyspecdata.find_file`.
"""

import re

import matplotlib.pyplot as plt
import pyspecdata as psd

Bname = r"$B_0$"
filename = "384K_1mM_TEMPO_SFO5_pt5G_MAmp_25dB_100G_22scans_2D.epr"

d = psd.find_file(re.escape(filename), exp_type="ciqtek")
g_axis = d.get_prop("ciqtek_g_axis")


# TODO ☐: don't have real and imag as two figures -- put on a single subplot, and bind the scale of the subplots together (autoscale to the larger data -- there should be a matplotlib way of doing this automatically)
plt.figure(1)
# TODO ☐: DO NOT use plt.plot!!! use psd.plot!!!
psd.plot(d.real, alpha=0.7)
plt.title("CIQTEK real signal")

plt.figure(2)
psd.plot(d.imag, alpha=0.7)
plt.title("CIQTEK imaginary signal")

plt.figure(3)
# TODO ☐: do not plot this as a separate plot -- rather use this to set up a twiny for the plots above that labels g factor along the top
plt.plot(d.getaxis(Bname), g_axis)
plt.xlabel(Bname)
plt.ylabel("g")
plt.gca().invert_yaxis()
plt.title("CIQTEK g-axis")

# TODO ☐: since this is 2D data, use psd.image to show the 2D data! put the field along the y axis (use reorder!)

plt.show()
