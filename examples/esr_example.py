"""
Bruker cw ESR Data
==================

Load cw Bruker ESR data, both 1D and 2D.

Check out the
`Simple UV <Cary_simple.html>`_
example to understand how
pySpecData locates the file here.
"""
from numpy import *
import pylab as plt
from pyspecdata import *
# %%
# Load some 1D ESR data with harmonic + phase info.
# The data is initial organized into two dimensions -- `harmonic` and `$B_0$`.
# 

d = find_file("S175R1a.*DHPC.*200304",
        exp_type='francklab_esr/Sam')
print(ndshape(d))
print("here, we see the harmonic axis contains both harmonic and phase info",repr(d.getaxis('harmonic')))
d.chunk_auto('harmonic','phase')

# %%
# `chunk_auto` breaks the `harmonic` dimensions since it was labeled with an axis that had 2 fields.

print(ndshape(d))

figure(1)
plot(d['phase',0], alpha=0.5)
plot(d['phase',1], ':', alpha=0.5)
title("1D Data with Multiple Harmonics")

# %%
# Next, let's load some power-saturation data

d = find_file("Power.*Sat.*200303",
        exp_type='francklab_esr/Sam')
d.chunk_auto('harmonic','phase')
figure(2)
image(d['harmonic',0]['phase',0].C.setaxis('Microwave Power','#').set_units('Microwave Power','scan #'))
title("2D Power Saturation")
plt.gca().set_aspect('auto')
plt.show()
