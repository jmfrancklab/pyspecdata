"""
Bruker cw ESR Data
==================

Load cw Bruker ESR data, both 1D and 2D.
"""
from numpy import *
import pylab as plt
from pyspecdata import *
# %
# load some 1D ESR data with harmonic + phase info

d = find_file("S175R1a.*DHPC.*200304",
        exp_type='francklab_esr/Sam')
print("here, we see the harmonic axis contains both harmonic and phase info",d.getaxis('harmonic'))
d.chunk_auto('harmonic','phase')
plot(d['phase',0], alpha=0.5)

# %
# Next, let's load some power-saturation data

d = find_file("Power.*Sat.*200303",
        exp_type='francklab_esr/Sam')
d.chunk_auto('harmonic','phase')
figure()
image(d['harmonic',0]['phase',0].C.setaxis('Microwave Power','#').set_units('Microwave Power','scan #'))
plt.show()
