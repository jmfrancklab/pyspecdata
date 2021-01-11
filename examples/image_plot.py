"""Plotting Routines
=================

A simple demonstration of a multi-line *vs.*
image plot with
`<https://en.wikipedia.org/wiki/Domain_coloring> domain coloring`_, nested
inside a figure list.

We specifically test a safety feature that doesn't allow image plots
to have unevenly spaced axes,
and show how to deal with this.
"""
from pyspecdata import *
from numpy import *
# let's make some fake inversion recovery data
vd = nddata(r_[0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0,2.0,3.0,5.0],'vd')
signal_amp = 1-2*exp(-vd/1.0)
t2 = nddata(r_[0:2:256j],'t2')
signal_amp *= exp(-t2/0.1)
signal_amp.add_noise(0.1)
signal_amp.set_units('s')
signal_amp.ft('t2', shift=True)
with figlist_var() as fl:
    fl.next('1D data')
    fl.plot(signal_amp)
    generated_error = False
    # If I were to just run the following command (not nested in try/except)
    fl.next('this figure intentionally blank!')
    try:
        fl.image(signal_amp)
    except:
        generated_error = True
    # (try it).
    # I would get an error that tells me to do this
    fl.next('image plot -- good axis')
    fl.image(signal_amp.C.setaxis('vd','#').set_units('vd','scan #'))
    assert generated_error
