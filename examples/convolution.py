'''Simple Convolution Example
==========================
Since we use convolution a bit for signal analysis, test it here.

Note how the default function preserves the integral, and therefore not the energy.
'''
from pylab import *
from pyspecdata import *
t = nddata(r_[0:4:1024j],'t')
signal = exp(-1j*2*pi*100*t-20*t/pi)
signal.add_noise(0.01)
with figlist_var() as fl:
    fl.next('Time domain')
    fl.plot(signal, label='original')
    fl.next('Fourier transform', legend=True)
    signal.ft('t', shift=True)
    fl.plot(signal, label='original')
    signal.convolve('t',5)
    fl.plot(signal, label='after convolve')
    fl.next('Time domain')
    signal.ift('t')
    fl.plot(signal, label='after convolve')
