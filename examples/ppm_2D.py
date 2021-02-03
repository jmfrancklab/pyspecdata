from pylab import *
from pyspecdata import *

# make peaks on a frequency axis
x = nddata(r_[0:5:1024j], 't2')
y = nddata(r_[0:5:1024j], 't1')

# generate time-domain signal
signal = 0*x*y # create an array of zeros that's the right shape
for nu1,nu2 = [....
    signal += exp(+1j*2*pi*nu2*x-x/2.+1j*2*pi*nu1*y-y/2)
signal.ft('t2',shift=True)
signal.ft('t1',shift=True)

signal.set_prop('acq',{'SFO1':400.})
signal.set_prop('procs',{'OFFSET':-10.})

with figlist_var() as fl:
    fl.next('the spectrum')
    fl.image(signal)
    signal.to_ppm('t2')
    signal.to_ppm('t1')
    fl.next('after converting to ppm')
    fl.image(signal)
