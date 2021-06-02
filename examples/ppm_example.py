from pylab import *
from pyspecdata import *

# make peaks on a frequency axis
x = nddata(r_[0:8:49152j], 't2')
x.set_units('t2','s')

# generate time-domain signal
signal = 0*x # create an array of zeros that's the right shape
SFO1 = 400.130438262389
for nu,R in [(25,0.08),(50,1.2),(20,0.45),(750,1.2),(10,0.08)]:
    nu = SFO1 + nu
    signal += exp(1j*2*pi*nu*x-x/R)
signal.set_units('t2','s')
signal.ft('t2',shift=True)

# Parameters needed by to_ppm to work
signal.set_prop('acq',{'SFO1':400.130438262389,'SW_h':4807.69230769231,'O1':438.262389})
signal.set_prop('proc',{'OFFSET':12.50096})

# Copy of signal to demonstrate truncation
signal_truncated = signal.C

with figlist_var() as fl:
    fl.next('full spectrum, Hz')
    fl.plot(signal)
    signal.to_ppm()
    fl.next('full spectrum, ppm')
    fl.plot(signal)
    signal_truncated = signal_truncated['t2':(0.2e3,1.4e3)]
    fl.next('truncated spectrum, Hz')
    fl.plot(signal_truncated)
    signal_truncated.to_ppm()
    fl.next('truncated spectrum, ppm')
    fl.plot(signal_truncated)
    fl.show();quit()
