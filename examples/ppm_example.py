from pylab import *
from pyspecdata import *

# make peaks on a frequency axis
x = nddata(r_[0:6.8157439:32768j], "t2")
x.set_units("t2", "s")

# generate time-domain signal
SW_H = 1./(x.getaxis('t2')[1] - x.getaxis('t2')[0])
signal = 0 * x  # create an array of zeros that's the right shape
SFO1 = 400.130438262389
TMS_shift = 12.5
for A, nu, R in [
    (0.1, TMS_shift, 0.04),
    (1, 25, 0.08),
    (1, 50, 1.2),
    (1, 20, 0.45),
    (1, 750, 1.2),
    (1, 10, 0.08),
]:
    #nu = SFO1 + nu # it's really unclear why this is done!
    signal += A * exp(1j * 2 * pi * nu * x - x / R)
signal.set_units("t2", "s")
signal.ft("t2", shift=True)

# Parameters needed by to_ppm to work
# pulled from BV_polymer_apr27_2021
signal.set_prop('acq',{'SFO1':SFO1,'SW_h':SW_H,'O1':438.262389})
signal.set_prop('proc',{'OFFSET':12.50096})

# Copy of signal to demonstrate truncation
signal_truncated = signal.C

with figlist_var() as fl:
    fl.next("full spectrum, Hz")
    fl.plot(signal)
    signal.to_ppm()
    fl.next("full spectrum, ppm")
    fl.plot(signal)
    signal_truncated = signal_truncated["t2":(0.2e3, 1.4e3)]
    fl.next("truncated spectrum, Hz")
    fl.plot(signal_truncated)
    signal_truncated.to_ppm()
    fl.next("truncated spectrum, ppm")
    fl.plot(signal_truncated)
