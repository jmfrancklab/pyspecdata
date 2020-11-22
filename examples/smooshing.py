"""Break a single dimension into may (chunk) and collapse multiple into one (smoosh)

In order to demonstrate this, we create some CPMG-like data (comprising many echos that are subject to a decay)
"""
from pylab import *
from pyspecdata import *
t_echo = nddata(r_[0:1:16j],'t_echo')
t2 = nddata(r_[-1/20:1/20:256j],'t2')
T2 = 0.4
T2star = 1e-2
fake_data = exp(-t_echo/T2)*exp(-abs(t2)/T2star) + 0j
fake_data.add_noise(0.2)
with figlist_var() as fl:
    fl.next("fake data")
    fl.image(fake_data)
    fake_data.smoosh(['t_echo','t2'],'t2')
    fl.next("visualize CPMG as 1D data")
    # smoosh stores the data in a structured array format in such a way that it
    # remembers how to automatically chunk the data again, but this also means
    # we need to relabel the axis if we want to plot it nicely
    fl.plot(fake_data.C.setaxis('t2','#')) # just number the time points
    fake_data.chunk_auto('t2')
    fl.next("restored fake data")
    fl.image(fake_data)
