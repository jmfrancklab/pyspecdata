from pyspecdata import *
tau = nddata(r_[0:2:100j], 'tau')
fake_data = 102*(1-2*exp(-tau*6.0))
fake_data.add_noise(5.0)
with figlist_var() as fl:
    fl.next('fake data')
    fl.plot(fake_data)
