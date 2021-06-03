from pylab import *
from pyspecdata import * 
from pyspecdata.plot_funcs.DCCT import DCCT
fl=figlist_var()
t_axis = nddata(r_[0:2:2048j],'t2')
s = exp(1j*2*pi*5*t_axis - t_axis/800e-3 )
s += exp(1j*2*pi*-30*t_axis - t_axis/800e-3)
ph1 = nddata(r_[0:4]/4.,'ph1')
ph2 = nddata(r_[0,2]/4.,'ph2')
s.add_noise(0.3)
s *= exp(1j*2*pi*ph1)
s *= exp(1j*2*pi*ph2)
s.set_units('t2','s')
s.reorder(['ph1','ph2','t2'])
print(ndshape(s))

DCCT(this_nddata=s,this_fig_obj=figure())
fl.show()
s.ft(['ph1','ph2'])
quit()
