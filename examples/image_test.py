from pylab import *
from pyspecdata import *

fl = figlist_var()

# Generate fake data
t_axis = nddata(r_[0:2:2048j],'t2')
s = exp(1j*2*pi*5*t_axis - t_axis/800e-3 )
s += exp(1j*2*pi*-30*t_axis - t_axis/800e-3)
ph1 = nddata(r_[0:4]/4.,'ph1')
ph2 = nddata(r_[0,2]/4.,'ph2')
s *= exp(1j*2*pi*ph1)
s *= exp(1j*2*pi*ph2)
s['t2',0] *= 0.5
s.ft('t2',shift=True)
s.reorder('t2',first=False)
s.ft(['ph1','ph2'])
fl.next('Time domain')
fl.image(s)
fl.next('F domain')
fl.image(s)
fl.show();quit()

#s = nddata(['t2','ph1','ph2'],[2048,4,2])
print(ndshape(s))
quit()
