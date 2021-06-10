#from pylab import *
from pyspecdata import * 
import numpy as np
from pyspecdata.plot_funcs.DCCT import DCCT
#from pyspecdata.core import *
fl=figlist_var()
t_axis = nddata(r_[0:2:256j],'t2')
vd = nddata(r_[0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0,2.0,3.0,5.0],'vd')
s = 1-2*np.exp(-vd/1.0)
s *= np.exp(-t_axis/0.1)
#s = np.exp(1j*2*pi*5*t_axis - t_axis/800e-3 )
#s += np.exp(1j*2*pi*-30*t_axis - t_axis/800e-3)
ph1 = nddata(r_[0:4]/4.,'ph1')
ph2 = nddata(r_[0,2]/4.,'ph2')
s.add_noise(0.3)
s *= np.exp(1j*2*pi*ph1)
s *= np.exp(1j*2*pi*ph2)
s.set_units('t2','s')
s.reorder(['ph1','ph2','vd','t2'])
print(ndshape(s))
quit()
DCCT(this_nddata=s,this_fig_obj=figure())
fl.show()
s.ft(['ph1','ph2'])
quit()
