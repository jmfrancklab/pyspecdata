from pyspecdata import *
from pyspecdata import fitdata 
from sympy import symbols
import numpy as np
#{{{creating a fake data recovery curve
tau = nddata(r_[0:2:100j], 'tau')
fake_data = 102*(1-2*exp(-tau*6.0))
fake_data.add_noise(5.0)
#}}}
#{{{ fitting data
f = fitdata(fake_data)
M0,Mi,R1,vd = sympy.symbols("M_0 M_inf R_1 tau",real=True)
f.functional_form = Mi + (M0-Mi)*sympy.exp(-vd*R1)
logger.info(strm("Functional Form", f.functional_form))
logger.info(strm("Functional Form", f.functional_form))
f.set_guess({M0:-500, Mi:500, R1:2})
f.settoguess()
save_guess = f.eval(100)
f.fit()
print("output:",f.output())
print("latex:",f.latex())
with figlist_var() as fl: 
    fl.next('fit with guess')
    fl.plot(f,'o',label='data')
    fl.plot(save_guess,label='guess')

