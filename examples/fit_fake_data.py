"""Fitting Data (Nonlinear + Symbolic)
======================================

This example creates fake data with noise
then fits the exponential with the fitdata
function."""  
from pyspecdata import *
import sympy as sp
# {{{ this is the contents of pylab.py -- works
# need to go through and figure out which lines
# are actually needed and which are not
# -- I have already stripped out some
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
from lmfitdata import lmfitdata
from matplotlib.pyplot import figure, subplot, show, xlim, ylim, plot, gca
from numpy import * # I think it wasn't importing from numpy b/c it seems we're inside sphinx
# }}}
fl=figlist_var()
#{{{creating a fake data recovery curve
tau = nddata(r_[0:2:100j], 'tau')
fake_data = 102*(1-2*exp(-tau*6.0))
fake_data.add_noise(5.0)
empty_data = nddata(r_[0:2:100j],'tau')
#}}}
#{{{ fitting data
f = fitdata(fake_data)
fM0, fMi, fR1, fvd = sp.symbols("M_0 M_inf R_1 tau",real=True)
f.functional_form = fMi + (fM0-fMi)*sp.exp(-fvd*fR1)
f.set_guess({fM0:-500, fMi:500, fR1:2})
f.settoguess()
fitdataguess = f.eval(100)
f.fit()
# {{{ this is just to show all the parameters
list_symbs = []
for j,k in f.output().items():
    s_repr = sp.latex(sp.Symbol(j))
    list_symbs.append(f'${s_repr} = {k:0.5g}$')
list_symbs = '\n'.join(list_symbs)
# }}}
T1 = 1./f.output('R_1')
# }}}
#}}}
#{{{lmfitdata method
M0,Mi,R1,vd = sp.symbols("M0 Mi R1 tau")
thisfit = lmfitdata(empty_data)
newfit = lmfitdata(fake_data)
newfit.functional_form = (Mi + (M0-Mi)*sp.exp(-(vd*R1)))
thisfit.functional_form = newfit.functional_form
newfit.set_guess(
        M0 = dict(value=-500, max=0, min=-501),
        Mi=dict(value=500, max = 501, min=0), 
        R1=dict(value=5, max = 6, min = 1))
newfit.settoguess()
newguess = newfit.eval(100)
newfit.fit()
out = minimize(
        thisfit.residual, newfit.pars, kws={"data":fake_data.data})
newerfit = empty_data.copy(data=False)
newerfit.data = newfit.residual(out.params)
#}}}
with figlist_var() as fl: 
    fl.next('fit with guess')
    fl.plot(fitdataguess,label='fitdata guess')
    fl.plot(newguess,label='lmfitdata guess')
    fl.plot(fake_data,'o',label='fake data')
    thisline = fl.plot(f.eval(100),label='fit data fit')
    thatline = fl.plot(newerfit,':',label='lmfitdata fit')

