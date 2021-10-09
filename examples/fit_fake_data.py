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
#}}}
#{{{ fitting data
f = fitdata(fake_data)
M0,Mi,R1,vd = sp.symbols("M_0 M_inf R_1 tau",real=True)
f.functional_form = Mi + (M0-Mi)*sp.exp(-vd*R1)
logger.info(strm("Functional Form", f.functional_form))
logger.info(strm("Functional Form", f.functional_form))
f.set_guess({M0:-500, Mi:500, R1:2})
f.settoguess()
guess = f.eval(100)
f.fit()
print("output:",f.output())
print("latex:",f.latex())
#}}}
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
newfit = lmfitdata(fake_data)
newfit.functional_form = (Mi + (M0-Mi)*sp.exp(-(vd*R1)))
newfit.set_guess(
        M_0 = dict(value=-500, max=0, min=-501),
        M_inf=dict(value=500, max = 501, min=0), 
        R_1=dict(value=5, max = 6, min = 1))
# the following is failing for me
newfit.settoguess()
newguess = newfit.eval(100)
newfit.fit()
# I don't understand -- minimize should be called as part of the fit
# method, above -- why are you calling it explicitly here?
out = minimize(
        thisfit.residual, newfit.pars, kws={"data":fake_data.data})
newerfit.data = newfit.residual(out.params)
#}}}
with figlist_var() as fl: 
    fl.next('fit with guess')
    fl.plot(fitdataguess, label='fitdata guess')
    fl.plot(newguess, label='lmfitdata guess')
    fl.plot(fake_data,'o',label='fake data')
    thisline = fl.plot(f.eval(100),label='fit data fit')
    thatline = fl.plot(newerfit,':',label='lmfitdata fit')
    # {{{ just put the text
    ax = gca()
    text(0.5,0.5,f.latex(),
            ha='center',va='center',
            color=thisline[0].get_color(),
            transform = ax.transAxes)
    text(0.5,0.5,(3*'\n')+list_symbs,
            ha='center',va='top',
            size=10,
            color=thisline[0].get_color(),
            transform = ax.transAxes)
    # }}}

