"""Fitting Data
==============

This example creates fake data with noise
then fits the exponential with the fitdata
function."""  
from pyspecdata import *
from sympy import symbols, latex, Symbol
from sympy import exp as s_exp
# {{{ this is the contents of pylab.py -- works
# need to go through and figure out which lines
# are actually needed and which are not
# -- I have already stripped out some
from matplotlib.pyplot import figure, subplot, show, xlim, ylim, plot, gca
from numpy import * # I think it wasn't importing from numpy b/c it seems we're inside sphinx
# }}}
#{{{creating a fake data recovery curve
tau = nddata(r_[0:2:100j], 'tau')
fake_data = 102*(1-2*exp(-tau*6.0))
fake_data.add_noise(5.0)
#}}}
#{{{ fitting data
f = fitdata(fake_data)
M0,Mi,R1,vd = symbols("M_0 M_inf R_1 tau",real=True)
f.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
logger.info(strm("Functional Form", f.functional_form))
logger.info(strm("Functional Form", f.functional_form))
f.set_guess({M0:-200, Mi:200, R1:2})
f.settoguess()
guess = f.eval(100)
f.fit()
print("output:",f.output())
print("latex:",f.latex())
#}}}
# {{{ this is just to show all the parameters
list_symbs = []
for j,k in f.output().items():
    s_repr = latex(Symbol(j))
    list_symbs.append(f'${s_repr} = {k:0.5g}$')
list_symbs = '\n'.join(list_symbs)
# }}}
T1 = 1./f.output('R_1')
# }}}
with figlist_var() as fl: 
    fl.next('fit with guess')
    fl.plot(guess,label='guess')
    fl.plot(fake_data,'o',label='fake data')
    thisline = fl.plot(f.eval(100),label='fit')
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
