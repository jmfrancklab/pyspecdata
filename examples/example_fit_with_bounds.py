"""
Fit Using Bounds
================

Adapt one of the examples from lmfit to use sympy and pyspecdata.
Eventually, we want to use parts of this in the base library, so we don't have
to retype them every time.

"""
import matplotlib.pyplot as plt
from numpy import exp, linspace, pi, random, sign, sin
import sympy as sp
from lmfit import Parameters, minimize
from pyspecdata import *
from lmfit.printfuncs import report_fit
import numpy as np
#{{{ creating data
p_true = Parameters()
p_true.add('amp', value=14.0)
p_true.add('period', value=5.4321)
p_true.add('shift', value=0.12345)
p_true.add('decay', value=0.01000)
def make_data(pars, x,data=None):
    argu = (x * pars['decay'])**2
    shift = pars['shift']
    if abs(shift) > pi/2:
        shift = shift - sign(shift)*pi
    model = pars['amp'] * sin(shift + x/pars['period']) * exp(-argu)
    if data is None:
        return model
    return model - data
random.seed(0)
x = np.linspace(0, 250, 32)
noise = random.normal(scale=2.80, size=x.size)
data = make_data(p_true, x) + noise
mydata = nddata(data,[-1],['x']).setaxis('x',x)
print(mydata)
#}}}
x_axis = nddata(np.linspace(0,250,1500),'x_axis')
A, shift, period, decay, x = sp.symbols('A shift period decay x')
expr = A*sp.sin(shift+x/period)*sp.exp(-(x*decay)**2)
print(expr.atoms(sp.Symbol))
fit_params = Parameters()
axis_names = set(mydata.dimlabels)
variable_names = axis_names & expr.atoms(sp.Symbol)
parameter_names = expr.atoms(sp.Symbol) - variable_names
for this_symbol in expr.atoms(sp.Symbol):
    fit_params.add('%s'%this_symbol)
for j in fit_params:
    print("fit param ---",j)
def make_fn(pars,x,data=None):
    x_axis =nddata(np.linspace(0,250,1500),'x_axis') 
    fn = lambdify([x_axis],expr.subs({A:pars['A'],shift:pars['shift'],
        period:pars['period'],decay:pars['decay']}),
        modules=[{'ImmutableMatrix':np.ndarray},'numpy','scipy'])
    l_fn = fn(x_axis)
    fn_nddata = nddata(l_fn,[-1],['x']).setaxis('x',x_axis)
    model = fn_nddata
    if data is None:
        return model
    return model - data
fit_params = Parameters()
fit_params.add('A', value=13.0, max=20, min=0.0)
fit_params.add('period', value=2, max=10)
fit_params.add('shift', value=0.0, max=pi/2., min=-pi/2.)
fit_params.add('decay', value=0.02, max=0.10, min=0.00)
out = minimize(make_fn, fit_params, args=(x_axis,), kws={'data': mydata})
fit = residual(out.params, x)

###############################################################################
# This gives the following fitting results:

report_fit(out, show_correl=True, modelpars=p_true)

###############################################################################
# and shows the plot below:
#
plt.plot(x, data, 'ro')
plt.plot(x, fit, 'b')
plt.show()
