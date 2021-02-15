"""
Fit Using Bounds
================

Adapt one of the examples from lmfit to use sympy and pyspecdata.
Eventually, we want to use parts of this in the base library, so we don't have
to retype them every time.

"""
import matplotlib.pyplot as plt
from numpy import exp, linspace, pi, random, sign, sin
from pyspecdata import *
from symfit import Parameter, Variable
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import sympy as sp
import numpy as np
import symfit as sf
from pylab import ndarray
#{{{make sympy expression with symfit and lambdify
A = Parameter('A', value = 14.0)
period = Parameter('period', value = 5.4321)
shift = Parameter('shift', value = 0.12345)
decay = Parameter('decay', value = 0.01000)
x = Variable('x')
expression = A*sp.sin(shift + x/period)*sp.exp(-(x*decay)**2).subs({A:A, shift:shift, 
    period:period,
    decay:decay})
#}}}
#{{{parameters for actual data
p_true = Parameters()
p_true.add('A', value=14.0)
p_true.add('period', value=5.4321)
p_true.add('shift', value=0.12345)
p_true.add('decay', value=0.01000)
#}}}
fit_function = sf.lambdify([x],expression.subs({A:A.value,
        period:period.value,
        shift:shift.value,
        decay:decay.value}),
        modules=[{'ImmutableMatrix':ndarray}, 'numpy','scipy'])

#{{{ define the equation that we are fitting to and residual
def residual(pars, x, data=None):
    shift = pars['shift']
    if abs(shift) > pi/2:
        shift = shift - sign(shift)*pi
    model = fit_function        
    if data is None:
        return model
    return model(x) - data(x)
#}}}
#{{{making data taht is to be fit
random.seed(0)
x = nddata(linspace(0, 250, 1500),'x')
noise = random.normal(scale=2.80)#, size=x.size)
data = residual(p_true, x)
#}}}
#{{{fitting parameters for lmfit
fit_params = Parameters()
fit_params.add('amp', value=13.0, max=20, min=0.0)
fit_params.add('period', value=2, max=10)
fit_params.add('shift', value=0.0, max=pi/2., min=-pi/2.)
fit_params.add('decay', value=0.02, max=0.10, min=0.00)
#}}}
#{{{minimizing lmfit function
out = minimize(residual, fit_params, args=(x,), kws={'data': data})
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
