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
from lmfit.printfuncs import report_fit
import numpy as np
from pyspecdata import *
#{{{ creating data
true_values = {'amp':14.0,
        'period':5.4321,
        'shift':0.12345,
        'decay':0.01000}
p_true = Parameters()
for k,v in true_values.items():
    p_true.add(k,value=v)
def residual_orig(pars, x, data=None):
    argu = (x * pars['decay'])**2
    shift = pars['shift']
    if abs(shift) > pi/2:
        shift = shift - sign(shift)*pi
    model = pars['amp'] * sin(shift + x/pars['period']) * exp(-argu)
    if data is None:
        return model
    return model - data
random.seed(0)
x = linspace(0, 250, 1500)
noise = random.normal(scale=2.80, size=x.size)
data = residual_orig(p_true, x) + noise
mydata = nddata(data,[-1],['x']).setaxis('x',x)
#}}}
A, shift, period, decay, x = sp.symbols('A shift period decay x')
expr = A*sp.sin(shift+x/period)*sp.exp(-(x*decay)**2)
print(expr.atoms(sp.Symbol))
fit_params = Parameters()
# {{{ decide which symbols are parameters vs. variables
all_symbols = expr.atoms(sp.Symbol)
axis_names = set([sp.Symbol(j) for j in mydata.dimlabels])
variable_symbols = axis_names & all_symbols
parameter_symbols = all_symbols - variable_symbols
variable_symbols = tuple(variable_symbols)
variable_names = tuple([str(j) for j in variable_symbols])
parameter_symbols = tuple(parameter_symbols)
parameter_names = tuple([str(j) for j in parameter_symbols])
print("all symbols are",all_symbols,"axis names are",axis_names,"variable names are",variable_names,"parameter names are",parameter_names)
# }}}
guesses = {'A':dict(value=13.0, max=20, min=0.0),
        'period':dict(value=2, max=10),
        'shift':dict(value=0.0, max=pi/2., min=-pi/2.),
        'decay':dict(value=0.02, max=0.10, min=0.00),}
for this_symbol in parameter_names:
    kwargs = {}
    if this_symbol in guesses.keys():
        print("applying bounds for",this_symbol)
        kwargs.update(guesses[str(this_symbol)])
    fit_params.add(this_symbol, **kwargs)
for j in fit_params:
    print("fit param ---",j)
print(fit_params)
fn = lambdify(variable_symbols + parameter_symbols,
        expr,
        modules=[{'ImmutableMatrix':np.ndarray},'numpy','scipy'])
def residual(pars, x, data=None):
    parlist = [fit_params[j] for j in parameter_names]
    model = fn(x, *parlist)
    if data is None:
        return model
    return model - data
out = minimize(residual, fit_params, args=(mydata.getaxis('x'),), kws={'data': mydata.data})
fit = mydata.copy(data=False)
guess = mydata.copy(data=False)
guess.data = residual(fit_params, mydata.getaxis('x'))
fit.data = residual(out.params, mydata.getaxis('x'))

###############################################################################
# This gives the following fitting results:

report_fit(out, show_correl=True, modelpars=p_true)

###############################################################################
# and shows the plot below:
#
plot(mydata, 'ro')
plot(fit, 'b', alpha=0.5)
plot(guess, 'g--')
plt.show()
