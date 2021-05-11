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
# {{{ helper function(s)
def gen_from_expr(expr, guesses={}):
    """generate parameter descriptions and a numpy (lambda) function from a sympy expresssion

    Parameters
    ==========
    expr: sympy expression
    guesses: dict
        dictionary of keyword arguments for guesses (value) or constraints
        (min/max)

    Returns
    =======
    pars: lmfit.Parameters
    parameter_names: tuple
        ordered list of the names of the arguments  to fn
    fn: function
        the fit function
    """
    # {{{ decide which symbols are parameters vs. variables
    all_symbols = expr.atoms(sp.Symbol)
    axis_names = set([sp.Symbol(j) for j in empty_data.dimlabels])
    variable_symbols = axis_names & all_symbols
    parameter_symbols = all_symbols - variable_symbols
    variable_symbols = tuple(variable_symbols)
    variable_names = tuple([str(j) for j in variable_symbols])
    parameter_symbols = tuple(parameter_symbols)
    parameter_names = tuple([str(j) for j in parameter_symbols])
    logger.info(strm("all symbols are", all_symbols, "axis names are", axis_names,
            "variable names are", variable_names, "parameter names are", parameter_names))
    # }}}
    pars = Parameters()
    for this_name in parameter_names:
        kwargs = {}
        if this_name in guesses.keys():
            logger.info(strm("applying bounds for",this_name))
            kwargs.update(guesses[str(this_name)])
        pars.add(this_name, **kwargs)
    for j in pars:
        logger.info(strm("fit param ---",j))
    logger.info(strm(pars))
    fn = lambdify(variable_symbols + parameter_symbols,
            expr,
            modules=[{'ImmutableMatrix':np.ndarray},'numpy','scipy'])
    return pars, parameter_names, fn
# }}}
#{{{ creating fake data
true_values = {'A':14.0,
        'period':5.4321,
        'shift':0.12345,
        'decay':0.01000}
p_true = Parameters()
for k,v in true_values.items():
    p_true.add(k,value=v)
logger.info(strm(p_true))
random.seed(0)
x_vals = linspace(0, 250, 1500)
empty_data = nddata(x_vals,'x').copy(data=False)
#}}}
#{{{making sympy expression
A, shift, period, decay, x = sp.symbols('A shift period decay x')
expr = A*sp.sin(shift+x/period)*sp.exp(-(x*decay)**2)
# seems likely that Parameters is an ordered list, in which case, we don't need
# parameter_names -- **however** we need to check the documentation to see that
# this is true
fit_params, parameter_names, fn = gen_from_expr(expr, {'A':dict(value=13.0, max=20, min=0.0),
            'period':dict(value=2, max=10),
            'shift':dict(value=0.0, max=pi/2., min=-pi/2.),
            'decay':dict(value=0.02, max=0.10, min=0.00),})
#}}}
def residual(pars, x, data=None):
    "calculate the residual OR if data is None, return fake data"
    print("PARAMETER NAMES ARE:",parameter_names)
    parlist = [pars[j] for j in parameter_names]
    print("PARLIST IS",parlist)
    logger.info(strm("parlist",parlist))
    shift = pars['shift']
    if abs(shift) > pi/2:
        shift = shift - sign(shift)*pi
    print("THIS IS PARLIST",*parlist)
    model = fn(x, *parlist)
    if data is None:
        return model
    return model - data
# {{{ nddata to generate the fake data
mydata = empty_data.copy(data=False)
mydata.data = residual(p_true, mydata.getaxis('x'))
mydata.add_noise(2.8)
# }}}
# {{{ nddata of the guess
guess = empty_data.copy(data=False)
guess.data = residual(fit_params, empty_data.getaxis('x'))
# }}}
# {{{ run the fit and generate nddata
out = minimize(residual, fit_params, args=(mydata.getaxis('x'),), kws={'data': mydata.data})
fit = empty_data.copy(data=False)
fit.data = residual(out.params, empty_data.getaxis('x'))
# }}}

# {{{ report the fit and generate the plot
report_fit(out, show_correl=True, modelpars=p_true)
plot(mydata, 'ro', label='data')
plot(fit, 'b', alpha=0.5, label='fit')
plot(guess, 'g--', label='guess')
# }}}
plt.legend()
plt.show()
