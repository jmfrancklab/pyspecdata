"""
Fit Multiple Data Sets
======================

Adapts the example from lmfit to use sympy and pyspecdata to
Fit multiple (simulated) Gaussian data sets simultaneously.

All minimizers require the residual array to be one-dimensional. Therefore, in
the ``residual`` we need to ```flatten``` the array before returning it.

TODO: this should be using the Model interface / built-in models!

"""
import matplotlib.pyplot as plt
from numpy import exp, linspace, pi, random
import numpy as np
import sympy as sp
from pyspecdata import *
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
from collections import ChainMap

#np.random.seed(15816)
# {{{ helper function(s)
def gen_from_expr(expr, global_params={}, n_datasets=3, guesses={}):
    """generate parameter descriptions and a numpy (lambda) function from a sympy expresssion

    Parameters
    ==========
    expr: sympy expression for the 1D function that describes every dataset.
        Note that the parameters for different datasets can be related to each
        other by one or more expressions that depend on global parameters
        (see `global_params`)
    guesses: dict
        dictionary of keyword arguments for guesses (value) or constraints
        (min/max)
    global_params: dict
        if any of the variables is a global variable,
        add to this dictionary as a key,
        whose value gives a new sympy expression,
        where the variables of the new expression
        are assumed to be global variables

        note that if a given symbol in expr is not found in this dictionary,
        then that means that a new, numbered global variable will
        be created for every dataset
    n_datasets: integer
        How many 1D datasets?

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
    global_symbols = set([sp.Symbol(j) for j in (global_params.keys())])
    variable_symbols = axis_names & all_symbols
    parameter_symbols = all_symbols - variable_symbols - global_symbols
    variable_symbols = tuple(variable_symbols)
    variable_names = tuple([str(j) for j in variable_symbols])
    parameter_symbols = tuple(parameter_symbols)
    parameter_names = tuple([str(j) for j in parameter_symbols])
    parameter_names = [
        "%s_%d" % (p, j) for j in range(n_datasets) for p in parameter_names
    ]  # organizes datasets together
    global_symbols = tuple(global_symbols)
    global_names = tuple([str(j) for j in global_symbols])
    print("all symbols are",
            all_symbols,
            "axis names are",
            axis_names,
            "variable names are",
            variable_names,
            "parameter names are",
            parameter_names,
            "global names are",
            global_names,
        )
    # }}}
    pars = Parameters()
    for this_name in parameter_names:
        kwargs = {}
        if this_name in guesses.keys():
            logger.info(strm("applying bounds for", this_name))
            kwargs.update(guesses[str(this_name)])
        pars.add(this_name, **kwargs)
    for j in pars:
        logger.info(strm("fit param ---", j))
    for k,v in global_params.items():    
        g_expr = v 
    g_symbols = g_expr.atoms(sp.Symbol)
    global_fn = lambdify(
            g_symbols,
            g_expr,
            modules = [{"ImmutableMatrix": np.ndarray},"numpy","scipy"],
            )
    local_fn = lambdify(
        variable_symbols + parameter_symbols + global_symbols,
        expr,
        modules=[{"ImmutableMatrix": np.ndarray}, "numpy", "scipy"],
    )
    n_vars = len(variable_names)
    n_fn_params = len(parameter_names)
    print(n_vars)
    print(n_fn_params)
    quit()
    g_vars = 1
    g_fn_params = 1
    def outer_fn(*args):
        var_args = args[0:n_vars]
        par_args = args[n_vars:]
        g_var_args = args[0:g_vars]
        g_par_args = args[g_vars:]
        data = np.empty((n_datasets, 151))
        g_data = np.empty((n_datasets,151))
        for j in range(n_datasets):    
            these_pars = par_args[j * n_fn_params : (j + 1) * n_fn_params]
            g_pars = g_par_args[j * g_fn_params : (j + 1) * g_fn_params]
            g_data = global_fn(*tuple(g_var_args + g_pars))
            data[j, :] = local_fn(*tuple(var_args + these_pars))
        return data

    return pars, parameter_names, outer_fn


# }}}
# {{{ creating fake data
p_true = Parameters()
for j in np.arange(3):
    values = {
        "cen_%d" % (j): -0.20 + 1.2 * np.random.rand(),
        "sig_%d" % (j): 0.25 + 0.03 * np.random.rand(),
        "amp_%d" % (j): 0.6 + 9.5 * np.random.rand(),
    }
    for k, v in values.items():
        p_true.add(k, value=v)
logger.info(strm("p_true is:", p_true))
x_vals = linspace(-1, 2, 151)
empty_data = ndshape([151, 3], ["x", "data"]).alloc(format=None)
# }}}
# {{{ making sympy expression
amp, cen, sig, x = sp.symbols("amp cen sig x")
B, y = sp.symbols("B y")
local_expr = (amp / sig) * sp.exp(
    -((x - cen) ** 2) / (2.0 * sig ** 2)
)  # preserves integral under curve
# seems likely that Parameters is an ordered list, in which case, we don't need
# parameter_names -- **however** we need to check the documentation to see that
# this is true

guess_dict = {("amp_%i" % j): dict(value=15, min=0.0, max=200) for j in range (3)}
guess_dict.update(
    {("cen_%i" % j): dict(value=0.5, min=-2.0, max=2.0) for j in range(3)}
)
guess_dict.update(
    {("sig_%i" % j): dict(value=0.3, min=0.01, max=3.0) for j in range(3)}
)
global_params = {'amp':B*y}
fit_params, parameter_names, myfunc = gen_from_expr(local_expr, global_params = global_params,guesses=guess_dict)
# }}}
def residual(pars, x, data=None):
    "calculate the residual OR if data is None, return fake data"
    logger.info(strm("PARAMETER NAMES ARE:", parameter_names))
    parlist = [pars[j].value for j in parameter_names]
    logger.info(strm("parlist", parlist))
    model = myfunc(x, *parlist)
    if data is None:
        return model
    return (data - model).ravel()


def make_nddata(data):
    return nddata(data, ["datasets", "x"]).setaxis("x", x_vals)


# {{{ nddata to generate fake data
#    (simulated gaussian datasets)
mydata = make_nddata(residual(p_true, x_vals))
mydata.add_noise(0.3)
# }}}
# {{{ nddata of the guess
guess = make_nddata(residual(fit_params, x_vals))
# }}}
# {{{ Run the global fit and generate nddata
out = minimize(
    residual, fit_params, args=(mydata.getaxis("x"),), kws={"data": mydata.data}
)
fits = make_nddata(residual(out.params, x_vals))
# }}}
# {{{ report the fit and generate the plot
report_fit(out, show_correl=True, modelpars=p_true)
for j in range(3):
    plot(x_vals, mydata["datasets", j].data, "o", label="data%d" % j)
    plot(x_vals, fits["datasets", j], "-", label="fitting%d" % j)
    plot(x_vals, guess["datasets", j], "--", label="guess%d" % j)
# }}}
plt.legend()
plt.show()
