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
np.random.seed(15816)
# {{{ helper function(s)
def gen_from_expr(expr, global_params=1, n_datasets=3, guesses={}):
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
    variable_symbols = axis_names & all_symbols
    parameter_symbols = all_symbols - variable_symbols
    variable_symbols = tuple(variable_symbols)
    variable_names = tuple([str(j) for j in variable_symbols])
    parameter_symbols = tuple(parameter_symbols)
    parameter_names = tuple([str(j) for j in parameter_symbols])
    parameter_names = ['%s_%d'%(p,j)
            for j in range(n_datasets)
            for p in parameter_names] # organizes datasets together
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
    fn = lambdify(variable_symbols + parameter_symbols,
            expr,
            modules=[{'ImmutableMatrix':np.ndarray},'numpy','scipy'])
    def outer_fun(*args,n_vars=1,n_fn_params=3):
        # outer function goes here
        var_args = args[0:n_vars]
        par_args = args[n_vars:]
        data = np.empty((n_datasets,151)) 
        for j in range(n_datasets):
            these_pars = par_args[j*n_fn_params:(j+1)*n_fn_params]
            data[j,:] = fn(*tuple(var_args+these_pars))
        return data    
    return pars, parameter_names, outer_fun
# }}}
def residual(pars, x, data=None):
    """Calculate total residual for fits of Gaussians to several data sets."""
    logger.info(strm("PARAMETER NAMES ARE:", parameter_names))
    parlist = [pars[j].value for j in parameter_names]
    logger.info(strm("parlist",parlist))
    model = function(x, *parlist, n_vars = 1, n_fn_params=3)
    if data is None:
        return model
    ndata = data.shape
    resid = 0.0*data[:]
    for i in range(3):
        resid[i, :] = data[i, :] - model[i,:]
        logger.info(strm("RESID",resid))
    # now flatten this to a 1D array, as minimize() needs
        np.concatenate(resid)
    for j in range(len(resid)):
        return resid.flatten()
#{{{making sympy expression
p_true = Parameters()
for j in np.arange(3):
    values = {'cen_%d'%(j):-0.20 + 1.2*np.random.rand(),
            'sig_%d'%(j):0.25 + 0.03*np.random.rand(),
            'amp_%d'%(j):0.6 + 9.5*np.random.rand()}
    for k,v in values.items():
            p_true.add(k,value=v)
logger.info(strm("p_true is:",p_true))            
x_vals = linspace(-1,2,151)
empty_data = ndshape([151,3],['x','data']).alloc(format=None)
amp,cen,sig,x=sp.symbols('amp cen sig x')
expression = (amp/sig) * sp.exp(-(x-cen)**2 / (2.*sig**2)) #preserves integral under curve
#seems likely that Parameters is an ordered list, in which case, we don't need
#parameter names -- **however** we need to check the documentation to see that
#this is true
fit_params, parameter_names, function = gen_from_expr(expression, guesses = {'amp_0':dict(value=15, min=0.0,max=200),
            'cen_0':dict(value=0.5,min=-2.0,max=2.0),
            'sig_0':dict(value=0.3,min=0.01,max=3.0),
            'amp_1':dict(value=15,min=0.0,max=200),
            'cen_1':dict(value=0.5,min=-2.0,max=2.0),
            'sig_1':dict(value=0.3,min=0.01,max=3.0),
            'amp_2':dict(value=15,min=0.0,max=200),
            'cen_2':dict(value=0.4,min=-2.0,max=2.0),
            'sig_2':dict(value=0.3,min=0.01,max=3.0)})
#}}}
#{{{ creating fake data
#    (simulated gaussian datasets)
dat = residual(p_true,x_vals)
mydata = nddata(dat.data,['datasets','x']).setaxis('x',x_vals)
mydata.add_noise(1.8)
#}}}
#{{{nddata of the guess
gues = residual(fit_params,x_vals)
guess = nddata(gues.data,['datasets','x']).setaxis('x',x_vals)
#}}}    
#{{{ Run the global fit and generate nddata
#raise ValueError("the way that you are calling minimize is absolutely wrong."
#        "In their example, the only call minimize once! "
#        "You are just doing a 1D minimization on each separate dataset -- this is not "
#        "a global fitting.")
out = minimize(residual,fit_params, args=(mydata.getaxis('x'),), kws={'data':mydata.data})
fit = residual(out.params, x_vals)
fits = nddata(fit.data,['datasets','x']).setaxis('x',x_vals)
#}}}
#{{{report the fit and generate the plot
report_fit(out,show_correl=True,modelpars=p_true)
for j in range(3):    
    plot(x_vals,mydata['datasets',j].data,'o',label='data%d'%j)
    plot(x_vals,fits['datasets',j],'-',label='fitting%d'%j)
    plot(x_vals,guess['datasets',j],'--',label='guess%d'%j)
#}}}    
plt.legend()
plt.show()
