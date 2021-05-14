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
def gen_from_expr(expr, global_params, n_datasets=3, guesses={}):
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
    logger.info(strm(pars))
    fn = lambdify(variable_symbols + parameter_symbols,
            expr,
            modules=[{'ImmutableMatrix':np.ndarray},'numpy','scipy'])
    n_vars = 1
    n_fn_params = 3
    def outer_fun(*args):
        # outer function goes here
        var_args = args[0:n_vars]
        par_args = args[n_vars:]
        for j in range(n_datasets):
            data = np.empty((n_datasets,151)) 
            these_pars = par_args[j*n_fn_params:(j+1)*n_fn_params]
            data[j,:] = fn(*tuple(var_args+these_pars))
        return data    
    return pars, parameter_names, outer_fun
# }}}
def residual(pars, x, data=None):
    """Calculate total residual for fits of Gaussians to several data sets."""
    #parameter_name = parameter_names[k]
    #print("PARAMETER NAMES ARE:", parameter_name)
    #print("pars are",pars)
    parlist = [pars[j] for j in parameter_names]
    #print("parlist",parlist)
    model = fn(x, *parlist)
    if data is None:
        return model
    ndata = data.shape
    resid = 0.0*data[:]
    for i in range(ndata):
        resid[i, :] = data[i, :] - model[i,:]
        print("RESID",resid)
        #resid.flatten()
    # now flatten this to a 1D array, as minimize() needs
        np.concatenate(resid)
    for j in range(len(resid)):
        return resid.flatten()
#{{{making sympy expression
p_true = Parameters()
for j in np.arange(3):
    values = {'amp_%d'%(j):20 + 2*np.random.rand(),
            'cen_%d'%(j):-0.20 + 3.0*np.random.rand(),
            'sig_%d'%(j):0.25 + 0.03*np.random.rand()}
    for k,v in values.items():
            p_true.add(k,value=v)
x_vals = linspace(-5,5,151)
empty_data = ndshape([151,3],['x','data']).alloc(format=None)
amp,cen,sig,x=sp.symbols('amp cen sig x')
expression = (amp) * sp.exp(-(x-cen)**2 / (2.*sig**2)) #preserves integral under curve
#seems likely that Parameters is an ordered list, in which case, we don't need
#parameter names -- **however** we need to check the documentation to see that
#this is true
fit_params, parameter_names, fn = gen_from_expr(expression, {'amp_%i'%(j+1):dict(value=21.0, min=0.0,max=200),
        'cen_%i'%(j+1):dict(value=0.5,min=-5.0,max=5.0),
        'sig_%i'%(j+1):dict(value=0.25,min=0.01,max=5.0),})
#}}}

#{{{ creating fake data
#    (simulated gaussian datasets)
mydata = ndshape(empty_data)
mydata.data = residual(p_true,x_vals)
#mydata.add_noise(0.8)
mydata = np.array(mydata)
print(mydata)
#}}}
#{{{nddata of the guess
guess = empty_data.copy(data=False)
guess.data = residual(fit_params, x_vals)
guess = np.array(mydata,ndmin=3)
#}}}    
#{{{ Run the global fit and generate nddata
#raise ValueError("the way that you are calling minimize is absolutely wrong."
#        "In their example, the only call minimize once! "
#        "You are just doing a 1D minimization on each separate dataset -- this is not "
#        "a global fitting.")
print("FIT PARAMS ARE",fit_params)
#out = minimize(residual,fit_params, args=(mydata.getaxis('x'),), kws={'data':mydata.data})
#fit = empty_data.copy(data=False)
#fit.data = residual(out.params, empty_data.getaxis('x'))
#quit()
#}}}
print(nddata(mydata['data':1]))
#{{{report the fit and generate the plot
#report_fit(out,show_correl=True,modelpars=mydata_params)
plot(mydata,label='data')
  #plot(fitting[j],'-',label='fitting%d'%j)
#plot(guess,'--',label='guess%d'%j)
#}}}    
plt.legend()
plt.show()
