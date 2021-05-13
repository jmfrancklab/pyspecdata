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
def gen_from_expr(expr, global_params, n_datasets, guesses={}):
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
    axis_names = set([sp.Symbol(j) for j in empty_data[0].dimlabels])
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
    def outer_fun(*args):
        # outer function goes here
    return pars, parameter_names, outer_fn
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
    print(ndata)
    print("MODEL IS",model)
    quit()
    for i in range(ndata):
        resid[i, :] = data[i, :] - model
    # now flatten this to a 1D array, as minimize() needs
    np.concatenate(resids)
    for j in range(len(resids)):
        return resids[j].flatten()
#{{{making sympy expression
expr = []
for j in np.arange(3):
    amp = [sp.symbols('amp_%d' %(j+1))]
    cen = [sp.symbols('cen_%d' %(j+1))]
    sig = [sp.symbols('sig_%d' %(j+1))]
    x = sp.symbols('x')
    expression = (amp[j]) * sp.exp(-(x-cen[j])**2 / (2.*sig[j]**2)) #preserves integral under curve
    #expr.append(expression)
#seems likely that Parameters is an ordered list, in which case, we don't need
#parameter names -- **however** we need to check the documentation to see that
#this is true
#parameter_names=[]
#fn=[]
#guess_params = []
#for j in np.arange(3):
    parameters, param_names, function = gen_from_expr(expr[j], {'amp_%i'%(j+1):dict(value=21.0, min=0.0,max=200),
        'cen_%i'%(j+1):dict(value=0.5,min=-5.0,max=5.0),
        'sig_%i'%(j+1):dict(value=0.25,min=0.01,max=5.0),})
    #guess_params.append(parameters)
    #parameter_names.append(param_names)
    #fn.append(function)

#}}}

#{{{ creating fake data
#    (simulated gaussian datasets)
#true_values = []
#mydata_params = []
#for j in np.arange(3):
    values = {'amp_%d'%(j+1):20 + 2*np.random.rand(),
            'cen_%d'%(j+1):-0.20 + 3.0*np.random.rand(),
            'sig_%d'%(j+1):0.25 + 0.03*np.random.rand()}
    #true_values.append(values)
    #p_true = [Parameters(),Parameters(),Parameters()]
    for k,v in values.items():
            p_true.add(k,value=v)
    #mydata_params.append(p_true[j])
    random.seed(0)
    x_vals = linspace(-5.0,5.0,501)
    #empty_data = []    
    empty_data = nddata(x_vals,'x').copy(data=False)
    #empty_data.append(edata)
#}}}
#{{{nddata to generate the fake data
#mydata = []
#print("mydata_params are",mydata_params)
#for j in np.arange(3):
    mydata = empty_data[j].copy(data=False)
    mydata.data = residual(p_true,mydata.getaxis('x'),data=None)
    dat.add_noise(0.8)
    mydata.append(dat)
mydata = np.array(mydata)
    
fit_params = Parameters()
for iy,y in enumerate(mydata):
    fit_params.add('amp_%i'%(iy+1),value=0.5,min=0.0,max=200)
    fit_params.add('cen_%i'%(iy+1),value=0.4,min=-2.0,max=2.0)
    fit_params.add('sig_%i'%(iy+1),value=0.3,min=0.01,max=3.0)
#}}}
#{{{nddata of the guess
guess = []
for j in np.arange(3):
    fit = empty_data[j].copy(data= False)
    fit.data = residual(guess_params[j], empty_data[j].getaxis('x'),k=j)
    guess.append(fit)
#}}}    
#{{{ Run the global fit and generate nddata
fitting = []
#raise ValueError("the way that you are calling minimize is absolutely wrong."
#        "In their example, the only call minimize once! "
#        "You are just doing a 1D minimization on each separate dataset -- this is not "
#        "a global fitting.")
out = minimize(objective,fit_params, args=(x_vals,), kws={'data':mydata.data})
#for j in np.arange(5):    
#    fits=empty_data[j].copy(data=False)
#    fits.data = residual(out.params,fits.getaxis('x'),k=j)
#    fitting.append(fits)
#}}}

#{{{report the fit and generate the plot
#report_fit(out,show_correl=True,modelpars=mydata_params)
for j in np.arange(3):
    plot(mydata[j],'o',label='data%d'%j)
    #plot(fitting[j],'-',label='fitting%d'%j)
    plot(guess[j],'--',label='guess%d'%j)
#}}}    
plt.legend()
plt.show()
