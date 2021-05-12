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
np.random.seed(15816)
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
    axis_names = set([sp.Symbol(j) for j in empty_data[0].dimlabels])
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
#{{{ creating fake data(simulated gaussian datasets)
true_values = []
for j in np.arange(3):
    values = {'amp_%d'%(j+1):20 + 2*np.random.rand(),
            'cen_%d'%(j+1):-0.20 + 3.0*np.random.rand(),
            'sig_%d'%(j+1):0.25 + 0.03*np.random.rand()}
    true_values.append(values)
mydata_params = []
p_true = [Parameters(),Parameters(),Parameters(),Parameters(),Parameters()]
for j in np.arange(3):
    for k,v in true_values[j].items():
            p_true[j].add(k,value=v)
            mydata_params.append(p_true[j])
logger.info(strm("p_true",p_true))
random.seed(0)
x_vals = linspace(-5.0,5.0,501)
empty_data = []    
for _ in np.arange(3):
    edata = nddata(x_vals,'x').copy(data=False)
    empty_data.append(edata)
#}}}
#{{{making sympy expression
amp = [sp.symbols('amp_%d' %(i+1)) for i in np.arange(3)]
cen = [sp.symbols('cen_%d' %(i+1)) for i in np.arange(3)]
sig = [sp.symbols('sig_%d' %(i+1)) for i in np.arange(3)]
x = sp.symbols('x')
expr = []
for j in np.arange(3):
    expression = (amp[j]) * sp.exp(-(x-cen[j])**2 / (2.*sig[j]**2)) #preserves integral under curve
    expr.append(expression)
#seems likely that Parameters is an ordered list, in which case, we don't need
#parameter names -- **however** we need to check the documentation to see that
#this is true
fit_params=[]
parameter_names=[]
fn=[]
for j in np.arange(3):
    parameters, param_names, function = gen_from_expr(expr[j], {'amp_%i'%(j+1):dict(value=21.0, min=0.0,max=200),
        'cen_%i'%(j+1):dict(value=1.0,min=-1.0,max=5.0),
        'sig_%i'%(j+1):dict(value=0.25,min=0.01,max=5.0),})    
    fit_params.append(parameters)
    parameter_names.append(param_names)
    fn.append(function)
#}}}
def residual(pars, x, k=2, data=None):
    """Calculate total residual for fits of Gaussians to several data sets."""
    parameter_name = parameter_names[k]
    logger.info(strm("PARAMETER NAMES ARE:", parameter_name))
    parlist = [pars[j] for j in parameter_name]
    logger.info(strm("parlist",parlist))
    models =[]
    model = fn[k](x, *parlist)
    if data is None:
        return model
    ndata = data.shape
    resids = []
    resid = 0.0*data[:]
    for i in range(len(ndata)):
        resid[i:] = data[i:] - model
        resids.append(resid)
    # now flatten this to a 1D array, as minimize() needs
    np.concatenate(resids)
    for j in range(len(resids)):
        return resids[j].flatten()
#{{{needed as for some reason there are repeats in mydata_params
for j in np.arange(3):
    mydata_params.pop(j)
    mydata_params.pop(j)
#}}}  
#{{{nddata to generate the fake data
mydata = []
for j in np.arange(3):
    dat = empty_data[j].copy(data=False)
    dat.data = residual(mydata_params[j],dat.getaxis('x'),k=j,data=None)
    dat.add_noise(0.8)
    mydata.append(dat)
#}}}
#{{{nddata of the guess
guess = []
for j in np.arange(3):
    fit = empty_data[j].copy(data= False)
    fit.data = residual(fit_params[j], empty_data[j].getaxis('x'),k=j)
    guess.append(fit)
#}}}    
#{{{ Run the global fit and generate nddata
fitting = []
out = []
#raise ValueError("the way that you are calling minimize is absolutely wrong."
#        "In their example, the only call minimize once! "
#        "You are just doing a 1D minimization on each separate dataset -- this is not "
#        "a global fitting.")
for j in np.arange(3):
    out = minimize(residual,fit_params[j], args=(mydata[j].getaxis('x'),j,), kws={'data':mydata[j].data})
    fits=empty_data[j].copy(data=False)
    fits.data = residual(out.params,fits.getaxis('x'),k=j)
    fitting.append(fits)
#}}}

#{{{report the fit and generate the plot
report_fit(out,show_correl=True,modelpars=mydata_params)
for j in np.arange(3):
    plot(mydata[j],'o',label='data%d'%j)
    plot(fitting[j],'-',label='fitting%d'%j)
    plot(guess[j],'--',label='guess%d'%j)
#}}}    
plt.legend()
plt.show()
