"""
Fit Multiple Data Sets
======================

Fitting multiple (simulated) Gaussian data sets simultaneously.

All minimizers require the residual array to be one-dimensional. Therefore, in
the ``residual`` we need to ```flatten``` the array before returning it.

TODO: this should be using the Model interface / built-in models!

"""
import matplotlib.pyplot as plt
from numpy import linspace, random
import numpy as np
import sympy as sp
from pyspecdata import *
from lmfit import Parameters, minimize, report_fit
logger = init_logging(level='debug')
np.random.seed(15816)
logger.debug(strm("first value",np.random.rand()))
fl=figlist_var()
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

#{{{ creating fake data
#    (five simulated gaussian datasets)
true_values = []
logger.debug(strm("second value",np.random.rand()))
for j in np.arange(5):
    values = {'amp_%d'%(j+1):0.60 + 9.50,
            'cen_%d'%(j+1):-0.20 + 1.20,
            'sig_%d'%(j+1):0.25 + 0.03}
    true_values.append(values)
logger.debug(strm("third value",np.random.rand()))
logger.debug(strm("values for true function",true_values))
mydata_params = []
p_true = [Parameters(),Parameters(),Parameters(),Parameters(),Parameters()]
for j in np.arange(5):
    for k,v in true_values[j].items():
            p_true[j].add(k,value=v)
            mydata_params.append(p_true[j])
logger.info(strm("p_true",p_true))
random.seed(0)
# }}}
x_vals = linspace(-5.0,5.0,151)
empty_data = []    
for _ in np.arange(5):
    edata = nddata(x_vals,'x').copy(data=False)
    empty_data.append(edata)
#empty_data=np.array(empty_data)
#{{{making sympy expression
amp = [sp.symbols('amp_%d' %(i+1)) for i in np.arange(5)]
cen = [sp.symbols('cen_%d' %(i+1)) for i in np.arange(5)]
sig = [sp.symbols('sig_%d' %(i+1)) for i in np.arange(5)]
x = sp.symbols('x')
expr = []
for j in np.arange(5):
    expression = (amp[j]/sig[j]) * sp.exp(-(x-cen[j])**2 / (2.*sig[j]**2)) #preserves integral under curve
    expr.append(expression)
fit_params=[]
parameter_names=[]
fn=[]
for j in np.arange(5):
    parameters, param_names, function = gen_from_expr(expr[j], {'amp_%i'%(j+1):dict(value=10, min=0.0,max=200),
        'cen_%i'%(j+1):dict(value=1.0,min=-2.0,max=2.0),
        'sig_%i'%(j+1):dict(value=0.3,min=0.01,max=3.0),})    
    fit_params.append(parameters)
    parameter_names.append(param_names)
    fn.append(function)
#}}}
def residual(pars, x, k=2, data=None):
    """Calculate total residual for fits of Gaussians to several data sets."""
    parameter_name = parameter_names[k]
    for j in np.arange(3):
        print("PARAMETER NAME IS", parameter_name)
        parlist = [pars[j] for j in parameter_name]
        logger.info(strm("parlist",parlist))
    models =[]
    model = fn[j](x, *parlist)
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
for j in np.arange(5):
    mydata_params.pop(j)
    mydata_params.pop(j)
#}}}    
mydata = []
for j in np.arange(5):
    dat = empty_data[j].copy(data=False)
    dat.data = residual(mydata_params[j],dat.getaxis('x'),k=j)
    dat.add_noise(2.8)
    mydata.append(dat)
guess = []
for j in np.arange(5):
    fit = empty_data[j].copy(data= False)
    fit.data = residual(fit_params[j], empty_data[j].getaxis('x'),k=j)
    guess.append(fit)
###############################################################################
# Run the global fit and show the fitting result
fitting = []
out = []
#raise ValueError("the way that you are calling minimize is absolutely wrong."
#        "In their example, the only call minimize once! "
#        "You are just doing a 1D minimization on each separate dataset -- this is not "
#        "a global fitting.")
for j in np.arange(5):
    out = minimize(residual,fit_params[j], args=(mydata[j].getaxis('x'),j,), kws={'data':mydata[j].data})
    fits=empty_data[j].copy(data=False)
    fits.data = residual(out.params,fits.getaxis('x'),k=j)
    fitting.append(fits)
###############################################################################
# Plot the data sets and fits
plt.figure()
for j in np.arange(5):
    plot(mydata[j],'o',label='data%d'%j)
    plot(fitting[j],'-',label='fitting%d'%j,alpha=0.5)
    plot(guess[j],'--',label='guess%d'%j)
plt.legend()
plt.show()
