"""
Fit Multiple Data Sets
======================

Fitting multiple (simulated) Gaussian data sets simultaneously.

All minimizers require the residual array to be one-dimensional. Therefore, in
the ``objective`` we need to ```flatten``` the array before returning it.

TODO: this should be using the Model interface / built-in models!

"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from pyspecdata import *
from lmfit import Parameters, minimize, report_fit

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

###############################################################################
# Create five simulated Gaussian data sets
x = np.linspace(-1,2,151)
true_values = []
for j in np.arange(5):
    values = {'amp_%d'%(j+1):0.60 + 9.50*np.random.rand(),
            'cen_%d'%(j+1):-0.20 + 1.20*np.random.rand(),
            'sig_%d'%(j+1):0.25 + 0.03*np.random.rand()}
    true_values.append(values)
print("TRUE VALUES ARE:",true_values)      
p_true = Parameters()
mydata_params = []
for j in np.arange(5):
    for k,v in true_values[j].items():
        p_true.add(k,value=v)
    mydata_params.append(p_true)    
  
empty_data = []    
for _ in np.arange(5):
    edata = nddata(x,'x').copy(data=False)
    empty_data.append(edata)
empty_data=np.array(empty_data)
fit_params = []
amp = [sp.symbols('amp_%d' %(i+1)) for i in np.arange(5)]
cen = [sp.symbols('cen_%d' %(i+1)) for i in np.arange(5)]
sig = [sp.symbols('sig_%d' %(i+1)) for i in np.arange(5)]
expr = []
for j in np.arange(5):
    expression = amp[j] * sp.exp(-(x-cen[j])**2 / (2.*sig[j]**2))
    expr.append(expression)
fit_params=[]
parameter_names=[]
fn=[]
for j in np.arange(5):
    parameters, param_names, function = gen_from_expr(expr[j], {'amp_%i'%(j+1):dict(value=0.5, min=0.0,max=200),
        'cen_%i'%(j+1):dict(value=0.4,min=-2.0,max=2.0),
        'sig_%i'%(j+1):dict(value=0.3,min=0.01,max=3.0)})
    fit_params.append(parameters)
    parameter_names.append(param_names)
    fn.append(function)
print("PARAMETER NAMES ARE",parameter_names)    
def objective(pars, x, data=None):
    """Calculate total residual for fits of Gaussians to several data sets."""
    parlist =[]
    for j in np.arange(5):
        params = [pars[j] for j in parameter_names]
        parlist.append(params)
    print("PARLIST:",parlist)    
    model = fn(x, *parlist)
    if data is None:
        return model
    ndata, _ = data.shape
    resid = 0.0*data[:]
     # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - model
    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()
mydata = []
for j in np.arange(5):
    dat = empty_data[j].copy(data=False)
    dat.data = objective(p_true,dat.getaxis('x'))
    mydata.append(dat)
print(mydata)
quit()
guess = []
for _ in np.arange(5):
    fit = empty_data.copy(data= False)
    fit.data = objective(fit_params, empty_data.getaxis('x'))
    guess.append(fit)
quit()    




data = []
for _ in np.arange(5):
    params = Parameters()
    amp = 0.60 + 9.50*np.random.rand()
    cen = -0.20 + 1.20*np.random.rand()
    sig = 0.25 + 0.03*np.random.rand()
    dat =gauss(x, amp, cen, sig) + np.random.normal(size=x.size, scale=0.1)
    data.append(dat)
print(data)
data = np.array(data)
print("NPARRAY OF DATA:",data)
fit_params = Parameters()
###############################################################################
# Run the global fit and show the fitting result

out = minimize(objective, fit_params, args=(x, data))
report_fit(out.params)

###############################################################################
# Plot the data sets and fits

plt.figure()
for i in range(5):
    y_fit = gauss_dataset(out.params, i, x)
    plt.plot(x, data[i, :], 'o', x, y_fit, '-')
plt.show()
