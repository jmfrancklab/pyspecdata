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
mydata_params = []
print("TRUEVALUE[1]=",true_values[1])
p_true = Parameters()
for k,v in true_values[0].items():
        p_true.add(k,value=v)
        mydata_params.append(p_true)
p_true1 = Parameters()
for k,v in true_values[1].items():
        p_true1.add(k,value=v)
        mydata_params.append(p_true1)
p_true2 = Parameters()        
for k,v in true_values[2].items():
        p_true2.add(k,value=v)
        mydata_params.append(p_true2)
p_true3 = Parameters()        
for k,v in true_values[3].items():
        p_true3.add(k,value=v)
        mydata_params.append(p_true3)
p_true4 = Parameters()        
for k,v in true_values[4].items():
        p_true4.add(k,value=v)
        mydata_params.append(p_true4) 
mydata_params.pop(1)
mydata_params.pop(1)
mydata_params.pop(2)
mydata_params.pop(2)
mydata_params.pop(3)
mydata_params.pop(3)
mydata_params.pop(4)
mydata_params.pop(4)
mydata_params.pop(5)
mydata_params.pop(5)
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
fit_pars = Parameters()
fits = []
for j in np.arange(5):
    for k,v in fit_params[j].items():
        fit_pars.add(k,value=v)
    fits.append(fit_pars)    
def objective(pars, x, k=1,data=None):
    """Calculate total residual for fits of Gaussians to several data sets."""
    print("PARAMETER NAMES ARE", parameter_names)
    parameter_name = parameter_names[k]
    print("PARAMETER NAME IS",parameter_name)
    for j in np.arange(3):
        parlist = [pars[j] for j in parameter_name]
        logger.info(strm("parlist",parlist))
    models =[]
    #for j in enumerate(parlist):
    print("PARLIST IS", *parlist)
    model = fn[j](*parlist)
    models.append(model)
    if data is None:
        print("DATA IS NOT NONE")
        return models
    ndata = data.shape
    resids = []
    resid = 0.0*data[:]
     # make residual per data set
    for i in range(len(ndata)):
        resid[i:] = data[i:] - models[i:]
        resids.append(resid)
    # now flatten this to a 1D array, as minimize() needs
    return resids.flatten()
mydata = []
for j in np.arange(5):
    dat = empty_data[j].copy(data=False)
    print("MYDATA_PARAMS ARE",mydata_params[j])
    dat.data = objective(mydata_params[j],dat.getaxis('x'),k=j)
    mydata.append(dat)
    print("MYDATA IS",mydata)
#for j in np.arange(5):
#    fit = empty_data[j].copy(data= False)
#    fit.data = objective(fit_pars, fit.getaxis('x'))
#    guess.append(fit)
###############################################################################
# Run the global fit and show the fitting result
fitting = []
#for j in np.arange(5):
#    out = minimize(objective, fit_pars, args=(mydata[j].getaxis('x'),),kws={'data':mydata[j].data})
#    fit=empty_data[j].copy(data=False)
#    fit.data = objective(out.params,empty_data[j].getaxis('x'))
#    report_fit(out, show_correl=True,modelpars=p_true)
#    fitting.append(fit)
###############################################################################
# Plot the data sets and fits
plt.figure()
print("MYDATA[1]",mydata[1])
print(type(mydata[1]))

for j in np.arange(5):
    #print('fit result',fitting[j])
    plot(mydata[j],'-',label='mydata')
    #plot(fitting[j],'o',label='fitting')
    #plot(guess[j],'--',label='guess')
plt.legend()
plt.show()
