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
    print("all symbols are", all_symbols, "axis names are", axis_names,
            "variable names are", variable_names, "parameter names are", parameter_names)
    # }}}
    pars = Parameters()
    for this_name in parameter_names:
        kwargs = {}
        if this_name in guesses.keys():
            print("applying bounds for",this_name)
            kwargs.update(guesses[str(this_name)])
        pars.add(this_name, **kwargs)
    for j in pars:
        print("fit param ---",j)
    print(pars)
    fn = lambdify(variable_symbols + parameter_symbols,
            expr,
            modules=[{'ImmutableMatrix':np.ndarray},'numpy','scipy'])
    return pars, parameter_names, fn
# }}}
#{{{creating true values for data 
true_p = []
for j in np.arange(5):    
    true_values = {'amp%d'%j:0.60 + 9.50*np.random.rand(),
            'cen%d'%j:-0.20 + 1.20*np.random.rand(),
            'sigma%d'%j:0.25 + 0.03*np.random.rand()}
    p_true = Parameters()
    for k,v in true_values.items():
        p_true.add(k,value=v)
    true_p.append(p_true)   
x_vals = np.linspace(-1, 2, 151)
empty_data = nddata(x_vals,'x').copy(data=False)
#}}}
#{{{sympy expression
amp, cen, sigma, x = sp.symbols('amp cen sigma x')
expr = amp*sp.exp(-(x-cen)**2 / (2*sigma**2))
#}}}
for j in np.arange(5):
    fit_params, parameter_names, fn = gen_from_expr(expr, {'amp%d'%j: dict(value=0.5, max=200, min=0.5),
        'cen%d'%j:dict(value=0.4, max=2.0, min=-2.0),
        'sigma%d'%j:dict(value=0.3, max=3.0,min=0.01)})
def gauss(pars, x, data = None):
    """Gaussian lineshape."""
    parlist = [pars[j] for j in parameter_names]
    print("parlist",parlist)
    model = fn(x, *parlist)
    if data is None:
        return model
    return model - data


def gauss_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    parlist = [params[j] for j in parameter_names]
    print("parlist",parlist)
    amp = params['amp_%i' % (i+1)]
    cen = params['cen_%i' % (i+1)]
    sigma = params['sigma_%i' % (i+1)]
    return gauss(x, amp, cen, sigma)


def objective(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]
     # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - gauss_dataset(params, i, x)

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()


###############################################################################
# Create five simulated Gaussian data sets
for iy, y in enumerate(mydata):
    fit_params.add('amp%i'%(iy+1),value=0.5,min=0.0, max=200)
    fit_params.add('cen%i'(iy+1),value=0.4,min=-2.0,max=2.0)
    fit_params.add('sig_%i'%(iy+1),value=0.3,min=0.01,max=3.0)
mydata= []
for j in np.arange(5):
    dat = empty_data.copy(data=False)
    dat = gen_from_expr(expr, {'amp%d'%j:dict(value)
    mydata.append(dat)
    mydata[j] = nddata(mydata[j],[-1],['x'])
quit()    
guess = []
for j in np.arange(5):
    guess_data = empty_data.copy(data=False)
    guess_data = gauss(fit_params, empty_data.getaxis('x'))
    guess.append(guess_data)
    guess[j] = nddata(guess[j],[-1],['x'])
# Constrain the values of sigma to be the same for all peaks by assigning
# sig_2, ..., sig_5 to be equal to sig_1.

for iy in (2, 3, 4, 5):
    fit_params['sig_%i' % iy].expr = 'sig_1'

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
