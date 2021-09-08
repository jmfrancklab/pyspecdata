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

init_logging(level="debug")
np.random.seed(15816)
# {{{ helper function(s)
class myfitclass (nddata):
    def __init__(self,*args,**kwargs):
        # copied from fitdata
        self.expression = None
        if isinstance(args[0],nddata):
            # move nddata attributes into the current instance
            myattrs = normal_attrs(args[0])
            for j in range(0,len(myattrs)):
                self.__setattr__(myattrs[j],args[0].__getattribute__(myattrs[j]))
        else:
            nddata.__init__(self,*args,**kwargs)
        return
    @property
    def functional_form(self):
        r'''A property of the fitdata class which is set by the user,
        takes as input a sympy expression of the desired fit
        expression'''
        print("Getting symbolic function")
        return self.symbolic_expr
    @functional_form.setter
    def functional_form(self,this_expr):
        """generate parameter descriptions and a numpy (lambda) function from a sympy expresssion

        Parameters
        ==========
        symbolic_expr: sympy expression
        guesses: dict
            dictionary of keyword arguments for guesses (value) or constraints
            (min/max)

        Returns
        =======
        pars: lmfit.Parameters
        fn: function
            the fit function
        """
        self.expression = this_expr
        # {{{ decide which symbols are parameters vs. variables
        if self.expression is None:
            raise ValueError("what expression are you fitting with??")
        all_symbols = self.expression.atoms(sp.Symbol)
        axis_names = set([sp.Symbol(j) for j in self.dimlabels])
        variable_symbols = axis_names & all_symbols
        parameter_symbols = all_symbols - variable_symbols
        variable_symbols = tuple(variable_symbols)
        self.variable_names = tuple([str(j) for j in variable_symbols])
        parameter_symbols = tuple(parameter_symbols)
        parameter_names = tuple([str(j) for j in parameter_symbols])
        logger.debug(
            strm(
                "all symbols are",
                all_symbols,
                "axis names are",
                axis_names,
                "variable names are",
                self.variable_names,
                "parameter names are",
                parameter_names,
            )
        )
        # }}}
        self.fn = lambdify(
            variable_symbols + parameter_symbols,
            self.expression,
            modules=[{"ImmutableMatrix": np.ndarray}, "numpy", "scipy"],
        )
        self.pars = Parameters()
        for this_name in parameter_names:
            self.pars.add(this_name)
    def set_guess(self,*args,**kwargs):
        """set both the guess and the bounds

        Parameters
        ==========
        guesses: dict of dicts
            each dict has a keyword giving the parameter and a value
            that comprises a dict with guesses (value) and/or constraints
            (min/max)

            Can be passed either as the only argument, or a kwarg called
            guesses, or as the kwargs themselves.
        """
        if len(args) == 1 and type(args[0]) == dict:
            guesses = args[0]
        elif len(kwargs) == 1 and 'guesses' in kwargs.keys():
            guesses = kwargs['guesses']
        else:
            guesses = kwargs
        for this_name in self.pars.keys():
            if this_name in guesses.keys():
                for k,v in guesses[this_name].items():
                    setattr(self.pars[this_name],k,v)
        for j in self.pars:
            logger.info(strm("fit param ---", j))
        logger.info(strm(self.pars))
        return
    def run_lambda(self,pars):
        """actually run the lambda function we separate this in case we want
        our function to involve something else, as well (e.g. taking a Fourier
        transform)"""
        return self.fn(*(self.getaxis(j) for j in self.variable_names), **pars.valuesdict())

    def residual(self, pars, data=None):
        "calculate the residual OR if data is None, return fake data"
        model = self.run_lambda(pars)
        if data is None:
            return model
        return model - data
    @property
    def function_string(self):
        r'''A property of the fitdata class which stores a string
        output of the functional form of the desired fit expression
        provided in func:`functional_form` in LaTeX format'''
        retval = sympy_latex(self.symbolic_expr).replace('$','')
        return r'$f(%s)='%(sympy_latex(sympy_symbol(self.fit_axis))) + retval + r'$'
    @function_string.setter
    def function_string(self):
        raise ValueError("You cannot set the string directly -- change the functional_form property instead!")


# }}}
# {{{ a lot of what's below depends on knowing what the shape and dimension labels of my data are, so define that here
x_vals = linspace(0, 250, 1500)
empty_data = nddata(x_vals, "x").copy(data=False)
# }}}
thisfit = myfitclass(empty_data)
# {{{making sympy expression
A, shift, period, decay, x = sp.symbols("A shift period decay x")
thisfit.functional_form = (A * sp.sin(shift + x / period) * sp.exp(-((x * decay) ** 2)))
thisfit.set_guess(
        A = dict(value=13.0, max=20, min=0.0),
        shift = dict(value=0.0, max=pi / 2.0, min=-pi / 2.0),
        period = dict(value=2, max=10),
        decay = dict(value=0.02, max=0.10, min=0.00),
)
# }}}


# {{{ nddata to generate the fake data
# {{{ create the "true" parameters for the fake data
true_values = {"A": 14.0, "period": 5.4321, "shift": 0.12345, "decay": 0.01000}
p_true = Parameters()
for k, v in true_values.items():
    p_true.add(k, value=v)
logger.info(strm("p_true is:", p_true))
# }}}
mydata = empty_data.copy(data=False)
mydata.data = thisfit.residual(p_true)
mydata.add_noise(2.8)
# }}}
# {{{ nddata of the guess
guess = empty_data.copy(data=False)
guess.data = thisfit.residual(thisfit.pars)
# }}}
# {{{ run the fit and generate nddata
out = minimize(
    thisfit.residual, thisfit.pars, kws={"data": mydata.data}
)
fit = empty_data.copy(data=False)
fit.data = thisfit.residual(out.params)
# }}}

# {{{ report the fit and generate the plot
report_fit(out, show_correl=True, modelpars=p_true)
plot(mydata, "ro", label="data")
plot(fit, "b", alpha=0.5, label="fit")
plot(guess, "g--", label="guess")
# }}}
plt.legend()
plt.show()
