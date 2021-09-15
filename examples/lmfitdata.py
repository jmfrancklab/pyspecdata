import matplotlib.pyplot as plt
from numpy import exp, linspace, pi, random, sign, sin
import sympy as sp
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import numpy as np
from pyspecdata import *
class lmfitdata (nddata):
    r'''Inherits from an nddata and enables curve fitting through use of a sympy expression.

    The user creates a lmfitdata class object from an existing nddata
    class object, and on this lmfitdata object can define the
    :func:`functional_form` of the curve it would like to fit to the
    data of the original nddata.
    This functional form must be provided as a sympy expression, with
    one of its variables matching the name of the dimension that the
    user would like to fit to.
    '''
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
        r'''A property of the myfitclass class which is set by the user,
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
        r'''A property of the myfitclass class which stores a string
        output of the functional form of the desired fit expression
        provided in func:`functional_form` in LaTeX format'''
        retval = sympy_latex(self.symbolic_expr).replace('$','')
        return r'$f(%s)='%(sympy_latex(sympy_symbol(self.fit_axis))) + retval + r'$'
    @function_string.setter
    def function_string(self):
        raise ValueError("You cannot set the string directly -- change the functional_form property instead!")
