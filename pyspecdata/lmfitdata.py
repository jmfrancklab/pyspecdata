# just put this in the package
import sympy as sp
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import numpy as np
from .core import nddata, normal_attrs, issympy
from .general_functions import strm
import logging
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
        fit_axis = None
        if 'fit_axis' in list(kwargs.keys()):
            fit_axis = kwargs.pop('fit_axis')
        if isinstance(args[0],nddata):
            # move nddata attributes into the current instance
            myattrs = normal_attrs(args[0])
            for j in range(0,len(myattrs)):
                self.__setattr__(myattrs[j],args[0].__getattribute__(myattrs[j]))
        else:
            nddata.__init__(self,*args,**kwargs)
        if fit_axis is None:
            if len(self.dimlabels) == 1:
                fit_axis = self.dimlabels[0]
            else:
                raise IndexError("Right now, we can only auto-determine the fit axis if there is a single axis")
        
        self.fit_axis = fit_axis
        self.set_to = None
        self.set_indices = None
        self.active_indices = None
        self.expression = None
        return
    
    @property
    def functional_form(self):
        r'''A property of the myfitclass class which is set by the user,
        takes as input a sympy expression of the desired fit
        expression'''
        print("Getting symbolic function")
        return self.expression
    @functional_form.setter
    def functional_form(self,this_expr):
        """generate parameter descriptions and a numpy (lambda) function from a sympy expresssion

        Parameters
        ==========
        symbolic_expr: sympy expression
        """
        assert issympy(this_expr), "for now, the functional form must be a sympy expression"
        self.expression = this_expr
        # {{{ decide which symbols are parameters vs. variables
        if self.expression is None:
            raise ValueError("what expression are you fitting with??")
        all_symbols = self.expression.atoms(sp.Symbol)
        axis_names = set([sp.Symbol(j) for j in self.dimlabels])
        variable_symbols = (axis_names & all_symbols)
        self.parameter_symbols = all_symbols - variable_symbols
        this_axis = variable_symbols
        variable_symbols = tuple(variable_symbols)
        self.variable_names = tuple([str(j) for j in variable_symbols])
        parameter_symbols = tuple(self.parameter_symbols)
        self.parameter_names = tuple([str(j) for j in self.parameter_symbols])
        self.fit_axis = set(self.dimlabels)
        self.symbol_list = [str(j) for j in parameter_symbols]
        logging.debug(
            strm(
                "all symbols are",
                all_symbols,
                "axis names are",
                axis_names,
                "variable names are",
                self.variable_names,
                "parameter names are",
                self.parameter_names,
            )
        )
        print( "all symbols are",
                all_symbols,
                "axis names are",
                axis_names,
                "variable names are",
                self.variable_names,
                "parameter names are",
                self.parameter_names,
            )
        self.symbolic_vars = all_symbols - axis_names
        self.fit_axis = list(self.fit_axis)[0]
        # }}}
        self.symbolic_vars = list(self.symbolic_vars)
        args = self.symbolic_vars + [str(*this_axis)]
        self.fitfunc_multiarg = sp.lambdify(
            args,
            self.expression,
            modules=[{"ImmutableMatrix": np.ndarray}, "numpy", "scipy"],
        )
        self.fitfunc_multiarg_v2 = sp.lambdify(
            variable_symbols + parameter_symbols,
            self.expression,
            modules=[{"ImmutableMatrix": np.ndarray}, "numpy", "scipy"],
        )
        def fn(p,x):
            p = self.add_inactive_p(p)
            assert len(p)==len(self.parameter_names),"length of parameter passed to fitfunc doesnt match number of symbolic parameters"
            return self.fitfunc_multiarg(*tuple(list(p) + [x]))
        self.fitfunc = fn
        self.pars = Parameters()
        for this_name in self.parameter_names:
            self.pars.add(this_name)
    def add_inactive_p(self,p):
        if self.set_indices is not None:
            #{{{uncollapse the function
            temp = p.copy()
            p = np.zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            #}}}
            p[self.set_indices] = self.set_to
        return p    
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
            logging.info(strm("fit param ---", j))
        logging.info(strm(self.pars))
        self.guess_dict = guesses
        return
    def guess(self):
        r'''Old code that we are preserving here -- provide the guess for our parameters; by default, based on pseudoinverse'''
        if hasattr(self,'guess_dict'):
            self.guess_dictionary = {k:self.guess_dict[k]['value'] for k in self.guess_dict.keys()}
            print(self.parameter_names)
            print(self.guess_dictionary)
            return [self.guess_dictionary[k] for k in self.parameter_names]
        else:
            return [1.0]*len(self.variable_names)
    def settoguess(self):
        'a debugging function, to easily plot the initial guess'
        self.fit_coeff = np.real(self.guess())
        return self
    def _taxis(self,taxis):
        r'You can enter None, to get the fit along the same range as the data, an integer to give the number of points, or a range of data, which will return 300 points'
        if taxis is None:
            taxis = self.getaxis(self.fit_axis).copy()
        elif isinstance(taxis, int):
            taxis = np.linspace(self.getaxis(self.fit_axis).min(),
                    self.getaxis(self.fit_axis).max(),
                    taxis)
        elif not np.isscalar(taxis) and len(taxis) == 2:
            taxis = np.linspace(taxis[0],taxis[1],300)
        return taxis    
    def eval(self,taxis,set_what = None, set_to = None):
        """Calculate the fit function along the axis taxis.

        Parameters
        ----------
        taxis: ndarray, int
            :if ndarray: the new axis coordinates along which we want to calculate the fit.
            :if int: number of evenly spaced points along the t-axis along the fit
        set_what: 'str', optional
            forcibly sets a specific symbol
        set_to: double, optional
            the specific value(int) you are assigning the symbol you included
        
        Returns
        -------
        self: nddata
            the fit function evaluated along the axis coordinates that were passed
        """
        if isinstance(set_what, dict):
            set_to = list(set_what.values())
            set_what = list(set_what.keys())
        taxis = self._taxis(taxis)
        if hasattr(self,'fit_coeff') and self.fit_coeff is not None:
            p = self.fit_coeff.copy()
        else:
            p = np.array([NaN]*len(self.variable_names))
        #{{{LOCALLY apply any forced values
        if set_what is not None:
            if self.set_indices is not None:
                raise ValueError("You're trying to set indices in an eval"
                        " function for a function that was fit constrained; this"
                        " is not currently supported")
            set_indices,set_to,active_mask = self.gen_indices(set_what,set_to)
            p[set_indices] = set_to
        #}}}
        #{{{ make a new blank np.array with the fit axis expanded to fit taxis
        newdata = ndshape(self)
        newdata[self.fit_axis] = np.size(taxis)
        newdata = newdata.alloc()
        newdata.set_plot_color(self.get_plot_color())
        #}}}
        #{{{keep all axis labels the same, except the expanded one
        newdata.axis_coords = list(newdata.axis_coords)
        newdata.labels([self.fit_axis],list([taxis]))
        #}}}
        newdata.data[:] = self.fitfunc(p,taxis).flatten()
        return newdata
    def thatresidual(self,p,x,y,sigma):
        fit = self.fitfunc(p,x)
        normalization = np.sum(1.0/sigma)
        sigma[sigma == 0.0] = 1
        try:
            # as noted here: https://stackoverflow.com/questions/6949370/scipy-leastsq-dfun-usage
            # this needs to be fit - y, not vice versa
            retval = (fit-y)/sigma #* normalization
        except ValueError as e:
            raise ValueError(strm('your error (',np.shape(sigma),
                    ') probably doesn\'t match y (',
                    np.shape(y),') and fit (',np.shape(fit),')')
                    + explain_error(e))
        return retval

    def fit(self, newfit, actualdata):
        r'''actually run the fit
        Parameters
        ==========
        thisfit: lmfitdata
                    original true data
        newfit:  lmfitdata
                    the lmfitdata instance of the actual data
        data:    nddata
                    data being fit
                    '''
        out = minimize(self.residual, newfit.pars, kws={"data":actualdata.data})
        fit = self.C
        fit.data = newfit.residual(out.params)
        report_fit(out,show_correl=True)
        self.fit_coeff =[]
        for name, param in out.params.items():
            self.fit_coeff.append(param.value)
        return fit

    def run_lambda(self,pars):
        """actually run the lambda function we separate this in case we want
        our function to involve something else, as well (e.g. taking a Fourier
        transform)"""
        print(self.getaxis(j) for j in self.variable_names)
        return self.fitfunc_multiarg_v2(*(self.getaxis(j) for j in self.variable_names), **pars.valuesdict())

    def residual(self, pars, data=None):
        "calculate the residual OR if data is None, return fake data"
        model = self.run_lambda(pars)
        if data is None:
            return model
        return model - data
    def copy(self): 
        namelist = []
        vallist = []
        for j in dir(self):
            if self._contains_symbolic(j):
                namelist.append(j)
                vallist.append(self.__getattribute__(j))
                self.__delattr__(j)
        new = deepcopy(self)
        for j in range(0,len(namelist)):
            new.__setattr__(namelist[j],vallist[j])
        for j in range(0,len(namelist)):
            self.__setattr__(namelist[j],vallist[j])
        return new    
    def gen_indices(self,this_set,set_to):
        r'''pass this this_set and this_set\_to parameters, and it will return:
        indices,values,mask
        indices --> gives the indices that are forced
        values --> the values they are forced to
        mask --> p[mask] are actually active in the fit'''
        if not isinstance(this_set, list):
            this_set = [this_set]
        if not isinstance(set_to, list):
            set_to = [set_to]
        if len(this_set) != len(set_to):
            raise ValueError(strm('length of this_set=',this_set,
                'and set_to', set_to, 'are not the same!'))
        logging.debug("*** *** *** *** *** ***")
        logging.debug(str(this_set))
        logging.debug("*** *** *** *** *** ***")
        set_indices = list(map(self.symbol_list.index,this_set))
        active_mask = np.ones(len(self.symbol_list),dtype=bool)
        active_mask[set_indices] = False
        return set_indices,set_to,active_mask
    def output(self, *name):
        r'''give the fit value of a particular symbol, or a dictionary of all values.

        Parameters
        -----------
        name: str (optional)
            name of the symbol.
            If no name is passed, then output returns a dictionary of the 
            resulting values.

        Returns
        -------
        retval: dict or float
            Either a dictionary of all the values, or the value itself
        '''
        if not hasattr(self,'fit_coeff') or self.fit_coeff is None:
            return None
        p = self.fit_coeff.copy()
        if self.set_indices is not None:
            temp = p.copy()
            p = np.zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            p[self.set_indices] = self.set_to
        if len(name) == 1:
            try:
                return p[self.symbol_list.index(name[0])]
            except:
                raise ValueError(strm("While running output: couldn't find",
                    name, "in", self.variable_names))
        elif len(name) == 0:
            return {self.variable_names[j]:p[j] for j in range(len(p))}
        else:
            raise ValueError(strm("You can't pass", len(name),"arguments to .output()"))
    def latex(self):
        r'''show the latex string for the function, with all the symbols substituted by their values'''
        # this should actually be generic to fitdata
        p = self.fit_coeff
        retval = self.function_string
        printfargs = []
        allsymb = []
        locations = []
        # {{{ I replace the symbols manually
        #     Note that I came back and tried to use sympy to do this,
        #     but then realize that sympy will automatically simplify,
        #     e.g. numbers in the denominator, so it ends up changing the
        #     way the function looks.  Though this is a pain, it's
        #     better.
        for j in range(0,len(self.symbol_list)):
            symbol = sympy_latex(self.symbolic_vars[j]).replace('$','')
            logger.debug(strm('DEBUG: replacing symbol "',symbol,'"'))
            location = retval.find(symbol)
            while location != -1:
                if retval[location-1] == '-':
                    newstring = retval[:location-1]+dp(-1*p[j])+retval[location+len(symbol):] # replace the symbol in the written function with the appropriate number
                else:
                    newstring = retval[:location]+dp(p[j])+retval[location+len(symbol):] # replace the symbol in the written function with the appropriate number
                logger.debug(strm(r"trying to replace",
                    retval[location:location+len(symbol)]))
                retval = newstring
                locations += [location]
                allsymb += [symbol]
                location = retval.find(symbol)
        # }}}
        logger.debug(strm(r"trying to generate",self.function_string,
            '\n',retval,'\n',[allsymb[x] for x in np.argsort(locations)],
            '\n',printfargs))
        return retval
    
    @property
    def function_string(self):
        r'''A property of the myfitclass class which stores a string
        output of the functional form of the desired fit expression
        provided in func:`functional_form` in LaTeX format'''
        retval = sympy_latex(self.expression).replace('$','')
        return r'$f(%s)='%(sympy_latex(sympy_symbol(self.fit_axis))) + retval + r'$'
    @function_string.setter
    def function_string(self):
        raise ValueError("You cannot set the string directly -- change the functional_form property instead!")
