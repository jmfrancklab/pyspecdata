# just put this in the package
import sympy as sp
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import numpy as np
from .core import nddata, normal_attrs, issympy, ndshape, sympy_latex, sympy_symbol, dp
from .general_functions import strm
import logging, warnings
from copy import deepcopy

class lmfitglobal(nddata):
    r"""Based on lmfitdata class, new class which
    will perform global fit on multitude of
    lmfitdata classes
    """
    def __init__(self, *args, **kwargs):
        self.global_expr = list(args[0].values())[0]
        self.global_func = list(args[0].keys())[0]

        self.datasets = []
        self.var_list = []

        self.global_vars_name = []
        self.global_vars_value = []
        self.local_vars_list = []
        self.global_func_list = []
        self.local_params_list = []
        self.global_params_list = []

        self.translation_list = []
        return

    def append(self, dataset, dataset_dict):
        # datasets is a list of each dataset
        self.datasets.append(dataset)

        temp_vars = dataset.variable_names
        temp_params = dataset.parameter_names
        print("Dataset's variables are",temp_vars)
        print("Dataset's parameters are",temp_params)


        if dataset_dict.keys() in self.global_vars_name:
            None
        else:    
            self.global_vars_name.append(list(dataset_dict.keys())[0])
        self.global_vars_value.append(list(dataset_dict.values())[0])
        for i,j in enumerate(temp_vars):
            if j in self.local_vars_list:
                None
            else:
                self.local_vars_list.append(temp_vars[i])
        for i,j in enumerate(temp_params):
            if j in self.global_func:
                    if j in self.global_func_list:
                        None
                    else:
                        self.global_func_list.append(j)
            else:
                    if j in self.local_params_list:
                        None
                    else:
                        self.local_params_list.append(j)

        print("Global vars",self.global_vars_name)
        print("Local vars",self.local_vars_list)
        print("Global func",self.global_func_list)
        print("Local params",self.local_params_list)
        
        # the order of the parameters used in run_lambda
        self.pars_order = list(dataset.pars.keys())

        # create elements of translation list
        translation_element = []
        for i,j in enumerate(self.global_func_list):
            translation_element.append( ('g',str(j))  )
        for i,j in enumerate(self.local_params_list):
            element_string = str(j+'_'+str(len(self.datasets)))
            if j in self.pars_order:
                element_position = self.pars_order.index(j)
            translation_element.append( ('l', element_string, str(j) ))

        ## add elements to translation list
        self.translation_list.append(translation_element)
        return

    def make_params(self):
        self.pars = Parameters()
        for i,j in enumerate(self.datasets):
            for k in (self.datasets[i].pars.keys()):
                for l in self.local_params_list:
                    if l == k:
                        temp = self.translation_list[i]
                        for temp_element in temp:
                            if temp_element[0] == 'g':
                                None
                            else:
                                temp_name = temp_element[1][:-2]
                                if temp_name == k:
                                    #new_elem = self.datasets[i].pars[k].name = temp_element[1]
                                    self.datasets[i].pars[k].name = temp_element[1]
                                    self.pars.add(deepcopy(self.datasets[i].pars[k]))
        return

    def inspect_model(self):
        # from functional_form.setter
        # {{{ decide which symbols are parameters vs. variables
        if self.global_expr is None:
            raise ValueError("Missing global expression.")
        all_symbols = self.global_expr.atoms(sp.Symbol)
        axis_names = set([sp.Symbol(j, real=True) for j in self.global_vars_name])
        variable_symbols = axis_names & all_symbols
        self.parameter_symbols = all_symbols - variable_symbols
        this_axis = variable_symbols
        variable_symbols = tuple(variable_symbols)
        self.variable_names = tuple([str(j) for j in variable_symbols])
        parameter_symbols = tuple(self.parameter_symbols)
        self.parameter_names = tuple([str(j) for j in self.parameter_symbols])
        self.fit_axis = set(self.global_vars_name)
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
        print(
            "all symbols are",
            all_symbols,
            "axis names are",
            axis_names,
            "variable names are",
            self.variable_names,
            "parameter names are",
            self.parameter_names,
        )
        self.global_params_list = self.parameter_names
        self.symbolic_vars = all_symbols - axis_names
        self.fit_axis = list(self.fit_axis)[0]

        self.symbolic_vars = list(self.symbolic_vars)
        args = self.symbolic_vars + [str(*this_axis)]
        print(variable_symbols)
        print(parameter_symbols)
        print(self.global_expr)

        self.fitfunc_multiarg_v2 = sp.lambdify(
            variable_symbols + parameter_symbols,
            self.global_expr,
            modules=[{"ImmutableMatrix": np.ndarray}, "numpy", "scipy"],
        )
        return
         #}}}
    def run_lambda(self, pars, variable_axis):
        """actually run the lambda function we separate this in case we want
        our function to involve something else, as well (e.g. taking a Fourier
        transform)"""
        self.setaxis = variable_axis
        return self.fitfunc_multiarg_v2(
            *(variable_axis.data,), **pars.valuesdict()
        )
    def make_model(self, pars, x, sigma=None):
        model_params = Parameters()
        temp_list = []
        for i,j in enumerate(self.pars.valuesdict()):
            for q in self.translation_list:
                for r in q:
                    if j == r[1]:
                        temp_list.append(j)
                    else:
                        None
        for i,j in enumerate(self.pars.valuesdict()):
            if j in temp_list:
                None
            else:
                model_params.add(self.pars[j])
        fit = self.run_lambda(model_params,x)
        return fit

    def member_model(self, member_idx, member_model_input):
        return self.datasets[member_idx].make_model(member_model_input)

    def residual(self, parameters):
        this_model = self.make_model(parameters,
                nddata(np.array(self.global_vars_value),str(self.fit_axis)))
        ndata = len(self.datasets)
        resid = 0.0*np.array(self.datasets[:])
        for i in range(ndata):
            resid[i] = (self.member_model(i, {'R1':this_model[i]}))[0] - self.datasets[i].data
        retval = np.concatenate(resid)
        retval = retval.real
        print(type(retval[0]))
        return retval

    def fit(self):
        x = self.global_vars_value
        this_model = self.make_model(self.pars,
                nddata(np.array(self.global_vars_value),str(self.fit_axis)))
        y = np.real(this_model.data)
        print(type(self.residual))
        print(type(self.pars))
        out = minimize(
                self.residual,
                self.pars,
                )
        return out

class lmfitdata(nddata):
    r"""Inherits from an nddata and enables curve fitting through use of a sympy expression.

    The user creates a lmfitdata class object from an existing nddata
    class object, and on this lmfitdata object can define the
    :func:`functional_form` of the curve it would like to fit to the
    data of the original nddata.
    This functional form must be provided as a sympy expression, with
    one of its variables matching the name of the dimension that the
    user would like to fit to.
    """

    def __init__(self, *args, **kwargs):
        # copied from fitdata
        fit_axis = None
        if "fit_axis" in list(kwargs.keys()):
            fit_axis = kwargs.pop("fit_axis")
        if isinstance(args[0], nddata):
            # move nddata attributes into the current instance
            myattrs = normal_attrs(args[0])
            for j in range(0, len(myattrs)):
                self.__setattr__(myattrs[j], args[0].__getattribute__(myattrs[j]))
        else:
            nddata.__init__(self, *args, **kwargs)
        if fit_axis is None:
            if len(self.dimlabels) == 1:
                fit_axis = self.dimlabels[0]
            else:
                raise IndexError(
                    "Right now, we can only auto-determine the fit axis if there is a single axis"
                )

        self.fit_axis = fit_axis
        self.set_to = None
        self.set_indices = None
        self.active_indices = None
        self.expression = None
        return

    @property
    def functional_form(self):
        r"""A property of the myfitclass class which is set by the user,
        takes as input a sympy expression of the desired fit
        expression"""
        print("Getting symbolic function")
        return self.expression

    @functional_form.setter
    def functional_form(self, this_expr):
        """generate parameter descriptions and a numpy (lambda) function from a sympy expresssion

        Parameters
        ==========
        this_expr: sympy expression
        """
        assert issympy(
            this_expr
        ), "for now, the functional form must be a sympy expression"
        self.expression = this_expr
        # {{{ decide which symbols are parameters vs. variables
        if self.expression is None:
            raise ValueError("what expression are you fitting with??")
        all_symbols = self.expression.atoms(sp.Symbol)
        axis_names = set([sp.Symbol(j, real=True) for j in self.dimlabels])
        variable_symbols = axis_names & all_symbols
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
        print(
            "all symbols are",
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
        #print(variable_symbols)
        #print(parameter_symbols)
        #print(self.expression)
        #quit()
        self.fitfunc_multiarg_v2 = sp.lambdify(
            variable_symbols + parameter_symbols,
            self.expression,
            modules=[{"ImmutableMatrix": np.ndarray}, "numpy", "scipy"],
        )

        def fn(p, x):
            p = self.add_inactive_p(p)
            assert len(p) == len(
                self.parameter_names
            ), "length of parameter passed to fitfunc doesnt match number of symbolic parameters"
            return self.fitfunc_multiarg(*tuple(list(p) + [x]))

        self.fitfunc = fn
        self.pars = Parameters()
        for this_name in self.parameter_names:
            self.pars.add(this_name)

    def add_inactive_p(self, p):
        if self.set_indices is not None:
            # {{{uncollapse the function
            temp = p.copy()
            p = np.zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            # }}}
            p[self.set_indices] = self.set_to
        return p

    def set_guess(self, *args, **kwargs):
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
        elif len(kwargs) == 1 and "guesses" in kwargs.keys():
            guesses = kwargs["guesses"]
        else:
            guesses = kwargs
        self.guess_dict = {}
        for this_name in self.pars.keys():
            if this_name in guesses.keys():
                if type(guesses[this_name]) is dict:
                    self.guess_dict[this_name] = {}
                    for k, v in guesses[this_name].items():
                        setattr(self.pars[this_name], k, v)
                        self.guess_dict[this_name][k] = v
                elif np.isscalar(guesses[this_name]):
                    self.pars[this_name].value = guesses[this_name]
                    self.guess_dict[this_name] = {"value":guesses[this_name]}
                else:
                    raise ValueError("what are the keys to your guesses???")
        for j in self.pars:
            logging.info(strm("fit param ---", j))
        logging.info(strm(self.pars))
        return

    def guess(self):
        r"""Old code that we are preserving here -- provide the guess for our
        parameters; by default, based on pseudoinverse"""
        if hasattr(self, "guess_dict"):
            self.guess_dictionary = {
                k: self.guess_dict[k]["value"] for k in self.guess_dict.keys()
            }
            return [self.guess_dictionary[k] for k in self.parameter_names]
        else:
            return [1.0] * len(self.variable_names)

    def settoguess(self):
        "a debugging function, to easily plot the initial guess"
        self.fit_coeff = np.real(self.guess())
        return self

    def _taxis(self, taxis):
        r"You can enter None, to get the fit along the same range as the data, an integer to give the number of points, or a range of data, which will return 300 points"
        if taxis is None:
            taxis = self.getaxis(self.fit_axis).copy()
        elif isinstance(taxis, int):
            taxis = np.linspace(
                self.getaxis(self.fit_axis).min(),
                self.getaxis(self.fit_axis).max(),
                taxis,
            )
        elif not np.isscalar(taxis) and len(taxis) == 2:
            taxis = np.linspace(taxis[0], taxis[1], 300)
        return taxis

    def eval(self, taxis=None, set_what=None, set_to=None):
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
        if taxis is None:
            taxis = self.getaxis(self.fit_axis)
        else:
            taxis = self._taxis(taxis)
        if hasattr(self, "fit_coeff") and self.fit_coeff is not None:
            p = self.fit_coeff.copy()
        else:
            p = np.array([NaN] * len(self.variable_names))
        # {{{LOCALLY apply any forced values
        if set_what is not None:
            if self.set_indices is not None:
                raise ValueError(
                    "You're trying to set indices in an eval"
                    " function for a function that was fit constrained; this"
                    " is not currently supported"
                )
            set_indices, set_to, active_mask = self.gen_indices(set_what, set_to)
            p[set_indices] = set_to
        # }}}
        # {{{ make a new blank np.array with the fit axis expanded to fit taxis
        newdata = ndshape(self)
        newdata[self.fit_axis] = np.size(taxis)
        newdata = newdata.alloc()
        newdata.set_plot_color(self.get_plot_color())
        # }}}
        # {{{keep all axis labels the same, except the expanded one
        newdata.axis_coords = list(newdata.axis_coords)
        newdata.labels([self.fit_axis], list([taxis]))
        # }}}
        newdata.data[:] = self.fitfunc(p, taxis).flatten()
        return newdata
    def fit(self):
        r"""actually run the fit"""
        # we can ignore set_what, since I think there's a mechanism in
        # lmfit to take care of that (it's for fixing parameters)
        # but the rest of what it's doing is to pull apart the
        # error, axis, etc, to be fed to minimize.
        #
        # It also automatically converts complex data to real data, and
        # does other things for error handling -- let's not just throw this out
        #
        # I think that a lot of this could be copied with little modification
        #
        # But you  should read through and see what the previous fit method is doing
        # and then copy over what you can
        x = self.getaxis(self.fit_axis)
        if np.iscomplex(self.data.flatten()[0]):
            logging.debug(strm("Warning, taking only real part of fitting data!"))
        y = np.real(self.data)
        sigma = self.get_error()
        out = minimize(
            self.residual,
            self.pars,
            args=(x, y, sigma),
        )
        # can you capture the following as a string? maybe return it?
        report_fit(out, show_correl=True)
        # {{{ capture the result for ouput, etc
        self.fit_coeff = [out.params[j].value for j in self.symbol_list]
        assert out.success
        self.covariance = out.covar
        # }}}
        return

    def run_lambda(self, pars):
        """actually run the lambda function we separate this in case we want
        our function to involve something else, as well (e.g. taking a Fourier
        transform)"""
        logging.info(strm(self.getaxis(j) for j in self.variable_names))
        print(self.getaxis(self.variable_names[0]))
        print(pars.valuesdict())
        return self.fitfunc_multiarg_v2(
            *(self.getaxis(j) for j in self.variable_names), **pars.valuesdict()
        )
    def make_model(self,model_input):
        """"called by member_model in global lmfitdata; provide values of each
        parameter and produce model output - right now just able to use one
        parameter (passed in :dict:`model_input`)"""
        if len(model_input) == 1:
            model_input_param = str(list(model_input.keys())[0])
            assert model_input_param in self.parameter_names,("parameter set in model_input is not found in this dataset's list of parameters; check the parameter names")
            self.pars[model_input_param].value = list(model_input.values())[0]
            # deal with the fact that Parameter name does not match
            # Parameters dictionary key for the member dataset
            for i,j in enumerate(self.pars):
                self.pars[j].name = j
            return self.run_lambda(self.pars),self.getaxis(self.variable_names[0])
        else:
            raise ValueError("only able to set one parameter with model_input at this time")
            return

    def residual(self, pars, x, y, sigma=None):
        "calculate the residual OR if data is None, return fake data"
        fit = self.run_lambda(pars)
        if sigma is not None:
            normalization = np.sum(1.0 / sigma[sigma != 0.0 and np.isfinite(sigma)])
            sigma[sigma == 0.0] = 1
            sigma[~np.isfinite(sigma)] = 1
        try:
            # as noted here: https://stackoverflow.com/questions/6949370/scipy-leastsq-dfun-usage
            # this needs to be fit - y, not vice versa
            if sigma is not None:
                retval = (fit - y) / sigma * normalization
            else:
                retval = fit - y
        except ValueError as e:
            raise ValueError(
                strm(
                    "your error (",
                    np.shape(sigma),
                    ") probably doesn't match y (",
                    np.shape(y),
                    ") and fit (",
                    np.shape(fit),
                    ")",
                )
                + explain_error(e)
            )
        return retval

    def copy(self):
        namelist = []
        vallist = []
        for j in dir(self):
            if self._contains_symbolic(j):
                namelist.append(j)
                vallist.append(self.__getattribute__(j))
                self.__delattr__(j)
        new = deepcopy(self)
        for j in range(0, len(namelist)):
            new.__setattr__(namelist[j], vallist[j])
        for j in range(0, len(namelist)):
            self.__setattr__(namelist[j], vallist[j])
        return new

    def gen_indices(self, this_set, set_to):
        r"""pass this this_set and this_set\_to parameters, and it will return:
        indices,values,mask
        indices --> gives the indices that are forced
        values --> the values they are forced to
        mask --> p[mask] are actually active in the fit"""
        if not isinstance(this_set, list):
            this_set = [this_set]
        if not isinstance(set_to, list):
            set_to = [set_to]
        if len(this_set) != len(set_to):
            raise ValueError(
                strm(
                    "length of this_set=",
                    this_set,
                    "and set_to",
                    set_to,
                    "are not the same!",
                )
            )
        logging.debug("*** *** *** *** *** ***")
        logging.debug(str(this_set))
        logging.debug("*** *** *** *** *** ***")
        set_indices = list(map(self.symbol_list.index, this_set))
        active_mask = np.ones(len(self.symbol_list), dtype=bool)
        active_mask[set_indices] = False
        return set_indices, set_to, active_mask

    def output(self, *name):
        r"""give the fit value of a particular symbol, or a dictionary of all values.

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
        """
        if not hasattr(self, "fit_coeff") or self.fit_coeff is None:
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
                raise ValueError(
                    strm(
                        "While running output: couldn't find",
                        name,
                        "in",
                        self.symbol_list,
                    )
                )
        elif len(name) == 0:
            return {self.symbol_list[j]: p[j] for j in range(len(p))}
        else:
            raise ValueError(
                strm("You can't pass", len(name), "arguments to .output()")
            )

    def latex(self):
        r"""show the latex string for the function, with all the symbols substituted by their values"""
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
        for j in range(0, len(self.symbol_list)):
            symbol = sympy_latex(self.symbolic_vars[j]).replace("$", "")
            logging.debug(strm('DEBUG: replacing symbol "', symbol, '"'))
            location = retval.find(symbol)
            while location != -1:
                if retval[location - 1] == "-":
                    newstring = (
                        retval[: location - 1]
                        + dp(-1 * p[j])
                        + retval[location + len(symbol) :]
                    )  # replace the symbol in the written function with the appropriate number
                else:
                    newstring = (
                        retval[:location] + dp(p[j]) + retval[location + len(symbol) :]
                    )  # replace the symbol in the written function with the appropriate number
                logging.debug(
                    strm(
                        r"trying to replace", retval[location : location + len(symbol)]
                    )
                )
                retval = newstring
                locations += [location]
                allsymb += [symbol]
                location = retval.find(symbol)
        # }}}
        logging.debug(
            strm(
                r"trying to generate",
                self.function_string,
                "\n",
                retval,
                "\n",
                [allsymb[x] for x in np.argsort(locations)],
                "\n",
                printfargs,
            )
        )
        return retval

    @property
    def function_string(self):
        r"""A property of the myfitclass class which stores a string
        output of the functional form of the desired fit expression
        provided in func:`functional_form` in LaTeX format"""
        retval = sympy_latex(self.expression).replace("$", "")
        return r"$f(%s)=" % (sympy_latex(sympy_symbol(self.fit_axis))) + retval + r"$"

    @function_string.setter
    def function_string(self):
        raise ValueError(
            "You cannot set the string directly -- change the functional_form property instead!"
        )

