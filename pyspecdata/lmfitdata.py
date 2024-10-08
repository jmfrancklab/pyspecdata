# just put this in the package
import sympy as sp
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import numpy as np
from .core import nddata, normal_attrs, issympy, ndshape, dp
from .general_functions import strm
import logging, warnings
from copy import deepcopy


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
                self.__setattr__(
                    myattrs[j], args[0].__getattribute__(myattrs[j])
                )
        else:
            nddata.__init__(self, *args, **kwargs)
        if fit_axis is None:
            if len(self.dimlabels) == 1:
                fit_axis = self.dimlabels[0]
            else:
                raise IndexError(
                    "Right now, we can only auto-determine the fit axis if"
                    " there is a single axis"
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
        logging.debug("Getting symbolic function")
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
        #     here, I discriminate "names" which are strings from "symbols"
        if self.expression is None:
            raise ValueError("what expression are you fitting with??")
        all_symbols = self.expression.atoms(sp.Symbol)
        all_symbol_names = set([str(j) for j in all_symbols])
        axis_names = set(self.dimlabels)
        variable_symbol_names = axis_names & all_symbol_names
        parameter_symbol_names = all_symbol_names - variable_symbol_names
        self.variable_names = tuple(variable_symbol_names)
        self.variable_symbols = [
            j for j in all_symbols if str(j) in variable_symbol_names
        ]
        self.parameter_symbols = [
            j for j in all_symbols if str(j) in parameter_symbol_names
        ]
        self.parameter_names = tuple([str(j) for j in self.parameter_symbols])
        self.fit_axis = set(self.dimlabels)
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
        self.fit_axis = list(self.fit_axis)[0]
        # }}}
        self.fitfunc_multiarg_v2 = sp.lambdify(
            self.variable_symbols + self.parameter_symbols,
            self.expression,
            modules=[{"ImmutableMatrix": np.ndarray}, "numpy", "scipy"],
        )

        self.guess_parameters = Parameters()
        for this_name in self.parameter_names:
            self.guess_parameters.add(this_name)

    def set_guess(self, *args, **kwargs):
        """set both the guess and the bounds

        Parameters
        ==========
        guesses: dict of dicts
            each dict has a keyword giving the parameter and a value
            that comprises a dict with guesses (value) and/or constraints
            (min/max)

            Can be passed either as the only argument, or a kwarg called
            guesses, or as the kwargs themselves (i.e. `.set_guesses(**guesses)`)

            For example (the last case:)

            >>> d.set_guess(
            >>>     param1=dict(value=1.0, min=0, max=10),
            >>>     param2=dict(value=2.0, min=0, max=10))
        """
        if len(args) == 1 and type(args[0]) == dict:
            guesses = args[0]
        elif len(kwargs) == 1 and "guesses" in kwargs.keys():
            guesses = kwargs["guesses"]
        else:
            guesses = kwargs
        self.guess_dict = {}
        logging.debug(strm("guesses", guesses))
        logging.debug(strm("guess_parameters", self.guess_parameters.keys()))
        for this_name in self.guess_parameters.keys():
            if this_name in guesses.keys():
                logging.debug(strm("adding", this_name))
                if type(guesses[this_name]) is dict:
                    self.guess_dict[this_name] = {}
                    for k, v in guesses[this_name].items():
                        setattr(self.guess_parameters[this_name], k, v)
                        self.guess_dict[this_name][k] = v
                elif np.isscalar(guesses[this_name]):
                    self.guess_parameters[this_name].value = guesses[this_name]
                    self.guess_dict[this_name] = {"value": guesses[this_name]}
                else:
                    raise ValueError("what are the keys to your guesses???")
                logging.debug(strm("now dict is", self.guess_dict))
        for j in self.guess_parameters:
            logging.debug(strm("fit param ---", j))
        logging.debug(strm(self.guess_parameters))
        logging.debug(strm(self.guess_dict))
        return self

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

    def eval(self, taxis=None):
        """Calculate the fit function along the axis taxis.

        Parameters
        ----------
        taxis: ndarray, int
            :if ndarray: the new axis coordinates along which we want to calculate the fit.
            :if int: number of evenly spaced points along the t-axis along the fit

        Returns
        -------
        self: nddata
            the fit function evaluated along the axis coordinates that were passed
        """
        if taxis is None:
            taxis = self.getaxis(self.fit_axis)
        else:
            taxis = self._taxis(taxis)
        if hasattr(self, "fit_coeff") and self.fit_coeff is not None:
            p = self.fit_coeff.copy()
            # here you see that fit_coeff stores the coefficients that
            # were previously fit, and these are stored in p, here
        else:
            p = np.array([NaN] * len(self.variable_names))
        # JF notes this is a copy of older code -- we should be able to
        # clean this up by using the newer copy functions -- currently I
        # see  the copy_props and copyaxes functions, but thought there
        # was a function to copy everything BUT data -- maybe this is on
        # the SVD branch or somesuch?
        # {{{ make a new blank np.array with the fit axis expanded to fit taxis
        newdata = ndshape(self)
        newdata[self.fit_axis] = np.size(taxis)
        newdata = newdata.alloc()
        newdata.set_plot_color(self.get_plot_color())
        # }}}
        # {{{keep all axis labels the same, except the expanded one
        newdata.axis_coords = list(newdata.axis_coords)
        newdata.labels([self.fit_axis], list([taxis]))
        if self.get_units(self.fit_axis) is not None:
            newdata.set_units(self.fit_axis, self.get_units(self.fit_axis))
        if self.get_error(self.fit_axis) is not None:
            newdata.set_error(self.fit_axis, self.get_error(self.fit_axis))
        # }}}
        variable_coords = {
            self.fit_axis: taxis
        }  # even though it's possible to
        #                                         combine this and the next line
        #                                         to make it more
        #                                         compact/efficient for one
        #                                         variable, we want to leave
        #                                         open the posisbility that we
        #                                         will be using more than one
        #                                         variable
        newdata.data[:] = self.fitfunc_multiarg_v2(*(
            tuple(
                (
                    variable_coords[j]
                    if j in variable_coords.keys()
                    else self.getaxis(j)
                )
                for j in self.variable_names
            )
            + tuple(self.fit_coeff)
        )).flatten()
        newdata.name(str(self.name()))
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
        logging.debug(
            strm(
                "I am attempting to actually run the fit",
                minimize.__module__,
                "here are the guess parameters",
                self.guess_parameters,
            )
        )
        x = self.getaxis(self.fit_axis)
        y = self.data
        sigma = self.get_error()
        out = minimize(
            self.residual,
            self.guess_parameters,
            args=(x, y, sigma),
        )
        logging.debug(
            strm(
                "here is the output success",
                out.success,
                "and parameters",
                out.params,
            )
        )
        # {{{ capture the result for ouput, etc
        self.fit_coeff = [out.params[j].value for j in self.parameter_names]
        assert out.success
        if hasattr(out, "covar"):
            self.covariance = out.covar
        # }}}
        return self

    def residual(self, pars, x, y, sigma=None):
        "calculate the residual OR if data is None, return fake data"
        fit = self.fitfunc_multiarg_v2(
            *(self.getaxis(j) for j in self.variable_names),
            **pars.valuesdict(),
        )
        if sigma is not None:
            normalization = np.sum(
                1.0 / sigma[np.logical_and(sigma != 0.0, np.isfinite(sigma))]
            )
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
        return retval.view(float)

    def copy(self, **kwargs):
        namelist = []
        vallist = []
        for j in dir(self):
            if self._contains_symbolic(j):
                namelist.append(j)
                vallist.append(self.__getattribute__(j))
                self.__delattr__(j)
        new = super().copy(**kwargs)
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
        set_indices = list(map(self.parameter_names.index, this_set))
        active_mask = np.ones(len(self.parameter_names), dtype=bool)
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
            p = np.zeros(len(self.parameter_names))
            p[self.active_mask] = temp
            p[self.set_indices] = self.set_to
        if len(name) == 1:
            try:
                return p[self.parameter_names.index(name[0])]
            except:
                raise ValueError(
                    strm(
                        "While running output: couldn't find",
                        name,
                        "in",
                        self.parameter_names,
                    )
                )
        elif len(name) == 0:
            return {self.parameter_names[j]: p[j] for j in range(len(p))}
        else:
            raise ValueError(
                strm("You can't pass", len(name), "arguments to .output()")
            )

    def latex(self):
        r"""show the latex string for the function, with all the symbols substituted by their values"""
        # this should actually be generic to fitdata
        thisdp = lambda x: "(" + dp(x) + ")"
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
        for j in range(0, len(self.parameter_names)):
            symbol = sp.printing.latex(self.parameter_symbols[j]).replace(
                "$", ""
            )
            logging.debug(strm('DEBUG: replacing symbol "', symbol, '"'))
            location = retval.find(symbol)
            while location != -1:
                if retval[location - 1] == "-":
                    newstring = (
                        retval[: location - 1]
                        + thisdp(-1 * p[j])
                        + retval[location + len(symbol) :]
                    )  # replace the symbol in the written function with the appropriate number
                else:
                    newstring = (
                        retval[:location]
                        + thisdp(p[j])
                        + retval[location + len(symbol) :]
                    )  # replace the symbol in the written function with the appropriate number
                logging.debug(
                    strm(
                        r"trying to replace",
                        retval[location : location + len(symbol)],
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
        retval = sp.printing.latex(self.expression).replace("$", "")
        return (
            r"$f(%s)=" % (sp.printing.latex(sp.core.Symbol(self.fit_axis)))
            + retval
            + r"$"
        )

    @function_string.setter
    def function_string(self):
        raise ValueError(
            "You cannot set the string directly -- change the functional_form"
            " property instead!"
        )
