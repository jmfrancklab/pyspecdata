fitdata Methods to Override
---------------------------

Mandatory
~~~~~~~~~

For easiest (thought not necessarily most efficient) usage, one should
take advantage of the symbolic math package (sympy), which can
automatically calculate the gradients based on the functional form.

.. todo::
    actually use the gradients

The calculation of these gradients makes some tedious programming – such
as generating an initial guess to be unnecessary. Therefore, one need
override only the following functions:

\__init__(self,*args,**kwargs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. currentmodule:: pyspecdata.fitdata
.. automethod:: __init__

One begins by defining a routine to initialize the class with the
correct variable names. This is done with a line of code like the
following, where one copies the code exactly, changing only the values
of ``symbol_list`` and the argument of ``gen_symbolic``, which give –
respectively – the names of the fit parameters (in the order used below)
and the name (i.e. :math:`y-`\ value) of the function.

::

   def __init__(self,*args,**kwargs):
       '''here, we give the particular latex representation and list of symbols for this particular child class'''
       fitdata.__init__(self,*args,**kwargs)
       self.symbol_list = [r'M(\infty)',r'M(0)',r'T_1']
       self.starting_guesses = map(double,[r_[1,1,1],r_[0,0,1],r_[-100,100,0.03],r_[0.001,0.001,0.001],r_[1,-1,4.0]])# a series of starting guesses used by the automatic guess routine
       self.guess_lb = r_[-inf,-inf,1e-4]# a lower bound applied when running pseudoinverses to generate the starting guesses
       self.guess_ub = r_[+inf,+inf,20.]# an upper bound for the same
       self.gen_symbolic(r'M(t)')
       return

Unfortunately, sympy imposes somewhat stringent restrictions on the
parameter names; while the parameters can be named as words or with a
word subscript, parameters named with multiple symbols or subscripts do
not appear to work correctly. In addition, none of the later parameter
names can contain one of the earlier parameter names as a substring.
Therefore, if unexpected errors occur, we recommend switching to simple
(i.e. single letter) parameter names.

If one is *not* using symbolic math (not recommended), one define the
attribute ``self.function_string`` and ``self.function_name`` by hand
(these are strings that give the functional format and name of the y
values, respectively).

If you want to see an example that allows “multiplicity” – i.e. a
biexponential, rather than an exponential – see the ``t2curve`` class in

fitfunc_raw(self,p,x) and fitfunc_raw_symb(self,p,x)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. todo::
    This should be changed.
    We should not have to overload these methods.

    Rather, we should
    use a property, which we set to a sympy expression,
    and let the setter generate whatever we need:

    We should generate the symbol list by using this:
    `this <https://stackoverflow.com/questions/30018977/how-can-i-get-a-list-of-the-symbols-in-a-sympy-expression>`__
    and we should determine the axis being fit by comparing the variables to dimlabels.

    If we actually need `fitfunc_raw_symb` in the form used before,
    we can generate `fitfunc_raw_symb` in a straightforward way.

    We can generate `fitfunc_raw` from `lambdify`.

    We can generate the gradient as part of teh setter.


These defines the functional form of the fit they are defined in terms
of the parameter vector, :math:`p` (this is a list of fit parameters
that are named in the in the \__init__() method), and the data vector
:math:`x` (which is the data along the fit dimension); fitfunc_raw uses
nddata-compatible functions, which have obvious names

::

   from numpy import *
   # this code written for a module with the above declaration
   def fitfunc_raw(self,p,x):
       '''just the actual fit function to
       return the array y as a function of p
       and x'''
       return p[0]+(p[1]-p[0])*exp(-x/p[2])

one must also generate a **mathematically identical** function that
rather uses functions from the sympy package. For instance, one must use
an “exp()” function here than can operate on symbolic variables to
generate an analytical expression.

::

   import sympy
   # this code written for a module with the above
   #declaration
   def fitfunc_raw_symb(self,p,x):
       '''if I'm using a named function, I have
       to define separately in terms of sympy
       rather than numpy functions'''
       return p[0]+(p[1]-p[0])*sympy.exp(-x/p[2])

It is **highly recommended** that after writing a new class, one first
checks that the two functions above are mathematically identical (by
inspection), and then checks that the parameter indices here line up in
the expected way with the parameter names given in the ``__init__()``
method with the ``function_string`` method, as used in the example below
(`1.2.7 <#sec:writeup_software_fitfunc_thingstouse_example>`__)

Non-mandatory
~~~~~~~~~~~~~

guess(self)
^^^^^^^^^^^

If desired, this function makes an initial guess for the parameter
vectors. If the fit does not use symbolic algebra, this step is
mandatory. For instance, one can guess the parameters for the
:math:`T_1` example here based on the initial slope and values near the
end of the recovery curve as follows:

::

   def guess(self):
       r'''provide the guess for our parameters, which is specific to the type of function'''
       x = self.getaxis(self.fit_axis)
       y = self.data
       testpoint = argmin(abs(x-x.max()/3)) # don't just pull 1/3 of the index, because it can be unevenly spaced
       initial_slope = (y[testpoint]-y[0])/(x[testpoint]-x[0])
       A = y[-1]
       B = y[testpoint]-x[testpoint]*initial_slope
       C = (A-B)/initial_slope
       if (C < 0):
           raise CustomError(maprep('Negative T1!!! A-B=',A-B,'initial_slope=',initial_slope,x,y))
       oldguess = r_[A,B,C/2.0] # guesses for the parameters, in the same order
       return oldguess

However, we note that the symbolic algebra version does work quite
consistently. It employs several steps with the “regularized
pseudoinverse” routine (i.e. Tikhonov regularized solution via SVD, see
`[sec:pinvr] <#sec:pinvr>`__). 

.. todo::
    While this seems to work very well,
    it is unclear why, since Levenberg-Marquardt, the algorithm at the core
    of scipy’s nonlinear fitting procedure, consists of a series of steps,
    many of which are mathematically identical to a regularized
    pseudoinverse. Maybe this is only when we are using the numerical
    derivative, and it will now happen when we use the analytical
    derivative.

linfunc(self,x,y,xerr = None,yerr = None)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the event that we want to use a “linearized format” of the fit
function, we can use “linfunc” to return this format. This routine
(which is not designed to be used directly), takes the inputs :math:`x`,
which are the values (i.e. labels) of the fit axis and :math:`y`, which
are the data values at the points given by :math:`x`.

For instance, in the case of a :math:`T_1` recovery curve, we may want
to check the fit by plotting the value :math:`\ln(M(t)-M(\infty))` as a
function of :math:`t`, which should be linear, since

.. math::

   \begin{aligned}
           \ln(M(t)-M(\infty)) = \ln(M(0)-M(\infty)) - \frac{t}{T_1}
           .
       \end{aligned}

\ This is coded as follows:

::

   def linfunc(self,x,y,xerr = None,yerr = None):
       '''just the actual fit function to return the pair of arrays x',y' that should be linear
       it accepts as inputs x and y, and it uses the output from the fit, where necessary
       also optionally propagates the error based on yerr and xerr, which can be passed in to it
       For the case of T1, we want to return ln(y-M(\infty)) = ln(M(0)-M(\infty)) - t/T_1
       '''
       temp = self.output(r'M(\infty)')-y # the argument for log
       rety = log(temp)
       if yerr != None:
           reterr = yerr/abs(temp)
       mask = isfinite(rety)
       retx = x # for instance, in Emax, this is not just x
       xname = self.fit_axis # same as the fit axis
       yname = r'$ln(M(\infty)-M(t))$'
       #{{{ this should be pretty standardized
       retval = nddata(rety,
               [size(rety),1],
               [xname,yname])
       retval.labels([self.fit_axis],
               [retx.copy()])
       if yerr != None:
           retval.set_error(reterr)
       #}}}
       return retval

Methods and attributes to Use
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because ``fitdata`` inherits from ``nddata``, all of the standard nddata
methods are available. In addition, the following methods are available.

Available before fit
~~~~~~~~~~~~~~~~~~~~

These are the functions available before the ``fit`` routine is called.

instance.function_string
^^^^^^^^^^^^^^^^^^^^^^^^

instance.fit(…)
^^^^^^^^^^^^^^^

.. _instance.fit-1:

instance.fit()
''''''''''''''

For various reasons, it is best to separate the actual fitting step from
initialization routine (i.e. function called when we create a new
instance). This actually fits the data to the curve format specified by
the particular class (i.e. ``t1curve, ksp``, etc.).

instance.fit(set = {‘p1’:1.0,‘p3’:2.0})
'''''''''''''''''''''''''''''''''''''''

This example will constrain parameter :math:`p1` to 1.0 and parameter
:math:`p3` to 2.0, and fit the remaining parameters. One can replace
“{‘p1’:1.0,‘p3’:2.0}” with any dictionary, where the keys must be the
names of fit parameters for this class.

instance.fit(set = [‘p1’,‘p3’], set_to = [1.0,2.0])
'''''''''''''''''''''''''''''''''''''''''''''''''''

This (older format) does the same thing as the previous example.

instance.guess()
^^^^^^^^^^^^^^^^

This evaluates the initial guess along the fit axis.

.. todo::
    how is the fit axis determined, and is it possible to fit if the
    data has more than one dimension?

instance.settoguess()
^^^^^^^^^^^^^^^^^^^^^

This function is for **debugging purposes only**. This works similar to
``instance.fit(…)``, except that it sets the “fit result” to the initial
guess, and does not take any fixed parameters.

After fitting
~~~~~~~~~~~~~

These are the functions available after the ``fit`` routine is called;
these are not supplied in order, but rather order of importance.

instance.latex()
^^^^^^^^^^^^^^^^

instance.output(…)
^^^^^^^^^^^^^^^^^^

instance.output(‘parametername’)
''''''''''''''''''''''''''''''''

Return the value of the parameter named :math:`parametername`.

.. _instance.output-1:

instance.output()
'''''''''''''''''

Output a numpy record array with all the symbols and their values. The
same result is obtained by calling ``instance.output(’parametername’)``
and ``myoutputs = instance.output(); myoutputs[’parametername’]``.

instance.eval(…)
^^^^^^^^^^^^^^^^

We may wish to evaluate the fit curve, thus generating the smooth curve
for purposes of either plotting or further data processing. Therefore,
this function evaluates the curve fit, and returns an nddata object with
the same plot color property (see above) as the original data. It can be
called in several formats.

instance.eval(None)
'''''''''''''''''''

This just evaluates along the time axis for the data.

instance.eval(100)
''''''''''''''''''

Returns an nddata with 100 points; 100 can be replaced by any integer.

instance.eval(r_[0:0.2:100])
''''''''''''''''''''''''''''

This will evaluate the function along the fit axis, from 0 to 100, with
a datapoint every 0.2. Here, ``r_[0:0.2:100]`` can be replaced by any
ndarray.

instance.eval(…,set = listordict,set_to = list)
'''''''''''''''''''''''''''''''''''''''''''''''

Sometimes, one may want to see how the evaluated fit would change if a
parameter were altered. For this reason, this function takes the same
``set`` and ``set_to`` keyword arguments as ``fit``, except that the
parameters are set on a one-time basis, just for the evaluation.

instance.covar(…)
^^^^^^^^^^^^^^^^^

instance.covar(‘p1’)
''''''''''''''''''''

Returns the covariance for the fit parameter :math:`p1` (i.e. the
expected :math:`\sigma^2` for this parameter) 

.. todo::
    see the theory
    section about fitting errors.

instance.covar(‘p1’,‘p2’)
'''''''''''''''''''''''''

Returns the covariance between the fit parameters :math:`p1` and
:math:`p2`.

instance.covarmat(‘p1’,‘p2’,…,‘pN’)
'''''''''''''''''''''''''''''''''''

Returns an ndarray containing the covariance matrix for parameters
:math:`p1`\ …\ :math:`pN`.

.. _instance.covar-1:

instance.covar()
''''''''''''''''

Returns an ndarray record array with a labeled covariance matrix. This
function is ideal for printing with ``lrecordarray``; note that this
includes a field of data called “labels,”

.. todo::
    are the labels actually implemented?

which label the various rows.

.. _instance.latex-1:

instance.latex()
^^^^^^^^^^^^^^^^

Shows the function string, with the results of the fit substituted in
for the appropriate parameters.

instance.linear()
^^^^^^^^^^^^^^^^^

instance.errfunc(…)
^^^^^^^^^^^^^^^^^^^

Internal
~~~~~~~~

The following are functions used internally by routines in the fitdata
class:

-  instance._pn(…)

-  instance._taxis(…)

-  instance.add_inactive_p(…)

-  instance.analytical_covariance(…)

-  instance.errfunc(…)

-  instance.fitfunc(p,x) all references to the fit function should be
   made with this method

-  instance.gen_symbolic(…) used above

-  instance.gen_indices(…)

-  instance.parameter_derivatives(…)

-  instance.parameter_gradient(…)

-  instance.remove_inactive_p(…)

-  instance.makereal() efficiently finds the real part of the data to be
   fit

.. _sec:writeup_software_fitfunc_thingstouse_example:

example
~~~~~~~

:math:`T_1`
^^^^^^^^^^^

::

   from pyspecdata.nmr import * # includes the t1curve class defined here
   obs('Moved the guess function into the base class for 870')
   fl = figlistl() # make a figure list designed for output in latex
   t = double(r_[0.:2.:10j]) # a numpy array running from 0 to 1 second
   print 't is',lsafen(t)
   d = nddata(1.-2.1*exp(-t),[-1],['t']).labels('t',t) # generate data for a ``fake'' t1 curve from 0 to 1 second
   d.name('example $T_1$ curve') # give it a name
   d = t1curve(d) # now we initialize a new t1curve object from this example data
   print 'The functional format: ',d.function_string,'\n\n' #verify that this has the correct functional format.
   print 'd is',lsafen(d)
   fl.next('t1test') # move to the next (here a new) figure named t1test
   plot(d,'o',label = d.name()) # really, it should automatically pull the label from the name
   # now go ahead and fit it
   d.fit()
   print 'I fit d to',d.latex(),'\n\n'
   # then, show the fit
   plot(d.eval(100),label = d.name()+' fit')
   autolegend()
   fl.show('t1test120131.pdf') # dump out all our figures.

.. _sec:random_testofksp:

test of
^^^^^^^

::

   from pyspecdata.nmr import * # includes the t1curve class defined here
   from nmrfit import * # includes the t1curve class defined here
   obs('having problems with pinv, rerun!')
   fl = figlistl() # make a figure list designed for output in latex
   p = r_[0.:1.:10j] # a numpy array running from 0 to 1 second
   phalf = 0.1
   print 'p is',lsafen(p)
   d = nddata(p/(phalf+p),[-1],['p']).labels('p',p) # generate data for a ``fake'' asymptote from 0 to 1 second
   d.name('example asymptote curve') # give it a name
   d = ksp(d) # now we initialize a new t1curve object from this example data
   print 'The functional format: ',d.function_string,'\n\n' #verify that this has the correct functional format.
   print 'd is',lsafen(d)
   fl.next('asymptotetest') # move to the next (here a new) figure named asymptotetest
   plot(d,'o',label = d.name()) # really, it should automatically pull the label from the name
   # now go ahead and fit it
   d.fit() # more dramatic guessing
   print 'I fit d to',d.latex(),'\n\n'
   # then, show the fit
   plot(d.eval(100),label = d.name()+' fit')
   autolegend()
   fl.show('asymptotetest120201.pdf') # dump out all our figures.

