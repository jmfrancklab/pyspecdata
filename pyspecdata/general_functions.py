r"""These are general functions that need to be accessible to everything inside
pyspecdata.core.  I can't just put these inside pyspecdata.core, because that
would lead to cyclic imports, and e.g. submodules of pyspecdata can't find
them."""

from pylab import *
import logging
import os
from paramset_pyspecdata import myparams

def process_kwargs(listoftuples, kwargs, pass_through=False, as_attr=False):
    '''This function allows dynamically processed (*i.e.* function definitions with `**kwargs`) kwargs (keyword arguments) to be dealt with in a fashion more like standard kwargs.
    The defaults set in `listoftuples` are used to process `kwargs`, which are then returned as a set of values (that are set to defaults as needed).

    Note that having `kwargs` as an explicit argument avoids errors where the user forgets to pass the `kwargs`.
    
    Parameters
    ==========
    kwargs : **dictionary

        The keyword arguments that you want to process.

    listoftuples : list of tuple pairs

        Tuple pairs, consisting of ``('param_name',param_value)``, that give the default values for the various parameters.

    pass_through : bool

        Defaults to False.  If it's true, then it's OK not to process all the kwargs here.
        In that case, the used kwargs are popped out of the dictionary, and you are expected to pass the unprocessed values (in the dictionary after the call) on to subsequent processing.
        Importantly, you should *always* end with a `pass_through`=`False` call of this function, or by passing **kwargs to a standard function in the standard way.
        Otherwise it's possible for the user to pass kwargs that are never processed!
    as_attr : bool, object

        Defaults to False. If not False, it must be an object whose attributes are set to the value of the respective kwargs.

    return : tuple

        It's expected that the output is assigned to variables with the **exact same** names as the string in the first half of the tuples, in the **exact same** order.
        These parameters will then be set to the appropriate values.
    '''
    kwargnames,kwargdefaultvals = zip(*listoftuples)
    output = []
    for j,val in enumerate(kwargnames):
        if val in kwargs.keys():
            output.append(kwargs.pop(val))
        else:
            output.append(kwargdefaultvals[j])
    if not pass_through and len(kwargs) > 0:
        raise ValueError("I didn't understand the kwargs:",repr(kwargs))
    return tuple(output)
def autostringconvert(arg):
    if isinstance(arg,basestring):
        return str(arg)
    else:
        return arg
def check_ascending_axis(u,tolerance = 1e-7,additional_message = []):
    r"""Check that the array `u` is ascending and equally spaced, and return the
    spacing, `du`.  This is a common check needed for FT functions, shears,
    etc.
    
    Parameters
    ----------

    tolerance : double
        The relative variation in `du` that is allowed.
        Defaults to 1e-7.

    additional_message : str
        So that the user can easily figure out where the assertion error is
        coming from, supply some extra text for the respective message.

    Returns
    -------

    du : double
        the spacing between the elements of u
    """
    if type(additional_message) is str:
        additional_message = [additional_message]
    du = (u[-1]-u[0])/(len(u)-1.) # the dwell gives the bandwidth, whether or not it has been zero padded -- I calculate this way for better accuracy
    thismsg = ', '.join(additional_message + ["the axis must be ascending (and equally spaced)"])
    assert du > 0, thismsg
    thismsg = ', '.join(additional_message + ["the axis must be equally spaced (and ascending)"])
    assert all(abs(diff(u) - du)/du < tolerance), thismsg# absolute
    #   tolerance can be large relative to a du of ns -- don't use
    #   allclose/isclose, since they are more recent numpy additions
    assert du > 0, thismsg
    return du

def init_logging(level=logging.INFO, filename='pyspecdata.log'):
    "Initialize logging on pyspecdata.log -- do NOT log if run from within a notebook (it's fair to assume that you will run first before embedding)"
    if myparams['figlist_type'] == 'figlistl':
        return
    FORMAT = "--> %(filename)s(%(lineno)s):%(name)s %(funcName)20s %(asctime)20s\n%(levelname)s: %(message)s"
    log_filename = os.path.join(os.path.expanduser('~'),filename)
    if os.path.exists(log_filename):
        # manually remove, and then use append -- otherwise, it won't write to
        # file immediately
        os.remove(log_filename)
    logging.basicConfig(format=FORMAT,
            filename=log_filename,
            filemode='a',
            level=level,
            )

def strm(*args):
    return ' '.join(map(str,args))
