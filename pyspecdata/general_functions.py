r"""These are general functions that need to be accessible to everything inside
pyspecdata.core.  I can't just put these inside pyspecdata.core, because that
would lead to cyclic imports, and e.g. submodules of pyspecdata can't find
them."""

import os
import sys
from matplotlib.pylab import gci
from numpy import pi
def inside_sphinx():
    if len(sys.argv) > 0:
        return os.path.basename(sys.argv[0]) == "sphinx-build"
    else:
        return False
import numpy as np
import logging
import re

def balance_clims():
    """works with matplotlib to generate a plot
    appropriate for positive and negative
    from here:
        https://stackoverflow.com/questions/13060450/how-to-get-current-plots-clim-in-matplotlib
    """
    thisi = gci()
    these_clims = thisi.get_clim()
    a = max(abs(array(these_clims)))
    thisi.set_clim((-a,a))
    return
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
    kwargnames,kwargdefaultvals = list(zip(*listoftuples))
    output = []
    for j,val in enumerate(kwargnames):
        if val in list(kwargs.keys()):
            output.append(kwargs.pop(val))
        else:
            output.append(kwargdefaultvals[j])
    if not pass_through and len(kwargs) > 0:
        raise ValueError("I didn't understand the kwargs:",repr(kwargs))
    if len(output) > 1:
        return tuple(output)
    elif len(output) == 1:
        return output[0]
def autostringconvert(arg):
    if isinstance(arg,str):
        return str(arg)
    else:
        return arg
def check_ascending_axis(u,tolerance = 1e-7,additional_message = [], allow_descending=False):
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
    if isinstance(additional_message, str):
        additional_message = [additional_message]
    du = (u[-1]-u[0])/(len(u)-1.) # the dwell gives the bandwidth, whether or not it has been zero padded -- I calculate this way for better accuracy
    thismsg = ', '.join(additional_message + ["the axis must be ascending (and equally spaced)"])
    assert du > 0, thismsg
    thismsg = ', '.join(additional_message + ["the axis must be equally spaced (and ascending)"])
    assert all(abs(np.diff(u) - du)/du < tolerance), thismsg# absolute
    #   tolerance can be large relative to a du of ns -- don't use
    #   allclose/isclose, since they are more recent numpy additions
    if not allow_descending:
        assert du > 0, thismsg
    return du

def level_str_to_int(level):
    if type(level) is str:
        if level.lower() == 'info':
            level=logging.INFO
        elif level.lower() == 'debug':
            level=logging.DEBUG
        else:
            raise ValueError("if you give me level as a string, give me 'info' or 'debug'")
    return level
def init_logging(level=logging.DEBUG, stdout_level=logging.INFO, filename='pyspecdata.%d.log', fileno=0):
    r"""Initialize a decent logging setup to log to `~/pyspecdata.log` (and `~/pyspecdata.XX.log` if that's taken).

    By default, everything above "debug" is logged to a
    file, while everything above "info" is printed to
    stdout.

    Do NOT log if run from within a notebook (it's fair to
    assume that you will run first before embedding)
    """
    FORMAT = "--> %(filename)s(%(lineno)s):%(name)s %(funcName)20s %(asctime)20s\n%(levelname)s: %(message)s"
    level = level_str_to_int(level)
    stdout_level = level_str_to_int(stdout_level)
    min_level = min([level,stdout_level])
    formatter = logging.Formatter(FORMAT)
    log_filename = os.path.join(os.path.expanduser('~'),filename%fileno)
    if os.path.exists(log_filename):
        # manually remove, and then use append -- otherwise, it won't write to
        # file immediately
        try:
            os.remove(log_filename)
        except:
            if fileno == 0:
                print(f"{log_filename} appears to be locked or otherwise inaccessible: I'm going to explore other options for fileno")
            if fileno > 20:
                raise ValueError("I'm not going to increase fileno above 20 -- that's crazy time!")
            return init_logging(level=level, filename=filename, fileno=fileno+1)
    print(f"logging output to {log_filename}")
    logger = logging.getLogger()
    logger.setLevel(min_level) # even if I set the handler level, it won't
    #                        print w/out this
    file_handler = logging.FileHandler(log_filename, mode='a')
    stdout_handler = logging.StreamHandler(sys.stdout)
    # can set levels independently with:
    stdout_handler.setLevel(stdout_level)
    stdout_handler.setFormatter(formatter)
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)
    logger.addHandler(file_handler)
    return logger

def strm(*args):
    return ' '.join(map(str,args))

exp_re = re.compile(r'(.*)e([+\-])0*([0-9]+)')
def reformat_exp(arg):
    "reformat scientific notation in a nice latex format -- used in both pdf and jupyter notebooks"
    m = exp_re.match(arg)
    if 'i' not in arg and float(arg) == 0:
        return ''
    if m:
        retstr,pm,fin_numb = m.groups()
        retstr += r'\times 10^{'
        retstr += pm
        #retstr += pm if pm == '-' else ''
        retstr += fin_numb
        retstr += '}'
        return retstr
    else:
        return arg
def complex_str(arg, fancy_format=False, format_code='%.4g'):
    "render a complex string -- leaving out imaginary if it's real"
    retval = [format_code%arg.real]
    if arg.imag != 0.0:
        retval.append((format_code+"i")%arg.imag)
    retval = [reformat_exp(j) for j in retval]
    if len(retval)>1 and retval[1][0] not in '+-':
        retval[1] = '+'+retval[1]
    return ''.join(retval)
def render_matrix(arg, format_code='%.4g'):
    "return latex string representing 2D matrix"
    math_str = r'\begin{bmatrix}'
    math_str += '\n'
    if hasattr(arg.dtype,'fields') and arg.dtype.fields is not None:
        math_str += '\\\\\n'.join([' & '.join([', '.join([r'\text{'+f[0]+r'}\!=\!\text{"'+elem[f[0]]+'"}'
                                                          if isinstance(elem[f[0]],str)
                                                          else r'\text{%s}\!=\!%g'%(f[0],elem[f[0]])
                                                          for f in arg.dtype.descr])# f[0] is the name (vs. size)
                                               for elem in arg[k,:]]) for k in range(arg.shape[0])])
    else:
        math_str += '\\\\\n'.join([' & '.join([complex_str(j, format_code=format_code) for j in arg[k,:]]) for k in range(arg.shape[0])])
    math_str += '\n'
    math_str += r'\end{bmatrix}'
    return math_str
def redim_F_to_C(a):
    r"""the following creates a C array, reversing the *apparent* order of dimensions, while preserving the order in memory"""
    return a.ravel(order='F').reshape(a.shape[::-1], order='C')# 'C' not required, but for clarity
def redim_C_to_F(a):
    "see redim_F_to_C"
    return a.ravel(order='C').reshape(a.shape[::-1], order='F')
