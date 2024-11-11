r"""These are general functions that need to be accessible to everything inside
pyspecdata.core.  I can't just put these inside pyspecdata.core, because that
would lead to cyclic imports, and e.g. submodules of pyspecdata can't find
them."""

import os
import sys
from matplotlib.pylab import gci
from numpy import pi
import numpy as np
import logging
import re
import pint
import textwrap

ureg = pint.UnitRegistry()
ureg.define("cyc = [cycle]")  # 'cycle' is a new dimension
ureg.define("Hz = cyc / s")  # Redefine 'Hz' to be cycles per second
ureg.define("rad = 2*pi*cyc")  # Redefine 'rad' in terms of cycles
Q_ = ureg.Quantity
if "√" in str(Q_("√W")):
    pass
else:
    print(
        "**Warning!** I'm hacking the sqrt behavior of pint.  Consider using"
        " the jmfranck/pint fork"
    )

    def Q_(*args):
        if len(args) == 1:
            b = args[0]
            a = 1
        elif len(args) == 2:
            a, b = args
        else:
            raise ValueError("I don't know what to do with more than 2 args!")
        m = re.match(r"(.*)√(\w+)(.*)", b)
        if m:
            g1, g2, g3 = m.groups()
            b = g1 + f" {g2}" + "^{0.5} " + g3
            print(b)
        return ureg.Quantity(a, b)


def inside_sphinx():
    if len(sys.argv) > 0:
        return os.path.basename(sys.argv[0]) == "sphinx-build"
    else:
        return False


class CustomError(Exception):
    def __init__(self, *value, **kwargs):
        raise NotImplementedError(
            "You should get rid of CustomError and use explain_error instead"
        )
        return


def emptytest(x):  # test is it is one of various forms of np.empty
    if type(x) in [list, np.array]:
        if len(x) == 0:
            return True
        elif x is np.array(None):
            return True
        elif len(x) > 0:
            return False
        # don't want the following, because then I may need to pop, etc
        # if type(x) is list and all(map(lambda x: x is None,x)): return True
    if np.size(x) == 1 and x is None:
        return True
    if np.size(x) == 0:
        return True
    return False


def balance_clims():
    """works with matplotlib to generate a plot
    appropriate for positive and negative
    from here:
        https://stackoverflow.com/questions/13060450/\
                how-to-get-current-plots-clim-in-matplotlib
    """
    thisi = gci()
    these_clims = thisi.get_clim()
    a = max(abs(np.array(these_clims)))
    thisi.set_clim((-a, a))
    return


def process_kwargs(listoftuples, kwargs, pass_through=False, as_attr=False):
    """This function allows dynamically processed (*i.e.* function definitions
    with `**kwargs`) kwargs (keyword arguments) to be dealt with in a fashion
    more like standard kwargs.
    The defaults set in `listoftuples` are used to process `kwargs`, which are
    then returned as a set of values (that are set to defaults as needed).

    Note that having `kwargs` as an explicit argument avoids errors where the
    user forgets to pass the `kwargs`.

    Parameters
    ==========
    kwargs : **dictionary

        The keyword arguments that you want to process.

    listoftuples : list of tuple pairs

        Tuple pairs, consisting of ``('param_name',param_value)``, that give
        the default values for the various parameters.

    pass_through : bool

        Defaults to False.  If it's true, then it's OK not to process all the
        kwargs here.
        In that case, the used kwargs are popped out of the dictionary, and you
        are expected to pass the unprocessed values (in the dictionary after
        the call) on to subsequent processing.
        Importantly, you should *always* end with a `pass_through`=`False` call
        of this function, or by passing **kwargs to a standard function in the
        standard way.
        Otherwise it's possible for the user to pass kwargs that are never
        processed!
    as_attr : bool, object

        Defaults to False. If not False, it must be an object whose attributes
        are set to the value of the respective kwargs.

    return : tuple

        It's expected that the output is assigned to variables with the **exact
        same** names as the string in the first half of the tuples, in the
        **exact same** order.
        These parameters will then be set to the appropriate values.
    """
    kwargnames, kwargdefaultvals = list(zip(*listoftuples))
    output = []
    for j, val in enumerate(kwargnames):
        if val in list(kwargs.keys()):
            output.append(kwargs.pop(val))
        else:
            output.append(kwargdefaultvals[j])
    if not pass_through and len(kwargs) > 0:
        raise ValueError("I didn't understand the kwargs:", repr(kwargs))
    if len(output) > 1:
        return tuple(output)
    elif len(output) == 1:
        return output[0]


def autostringconvert(arg):
    if isinstance(arg, str):
        return str(arg)
    else:
        return arg


def check_ascending_axis(
    u, tolerance=1e-7, additional_message=[], allow_descending=False
):
    r"""Check that the array `u` is ascending and equally spaced, and return
    the spacing, `du`.  This is a common check needed for FT functions, shears,
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
    du = (u[-1] - u[0]) / (
        len(u) - 1.0
    )  # the dwell gives the bandwidth, whether or not it has been zero padded
    #    -- I calculate this way for better accuracy
    thismsg = ", ".join(
        additional_message + ["the axis must be equally spaced"]
    )
    assert all(abs(np.diff(u) - du) / du < tolerance), thismsg  # absolute
    #   tolerance can be large relative to a du of ns -- don't use
    #   allclose/isclose, since they are more recent numpy additions
    if not allow_descending:
        thismsg = ", ".join(
            additional_message
            + ["the axis must be equally spaced (and ascending)"]
        )
        assert du > 0, thismsg
    return du


def level_str_to_int(level):
    if type(level) is str:
        if level.lower() == "info":
            level = logging.INFO
        elif level.lower() == "debug":
            level = logging.DEBUG
        else:
            raise ValueError(
                "if you give me level as a string, give me 'info' or 'debug'"
            )
    return level


if (
    "pyspecdata_figures" in os.environ
    and os.environ["pyspecdata_figures"] == "latex"
):
    print_log_info = False
else:
    print_log_info = True


def init_logging(
    level=logging.DEBUG,
    stdout_level=logging.INFO,
    filename="pyspecdata.%d.log",
    fileno=0,
):
    r"""Initialize a decent logging setup to log to `~/pyspecdata.log` (and
    `~/pyspecdata.XX.log` if that's taken).

    By default, everything above "debug" is logged to a
    file, while everything above "info" is printed to
    stdout.

    Do NOT log if run from within a notebook (it's fair to
    assume that you will run first before embedding)
    """
    FORMAT = (
        "--> %(filename)s(%(lineno)s):%(name)s %(funcName)20s"
        " %(asctime)20s\n%(levelname)s: %(message)s"
    )
    level = level_str_to_int(level)
    stdout_level = level_str_to_int(stdout_level)
    min_level = min([level, stdout_level])
    formatter = logging.Formatter(FORMAT)
    log_filename = os.path.join(os.path.expanduser("~"), filename % fileno)
    local_print = True
    if os.path.exists(log_filename):
        # manually remove, and then use append -- otherwise, it won't write to
        # file immediately
        try:
            os.remove(log_filename)
        except Exception:
            if fileno == 0:
                if print_log_info:
                    print(
                        f"{log_filename} appears to be locked or otherwise"
                        " inaccessible: I'm going to explore other options"
                        " for fileno"
                    )
            if fileno > 20:
                raise ValueError(
                    "I'm not going to increase fileno above 20 -- that's crazy"
                    " time!"
                )
            local_print = False
            return init_logging(
                level=level, filename=filename, fileno=fileno + 1
            )
    if print_log_info and local_print:
        print(
            "-" * 10
            + "  "
            + f"logging output to {log_filename}"
            + "  "
            + "-" * 10
        )
    logger = logging.getLogger()
    logger.setLevel(min_level)  # even if I set the handler level, it won't
    #                        print w/out this
    file_handler = logging.FileHandler(
        log_filename, mode="a", encoding="utf-8"
    )
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
    return " ".join(map(str, args))


exp_re = re.compile(r"(.*)e([+\-])0*([0-9]+)")


def reformat_exp(arg):
    """reformat scientific notation in a nice latex format -- used in both pdf
    and jupyter notebooks"""
    m = exp_re.match(arg)
    if "i" not in arg and float(arg) == 0:
        return ""
    if m:
        retstr, pm, fin_numb = m.groups()
        retstr += r"\times 10^{"
        retstr += pm
        # retstr += pm if pm == '-' else ''
        retstr += fin_numb
        retstr += "}"
        return retstr
    else:
        return arg


def complex_str(arg, fancy_format=False, format_code="%.4g"):
    "render a complex string -- leaving out imaginary if it's real"
    retval = [format_code % arg.real]
    if arg.imag != 0.0:
        retval.append((format_code + "i") % arg.imag)
    retval = [reformat_exp(j) for j in retval]
    if len(retval) > 1 and retval[1][0] not in "+-":
        retval[1] = "+" + retval[1]
    return "".join(retval)


def render_matrix(arg, format_code="%.4g"):
    "return latex string representing 2D matrix"
    math_str = r"\begin{bmatrix}"
    math_str += "\n"
    if hasattr(arg.dtype, "fields") and arg.dtype.fields is not None:
        math_str += "\\\\\n".join([
            " & ".join([
                ", ".join([
                    (
                        r"\text{" + f[0] + r'}\!=\!\text{"' + elem[f[0]] + '"}'
                        if isinstance(elem[f[0]], str)
                        else r"\text{%s}\!=\!%g" % (f[0], elem[f[0]])
                    )
                    for f in arg.dtype.descr
                ])  # f[0] is the name (vs. size)
                for elem in arg[k, :]
            ])
            for k in range(arg.shape[0])
        ])
    else:
        math_str += "\\\\\n".join([
            " & ".join(
                [complex_str(j, format_code=format_code) for j in arg[k, :]]
            )
            for k in range(arg.shape[0])
        ])
    math_str += "\n"
    math_str += r"\end{bmatrix}"
    return math_str


def redim_F_to_C(a):
    r"""the following creates a C array, reversing the *apparent* order of
    dimensions, while preserving the order in memory"""
    return a.ravel(order="F").reshape(
        a.shape[::-1], order="C"
    )  # 'C' not required, but for clarity


def redim_C_to_F(a):
    "see redim_F_to_C"
    return a.ravel(order="C").reshape(a.shape[::-1], order="F")


def fname_makenice(fname):
    fname = fname.replace(" ", "_")
    fname = fname.replace("-", "m")
    fname = fname.replace("+", "p")
    fname = fname.replace(",", "_")
    fname = fname.replace("\\", "_")
    fname = fname.replace("$", "")
    fname = fname.replace("(", "")
    fname = fname.replace(")", "")
    fname = fname.replace('"', "")
    fname = fname.replace("=", "_")
    fname = fname.replace("\n", "_")
    fname = fname.replace("*", "_star_")
    fname = fname.replace(":", "")
    fname = fname.replace("^", "")
    fname = fname.replace("}", "")
    fname = fname.replace("{", "")
    return fname


def lsafen(*string, **kwargs):
    "see lsafe, but with an added double newline"
    string = list(string)
    string += ["\n\n"]
    return lsafe(*tuple(string), **kwargs)


def lsafe(*string, **kwargs):
    "Output properly escaped for latex"
    if len(string) > 1:
        lsafewkargs = lambda x: lsafe(x, **kwargs)
        return " ".join(list(map(lsafewkargs, string)))
    else:
        string = string[0]
    # {{{ kwargs
    spaces = False
    if "spaces" in list(kwargs.keys()):
        spaces = kwargs.pop("spaces")
    if "wrap" in list(kwargs.keys()):
        wrap = kwargs.pop("wrap")
    else:
        wrap = None
    # }}}
    if not isinstance(string, str):
        string = str(string)
    if wrap is True:
        wrap = 60
    if wrap is not None:
        string = "\n".join(textwrap.wrap(string, wrap))
    string = string.replace("\\", "\\textbackslash ")
    if spaces:
        string = string.replace(" ", "\\ ")
    string = string.replace("\n\t", "\n\n\\quad ")
    string = string.replace("\t", "\\quad ")
    string = string.replace("_", r"\_")
    string = string.replace("{", r"\{")
    string = string.replace("}", r"\}")
    string = string.replace("$$", r"ACTUALDOUBLEDOLLAR")
    string = string.replace("]", r"$]$")
    string = string.replace("[", r"$[$")
    string = string.replace("<", r"$<$")
    string = string.replace(">", r"$>$")
    string = string.replace("$$", r"")
    string = string.replace("ACTUALDOUBLEDOLLAR", r"$$")
    string = string.replace("^", r"\^")
    string = string.replace("#", r"\#")
    string = string.replace("%", r"\%")
    string = string.replace("&", r"\&")
    string = string.replace("+/-", r"\ensuremath{\pm}")
    string = string.replace("|", r"$|$")
    return string


def copy_maybe_none(input):
    if input is None:
        return None
    else:
        if isinstance(input, list):
            return list(map(np.copy, input))
        else:
            return input.copy()


def whereblocks(a):
    """returns contiguous chunks where the condition is true
    but, see the "contiguous" method, which is more OO"""
    parselist = np.where(a)[0]
    jumps_at = np.where(np.diff(parselist) > 1)[0] + 1
    retlist = []
    lastjump = 0
    for jump in jumps_at:
        retlist += [parselist[lastjump:jump]]
        lastjump = jump
    retlist += [parselist[lastjump:]]
    return retlist


def box_muller(length, return_complex=True):
    r"""algorithm to generate normally distributed noise"""
    s1 = np.random.rand(length)
    s2 = np.random.rand(length)
    n1 = np.sqrt(-2 * np.log(s1)) * np.cos(2 * pi * s2)
    if return_complex:
        n2 = np.sqrt(-2 * np.log(s1)) * np.sin(2 * pi * s2)
        return (n1 + 1j * n2) * 0.5
    else:
        return (n1) * 0.5


def dp(number, decimalplaces=2, scientific=False, max_front=3):
    """format out to a certain decimal places, potentially in scientific
    notation

    Parameters
    ----------
    decimalplaces: int (optional, default 3)
        number of decimal places
    scientific: boolean (optional, default False)
        use scientific notation
    max_front: int (optional, default 3)
        at most this many places in front of the decimal before switching
        automatically to scientific notation.
    """
    if scientific:
        logging.debug(
            strm("trying to convert", number, "to scientific notation")
        )
        tenlog = int(np.floor(np.log10(abs(number))))
        number /= 10**tenlog
        fstring = "%0." + "%d" % decimalplaces + r"f\times 10^{%d}" % tenlog
    else:
        fstring = "%0." + "%d" % decimalplaces + "f"
        if len(fstring % number) > 1 + decimalplaces + max_front:
            return dp(number, decimalplaces=decimalplaces, scientific=True)
    return fstring % number


def fa(input, dtype="complex128"):  # make a fortran array
    return np.array(
        input, order="F", dtype=dtype
    )  # will need transpose reverses the dimensions, since the bracketing
    #    still works in C order (inner is last index), but F tells it to store
    #    it appropriately in memory


def ndgrid(*input):
    thissize = list([1])
    thissize = thissize * len(input)
    output = list()
    for j in range(0, len(input)):
        tempsize = np.copy(thissize)
        tempsize[j] = input[j].size
        output.append(input[j].reshape(tempsize))
    return output


def pinvr(C, alpha):
    U, S, V = np.linalg.svd(C, full_matrices=0)
    # print 'U S V shapes:'
    # print U.shape
    # print S.shape
    # print V.shape
    if np.any(~np.isfinite(U)):
        raise ValueError("pinvr error, U is not finite")
    if np.any(~np.isfinite(V)):
        raise ValueError("pinvr error, V is not finite")
    if np.any(~np.isfinite(S)):
        raise ValueError("pinvr error, S is not finite")
    S = np.diag(S / (S**2 + alpha**2))
    if np.any(~np.isfinite(S)):
        raise ValueError(
            "pinvr error, problem with S/(S^2+alpha^2) --> set your"
            " regularization higher"
        )
    return np.dot(
        np.conj(np.transpose(V)), np.dot(S, np.conj(np.transpose(U)))
    )


def sech(x):
    return 1.0 / np.cosh(x)


def myfilter(x, center=250e3, sigma=100e3):
    x = (x - center) ** 2
    x /= sigma**2
    return np.exp(-x)


def det_unit_prefactor(thisstr):
    "use pint to determine the prefactor of the string-formatted unit thisstr"
    return 3 * int(np.log10(ureg(thisstr).to_base_units().magnitude) // 3)
