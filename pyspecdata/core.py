r"""Provides the core components of pyspecdata.
Currently, this is a very large file that we will slowly break down into
separate modules or packages.

The classes :class:`nddata`, :class:`nddata_hdf`, :class:`ndshape`, the
function :func:`plot`, and the class :class:`fitdata`
are the core components of the N-Dimensional processing routines.
Start by familiarizing yourself with those.

The :class:`figlist` is the base class for "Figure lists."
Figure lists allows you to organize plots and text and to refer to plots
by name, rather than number.
They are designed so that same code can be used seamlessly from within
ipython, jupyter, a python script, or a python environment within latex
(JMF can also distribute latex code for this -- nice python based
installer is planned).
The user does not initialize the figlist class directly,
but rather initializes ``figlist_var``.
At the end of this file,
there is a snippet of code that sets
``figlist_var`` to choice that's appropriate for the working environment
(*i.e.*, python, latex environment, *etc.)

There are many helper and utility functions that need to be sorted an
documented by JMF,
and can be ignored.
These are somewhat wide-ranging in nature.
For example, :func:`box_muller` is a helper function (based on numerical
recipes) used by :func:`nddata.add_noise`,
while h5 functions are helper functions for using pytables in a fashion that
will hopefull be intuitive to those familiar with SQL, etc.
"""

from .datadir import pyspec_config
from .matrix_math.svd import svd as MM_svd
from .matrix_math.dot import dot as MM_dot
from .matrix_math.dot import matmul as MM_matmul
from .matrix_math.dot import along as MM_along
from .matrix_math.nnls import nnls as MM_nnls
from os import environ
import numpy as np
import sympy as sp
from numpy import r_, c_, nan, inf, pi
from mpl_toolkits.mplot3d import axes3d
import textwrap, atexit, scipy, warnings, inspect, re
from copy import deepcopy
from sympy.functions.elementary.miscellaneous import sqrt as sympy_sqrt
import scipy.sparse as sparse
import numpy.lib.recfunctions as recf
from scipy.interpolate import interp1d
import logging
from .ndshape import ndshape_base
from . import fourier as this_fourier
from . import axis_manipulation
from . import plot_funcs as this_plotting
from .general_functions import lsafe as orig_lsafe
from .general_functions import (
    emptytest,
    process_kwargs,
    autostringconvert,
    strm,
    lsafen,
    dp,
    pinvr,
    Q_,
)
from .hdf_utils import (
    h5loaddict,
    h5child,
    h5table,
    h5nodebypath,
    h5attachattributes,
)
from .mpl_utils import (
    plot_label_points,
    default_cycler,
)

# {{{ determine the figure style, and load the appropriate modules
_figure_mode_setting = pyspec_config.get_setting(
    "figures", section="mode", environ="pyspecdata_figures"
)
if _figure_mode_setting is None:
    print(
        "Warning!  Figure mode is not set, so I'm going to set it to standard"
        " by default!!!"
    )
    _figure_mode_setting = "standard"
    pyspec_config.set_setting("mode", "figures", "standard")
    import matplotlib.pyplot as plt
elif _figure_mode_setting == "latex":
    environ["ETS_TOOLKIT"] = "qt4"
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt
# }}} -- continued below

# rc('image',aspect='auto',interpolation='bilinear') # don't use this, because
# it gives weird figures in the pdf
plt.rc("image", aspect="auto", interpolation="nearest")
# rcParams['text.usetex'] = True
plt.rc("font", family="Arial")  # I need this to render unicode
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
# rcParams['ytick.major.size'] = 12
# rcParams['ytick.minor.size'] = 6
# rcParams['lines.linewidth'] = 3.0
# rcParams['legend.fontsize'] = 12
# rcParams['font.size'] = 6
plt.rcParams["axes.grid"] = False
plt.rcParams["image.cmap"] = "jet"
plt.rcParams["figure.figsize"] = (7 * (1 + np.sqrt(5)) / 2, 7)
mat2array = [
    {"ImmutableMatrix": np.array},
    "numpy",
]  # for sympy returns arrays rather than the stupid matrix class
logger = logging.getLogger("pyspecdata.core")
# {{{ constants
k_B = 1.380648813e-23
mu_0 = 4e-7 * pi
mu_B = 9.27400968e-24  # Bohr magneton
epsilon_0 = 8.854187817e-12
hbar = 6.6260695729e-34 / 2.0 / pi
N_A = 6.02214179e23
gammabar_H = 4.258e7
gammabar_D = gammabar_H * 61.422391 / 400.13  # ratio from Bruker BF
gammabar_e = 2.807e10  # this is for a nitroxide


# }}}
def det_oom(data_to_test):
    """determine the average order of magnitude -- for prefixing units

    Parameters
    ==========
    data_to_test: ndarray
        a numpy array (e.g. the result of a .getaxis( call
    Returns
    =======
    average_oom: int
        the average order of magnitude, rounded to the nearest multiple of 3
    """
    try:
        data_to_test = data_to_test[np.isfinite(data_to_test)]
    except Exception:
        raise ValueError(
            strm(
                "data_to_test is",
                data_to_test,
                "isfinite is",
                np.isfinite(data_to_test),
            )
        )
    if len(data_to_test) == 0:
        raise ValueError("this axis doesn't seem to have any sensible values!")
    # {{{ find the average order of magnitude, rounded down to the nearest
    #     power of 3
    average_oom = abs(data_to_test)
    average_oom = average_oom[average_oom != 0]
    average_oom = np.log10(average_oom) / 3.0
    logger.debug(strm("dtype", data_to_test.dtype))
    logger.debug(strm("oom:", average_oom))
    average_oom = average_oom[np.isfinite(average_oom)].mean()
    # }}}
    logger.debug(strm("the average oom is", average_oom * 3))
    average_oom = 3 * np.floor(average_oom)
    logger.debug(strm("I round this to", average_oom))
    return average_oom


def apply_oom(average_oom, numbers, prev_label=""):
    """scale numbers by the order of magnitude average_oom and change the
    name of the units by adding the appropriate SI prefix

    Parameters
    ----------
    average_oom: int or float
        the average order of magnitude to use
    numbers: np.ndarray
        The numbers to be scaled by average_oom.
        The np.array is modified in-place.
    prev_label: str
        a string representing the units

    Returns
    -------
    new_label: str
        prev_label is prefixed by the appropriate SI prefix
    """
    oom_names = ["T", "G", "M", "k", "", "m", "μ", "n", "p"]
    oom_values = r_[12, 9, 6, 3, 0, -3, -6, -9, -12]
    eq = oom_values == average_oom
    if not np.any(eq):
        if all(average_oom < oom_values):
            oom_index = len(oom_values) - 1
        elif all(average_oom > oom_values):
            oom_index = 0
        else:
            raise ValueError(
                strm(
                    "you passed",
                    average_oom,
                    "which I can't find a prefix for",
                )
            )
    else:
        oom_index = np.nonzero(eq)[0][0]
    if numbers.dtype in ["int32", "int64"]:
        # this is not necessary unless we have an integer type
        logger.warning(
            "you are trying to determine the SI prefix of a"
            "set of numbers that are described by integers.  This is"
            "probably not a good idea!!"
        )
        new_values = numbers / 10.0 ** oom_values[oom_index]
        numbers[:] = new_values.astype(numbers.dtype)
    else:
        numbers[:] /= 10.0 ** oom_values[oom_index]
    return oom_names[oom_index] + prev_label


def issympy(x):
    "tests if something is sympy (based on the module name)"
    return isinstance(x, sp.core.Expr)


# {{{ function trickery
def mydiff(data, axis=-1):
    """this will replace np.diff with a version that has the same number of
    indices, with the last being the copy of the first"""
    newdata = np.empty(np.shape(data), dtype=data.dtype)
    indices = [slice(None, None, None)] * len(data.shape)
    indices[axis] = slice(None, -1, None)
    newdata[tuple(indices)] = np.diff(data, axis=axis)
    # setfrom = list(indices)
    # indices[axis] = -1
    # setfrom[axis] = 0
    # newdata[indices] = newdata[setfrom]
    return newdata


# }}}
def normal_attrs(obj):
    myattrs = [
        x for x in dir(obj) if not inspect.ismethod(obj.__getattribute__(x))
    ]
    myattrs = [x for x in myattrs if not x[0:2] == "__"]
    # next line filters out properties
    myattrs = [
        x for x in myattrs if x not in ["C", "angle", "imag", "real", "shape"]
    ]
    return myattrs


def showtype(x):
    if isinstance(x, np.ndarray):
        return np.ndarray, x.dtype
    else:
        return type(x)


def emptyfunction():
    pass


# {{{ structured np.array helper functions
def make_bar_graph_indices(
    mystructarray, list_of_text_fields, recursion_depth=0, spacing=0.1
):
    """This is a recursive function that is used as part of textlabel_bargraph;
    it does NOT work without the sorting given at the beginning of that
    function"""
    # {{{ if there are still text fields left, then break down the np.array
    #     further, otherwise, just return the indices for this subarray
    if len(list_of_text_fields) > 0:
        unique_values = np.unique(
            mystructarray[list_of_text_fields[0]]
        )  # the return_index argument doesn't do what it's supposed to all the
        #    time, so I have to manually find the start indices, as given in
        #    the following line
        start_indices = [
            np.nonzero(mystructarray[list_of_text_fields[0]] == val)[0][0]
            for val in unique_values
        ]
        # find the structured np.array for the unique value
        index_values = []
        label_values = []
        start_indices = r_[
            start_indices, len(mystructarray)
        ]  # I add this so I can do the next step
        logger.debug(
            strm(
                "recursion depth is",
                recursion_depth,
                "and I am analyzing",
                list_of_text_fields[0],
                ": ",
            )
        )
        logger.debug(
            strm(
                "I found these unique values:",
                unique_values,
                "at these start indices:",
                start_indices[:-1],
            )
        )
        for k in range(0, len(start_indices) - 1):
            logger.debug(
                strm(
                    "recursion depth is",
                    recursion_depth,
                    "and I am analyzing",
                    list_of_text_fields[0],
                    ": ",
                )
            )
            logger.debug(
                strm(
                    "trying to extract unique value",
                    unique_values[k],
                    "using the range",
                    start_indices[k],
                    start_indices[k + 1],
                )
            )
            logger.debug(strm("which has this data"))
            indiv_struct_array = mystructarray[
                start_indices[k] : start_indices[k + 1]
            ]
            logger.debug(strm(lsafen(indiv_struct_array)))
            these_index_values, these_labels = make_bar_graph_indices(
                indiv_struct_array,
                list_of_text_fields[1:],
                recursion_depth=recursion_depth + 1,
            )
            index_values.append(these_index_values)
            label_values.append(
                [str(unique_values[k]) + "," + j for j in these_labels]
            )
        # {{{ scale the result of each call down to the equal size (regardless
        #     of number of elements), shift by the position in this np.array,
        #     and return
        logger.debug(
            strm(
                "recursion depth is",
                recursion_depth,
                "and I just COMPLETED THE LOOP, which gives a list of index"
                " values like this",
                index_values,
            )
        )
        max_indices = max(
            np.array(list(map(len, index_values)), dtype="double")
        )  # the maximum width of the np.array inside
        index_values = [
            x + (max_indices - len(x)) / 2.0 for x in index_values
        ]  # if the bar is less than max indices, shift it over, so it's still
        #    in the center
        logger.debug(
            strm(
                "recursion depth is",
                recursion_depth,
                "and I centered each set like this",
                index_values,
            )
        )
        index_values = [
            x / max_indices * (1 - spacing) + (1 - spacing) / 2
            for x in index_values
        ]  # scale down, so the width from left edge of bar to right edge of
        #    largest bar runs 0--> 1
        logger.debug(
            strm(
                "recursion depth is",
                recursion_depth,
                "and I scaled down so each runs zero to one*(1-spacing)"
                " (centered) like this",
                index_values,
            )
        )
        # this adds an index value, and also collapses down to a single
        # dimension list
        retval_indices = [
            x + num for num, val in enumerate(index_values) for x in val
        ]
        # now collapse labels down to a single dimension
        retval_labels = [k for j in label_values for k in j]
        logger.debug(
            strm(
                "recursion depth is",
                recursion_depth,
                "and I am passing up indices",
                retval_indices,
                "and labels",
                retval_labels,
            )
        )
        return retval_indices, retval_labels
        # }}}
    else:
        logger.debug(strm("recursion depth is", recursion_depth))
        N = len(mystructarray)
        logger.debug(
            strm(
                "hit innermost (no text labels left) and passing up a list of"
                " indices that looks like this:",
                r_[0:N],
            )
        )
        return r_[0:N], [""] * N
    # }}}


def textlabel_bargraph(
    mystructarray, othersort=None, spacing=0.1, ax=None, tickfontsize=8
):
    if ax is None:
        thisfig = plt.gcf()
        ax = thisfig.add_axes([0.2, 0.5, 0.8, 0.5])
        try:
            ax.tick_params(axis="both", which="major", labelsize=tickfontsize)
            ax.tick_params(axis="both", which="minor", labelsize=tickfontsize)
        except Exception:
            print(
                "Warning, in this version I can't set the tick params method"
                " for the axis"
            )
    # {{{ find the text fields, put them first, and sort by them
    mystructarray = mystructarray.copy()
    list_of_text_fields = [
        str(j[0]) for j in mystructarray.dtype.descr if j[1][0:2] == "|S"
    ]
    mystructarray = sorted(
        mystructarray[
            list_of_text_fields
            + [
                x[0]
                for x in mystructarray.dtype.descr
                if x[0] not in list_of_text_fields
            ]
        ]
    )
    logger.debug(
        strm("test --> now, it has this form:", lsafen(mystructarray))
    )
    # }}}
    error_fields = [
        str(j) for j in mystructarray.dtype.names if j[-6:] == "_ERROR"
    ]
    if len(error_fields) > 0:
        mystructarray_errors = mystructarray[error_fields]
        logger.debug(strm("found error fields:", mystructarray_errors))
    mystructarray = mystructarray[
        [str(j) for j in mystructarray.dtype.names if j not in error_fields]
    ]
    if othersort is not None:
        list_of_text_fields.append(othersort)
    logger.debug(strm("list of text fields is", lsafen(list_of_text_fields)))
    indices, labels = make_bar_graph_indices(
        mystructarray, list_of_text_fields, spacing=spacing
    )
    temp = list(zip(indices, labels))
    logger.debug(strm("(indices,labels) (len %d):" % len(temp), lsafen(temp)))
    logger.debug(
        strm(
            "I get these labels (len %d):" % len(labels),
            labels,
            "for the data (len %d)" % len(mystructarray),
            lsafen(mystructarray),
        )
    )
    indices = np.array(indices)
    indiv_width = min(np.diff(indices)) * (1 - spacing)
    remaining_fields = [
        x for x in mystructarray.dtype.names if x not in list_of_text_fields
    ]  # so they are in the right order, since set does not preserve order
    logger.debug(
        strm(
            "The list of remaining (i.e. non-text) fields is",
            lsafen(remaining_fields),
        )
    )
    colors = ["b", "g", "r", "c", "m", "k"]
    rects = []
    for j, thisfield in enumerate(remaining_fields):
        field_bar_width = indiv_width / len(remaining_fields)
        thiserror = None
        if thisfield + "_ERROR" in error_fields:
            thiserror = mystructarray_errors[thisfield + "_ERROR"]
        try:
            rects.append(
                ax.bar(
                    indices + j * field_bar_width,
                    mystructarray[thisfield],
                    field_bar_width,
                    color=colors[j],
                    yerr=thiserror,  # just to test
                    ecolor="k",
                    label="$%s$" % thisfield,
                )
            )
        except Exception:
            raise RuntimeError(
                strm(
                    "Problem with bar graph: there are %d indices, but %d"
                    " pieces of data"
                    % (len(indices), len(mystructarray[thisfield])),
                    "indices:",
                    indices,
                    "data",
                    mystructarray[thisfield],
                )
            )
    ax.set_xticks(indices + indiv_width / 2)
    ax.set_xticklabels(labels)
    ax.legend(
        [j[0] for j in rects],
        ["$%s$" % j for j in remaining_fields],
        loc="best",
    )
    return


def lookup_rec(A, B, indexpair):
    r"""look up information about A in table B (i.e. chemical by index, etc)
    indexpair is either the name of the index
    or -- if it's differently named -- the pair of indices
    given in (A,B) respectively

    This will just drop any fields in B that are also in A,
    and the output uses the first indexname

    note that it it seems like the join_rec function above may be more
    efficient!!
    """
    raise RuntimeError("You should now use decorate_rec!!")
    if type(indexpair) not in [tuple, list]:
        indexpair = (indexpair, indexpair)
    B = recf.drop_fields(
        B, (set(B.dtype.names) & set(A.dtype.names)) - {indexpair[1]}
    )  # indexpair for B gets dropped later anyways
    joined = []
    for j in A:
        matchedrows = B[B[indexpair[1]] == j[indexpair[0]]]
        for matchedrow in matchedrows:
            joined.append((j, matchedrow))
    if len(joined) == 0:
        raise IndexError(
            strm(
                "Unable to find any matches between",
                A[indexpair[0]],
                "and",
                B[indexpair[1]],
                "!",
            )
        )
    whichisindex = joined[0][1].dtype.names.index(indexpair[1])
    allbutindex = (
        lambda x: list(x)[0:whichisindex] + list(x)[whichisindex + 1 :]
    )
    joined = np.concatenate([
        np.array(
            tuple(list(j[0]) + allbutindex(j[1])),
            dtype=np.dtype(j[0].dtype.descr + allbutindex(j[1].dtype.descr)),
        ).reshape(1)
        for j in joined
    ])
    return joined


def reorder_rec(myarray, listofnames, first=True):
    try:
        indices_to_move = [myarray.dtype.names.index(j) for j in listofnames]
    except Exception:
        stuff_not_found = [
            j for j in listofnames if j not in myarray.dtype.names
        ]
        if len(stuff_not_found) > 0:
            raise IndexError(
                strm(
                    stuff_not_found,
                    "is/are in the list you passed,",
                    "but not one of the fields, which are",
                    myarray.dtype.names,
                )
            )
        else:
            raise RuntimeError("unknown problem")
    old_type = list(myarray.dtype.descr)
    new_type = [old_type[j] for j in indices_to_move] + [
        old_type[j]
        for j in range(0, len(old_type))
        if j not in indices_to_move
    ]
    new_list_of_data = [myarray[j[0]] for j in new_type]
    return np.core.rec.fromarrays(new_list_of_data, dtype=new_type)


def lambda_rec(myarray, myname, myfunction, *varargs):
    r"""make a new field "myname" which consists of "myfunction" evaluated with
    the fields given by "myargs" as arguments
    the new field is always placed after the last argument name
    if myname is in myargs, the original row is popped"""
    if len(varargs) == 1:
        myargs = varargs[0]
    elif len(varargs) == 0:
        myargs = [myname]
    else:
        raise IndexError(
            "For the fourth argument, you must pass either a list"
            " with the names of the arguments, or nothing (to use the field"
            " itself as an argument)"
        )
    myargs = autostringconvert(myargs)
    if isinstance(myargs, str):
        myargs = (myargs,)
    if not isinstance(myargs, tuple):
        myargs = tuple(myargs)
    argdata = list(map((lambda x: myarray[x]), myargs))
    try:
        newrow = myfunction(*tuple(argdata))
    except TypeError:
        newrow = np.array([
            myfunction(*tuple([x[rownumber] for x in argdata]))
            for rownumber in range(0, len(argdata[0]))
        ])
    if isinstance(newrow, list) and isinstance(newrow[0], str):
        newrow = np.array(newrow, dtype="|S100")
    try:
        new_field_type = list(newrow.dtype.descr[0])
    except AttributeError:
        raise IndexError(
            strm(
                "evaluated function on",
                argdata,
                "and got back",
                newrow,
                "which appears not to be a numpy np.array",
            )
        )
    new_field_type[0] = myname
    starting_names = myarray.dtype.names
    # {{{ make the dtype
    new_dtype = list(myarray.dtype.descr)
    # {{{ determine if I need to pop one of the existing rows due to a name
    #     conflict
    eliminate = None
    if myname in myargs:
        eliminate = myname
        insert_before = starting_names.index(
            myname
        )  # if we are replacing, we want it in the same place
        new_dtype = [j for j in new_dtype if j[0] != eliminate]
    # }}}
    # if we haven't already eliminated, determine where to put it
    if eliminate is None:
        insert_before = starting_names.index(myargs[-1]) + 1
    # insert the new field where it goes
    new_dtype.insert(insert_before, tuple(new_field_type))
    # }}}
    # {{{ separate starting_names and ending_names
    if eliminate is None:
        ending_names = starting_names[insert_before:]
        starting_names = starting_names[:insert_before]
    else:  # if I'm eliminating, I don't want to include the eliminated one
        ending_names = starting_names[insert_before + 1 :]
        starting_names = starting_names[:insert_before]
    # }}}
    return np.core.rec.fromarrays(
        [myarray[x] for x in starting_names if x != eliminate]
        + [newrow]
        + [myarray[x] for x in ending_names if x != eliminate],
        dtype=new_dtype,
    )


def join_rec(xxx_todo_changeme, xxx_todo_changeme1):
    (A, a_ind) = xxx_todo_changeme
    (B, b_ind) = xxx_todo_changeme1
    raise RuntimeError("You should now use decorate_rec!!")


def decorate_rec(xxx_todo_changeme2, xxx_todo_changeme3, drop_rows=False):
    r"""Decorate the rows in A with information in B --> if names overlap,
    keep the np.ones in A
    b_ind and a_ind can be either a single key, or a list of keys;
    if more than one element in B matches that in A, include both options!!"""
    (A, a_ind) = xxx_todo_changeme2
    (B, b_ind) = xxx_todo_changeme3
    dropped_rows = None
    # first find the list of indices that give us the data we want
    # {{{ process the arguments
    if (isinstance(b_ind, str)) and (isinstance(a_ind, str)):
        b_ind = [b_ind]
        a_ind = [a_ind]
    if ((isinstance(b_ind, list)) and (isinstance(a_ind, list))) and (
        len(b_ind) == len(a_ind)
    ):
        pass
    else:
        raise ValueError(
            "If you call a list for b_ind and/or a_ind, they must match in"
            " length!!!"
        )
    if np.any([x not in B.dtype.names for x in b_ind]):
        problem_index = [x for x in b_ind if x not in B.dtype.names]
        raise ValueError(
            repr(problem_index)
            + " not in second argument, which has fields"
            + repr(B.dtype.names)
            + "!!!"
        )
    if np.any([x not in A.dtype.names for x in a_ind]):
        problem_index = [x for x in a_ind if x not in A.dtype.names]
        raise ValueError(
            repr(problem_index)
            + " not in first argument, which has fields"
            + repr(A.dtype.names)
            + "!!!"
        )
    # }}}
    B_reduced = B[b_ind]  # a version of B reduced to only include the keys
    B_reduced = reorder_rec(
        B_reduced, b_ind
    )  # again, because it doesn't do this just based on the indexing
    A_reduced = A[a_ind]  # same for A
    A_reduced = reorder_rec(
        A_reduced, a_ind
    )  # again, because it doesn't do this just based on the indexing
    # now, I need to generate a mapping from the b_ind to a_ind
    field_mapping = dict(list(zip(b_ind, a_ind)))
    # now I change the names so they match and I can compare them
    B_reduced.dtype.names = tuple(
        [field_mapping[x] for x in B_reduced.dtype.names]
    )
    # {{{ now find the list of indices for B that match each value of A
    old_B_reduced_names, old_B_reduced_types = tuple(
        zip(*tuple(B_reduced.dtype.descr))
    )
    B_reduced.dtype = np.dtype(
        list(zip(A_reduced.dtype.names, old_B_reduced_types))
    )
    if A_reduced.dtype != B_reduced.dtype:
        B_reduced.dtype = np.dtype(
            list(zip(old_B_reduced_names, old_B_reduced_types))
        )
        raise TypeError(
            strm(
                "The datatype of A_reduced=",
                A_reduced.dtype,
                "and B_reduced=",
                B_reduced.dtype,
                "are not the same,  which is going to create problems!",
            )
        )
    try:
        list_of_matching = [np.nonzero(B_reduced == j)[0] for j in A_reduced]
    except Exception:
        raise RuntimeError(
            strm(
                "When trying to decorate,  A_reduced=",
                A_reduced,
                "with B_reduced=",
                B_reduced,
                "one or more of the following is an np.empty tuple,  which is"
                " wrong!:",
                [np.nonzero(B_reduced == j) for j in A_reduced],
            )
        )
    logger.debug(
        strm("(decorate\\_rec):: original list of matching", list_of_matching)
    )
    length_of_matching = np.array([len(j) for j in list_of_matching])
    logger.debug(
        strm("(decorate\\_rec):: length of matching is", length_of_matching)
    )
    if np.any(length_of_matching == 0):
        if drop_rows:
            if drop_rows == "return":
                dropped_rows = A[length_of_matching == 0].copy()
            else:
                dropped_rows = A_reduced[length_of_matching == 0]
                print(
                    r"{\color{red}Warning! decorate\_rec dropped fields"
                    r" in the first argument",
                    lsafen(
                        repr(
                            list(
                                zip(
                                    A_reduced.dtype.names * len(dropped_rows),
                                    dropped_rows.tolist(),
                                )
                            )
                        )
                    ),
                    r"}",
                )
            # {{{ now, remove all trace of the dropped fields
            A = A[length_of_matching != 0]
            list_of_matching = [j for j in list_of_matching if len(j) > 0]
            length_of_matching = [len(j) for j in list_of_matching]
            # }}}
        else:
            raise ValueError(
                strm(
                    "There is no data in the second argument that has",
                    b_ind,
                    "fields to match the",
                    a_ind,
                    "fields of the first argument for the following records:",
                    A_reduced[length_of_matching == 0],
                    "if this is correct, you can set the drop_rows = True",
                    "keyword argument to drop these fields",
                )
            )
    # now, do a neat trick of stackoverflow to collapse a nested list
    # this gives just the indices in B that match the values of A
    list_of_matching = [j for i in list_of_matching for j in i]
    # }}}
    logger.debug(
        strm("(decorate\\_rec):: list of matching is", list_of_matching)
    )
    # now grab the data for these rows
    add_data = B[list_of_matching]
    # {{{ finally, smoosh the two sets of data together
    # {{{ Now, I need to replicate the rows that have multiple matchesjk
    if np.any(length_of_matching > 1):
        index_replication_vector = [
            k
            for j in range(0, len(length_of_matching))
            for k in [j] * length_of_matching[j]
        ]
        retval = A[index_replication_vector]
    else:
        retval = A.copy()
    # }}}
    # {{{ add the new fields
    new_dtypes = [j for j in B.dtype.descr if j[0] not in A.dtype.names]
    logger.debug(strm("(decorate\\_rec):: new dtypes:", repr(new_dtypes)))
    try:
        retval = newcol_rec(retval, new_dtypes)
    except Exception:
        raise ValueError(
            strm(
                "Problem trying to add new columns with the dtypes", new_dtypes
            )
        )
    # }}}
    logger.debug(strm("(decorate\\_rec):: add data:", repr(add_data)))
    for name in np.dtype(new_dtypes).names:
        logger.debug(
            strm(
                "(decorate\\_rec):: trying to add data for",
                name,
                ":",
                add_data[name][:],
            )
        )
        retval[name][:] = add_data[name][:]
    # }}}
    if drop_rows == "return":
        return retval, dropped_rows
    else:
        return retval


def newcol_rec(A, new_dtypes):
    r"""add new, np.empty (i.e. random numbers) fields to A, as given by
    new_dtypes --> note that there are deeply nested numpy functions to do
    this, but the options are confusing, and I think the way these work is
    efficient
    """
    if isinstance(new_dtypes, np.dtype):
        new_dtypes = new_dtypes.descr
    elif isinstance(new_dtypes, tuple):
        new_dtypes = [new_dtypes]
    elif isinstance(new_dtypes, list):
        if not isinstance(new_dtypes[0], tuple):
            new_dtypes = [tuple(new_dtypes)]
    retval = np.empty(A.shape, dtype=A.dtype.descr + new_dtypes)
    for name in A.dtype.names:
        retval[name][:] = A[name][:]
    return retval


def applyto_rec(myfunc, myarray, mylist):
    r"""apply myfunc to myarray with the intention of collapsing it to a
    smaller number of values"""
    if not isinstance(mylist, list) and isinstance(mylist, str):
        mylist = [mylist]
    combined = []
    j = 0
    # {{{ make the list "combined", which I later concatenate
    while len(myarray) > 0:
        thisitem = myarray[0]  # always grab the first row of what's left
        # {{{ initialize the np.empty new row
        if j == 0:
            newrow = thisitem.reshape(1)
        newrow = newrow.copy()
        # }}}
        # {{{ make a mask for all items that are identified as the same data
        # and copy the identical data to newrow
        mask = myarray[mylist[0]] == thisitem[mylist[0]]
        newrow[mylist[0]] = thisitem[mylist[0]]
        for k in range(1, len(mylist)):
            mask &= myarray[mylist[k]] == thisitem[mylist[k]]
            newrow[mylist[k]] = thisitem[mylist[k]]
        # }}}
        logger.debug(
            strm(
                lsafen(
                    "(applyto np.core.rec): for row %d, I select these:" % j
                )
            )
        )
        myarray_subset = myarray[mask]
        logger.debug(
            strm(lsafen("(applyto np.core.rec): ", repr(myarray_subset)))
        )
        other_fields = set(mylist) ^ set(thisitem.dtype.names)
        logger.debug(
            strm(
                lsafen(
                    "(applyto np.core.rec): other fields are:", other_fields
                )
            )
        )
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = myfunc(myarray_subset[thisfield])
            except Exception:
                raise ValueError(
                    strm(
                        "error in applyto_rec:  You usually get this",
                        "when one of the fields that you have NOT passed"
                        " in the",
                        "second argument is a string.  The fields and types",
                        "are:",
                        repr(myarray_subset.dtype.descr),
                    )
                )
        logger.debug(
            strm(
                lsafen(
                    "(applyto np.core.rec): for row %d, I get this as a"
                    " result:" % j,
                    newrow,
                )
            )
        )
        combined.append(newrow)  # add this row to the list
        myarray = myarray[
            ~mask
        ]  # mask out everything I have used from the original matrix
        logger.debug(
            strm(
                lsafen(
                    "(applyto np.core.rec): the np.array is now", repr(myarray)
                )
            )
        )
        j += 1
    # }}}
    combined = np.concatenate(combined)
    logger.debug(
        strm(
            lsafen(
                "(applyto np.core.rec): final result",
                repr(combined),
                "has length",
                len(combined),
            )
        )
    )
    return combined


def meanstd_rec(myarray, mylist, standard_error=False):
    r"""this is something like applyto_rec, except that it applies the mean and
    creates new rows for the "error," where it puts the standard deviation"""
    if not isinstance(mylist, list) and isinstance(mylist, str):
        mylist = [mylist]
    combined = []
    other_fields = set(mylist) ^ set(myarray.dtype.names)
    logger.debug(strm("(meanstd_rec): other fields are", lsafen(other_fields)))
    newrow_dtype = [
        [j, ("%s_ERROR" % j[0],) + j[1:]] if j[0] in other_fields else [j]
        for j in myarray.dtype.descr
    ]
    newrow_dtype = [k for j in newrow_dtype for k in j]
    logger.debug(
        strm(lsafen("(meanstd np.core.rec): other fields are:", other_fields))
    )
    # {{{ make the list "combined", which I later concatenate
    j = 0
    while len(myarray) > 0:
        thisitem = myarray[0]  # always grab the first row of what's left
        # {{{ initialize the np.empty new row
        newrow = np.zeros(1, dtype=newrow_dtype)
        # }}}
        # {{{ make a mask for all items that are identified as the same data
        # and copy the identical data to newrow
        mask = myarray[mylist[0]] == thisitem[mylist[0]]
        newrow[mylist[0]] = thisitem[mylist[0]]
        for k in range(1, len(mylist)):
            mask &= myarray[mylist[k]] == thisitem[mylist[k]]
            newrow[mylist[k]] = thisitem[mylist[k]]
        # }}}
        logger.debug(
            strm(
                lsafen(
                    "(meanstd np.core.rec): for row %d, I select these:" % j
                )
            )
        )
        myarray_subset = myarray[mask]
        logger.debug(
            strm(lsafen("(meanstd np.core.rec): ", repr(myarray_subset)))
        )
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = np.mean(myarray_subset[thisfield])
                if standard_error:
                    newrow[thisfield + "_ERROR"] = np.std(
                        myarray_subset[thisfield]
                    ) / sqrt(len(myarray_subset[thisfield]))
                else:
                    newrow[thisfield + "_ERROR"] = np.std(
                        myarray_subset[thisfield]
                    )
            except Exception:
                raise RuntimeError(
                    "error in meanstd_rec:  You usually get this",
                    "when one of the fields that you have NOT passed in the",
                    "second argument is a string.  The fields and types",
                    "are:",
                    repr(myarray_subset.dtype.descr),
                )
        logger.debug(
            strm(
                lsafen(
                    "(meanstd np.core.rec): for row %d, I get this as a"
                    " result:" % j,
                    newrow,
                )
            )
        )
        combined.append(newrow)  # add this row to the list
        myarray = myarray[
            ~mask
        ]  # mask out everything I have used from the original matrix
        logger.debug(
            strm(
                lsafen(
                    "(meanstd np.core.rec): the np.array is now", repr(myarray)
                )
            )
        )
        j += 1
    # }}}
    combined = np.concatenate(combined)
    logger.debug(
        strm(
            lsafen(
                "(meanstd np.core.rec): final result",
                repr(combined),
                "has length",
                len(combined),
            )
        )
    )
    return combined


def make_rec(*args, **kwargs):
    r"""input,names or a single argument, which is a dictionary
    strlen = 100 gives length of the strings (which need to be specified in
    record arrays) you can also specify (especially useful with the dictionary
    format) the list order = [str1,str2,...] which orders the output records
    with the field containing str1 first, then the field containing str2, then
    any remaining fields"""
    strlen, order, zeros_like = process_kwargs(
        [("strlen", 100), ("order", None), ("zeros_like", False)], kwargs
    )
    if len(args) == 1 and (isinstance(args[0], dict)):
        names = list(args[0].keys())
        input = list(args[0].values())
    elif len(args) == 2:
        input = args[0]
        names = args[1]
    else:
        raise ValueError(
            strm(
                "I don't understand the arguments you passed to",
                "make_rec!!!\nshould be (list of values, list of field"
                " names),",
                "or a dictionary",
            )
        )
    # {{{ apply the order kwarg
    if order is not None:
        newindices = []
        for orderitem in order:
            newindices += [
                j
                for j, k in enumerate(names)
                if (k.find(orderitem) > -1 and j not in newindices)
            ]
        newindices += [j for j, k in enumerate(names) if j not in newindices]
        names = [names[j] for j in newindices]
        input = [input[j] for j in newindices]
    # }}}
    if not (isinstance(input, list) and isinstance(names, list)):
        raise TypeError("you must enter a list for both")
    types = list(map(type, input))
    shapes = list(map(np.shape, input))
    if all([j == shapes[0] for j in shapes]):
        if shapes[0] == ():  # if it's one dimensional
            equal_shapes = False
            shapes = [(1)] * len(shapes)
        else:
            equal_shapes = True
            shape_of_array = shapes[0]
            shapes = [()] * len(shapes)
    else:
        equal_shapes = False
    for j, k in enumerate(input):
        if isinstance(k, list) and equal_shapes:
            k = k[0]
        if isinstance(k, str):
            types[j] = "|S%d" % strlen
        if isinstance(k, np.ndarray):
            types[j] = k.dtype
    try:
        mydtype = np.dtype(list(zip(names, types, shapes)))
    except Exception:
        raise ValueError(
            strm(
                "problem trying to make names",
                names,
                " types",
                types,
                "shapes",
                shapes,
            )
        )
    if np.zeros_like:
        retval = np.zeros(np.zeros_like, dtype=mydtype)
        return retval
    if equal_shapes:
        retval = np.empty(shape_of_array, dtype=mydtype)
        for j, thisname in enumerate(names):
            try:
                retval[thisname][:] = input[j][:]
            except Exception:
                raise RuntimeError(
                    "error trying to load input for '"
                    + thisname
                    + "' of shape "
                    + repr(np.shape(input[j]))
                    + " into retval field of shape "
                    + repr(np.shape(retval[thisname]))
                )
        return retval
    else:
        try:
            return np.array([tuple(input)], dtype=mydtype)
        except Exception:
            raise ValueError(
                strm(
                    "problem trying to assign data of type",
                    list(map(type, input)),
                    "\nvalues",
                    input,
                    "\nonto",
                    mydtype,
                    "\ndtype made from tuple:",
                    list(zip(names, types, shapes)),
                )
            )


# }}}
def maprep(*mylist):
    mylist = list(mylist)
    for j in range(0, len(mylist)):
        if not isinstance(mylist[j], str):
            mylist[j] = mylist[j].__repr__()
    return " ".join(mylist)


# {{{ plot wrapper
global myplotfunc
myplotfunc = plt.plot


def plot(*args, **kwargs):
    """The base plotting function that wraps around matplotlib to do a couple
    convenient things.

    Parameters
    ----------
    label_format_string: str
        If supplied, it formats the values of the other dimension to turn them
        into a label string.
    human_units: bool
    """
    global myplotfunc
    has_labels = False
    # {{{ deal with axes and some other kwargs
    (
        ax,
        human_units,
        label_format_string,
        normalize,
        noerr,
        longest_is_x,
    ) = process_kwargs(
        [
            ("ax", plt.gca()),
            ("human_units", False),
            ("label_format_string", None),
            ("normalize", False),
            ("noerr", False),
            ("longest_is_x", True),
        ],
        kwargs,
        pass_through=True,
    )
    # }}}
    myplotfunc = ax.plot  # default
    # {{{ all possible properties
    myformat = None
    myxlabel = None
    myylabel = None
    myx = None
    myy = None
    # }}}
    # {{{assign all the possible combinations
    if len(args) == 1:
        myy = args[0]
    elif (len(args) == 2) and (isinstance(args[1], str)):
        myy = args[0]
        myformat = args[1]
    else:
        myx = args[0]
        myy = args[1]
    if len(args) == 3:
        myformat = args[2]
    if np.isscalar(myx):
        myx = np.array([myx])
    if np.isscalar(myy):
        myy = np.array([myy])
    # }}}
    x_inverted = False
    # {{{ parse nddata
    if isinstance(myy, nddata):
        myy = myy.copy()
        if myy.get_error() is not None:
            logging.debug(
                strm(
                    "shapes at top of function",
                    ndshape(myy),
                    myy.data.shape,
                    myy.data_error.shape,
                )
            )
        # {{{ automatically reduce any singleton dimensions
        if not len(myy.dimlabels) == 1:
            if np.any(np.array(myy.data.shape) == 1):
                for singleton_dim in [
                    lb
                    for j, lb in enumerate(myy.dimlabels)
                    if myy.data.shape[j] == 1
                ]:
                    myy = myy[singleton_dim, 0]
        # }}}
        if len(myy.data.shape) > 1 and longest_is_x:
            longest_dim = np.argmax(myy.data.shape)
            all_but_longest = set(range(len(myy.data.shape))) ^ set((
                longest_dim,
            ))
            if len(all_but_longest) > 0:
                last_not_longest = max(all_but_longest)
            else:
                last_not_longest = -1
            all_but_longest = list(
                all_but_longest
            )  # seems to be sorted by default
        else:
            longest_dim = 0  # treat first as x, like before
            last_not_longest = -1
            if len(myy.data.shape) > 1:
                all_but_longest = set(range(len(myy.data.shape))) ^ set((
                    longest_dim,
                ))
                all_but_longest = list(all_but_longest)
            else:
                all_but_longest = []
        if human_units:
            myy = myy.human_units()
        if myy.get_plot_color() is not None and "color" not in list(
            kwargs.keys()
        ):  # allow override
            kwargs.update({"color": myy.get_plot_color()})
        if myy.name() is not None:
            myylabel = myy.name()
        else:
            myylabel = "data"
        myylabel = myy.unitify_axis(myylabel, is_axis=False)
        if len(myy.dimlabels) > 0:
            myxlabel = myy.unitify_axis(longest_dim)
        if myx is None:
            try:
                myx = myy.getaxis(myy.dimlabels[longest_dim])
            except Exception:
                if len(myy.data.shape) == 0:
                    raise ValueError(
                        "I can't plot zero-dimensional data (typically arises"
                        " when you have a dataset with one point)"
                    )
                myx = r_[0 : myy.data.shape[longest_dim]]
        if (
            not noerr
            and isinstance(myy.data_error, np.ndarray)
            and len(myy.data_error) > 0
        ):  # then this should be an errorbar plot

            def thiserrbarplot(*tebargs, **tebkwargs):
                if "capsize" not in tebkwargs:
                    tebkwargs.update({"capsize": 6})
                if isinstance(tebargs[-1], str):
                    tebkwargs.update({"fmt": tebargs[-1]})
                    return ax.errorbar(*tebargs[:-1], **tebkwargs)
                else:
                    return ax.errorbar(*tebargs, **tebkwargs)

            myplotfunc = thiserrbarplot
            logger.debug("this is an errorbar plot")
            # {{{ pop any singleton dims
            myyerror = myy.get_error()
            myyerror = np.squeeze(myyerror)
            # }}}
            kwargs.update({"yerr": None})
            valueforxerr = myy.get_error(myy.dimlabels[longest_dim])
            if valueforxerr is not None:  # if we have x errorbars too
                # print "DEBUG decided to assign to xerr:",valueforxerr
                kwargs.update({"xerr": valueforxerr})
            logging.debug(
                strm(
                    "shapes after splitting nddata",
                    myy.data.shape,
                    myyerror.shape,
                )
            )
        # {{{ deal with axis labels along y
        try:
            yaxislabels = myy.getaxis(myy.dimlabels[last_not_longest])
        except Exception:
            yaxislabels = None
        # at this point, if there is no axis label, it will break and go to
        # pass
        if yaxislabels is not None:
            if len(yaxislabels) > 0:
                if isinstance(yaxislabels[0], np.string_):
                    has_labels = True
                elif label_format_string is not None:
                    yaxislabels = [
                        label_format_string % j for j in yaxislabels
                    ]
                    has_labels = True
        # }}}
        # {{{ add label if name is present, and squeeze -- could do this
        #     instead of ylabel, above
        if myy.get_prop("x_inverted"):
            x_inverted = True
        # myy_name = myy.name()
        if len(myy.data.shape) == 1:
            myy = myy.data
        else:
            myy = np.squeeze(
                myy.data.transpose([longest_dim] + all_but_longest)
            )
            if "yerr" in kwargs.keys():
                myyerror = np.squeeze(
                    myyerror.transpose([longest_dim] + all_but_longest)
                )
        if "yerr" in kwargs.keys():
            logging.debug(
                strm("halfway checkpoint", myy.shape, myyerror.shape)
            )
        # if ((len(myy.data) == 1 and 'label' not in kwargs.keys())
        #     and (myy_name is not None)):
        #    kwargs.update('label',myy_name)
        # }}}
    # }}}
    # {{{ allow list arguments
    if type(myy) is list:
        myy = np.array(myy)
    if type(myx) is list:
        myx = np.array(myx)
    # }}}
    # {{{ semilog where appropriate
    if (
        (myx is not None) and (len(myx) > 1) and all(myx > 0.0)
    ):  # by doing this and making myplotfunc global, we preserve the plot
        # style if we want to tack on one point
        try:
            b = np.diff(np.log10(myx))
        except Exception:
            raise Exception(
                strm(
                    "likely a problem with the type of the x label, which is",
                    myx,
                )
            )
        if (
            (np.size(b) > 3)
            and all(abs((b - b[0]) / b[0]) < 1e-4)
            and not ("nosemilog" in list(kwargs.keys()))
        ):
            if "plottype" not in list(kwargs.keys()):
                myplotfunc = ax.semilogx
    if "nosemilog" in list(kwargs.keys()):
        # print 'this should pop nosemilog'
        kwargs.pop("nosemilog")
    if "yerr" in kwargs.keys():
        logging.debug(
            strm("halfway checkpoint", myy.data.shape, myyerror.shape)
        )
    if "plottype" in list(kwargs.keys()):
        if kwargs["plottype"] == "semilogy":
            myplotfunc = ax.semilogy
        elif kwargs["plottype"] == "semilogx":
            myplotfunc = ax.semilogx
        elif kwargs["plottype"] == "loglog":
            myplotfunc = ax.loglog
        elif kwargs["plottype"] == "linear":
            myplotfunc = ax.plot
        else:
            raise ValueError(
                strm("plot type", kwargs["plottype"], "not allowed!")
            )
        kwargs.pop("plottype")
    # }}}
    # {{{ take care of manual colors
    if myformat is not None:
        colorpos = myformat.find("#")
        if colorpos > -1:
            kwargs.update({"color": myformat[colorpos : colorpos + 7]})
            myformat = myformat[0:colorpos] + myformat[colorpos + 7 :]
        ##kwargs.update({'fmt':myformat})
        linematched = False
        for linestyle in ["-", "--", "-.", ":", "None", "  "]:
            if linestyle in myformat:
                linematched = True
                myformat = myformat.replace(linestyle, "")
                kwargs.update({"linestyle": linestyle})
        for markerlabel in ["o", ".", "d"]:
            if markerlabel in myformat:
                if not linematched:
                    kwargs.update({"linestyle": ""})
                myformat = myformat.replace(markerlabel, "")
                kwargs.update({"marker": markerlabel})
        if len(myformat) == 0:
            myformat = None
    # }}}
    if normalize is not None and normalize:
        myy /= myy.max()
    # {{{ hsv plots when we have multiple lines
    if len(np.shape(myy.squeeze())) > 1 and np.sum(
        np.array(np.shape(myy)) > 1
    ):
        # {{{ hsv plots
        retval = []
        if "yerr" in kwargs.keys():
            logger.debug("I see error")
        else:
            logger.debug("I do not see error")
        for j in range(0, myy.shape[1]):
            # {{{ this is the way to assign plot arguments
            plotargs = [k for k in (myx, myy[:, j], myformat) if k is not None]
            # }}}
            if "yerr" in kwargs.keys():
                kwargs["yerr"] = myyerror[:, j]
            # {{{ here, i update the kwargs to include the specific color for
            #     this line
            newkwargs = kwargs.copy()  # kwargs is a dict
            newkwargs.update(
                {"color": plt.cm.hsv(np.double(j) / np.double(myy.shape[1]))}
            )
            # }}}
            # {{{ here, I update to use the labels
            if has_labels:
                newkwargs.update({"label": yaxislabels[j]})
            # }}}
            if "yerr" in newkwargs.keys():
                logging.debug(
                    strm(
                        "shapes before plot",
                        plotargs[0].shape,
                        plotargs[1].shape,
                        newkwargs["yerr"].shape,
                    )
                )
                # myplotfunc = ax.plot
                # newkwargs.pop('yerr')
            elif len(plotargs) > 1 and isinstance(plotargs[1], np.ndarray):
                logging.debug(
                    strm(
                        "shapes before plot",
                        plotargs[0].shape,
                        plotargs[1].shape,
                    )
                )
            if np.any(np.isinf(myy)):
                myy[np.isinf(myy)] = (
                    np.nan
                )  # added this to prevent an overflow error
            try:
                retval += [myplotfunc(*tuple(plotargs), **newkwargs)]
            except Exception:
                raise RuntimeError(
                    strm(
                        "Error trying to plot using function",
                        myplotfunc,
                        "\nwith",
                        len(plotargs),
                        "arguments",
                        "\nwhich were\n",
                        plotargs,
                        "\nand had len\n",
                        list(map(len, plotargs)),
                        "and",
                        len(newkwargs),
                        "\noptions",
                        newkwargs,
                        "of len",
                        ", ".join([
                            (
                                str(type(j)) + " " + str(j)
                                if np.isscalar(j)
                                else str(len(j))
                            )
                            for j in list(newkwargs.values())
                        ]),
                    )
                )
            if x_inverted:
                these_xlims = ax.get_xlim()
                ax.set_xlim((max(these_xlims), min(these_xlims)))
        # }}}
        # }}}
    else:
        logger.debug(strm("here are the kwargs", kwargs))
        if "yerr" in kwargs.keys() and kwargs["yerr"] is None:
            kwargs["yerr"] = myyerror
        plotargs = [j for j in [myx, np.real(myy), myformat] if j is not None]
        try:
            logger.debug(
                strm(
                    "plotting with args",
                    plotargs,
                    "(myformat is",
                    myformat,
                    ") and kwargs",
                    kwargs,
                )
            )
            retval = myplotfunc(*plotargs, **kwargs)
        except Exception:
            raise RuntimeError(
                strm(
                    "error trying to plot",
                    type(myplotfunc),
                    "with value",
                    myplotfunc,
                    "\nlength of the np.ndarray arguments:",
                    [
                        (
                            "shape:" + str(np.shape(j))
                            if isinstance(j, np.ndarray)
                            else j
                        )
                        for j in plotargs
                    ],
                    "\nsizes of np.ndarray kwargs",
                    dict([
                        (
                            (j, np.shape(kwargs[j]))
                            if isinstance(kwargs[j], np.ndarray)
                            else (j, kwargs[j])
                        )
                        for j in list(kwargs.keys())
                    ]),
                    "\narguments = ",
                    plotargs,
                    "\nkwargs =",
                    kwargs,
                )
            )
        if x_inverted:
            these_xlims = ax.get_xlim()
            ax.set_xlim((max(these_xlims), min(these_xlims)))
    # {{{ attach labels and such
    if myxlabel is not None:
        ax.set_xlabel(myxlabel)
    if myylabel is not None:
        ax.set_ylabel(myylabel)
    try:
        ax.axis("tight")
    except Exception:
        raise Exception(
            strm(
                "error trying to set axis tight after plot",
                myplotfunc,
                "with arguments",
                plotargs,
                "and kwargs",
                kwargs,
                "\nsizes of arguments:",
                [np.shape(j) for j in plotargs],
                "\nsizes of np.ndarray kwargs:",
                dict([
                    (j, np.shape(kwargs[j]))
                    for j in list(kwargs.keys())
                    if isinstance(kwargs[j], np.ndarray)
                ]),
            )
        )
    # plt.grid(True)
    # }}}
    return retval


# }}}


# {{{ concatenate datalist along dimname
def concat(datalist, dimname, chop=False):
    """concatenate multiple datasets together along a new dimension.

    Parameters
    ==========
    datalist: list of nddata
        the data you want to concatenate -- they must have the same ndshape!
    dimname: str
        name of the new dimension
    """
    # {{{ allocate a new datalist structure
    newdimsize = 0
    logging.debug("type(datalist) " + str(type(datalist)))
    try:
        shapes = list(map(ndshape, datalist))
    except Exception:
        if not isinstance(datalist, list):
            raise TypeError(
                strm("You didn't pass a list, you passed a", type(datalist))
            )
        raise RuntimeError(
            strm(
                "Problem with what you passed to concat, list of types,",
                list(map(type, datalist)),
            )
        )
    other_info_out = datalist[0].other_info
    if dimname in datalist[0].dimlabels:
        dim_idx = datalist[0].axn(dimname)
        assert all([
            datalist[j].dimlabels == datalist[0].dimlabels
            for j in range(len(datalist))
        ]), (
            "the dimlabels for all your datasets do no match and/or are not"
            " ordered the same way"
        )
        # {{{ check that all the shapes match, too, aside from the dim we
        #     concat along
        shape_check = [
            list(datalist[j].data.shape) for j in range(len(datalist))
        ]
        for j in shape_check:
            j.pop(dim_idx)
        assert all(
            [shape_check[j] == shape_check[0] for j in range(len(shape_check))]
        ), ("shapes " + str(shape_check) + " are not all equal")
        # }}}
        # {{{ concatenate the data ndarrays and dimname axis ndarrays
        concated_data = np.concatenate(
            tuple(datalist[j].data for j in range(len(datalist))),
            axis=dim_idx,
        )
        concated_ax_coords = np.concatenate(
            tuple(datalist[j][dimname] for j in range(len(datalist)))
        )
        # }}}
        retval = datalist[0].copy(data=False)
        retval.axis_coords[dim_idx] = concated_ax_coords
        retval.data = concated_data
    else:
        for j in range(0, len(datalist)):
            # {{{ make list for the shape to check, which contains the
            #     dimensions we are NOT concatting along
            newdimsize += 1
            shapetocheck = list(shapes[j].shape)
            # }}}
            if j == 0:
                shapetocheckagainst = shapetocheck
            else:
                if np.any(
                    ~(np.array(shapetocheck) == np.array(shapetocheckagainst))
                ):
                    if chop:
                        logger.debug(
                            strm(
                                repr(shapetocheck),
                                lsafen(repr(shapetocheckagainst)),
                            )
                        )
                        raise ValueError(
                            strm(
                                "For item ",
                                j,
                                "in concat, ",
                                shapetocheck,
                                "!=",
                                shapetocheckagainst,
                                "where all the shapes of the things",
                                "you're trying to concat are:",
                                shapes,
                            )
                        )
                    else:
                        raise ValueError(
                            strm(
                                "For item ",
                                j,
                                "in concat, ",
                                shapetocheck,
                                "!=",
                                shapetocheckagainst,
                                "where all the shapes of the things you're"
                                " trying to concat are:",
                                shapes,
                            )
                        )
        retval = ndshape(datalist[-1])
        retval += ([newdimsize], [dimname])
        try:
            retval = retval.alloc()
        except Exception:
            raise ValueError(
                strm(
                    "trying to alloc the retval",
                    retval,
                    "created a problem",
                )
            )
        if datalist[0].get_error() is not None:
            retval.set_error(np.zeros(np.shape(retval.data)))
        # }}}
        # {{{ actually contract the datalist
        newdimsize = 0  # now use it to track to position
        for j in range(0, len(datalist)):
            retval[dimname, newdimsize : newdimsize + 1] = datalist[j]
            newdimsize += 1
        # }}}
        # {{{ pull the axis labels from the last item in the list
        if len(datalist[-1].axis_coords) > 0:
            dimlabels = list(datalist[-1].dimlabels)
            axis_coords = list(datalist[-1].axis_coords)
            # print "axis_coords are",axis_coords,"for",dimlabels
            dimlabels += [dimname]
            axis_coords += [r_[0:newdimsize]]
            try:
                retval.labels(dimlabels, axis_coords)
            except Exception:
                raise ValueError(
                    strm(
                        "trying to attach axes of lengths",
                        list(map(len, axis_coords)),
                        "to",
                        dimlabels,
                    )
                )
        # }}}
    retval.other_info = other_info_out
    return retval


# }}}
class nddata(object):
    """This is the detailed API reference.
    For an introduction on how to use ND-Data, see the
    :ref:`Main ND-Data Documentation <nddata-summary-label>`.
    """

    want_to_prospa_decim_correct = False

    # {{{ initialization
    def __init__(self, *args, **kwargs):
        """initialize nddata -- several options.
        Depending on the information available, one of several formats can be
        used.

        3 arguments:
            ``nddata(inputarray, shape, dimlabels)``

            :inputarray:
                np.ndarray storing the data -- note that the size is ignored
                and the data is reshaped as needed
            :shape:
                a list (or np.array, *etc.*) giving the size of each dimension,
                in order
            :dimlabels:
                a list giving the names of each dimension, in order
        2 arguments:
            ``nddata(inputarray, dimlabels)``

            :inputarray:
                np.ndarray storing the data -- the data is *not* reshaped
            :dimlabels:
                a list giving the names of each dimension, in order
        2 arguments:
            ``nddata(inputarray, single_dimlabel)``

            :inputarray:
                np.ndarray storing the data -- must be 1D
                inputarray is *also* used to label the single axis
            :single_dimlabel:
                a list giving the name of the single axis
        1 argument:
            ``nddata(inputarray, shape, dimlabels)``

            :inputarray:
                np.ndarray storing the data -- reduced to 1D
                A single dimension, called "INDEX" is set.
                This suppresses the printing of axis labels.
                This is used to store numbers and arrays
                that might have error and units,
                but aren't gridded data.
        keyword args
            these can be used to set the labels, etc, and are passed to
            :func:`__my_init__`

        """
        logger.debug("called init")
        if len(args) > 1:
            logger.debug("more than one argument -- args: " + strm(args))
            if isinstance(args[0], nddata):
                raise ValueError(
                    "you can't initialize an nddata from another nddata!!!"
                )
            if len(args) == 2:
                if len(args[0].shape) == 1 and isinstance(args[1], str):
                    logger.debug("constructing 1D np.array")
                    self.__my_init__(args[0], [len(args[0])], [args[1]])
                    self.labels(
                        args[1], args[0].copy()
                    )  # needs to be a copy, or when we write data, we will
                    #    change the axis
                elif all([isinstance(j, str) for j in args[1]]):
                    logger.debug("passed only axis labels")
                    self.__my_init__(args[0], list(args[0].shape), args[1])
                else:
                    raise ValueError(
                        "You can pass two arguments only if you pass a 1d"
                        " np.ndarray and a name for the axis"
                    )
            elif len(args) == 3:
                self.__my_init__(args[0], args[1], args[2], **kwargs)
            else:
                raise ValueError(
                    strm(
                        "You passed",
                        len(args),
                        "to nddata.  I don't know what to do with this.",
                    )
                )
        else:
            logger.debug("only one argument")
            self.__my_init__(args[0], [-1], ["INDEX"], **kwargs)
        return

    def __my_init__(
        self,
        data,
        sizes,
        dimlabels,
        axis_coords=[],
        ft_start_time=None,
        data_error=None,
        axis_coords_error=None,
        axis_coords_units=None,
        data_units=None,
        other_info={},
    ):
        if ft_start_time is not None:
            raise ValueError(
                "ft_start_time is obsolete -- you will want to pass a float"
                " value to the shift keyword argument of either .ft() or"
                " .ift()"
            )
        self.genftpairs = False
        if not (isinstance(data, np.ndarray)):
            if (
                np.isscalar(data)
                or (isinstance(data, list))
                or (isinstance(data, tuple))
            ):
                data = np.array(data)
            else:
                raise TypeError(
                    strm("data is not an np.array, it's", type(data), "!")
                )
        if not (isinstance(dimlabels, list)):
            raise TypeError(
                strm(
                    "you provided a multi-dimensional np.ndarray but a set of"
                    " dimension labels of type",
                    type(dimlabels),
                    "if you want a 1D nddata, give a 1D np.array, or if you"
                    " want a ND nddata, give a list of dimensions",
                )
            )
        try:
            self.data = np.reshape(data, sizes)
        except Exception:
            try:
                error_string = strm(
                    "While initializing nddata, you are trying trying to"
                    " reshape a",
                    data.shape,
                    "array (",
                    data.size,
                    "data elements) with list of sizes",
                    list(zip(dimlabels, sizes)),
                    "(implying that there are ",
                    np.prod(sizes),
                    "data elements)",
                )
            except TypeError:
                error_string = strm(
                    "While initializing nddata, you are trying trying to"
                    " reshape a",
                    data.shape,
                    "array (",
                    data.size,
                    "data elements) with list of sizes",
                    sizes,
                )
            raise ValueError(error_string)
        self.dimlabels = dimlabels
        self.axis_coords = axis_coords
        self.data_error = data_error
        self.data_units = data_units
        self.other_info = deepcopy(other_info)
        if axis_coords_error is None:
            self.axis_coords_error = [None] * len(axis_coords)
        else:
            self.axis_coords_error = axis_coords_error
        if axis_coords_units is None:
            self.axis_coords_units = [None] * len(axis_coords)
        else:
            self.axis_coords_units = axis_coords_units
        return

    # }}}
    def _contains_symbolic(self, string):
        return string[:9] == "symbolic_" and hasattr(self, string)

    # {{{ for printing
    def __repr_pretty__(self, p, cycle):
        if cycle:
            p.text("...")
        else:
            p.text(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self):
        def show_array(x, indent=""):
            x = repr(x)
            if x.startswith("np.array("):
                x = x.split("\n")
                # need to remove the "np.array(" and aligning spaces
                return ("\n" + indent).join(j[6:-1] for j in x)
            else:
                return x

        retval = show_array(self.data)
        retval += "\n\t\t+/-"
        retval += show_array(self.get_error())
        if (
            len(self.dimlabels) > 1
            or len(self.dimlabels) == 0
            or self.dimlabels[0] != "INDEX"
        ):
            retval += "\n\tdimlabels="
            retval += repr(self.dimlabels)
            retval += "\n\taxes="

            def rep_this_dict(starting_indent, thisdict, errordict):
                dictrep = []
                for k, v in thisdict.items():
                    dictrep.append(
                        "`"
                        + k
                        + "':"
                        + show_array(v, starting_indent)
                        + starting_indent
                        + "\t\t+/-"
                        + repr(errordict[k])
                    )
                return (
                    "{" + ("," + starting_indent + "\t").join(dictrep) + "}"
                )  # separate with an extra comma, the existing indent, and a
                #    tab

            retval += rep_this_dict(
                "\n\t",
                self.mkd(self.axis_coords),
                self.mkd(self.axis_coords_error),
            )
        # retval += '\n\t\t+/-'
        # retval += rep_this_dict('\n\t\t',self.mkd(self.axis_coords_error))
        retval += "\n"
        return retval

    # }}}
    # {{{ for plotting
    def gnuplot_save(self, filename):
        x = self.getaxis(self.dimlabels[0])[:5]
        y = self.getaxis(self.dimlabels[1])[:5]
        z = self.data[:5, :5]
        print(
            "size of x",
            np.size(x),
            "size of y",
            np.size(y),
            "size of z",
            np.size(z),
        )
        print("x", x, "y", y, "z", z)
        data = np.empty((z.shape[0] + 1, z.shape[1] + 1))
        data[1:, 1:] = z[:]
        data[0, 0] = z.shape[1]
        data[0, 1:] = y.flatten()
        data[1:, 0] = x.flatten()
        print("data", data)
        fp = open("auto_figures/" + filename + ".dat", "w")
        fp.write(np.float32(data).tostring())
        fp.write("\n")
        fp.close()
        return

    # {{{ sort and shape the data for 3d plotting
    def sort_and_xy(self):
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if len(self.dimlabels) > 2:
            raise ValueError(
                "I don't know how to handle something with more than two"
                " dimensions for a surface plot!"
            )
        # {{{ shared to both
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        x_axis = self.retaxis(x_dim).data
        y_axis = self.retaxis(y_dim).data
        # }}}
        return x_axis, y_axis

    def matrices_3d(
        self,
        also1d=False,
        invert=False,
        max_dimsize=1024,
        downsample_self=False,
    ):
        """returns X,Y,Z,x_axis,y_axis
        matrices X,Y,Z, are suitable for a variety of mesh plotting, etc,
        routines x_axis and y_axis are the x and y axes
        """
        this_size = np.array(self.data.shape)
        sortedself = self.copy()
        if np.any(this_size > max_dimsize):
            print(
                lsafen(
                    "Warning! The data is big (%s), so I'm automatically"
                    " downsampling" % (ndshape(self))
                )
            )
            for j in np.where(this_size > max_dimsize)[0]:
                downsampling = int(
                    np.ceil(np.double(this_size[j]) / max_dimsize)
                )
                print("downsampling", self.dimlabels[j], "by", downsampling)
                sortedself = sortedself[self.dimlabels[j], 0::downsampling]
            print(
                lsafen(
                    "I reduced to a max of max_dimsize = %d so the data is"
                    " now %s" % (max_dimsize, ndshape(sortedself))
                )
            )

        x_axis, y_axis = sortedself.sort_and_xy()
        if invert:
            print("trying to invert meshplot-like data")
        X = x_axis * np.ones(np.shape(y_axis))
        Y = np.ones(np.shape(x_axis)) * y_axis
        Z = np.real(sortedself.data)
        if invert:
            X = X[:, ::-1]
            Y = Y[:, ::-1]
            Z = Z[:, ::-1]
        if downsample_self:
            self.data = sortedself.data
            self.setaxis(self.dimlabels[0], x_axis)
            self.setaxis(self.dimlabels[1], y_axis)
        if also1d:
            if invert:
                return X, Y, Z, x_axis[::-1], y_axis[::-1]
            else:
                return X, Y, Z, x_axis, y_axis
        else:
            return X, Y, Z

    # }}}
    def mayavi_surf(self):
        """use the mayavi surf function, assuming that we've already loaded
        mlab during initialization"""
        X, Y, Z = self.matrices_3d()
        s = self.mlab.surf(X, Y, Z)
        return s

    # {{{ 3D mesh plot
    def meshplot(
        self,
        stride=None,
        alpha=1.0,
        onlycolor=False,
        light=None,
        rotation=None,
        cmap=plt.cm.gray,
        ax=None,
        invert=False,
        **kwargs,
    ):
        r"""takes both rotation and light as elevation, azimuth
        only use the light kwarg to generate a black and white shading display
        """
        X, Y, Z = self.matrices_3d()
        if light is True:
            light = [
                0,
                0,
            ]  # I think this is 45 degrees up shining down from the left of
            #    the y axis
        if not onlycolor:
            if ax is None:
                ax = self._init_3d_axis(ax, rotation=rotation)
            else:
                if rotation is not None:
                    raise ValueError(
                        "you can only set the rotation once! (you tried"
                        + repr(rotation)
                        + ")"
                    )
        rstride = 1
        cstride = 1
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        if stride is not None:
            if x_dim in list(stride.keys()):
                rstride = stride[x_dim]
            if y_dim in list(stride.keys()):
                cstride = stride[y_dim]
        if light is not None:
            ls = mpl.colors.LightSource(azdeg=light[1], altdeg=light[0])
            if cmap is not None:
                rgb = ls.shade(Z, cmap)
        else:
            mask = np.isfinite(Z.flatten())
            for_rgb = Z - Z.flatten()[mask].min()
            for_rgb /= for_rgb.flatten()[mask].max()
            if cmap is not None:
                rgb = cmap(for_rgb)
        if onlycolor:
            plt.imshow(rgb)
        else:
            if light is None:
                if cmap is not None:
                    kwargs.update({"cmap": cmap})
                ax.plot_surface(
                    X,
                    Y,
                    Z,
                    rstride=rstride,
                    cstride=cstride,
                    shade=True,
                    **kwargs,
                )
            else:
                newkwargs = {}
                newkwargs["linewidth"] = 0.0
                newkwargs.update(kwargs)
                if cmap is not None:
                    newkwargs["facecolors"] = rgb
                ax.plot_surface(
                    X,
                    Y,
                    Z,
                    rstride=rstride,
                    cstride=cstride,
                    alpha=alpha,
                    shade=False,
                    **newkwargs,
                )
            ax.set_xlabel(x_dim)
            ax.set_ylabel(y_dim)
            ax.set_zlabel(self.name())
        if onlycolor:
            return
        else:
            return ax

    def contour(self, labels=True, **kwargs):
        """Contour plot -- kwargs are passed to the matplotlib
        `contour` function.

        See docstring of `figlist_var.image()` for an example

        Attributes
        ----------
        labels : boolean
            Whether or not the levels should be labeled.
            Defaults to True
        """
        x_axis, y_axis = self.dimlabels
        x = self.getaxis(x_axis)[:, None]
        y = self.getaxis(y_axis)[None, :]
        if "levels" not in list(kwargs.keys()):
            kwargs.update(
                {"levels": r_[self.data.min() : self.data.max() : 30j]}
            )
        cs = plt.contour(
            x * np.ones_like(y), np.ones_like(x) * y, self.data, **kwargs
        )
        if labels:
            plt.clabel(cs, inline=1)
        plt.xlabel(self.unitify_axis(x_axis))
        plt.ylabel(self.unitify_axis(y_axis))
        return cs

    pcolor = this_plotting.pcolormesh.pcolormesh

    def waterfall(
        self, alpha=0.3, ax=None, rotation=None, color="b", edgecolor="k"
    ):
        if ax is None:
            ax = self._init_3d_axis(ax, rotation=rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        if len(self.dimlabels) > 2:
            raise ValueError(
                "I don't know how to handle something with more than two"
                " dimensions for a surface plot!"
            )
        # {{{ shared to both
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        try:
            x_axis = self.retaxis(x_dim).data
        except Exception:
            raise ValueError(
                strm(
                    "trying to get the info on axis",
                    x_dim,
                    "which is",
                    self.getaxis(x_dim),
                )
            )
        y_axis = self.retaxis(y_dim).data
        # }}}
        ax.set_xlabel(self.unitify_axis(x_dim))
        ax.set_ylabel(self.unitify_axis(y_dim))
        ax.set_zlabel(self.unitify_axis(self.name(), is_axis=False))
        verts = []
        xs = x_axis.flatten()
        xs = r_[
            xs[0], xs, xs[-1]
        ]  # add points for the bottoms of the vertices
        ys = y_axis.flatten()
        for j in range(0, len(ys)):
            zs = self[y_dim, j].data.flatten()
            zs = r_[0, zs, 0]
            verts.append(list(zip(xs, zs)))  # one of the faces
        poly = mpl.collections.PolyCollection(
            verts, facecolors=[color] * len(verts), edgecolors=edgecolor
        )  # the individual facecolors would go here
        poly.set_alpha(alpha)
        ax.add_collection3d(poly, zs=ys, zdir="y")
        ax.set_zlim3d(self.data.min(), self.data.max())
        ax.set_xlim3d(xs.min(), xs.max())
        ax.set_ylim3d(ys.min(), ys.max())
        return ax

    def _init_3d_axis(self, ax, rotation=None):
        # other things that should work don't work correctly, so use this to
        # initialize the 3D axis
        # ax.view_init(elev = rotation[0],azim = rotation[1])
        if rotation is None:
            rotation = [0, 0]
        if ax is None:
            fig = plt.gcf()
            ax = axes3d.Axes3D(fig)
            print("I'm trying to rotate to", rotation)
            # ax.view_init(20,-120)
            # ax.view_init(elev = 20 + rotation[1],azim = -120 + rotation[0])
            ax.view_init(azim=rotation[0], elev=rotation[1])
        return ax

    def oldtimey(
        self,
        alpha=0.5,
        ax=None,
        linewidth=None,
        sclinewidth=20.0,
        light=True,
        rotation=None,
        invert=False,
        **kwargs,
    ):
        sortedself = self.copy()
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if invert:
            print("trying to invert oldtimey")
        if linewidth is None:
            linewidth = sclinewidth / sortedself.data.shape[1]
            print("setting linewidth to %0.1f" % linewidth)
        if ax is None:
            ax = sortedself._init_3d_axis(ax, rotation=rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        ax = sortedself.meshplot(
            linewidth=0, light=light, ax=ax, invert=invert
        )
        # return
        if len(sortedself.dimlabels) > 2:
            raise ValueError(
                "I don't know how to handle something with more than two"
                " dimensions for a surface plot!"
            )
        # {{{ shared to both
        x_dim = sortedself.dimlabels[0]
        y_dim = sortedself.dimlabels[1]
        x_axis = sortedself.retaxis(x_dim).data
        y_axis = sortedself.retaxis(y_dim).data
        # }}}
        xs = x_axis.flatten()
        ys = y_axis.flatten()  # this is the depth dimension
        if invert:
            ys = ys[::-1]
        for j in range(0, len(ys)):
            zs = sortedself[
                y_dim, j
            ].data.flatten()  # pulls the data (zs) for a specific y slice
            if invert:
                zs = zs[::-1]
            ax.plot(xs, np.ones(len(xs)) * ys[j], zs, "k", linewidth=linewidth)
        ax.set_zlim3d(sortedself.data.min(), sortedself.data.max())
        ax.set_xlim3d(xs.min(), xs.max())
        # if invert:
        #    ax.set_ylim3d(ys.max(),ys.min())
        # else:
        ax.set_ylim3d(ys.min(), ys.max())
        return ax

    # }}}
    # }}}
    # {{{ error-related functions
    def normalize(
        self, axis, first_figure=None
    ):  # ,whichpoint = slice(0,1,None)):
        x = self.data
        n = len(x)
        S = sparse.lil_matrix((n, n))
        S.setdiag((self.get_error()) ** 2)
        self.set_error(None)
        first_point = self[
            axis, 0:1
        ].copy()  # this makes another instance that contains just the first
        #           point, for error propagation
        B = sparse.lil_matrix((n, n))
        B.setdiag(1.0 / x)
        B[0, :] = -x / (
            x[0] ** 2
        )  # Sparse seems to only support row assignment, so make the transpose
        #    to give it what it wants
        B[0, 0] = 0.0
        B = B.T
        E = B * S * (B.T)  # verified that this is matrix multiplication
        self /= first_point  # this gives the experimentally measured E
        # {{{ now, chop out the first point, which is meaningless
        self = self[axis, 1:]
        E = E[1:, 1:]
        # }}}
        self.set_error(sqrt(E.diagonal()))
        E.setdiag(np.zeros(n - 1))
        self.data_covariance = E
        return self

    def get_covariance(self):
        """this returns the covariance matrix of the data"""
        if hasattr(self, "data_covariance"):
            E = self.data_covariance.copy()
        else:
            n = np.size(self.data)
            E = sparse.lil_matrix((n, n))
        try:
            E.setdiag(self.get_error() ** 2)
        except Exception:
            raise ValueError(
                strm(
                    "Problem getting covariance because error is",
                    self.get_error(),
                )
            )
        return E.toarray()

    # }}}
    # {{{ shortcuts for axes
    def axlen(self, axis):
        r"""return the size (length) of an axis, by name

        Parameters
        ----------

        axis: str
            name of the axis whos length you are interested in
        """
        return np.shape(self.data)[self.axn(axis)]

    def axn(self, axis):
        r"""Return the index number for the axis with the name "axis"

        This is used by many other methods.
        As a simple example,
        self.:func:`axlen`(axis) (the axis length) returns
        ``np.shape(self.data)[self.axn(axis)]``

        Parameters
        ----------

        axis: str
            name of the axis
        """
        try:
            return self.dimlabels.index(axis)
        except Exception:
            raise ValueError(
                " ".join(
                    map(
                        repr,
                        [
                            "there is no axis named",
                            axis,
                            "all axes are named",
                            self.dimlabels,
                        ],
                    )
                )
            )

    def indices(self, axis_name, values):
        """Return a string of indeces that most closely match the axis labels
        corresponding to values. Filter them to make sure they are unique."""
        x = self.getaxis(axis_name)
        retval = []
        for j in values:
            retval.append(np.argmin(abs(x - j)))
        retval = np.array(retval)
        return np.unique(retval)

    # }}}
    # {{{ dictionary functions -- these convert between two formats:
    #     dictionary -- stuff labeled according the dimension label.  list --
    #     same information, but it's assumed they are listed in the order given
    #     by "dimlabels"
    def mkd(self, *arg, **kwargs):
        "make dictionary format"
        give_None = process_kwargs([("give_None", True)], kwargs)
        if len(arg) == 1:
            input_list = arg[0]
            if emptytest(input_list):
                return dict(zip(self.dimlabels, [None] * len(self.dimlabels)))
            if len(input_list) != len(self.dimlabels):
                print(
                    r"{\color{red}WARNING! mkd error (John will fix this"
                    r" later):}"
                )
                print(
                    "When making a dictionary with mkd, you must pass a list"
                    " that has one element for each dimension!  dimlabels is "
                    + repr(self.dimlabels)
                    + " and you passed "
                    + repr(arg)
                    + "\n\n"
                )
                raise ValueError(
                    "When making a dictionary with mkd, you must pass a list"
                    " that has one element for each dimension!  dimlabels is "
                    + repr(self.dimlabels)
                    + " and you passed "
                    + repr(arg)
                )
            for i, v in enumerate(input_list):
                if isinstance(v, np.ndarray):
                    if v.shape == () and v.size == 0:
                        input_list[i] = None
                    if v.dtype.type in [np.str_, np.bytes_]:
                        input_list[i] = str(v)
            if give_None:
                return dict(zip(self.dimlabels, input_list))
            else:
                # {{{ don't return values for the things that are None
                mykeys = [
                    self.dimlabels[j]
                    for j in range(0, len(self.dimlabels))
                    if input_list[j] is not None
                ]
                myvals = [
                    input_list[j]
                    for j in range(0, len(self.dimlabels))
                    if input_list[j] is not None
                ]
                return dict(zip(mykeys, myvals))
                # }}}
        elif len(arg) == 0:
            if not give_None:
                raise ValueError(
                    "You can't tell me not to give none and then not pass me"
                    " anything!!"
                )
            return dict(zip(self.dimlabels, [None] * len(self.dimlabels)))
        else:
            raise ValueError(
                strm(
                    ".mkd() doesn't know what to do with %d arguments",
                    len(arg),
                )
            )

    def fld(self, dict_in, noscalar=False):
        "flatten dictionary -- return list"
        return [dict_in[x] for x in self.dimlabels]

    # }}}
    # {{{ set + get the error + units
    # {{{ set units
    def set_units(self, *args):
        if len(args) == 2:
            unitval = args[1]  # later, have some type of processing bojive
            if (
                self.axis_coords_units is None
                or len(self.axis_coords_units) == 0
            ):
                self.axis_coords_units = [None] * len(self.dimlabels)
            self.axis_coords_units[self.axn(args[0])] = unitval
        elif len(args) == 1:
            unitval = args[0]  # later, have some type of processing bojive
            self.data_units = unitval
        else:
            raise TypeError(
                ".set_units() takes data units or 'axis' and axis units"
            )
        return self

    def human_units(self):
        """This function attempts to choose "human-readable" units for axes or
        *y*-values of the data.
        (Terminology stolen from "human readable" file
        sizes when running shell commands.)
        This means that it looks at the axis or at the
        *y*-values and converts *e.g.* seconds to milliseconds where
        appropriate, also multiplying or dividing the data in an appropriate
        way.
        """
        prev_label = self.get_units()
        for thisaxis in self.dimlabels:
            prev_label = self.get_units(thisaxis)
            if prev_label is not None and len(prev_label) > 0:
                data_to_test = self.getaxis(thisaxis)
                if data_to_test is not None:
                    average_oom = det_oom(data_to_test)
                    x = self.getaxis(thisaxis)
                    result_label = apply_oom(
                        average_oom, x, prev_label=prev_label
                    )
                    self.set_units(thisaxis, result_label)
                else:
                    logger.debug(strm(thisaxis, "does not have an axis label"))
            else:
                logger.debug(strm(thisaxis, "does not have a unit label"))
        return self

    # }}}
    # {{{ get units
    def units_texsafe(self, *args):
        retval = self.get_units(*args)
        if retval is None:
            return None
        if retval.find("\\") > -1:
            retval = "$" + retval + "$"
        return retval

    def replicate_units(self, other):
        for thisaxis in self.dimlabels:
            if other.get_units(thisaxis) is not None:
                self.set_units(thisaxis, other.get_units(thisaxis))
        if other.get_units() is not None:
            self.set_units(other.get_units(thisaxis))
        return self

    def div_units(self, arg1, arg2=None):
        """divide units of the data (or axis)
        by the units that are given, and
        return the multiplier as a number.

        If the result is not dimensionless,
        an error will be generated.

        e.g. call as `d.divide_by("axisname","s")`
        to divide the axis units by seconds

        or `d.divide_by("s")` to divide the
        data units by seconds.
        """
        if arg2 is not None:
            denom_units = Q_(arg2)
            numer_units = Q_(self.get_units(arg1))
        else:
            denom_units = Q_(arg1)
            numer_units = Q_(self.get_units())
        retval = (numer_units / denom_units).to_base_units()
        assert retval.check(""), (
            "the quotient you're asking for is not unitless -- it has units"
            f" of {retval.dimensionality}"
        )
        return retval.magnitude

    def get_units(self, *args):
        if len(args) == 1:
            if self.axis_coords_units is None:
                return None
            if len(self.axis_coords_units) == 0:
                return None
            try:
                return self.axis_coords_units[self.axn(args[0])]
            except Exception:
                raise RuntimeError(
                    strm(
                        "problem getting units for",
                        args[0],
                        "dimension",
                        self.dimlabels,
                        self.axis_coords_units,
                    )
                )
        elif len(args) == 0:
            return self.data_units
        else:
            raise ValueError(".set_units() takes axis or nothing")

    # }}}
    # {{{ set error
    def set_error(self, *args):
        r"""set the errors: either

        `set_error('axisname',error_for_axis)` or `set_error(error_for_data)`

        `error_for_data` can be a scalar, in which case, **all** the data
        errors are set to `error_for_data`

        .. todo::
                several options below -- enumerate them in the documentation
        """
        if (len(args) == 1) and np.isscalar(args[0]):
            if args[0] == 0:
                args = (np.zeros_like(self.data),)
            else:
                args = (np.ones_like(self.data) * args[0],)
        if (len(args) == 1) and (isinstance(args[0], np.ndarray)):
            self.data_error = np.reshape(args[0], np.shape(self.data))
        elif (len(args) == 1) and (isinstance(args[0], list)):
            self.data_error = np.reshape(
                np.array(args[0]), np.shape(self.data)
            )
        elif (
            (len(args) == 2)
            and (isinstance(args[0], str))
            and (isinstance(args[1], np.ndarray))
        ):
            self.axis_coords_error[self.axn(args[0])] = args[1]
        elif (
            (len(args) == 2)
            and (isinstance(args[0], str))
            and (np.isscalar(args[1]))
        ):
            self.axis_coords_error[self.axn(args[0])] = args[1] * np.ones_like(
                self.getaxis(args[0])
            )
        elif (len(args) == 1) and args[0] is None:
            self.data_error = None
        else:
            raise TypeError(
                " ".join(
                    map(
                        repr,
                        [
                            "Not a valid argument to set_error:",
                            list(map(type, args)),
                        ],
                    )
                )
            )
        return self

    # }}}
    # {{{ random mask -- throw out points
    def random_mask(self, axisname, threshold=np.exp(-1.0), inversion=False):
        r"""generate a random mask with about 'threshold' of the points thrown
        out"""
        if inversion:
            threshold = threshold / (1.0 - threshold)
        myr = np.random.rand(
            self.data.shape[self.axn(axisname)]
        )  # random np.array same length as the axis
        return myr > threshold

    # }}}
    # {{{ get error
    def get_error(self, *args):
        """get a copy of the errors\neither
        set_error('axisname',error_for_axis) or set_error(error_for_data)"""
        if len(args) == 0:
            if self.data_error is None:
                return None
            else:
                return np.real(self.data_error)
        elif len(args) == 1:
            thearg = args[0]
            if isinstance(thearg, np.str_):
                thearg = str(
                    thearg
                )  # like in the other spot, this became necessary with some
                #    upgrade, though I'm not sure that I should maybe just
                #    change the error functions to treat the numpy string in
                #    the same way
            if isinstance(thearg, str):
                if len(self.axis_coords_error) == 0:
                    self.axis_coords_error = [None] * len(
                        self.dimlabels
                    )  # is we have an np.empty axis_coords_error, need to fill
                    #    with None's
                try:
                    errorforthisaxis = self.axis_coords_error[self.axn(thearg)]
                except Exception:
                    raise RuntimeError(
                        strm(
                            "Problem trying to load error",
                            self.axn(thearg),
                            "for axis",
                            thearg,
                            "out of",
                            self.axis_coords_error,
                        )
                    )
                if errorforthisaxis is None:
                    return None
                else:
                    x = self.axis_coords_error[self.axn(thearg)]
                    if isinstance(x, np.ndarray):
                        if x.shape == ():
                            return None
                        else:
                            return np.real(
                                self.axis_coords_error[self.axn(thearg)]
                            )
                    else:
                        return np.real(
                            self.axis_coords_error[self.axn(thearg)]
                        )
        else:
            raise ValueError(
                strm(
                    "Not a valid argument to get_error: *args=",
                    args,
                    "map(type,args)=",
                    list(map(type, args)),
                )
            )
        # }}}

    # }}}
    # {{{ match dims --
    def matchdims(self, other):
        r"add any dimensions to self that are not present in other"
        # print 'diagnose: matching',ndshape(self),'to',ndshape(other)
        addeddims = list(set(self.dimlabels) ^ set(other.dimlabels))
        newdims = addeddims + self.dimlabels
        newshape = [1] * len(addeddims) + list(self.data.shape)
        # print 'diagnose: newshape',newshape,'newdims',newdims
        # {{{ reshape to the new dimensions
        new_axis_coords = [r_[1]] * len(addeddims) + self.axis_coords
        self.data = self.data.reshape(newshape)
        self.dimlabels = newdims
        if len(self.axis_coords) > 0:
            self.axis_coords = new_axis_coords
        # }}}
        # {{{ if we are adding dimensions, we will need to reorder to match the
        #     order of the other
        if len(addeddims) > 0:
            self.reorder(other.dimlabels)
        # }}}
        return self

    # }}}
    # {{{ rename
    def rename(self, previous, new):
        self.dimlabels = list(self.dimlabels)  # otherwise, it weirdly
        # changes the names of copies/sources
        self.dimlabels[self.dimlabels.index(previous)] = new
        return self

    # }}}
    @property
    def shape(self):
        return ndshape(self)

    @shape.setter
    def shape(self):
        raise ValueError(
            "You can't set the shape property -- right now, it's just used to"
            " read the shape"
        )

    # {{{ display and other properties
    # {{{ set and get prop
    def unset_prop(self, arg):
        "remove a 'property'"
        self.other_info.pop(arg)
        if len(self.other_info) == 0:
            del self.other_info
        return self

    def set_prop(self, *args):
        r"""set a 'property' of the nddata
        This is where you can put all unstructured information (e.g.
        experimental parameters, etc)
        """
        if len(args) == 2:
            propname, val = args
            self.other_info.update({propname: val})
        elif len(args) == 1 and isinstance(args[0], dict):
            self.other_info.update(args[0])
        else:
            raise ValueError("I don't know what you're passing to set prop!!!")
        return self

    def copy_props(self, other):
        r"""Copy all properties (see :func:`get_prop`) from another nddata
        object -- note that these include properties pertaining the the FT
        status of various dimensions."""
        self.other_info.update(deepcopy(other.other_info))
        return self

    def get_prop(self, propname=None):
        r"""return arbitrary ND-data properties (typically acquisition
        parameters *etc.*) by name (`propname`)

        In order to allow ND-data to store acquisition parameters and other
        info that accompanies the data,
        but might not be structured in a gridded format, nddata instances
        always have a `other_info` dictionary attribute,
        which stores these properties by name.

        If the property doesn't exist, this returns `None`.

        Parameters
        ----------
        propname: str
            Name of the property that you're want returned.
            If this is left out or set to "None" (not given), the names of the
            available properties are returned.
            If no exact match is found, and propname contains a . or * or [,
            it's assumed to be a regular expression.
            If several such matches are found, the error message is
            informative.

            .. todo::
                have it recursively search dictionaries (e.g. bruker acq)

        Returns
        -------
        The value of the property (can by any type) or `None` if the property
        doesn't exist.
        """
        if propname is None:
            return self.other_info.keys()
        if propname not in self.other_info.keys():
            if "." in propname or "*" in propname or "[" in propname:
                propname_re = re.compile(propname)
                matches = [
                    j for j in self.other_info.keys() if propname_re.match(j)
                ]
                if len(matches) == 0:
                    return None
                assert (
                    len(matches) == 1
                ), "I found %d matches for regexp %s in properties: %s" % (
                    len(matches),
                    propname,
                    " ".join(matches),
                )
                return self.other_info[matches[0]]
            else:
                return None
        return self.other_info[propname]

    def name(self, *arg):
        r"""args:
        .name(newname) --> Name the object (for storage, etc)
        .name() --> Return the name"""
        if len(arg) == 1:
            self.set_prop("name", arg[0])
            return self
        elif len(arg) == 0:
            return self.get_prop("name")
        else:
            raise ValueError("invalid number of arguments")

    # }}}
    # {{{ set and get plot color
    def set_plot_color(self, thiscolor):
        if thiscolor is None:
            return
        if thiscolor is str:
            colordict = {
                "r": [1, 0, 0],
                "g": [0, 1, 0],
                "b": [0, 0, 1],
                "k": [0, 0, 0],
                "y": [0.5, 0.5, 0],
                "o": [0.75, 0.25, 0],
                "c": [0, 0.5, 0.5],
            }
            try:
                thiscolor = colordict[thiscolor]
            except Exception:
                raise ValueError(strm("Color", thiscolor, "not in dictionary"))
        self.other_info.update({"plot_color": thiscolor})
        return

    def set_plot_color_next(self):
        self.set_plot_color(next(default_cycler))

    def get_plot_color(self):
        if "plot_color" in self.get_prop():
            return self.other_info["plot_color"]
        else:
            return None

    # }}}
    # }}}
    svd = MM_svd
    # {{{ arithmetic
    along = MM_along
    dot = MM_dot

    def __add__(self, arg):
        if np.isscalar(arg):
            A = self.copy()
            if isinstance(arg, complex) and self.data.dtype not in [
                np.complex128,
                np.complex64,
            ]:
                A.data = np.complex128(A.data)
            A.data += arg
            # error does not change
            return A
        # {{{ shape and add
        A, B = self.aligndata(arg)
        logger.debug(
            strm("after alignment, right data looks like:", ndshape(B))
        )
        retval = A.copy()
        retval.data = A.data + B.data
        # }}}
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0
        if Aerr is not None:
            Rerr += (Aerr) ** 2
        if Berr is not None:
            Rerr += (Berr) ** 2
        Rerr = sqrt(np.real(Rerr))  # convert back to stdev
        if Aerr is None and Berr is None:
            Rerr = None
        retval.set_error(Rerr)
        return retval

    def __sub__(self, arg):
        return self.__add__(-1 * arg)

    def __lt__(self, arg):
        if isinstance(arg, np.ndarray):
            retval = self.copy()
            retval.data = retval.data < arg
            return retval
        elif isinstance(arg, nddata):
            retval, B = self.aligndata(arg)
            retval.data = retval.data < B.data
            return retval
        elif np.isscalar(arg):
            retval = self.copy()
            retval.data = retval.data < arg
            return retval
        else:
            raise ValueError(
                "I don't know what to do with an argument of type"
                + repr(type(arg))
            )

    def __gt__(self, arg):
        if isinstance(arg, np.ndarray):
            retval = self.copy()
            retval.data = retval.data > arg
            return retval
        elif isinstance(arg, nddata):
            retval, B = self.aligndata(arg)
            retval.data = retval.data > B.data
            return retval
        elif np.isscalar(arg):
            retval = self.copy()
            retval.data = retval.data > arg
            return retval
        else:
            raise ValueError(
                "I don't know what to do with an argument of type"
                + repr(type(arg))
            )

    def __le__(self, arg):
        if isinstance(arg, np.ndarray):
            retval = self.copy()
            retval.data = retval.data <= arg
            return retval
        elif isinstance(arg, nddata):
            retval, B = self.aligndata(arg)
            retval.data = retval.data <= B.data
            return retval
        elif np.isscalar(arg):
            retval = self.copy()
            retval.data = retval.data <= arg
            return retval
        else:
            raise ValueError(
                "I don't know what to do with an argument of type"
                + repr(type(arg))
            )

    def __ge__(self, arg):
        if isinstance(arg, np.ndarray):
            retval = self.copy()
            retval.data = retval.data >= arg
            return retval
        elif isinstance(arg, nddata):
            retval, B = self.aligndata(arg)
            retval.data = retval.data >= B.data
            return retval
        elif np.isscalar(arg):
            retval = self.copy()
            retval.data = retval.data >= arg
            return retval
        else:
            raise ValueError(
                "I don't know what to do with an argument of type"
                + repr(type(arg))
            )

    __matmul__ = MM_matmul

    def __mul__(self, arg):
        # {{{ do scalar multiplication
        if np.isscalar(arg):
            A = self.copy()
            if isinstance(arg, complex) and self.data.dtype not in [
                np.complex128,
                np.complex64,
            ]:
                A.data = np.complex128(A.data)
            A.data *= arg
            if A.get_error() is not None:
                error = A.get_error()
                error *= abs(arg)
            return A
        # }}}
        # {{{ shape and multiply
        A, B = self.aligndata(arg)
        retval = A.copy()
        retval.data = A.data * B.data
        # }}}
        # {{{ if we have error for both the sets of data, I should propagate
        #     that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0  # we can have error on one or both, so we're going to need
        #             to add up the variances
        if Aerr is not None:
            Rerr += (Aerr * B.data) ** 2
        if Berr is not None:
            Rerr += (Berr * A.data) ** 2
        Rerr = sqrt(np.real(Rerr))  # convert back to stdev
        if Aerr is None and Berr is None:
            Rerr = None
        # }}}
        retval.set_error(Rerr)
        return retval

    def __rpow__(self, arg):
        result = self.copy()
        result.set_error(None)
        logger.debug(
            "error propagation for right power not currently supported (do you"
            " need this, really?)"
        )
        assert np.isscalar(arg) or isinstance(arg, np.ndarray), (
            "currently right power only supported for np.ndarray and scalars"
            " -- do you really need something else??"
        )
        result.data = arg**self.data
        return result

    def __pow__(self, arg):
        if arg == -1:
            x = self.get_error()
            result = self.copy()
            result.data = 1.0 / result.data
            if x is not None:
                result.set_error(abs(x.copy() / (self.data**2)))
            return result
        elif arg == 2:
            return self * self
        else:
            if self.get_error() is not None:
                raise ValueError(
                    strm(
                        "nothing but -1 and 2 supported yet! (you tried to"
                        " raise to a power of "
                        + repr(arg)
                        + ")"
                    )
                )
            else:
                result = self.copy()
                result.data = result.data**arg
                return result

    def __truediv__(self, arg):
        return self.__div__(arg)

    def __rtruediv__(self, arg):
        return self.__rdiv__(arg)

    def __div__(self, arg):
        if np.isscalar(arg):
            A = self.copy()
            A.data /= arg
            if A.get_error() is not None:
                error = A.get_error()
                error /= abs(arg)
            return A
        A, B = self.aligndata(arg)
        retval = A.copy()
        retval.data = A.data / B.data
        # {{{ if we have error for both the sets of data, I should propagate
        #     that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0  # we can have error on one or both, so we're going to need
        #             to add up the variances
        dt128 = np.dtype("complex128")
        if Aerr is not None:
            if (A.data.dtype is dt128) or (
                B.data.dtype is dt128
            ):  # this should avoid the error that Ryan gets
                Rerr += (np.complex128(Aerr) / np.complex128(B.data)) ** 2
            else:
                Rerr += (Aerr / B.data) ** 2
        if Berr is not None:
            if (
                (A.data.dtype is dt128)
                or (Berr.dtype is dt128)
                or (B.data.dtype is dt128)
            ):  # this should avoid the error that Ryan gets
                Rerr += (
                    np.complex128(A.data)
                    * np.complex128(Berr)
                    / (np.complex128(B.data) ** 2)
                ) ** 2
            else:
                try:
                    Rerr += (A.data * Berr / (B.data**2)) ** 2
                except Exception:
                    raise ValueError(
                        strm(
                            "self was",
                            self,
                            "arg was",
                            arg,
                            "dtype of A.data",
                            A.data.dtype,
                            "dtype of Berr",
                            Berr.dtype,
                            "dtype of B.data",
                            Berr,
                        )
                    )
        try:
            Rerr = sqrt(
                np.real(Rerr)
            )  # convert back to stdev --> note that this has problems with
            #    complex numbers, hence the "abs" above
        except AttributeError:
            raise AttributeError(
                strm("Rerr gave an attribute error when you passed", Rerr)
            )
        if Aerr is None and Berr is None:
            Rerr = None
        # }}}
        retval.set_error(Rerr)
        return retval

    def __invert__(self):
        if self.data.dtype is np.dtype("bool"):
            self.data = ~self.data
            return self
        else:
            raise ValueError("invert only implemented for boolean now")

    def __abs__(self):
        return self.runcopy(abs)

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, arg):
        return -1 * (self - arg)

    def __neg__(self):
        return -1 * self

    def __rdiv__(self, arg):
        return arg * (self ** (-1))

    # def real(self):
    #    self.data = np.real(self.data)
    #    return self
    # }}}
    # {{{ align data
    def aligndata(self, arg):
        r"""This is a fundamental method used by all of the arithmetic
        operations.
        It uses the dimension labels of `self` (the current instance) and `arg`
        (an nddata passed to this method) to generate two corresponding output
        nddatas that I refer to here, respectively, as `A` and `B`.  `A` and
        `B` have dimensions that are "aligned" -- that is, they are identical
        except for singleton dimensions (note that numpy automatically tiles
        singleton dimensions).  Regardless of how the dimensions of `self.data`
        and `arg.data` (the underlying numpy data) were ordered, `A.data` and
        `B.data` are now ordered identically, where dimensions with the same
        label (`.dimlabel`) correspond to the same numpy index.  This allows
        you do do math.

        Note that, currently, both `A` and `B` are given a full set of axis
        labels, even for singleton dimensions.  This is because we're assuming
        you're going to do math with them, and that the singleton dimensions
        will be expanded.

        Parameters
        ==========
        arg : nddata or np.ndarray
            The nddata that you want to align to `self`.
            If arg is an np.ndarray, it will try to match dimensions to self
            based
            on the length of the dimension.
            **Note:** currently there is an issue where this will only really
            work for 1D data, since it first makes an nddata instance based on
            arg, which apparently collapses multi-D data to 1D data.

        Returns
        =======
        A : nddata
            realigned version of `self`
        B : nddata
            realigned version of `arg` (the argument)
        """
        # {{{ if zero dimensional, fake a singleton dimension and recurse
        # {{{ unless both are zero dimensional, in which case, just leave alone
        logger.debug(strm("starting aligndata"))
        if np.isscalar(arg) or isinstance(arg, np.ndarray):
            arg = nddata(arg)
            index_dims = [
                j
                for j in r_[0 : len(arg.dimlabels)]
                if arg.dimlabels[j] == "INDEX"
            ]
            for j in index_dims:  # find dimension of matching length
                match_dims = np.nonzero(
                    arg.data.shape[j] == np.array(self.data.shape)
                )[0]
                if len(match_dims) > 0:
                    arg.dimlabels[j] = self.dimlabels[match_dims[0]]
                if len(match_dims) != len(index_dims):
                    raise ValueError(
                        "you seem to by multiplying by something with an"
                        " 'INDEX' data and something that doesn't have that --"
                        " is this really what you want?  (this is commonly"
                        " produced by multiplying a mismatched np.ndarray by"
                        " an nddata)"
                    )
        if ndshape(self).zero_dimensional and ndshape(arg).zero_dimensional:
            logger.debug(strm("(1) yes, I found something zero dimensional"))
            return self.copy(), arg.copy()
        # }}}
        elif ndshape(self).zero_dimensional:
            logger.debug(strm("(2) yes, I found something zero dimensional"))
            logger.debug(strm("yes, I found something zero dimensional"))
            A = self.copy()
            A.dimlabels = [arg.dimlabels[0]]
            A.data = A.data.reshape(1)
            return A.aligndata(arg)
        elif ndshape(arg).zero_dimensional:
            logger.debug(strm("(3) yes, I found something zero dimensional"))
            logger.debug(strm("yes, I found something zero dimensional"))
            arg = arg.copy()
            arg.dimlabels = [self.dimlabels[0]]
            arg.data = arg.data.reshape(1)
            return self.aligndata(arg)
        # }}}
        selfout = self.copy()  # copy self
        assert len(selfout.data.shape) != 0 and len(arg.data.shape) != 0, (
            "neither self nor arg should be zero dimensional at this point"
            " (previous code should have taken care of that"
        )
        # {{{create newdims, consisting of dimlabels for self, followed by the
        # names of the dimensions in arg that are not also in self -- order for
        # both is important; then create a matching selfshape
        augmentdims = [
            x
            for x in arg.dimlabels
            if x in set(self.dimlabels) ^ set(arg.dimlabels)
        ]  # dims in arg
        #                   but not self, ordered as they were in arg
        newdims = self.dimlabels + augmentdims
        selfshape = list(selfout.data.shape) + list(
            np.ones(len(augmentdims), dtype=np.uint64)
        )  # there is no need to
        #       transpose self, since its order is preserved
        # }}}
        argout = arg.copy()
        # {{{ now create argshape for the reshaped argument
        #
        #  only the labels valid for arg, ordered as they are in newdims
        new_arg_labels = [x for x in newdims if x in arg.dimlabels]
        # following should be a better solution
        argshape = list(np.ones(len(newdims), dtype=np.int64))
        logger.debug(
            strm(
                "DEBUG 2: shape of self",
                ndshape(self),
                "self data shape",
                self.data.shape,
                "shape of arg",
                ndshape(arg),
                "arg data shape",
                arg.data.shape,
            )
        )
        logger.debug(
            strm(
                "DEBUG 3: shape of selfout",
                ndshape(selfout),
                "selfout data shape",
                selfout.data.shape,
                "shape of argout",
                ndshape(argout),
                "argout data shape",
                argout.data.shape,
            )
        )
        # {{{ wherever the dimension already exists in arg, pull the shape from
        #     arg
        for j, k in enumerate(newdims):
            if k in argout.dimlabels:
                try:
                    argshape[j] = argout.data.shape[argout.axn(k)]
                except Exception:
                    raise ValueError(
                        "There seems to be a problem because the"
                        + "shape of argout is now len:%d"
                        % len(argout.data.shape),
                        argout.data.shape,
                        "while the dimlabels is len:%d"
                        % len(argout.dimlabels),
                        argout.dimlabels,
                    )
        # }}}
        # }}}
        # {{{ transpose arg to match newshape
        argorder = list(map(argout.dimlabels.index, new_arg_labels))  # for
        #          each new dimension, determine the position of the
        #          original dimension
        selfout.data = selfout.data.reshape(np.int64(selfshape))  # and reshape
        #          to its new shape
        selfout.dimlabels = newdims
        try:
            argshape = np.int64(argshape)
            argout.data = argout.data.transpose(argorder).reshape(
                argshape
            )  # and reshape the data
        except ValueError as Argument:
            raise ValueError(
                "the shape of the data is "
                + repr(argout.data.shape)
                + " the transpose "
                + repr(argorder)
                + " and the new shape "
                + repr(argshape)
                + " original arg: "
                + repr(Argument)
            )
        argout.dimlabels = newdims
        # }}}
        # {{{ transpose the data errors appropriately
        if selfout.get_error() is not None:
            try:
                temp = selfout.get_error().copy().reshape(selfshape)
            except ValueError as Argument:
                raise ValueError(
                    "The instance (selfout) has a shape of "
                    + repr(selfout.data.shape)
                    + " but its error has a shape of"
                    + repr(selfout.get_error().shape)
                    + "!!!\n\n(original argument:\n"
                    + repr(Argument)
                    + "\n)"
                )
            selfout.set_error(temp)
        if argout.get_error() is not None:
            try:
                temp = (
                    argout.get_error()
                    .copy()
                    .transpose(argorder)
                    .reshape(argshape)
                )
            except ValueError as Argument:
                if (
                    argout.data.shape == (1,)
                    and argout.get_error().shape == ()
                ):
                    temp = (
                        np.array(argout.get_error())
                        .reshape((1,))
                        .copy()
                        .transpose(argorder)
                        .reshape(argshape)
                    )
                else:
                    raise ValueError(
                        "The argument (argout) has a shape of "
                        + repr(argout.data.shape)
                        + " but its error has a shape of"
                        + repr(argout.get_error().shape)
                        + "(it's "
                        + repr(argout.get_error())
                        + ")!!!\n\n(original argument:\n"
                        + repr(Argument)
                        + "\n)"
                    )
            argout.set_error(temp)
        # }}}
        if (len(selfout.axis_coords) > 0) or (len(argout.axis_coords) > 0):
            # {{{ transfer the errors and the axis labels
            # {{{ make dictionaries for both, and update with info from both,
            #     giving preference to self
            axesdict = selfout.mkd()
            # print "DEBUG 4: original mkd",axesdict
            errordict = selfout.mkd()

            # {{{ define a function that allows me to only update non-zero axes
            def non_empty_axes(input_data, ret_err=False):
                if ret_err:
                    temp_dict = input_data.mkd(input_data.axis_coords_error)
                else:
                    temp_dict = input_data.mkd(input_data.axis_coords)
                temp_dict = {
                    k: v
                    for k, v in temp_dict.items()
                    if v is not None and len(v) > 0
                }
                return temp_dict

            # }}}
            # {{{ add the axes and errors for B
            if isinstance(arg.axis_coords, list):
                if len(arg.axis_coords) > 0:
                    axesdict.update(non_empty_axes(arg))
            if isinstance(arg.axis_coords_error, list):
                if len(arg.axis_coords_error) > 0 and not all(
                    [x is None for x in arg.axis_coords_error]
                ):
                    errordict.update(arg.mkd(arg.axis_coords_error))
            # }}}
            # {{{ add the axes and errors for A
            if isinstance(self.axis_coords, list):
                if len(self.axis_coords) > 0:
                    axesdict.update(non_empty_axes(self))
            if isinstance(self.axis_coords_error, list):
                if len(self.axis_coords_error) > 0 and not all(
                    [x is None for x in self.axis_coords_error]
                ):
                    errordict.update(self.mkd(self.axis_coords_error))
            # }}}
            # }}}
            selfout.axis_coords_error = selfout.fld(errordict)
            argout.axis_coords_error = selfout.fld(errordict)
            selfout.axis_coords = selfout.fld(axesdict)
            argout.axis_coords = selfout.fld(axesdict)
            # }}}
            selfout.axis_coords_units = [None] * len(newdims)
            argout.axis_coords_units = [None] * len(newdims)
            for thisdim in newdims:
                if thisdim in self.dimlabels:
                    selfout.set_units(thisdim, self.get_units(thisdim))
                    argout.set_units(thisdim, self.get_units(thisdim))
                elif thisdim in arg.dimlabels:
                    selfout.set_units(thisdim, arg.get_units(thisdim))
                    argout.set_units(thisdim, arg.get_units(thisdim))
        return selfout, argout

    # }}}
    # {{{ integrate, differentiate, and sum
    def integrate(self, thisaxis, backwards=False, cumulative=False):
        r"""Performs an integration -- which is similar to a sum, except that
        it takes the axis into account, *i.e.*, it performs:
        :math:`\int f(x) dx`
        rather than
        :math:`\sum_i f(x_i)`

        Gaussian quadrature, etc, is planned for a future version.

        Parameters
        ==========
        thisaxis:
            The dimension that you want to integrate along
        cumulative: boolean (default False)
            Perform a cumulative integral (analogous to a cumulative sum)
            -- *e.g.* for ESR.
        backwards: boolean (default False)
            for cumulative integration -- perform the integration backwards
        """
        if backwards is True:
            self.data = self[thisaxis, ::-1].data
        t = None
        if len(self.axis_coords) > 0:
            t = self.getaxis(thisaxis)
            dt_array = np.diff(t)
            dt = dt_array[0]
            if np.allclose(dt, dt_array):
                simple_integral = True
            else:
                print("not a simple integral")
                simple_integral = False
                dt_array = (
                    0.5
                    * r_[
                        dt_array, dt_array[-1]
                    ]  # diff interval after current point
                    + 0.5 * r_[dt_array[0], dt_array]
                )  # diff interval before current point
                time_slices = self.fromaxis(thisaxis)
                time_slices.data[:] = dt_array
        if t is None:
            raise ValueError("You can't call integrate on an unlabeled axis")
        if not simple_integral:
            result = self * time_slices
            self.data = result.data
        if cumulative:
            self.run_nopop(np.cumsum, thisaxis)
            if backwards is True:
                self.data = self[thisaxis, ::-1].data
        else:
            self.run(np.sum, thisaxis)
        if simple_integral:
            self.data *= dt
        return self

    def phdiff(self, axis, return_error=True):
        """calculate the phase gradient (units: cyc/Δx) along axis,
        setting the error appropriately

        For example, if `axis` corresponds to a time
        axis, the result will have units of frequency
        (cyc/s=Hz).
        """
        if self.get_ft_prop(axis):
            dt = self.get_ft_prop(axis, "df")
        else:
            dt = self.get_ft_prop(axis, "dt")
        if dt is None:
            dt = np.diff(self[axis][r_[0, -1]]).item() / len(self[axis])
        A = self[axis, 1:]
        B = self[axis, :-1]
        if return_error:
            A_sigma = A.get_error()
            A_sigma = 1 if A_sigma is None else A_sigma
            B_sigma = B.get_error()
            B_sigma = 1 if B_sigma is None else B_sigma
        self.data = np.angle(A.data / B.data) / 2 / pi / dt
        self.setaxis(axis, A.getaxis(axis))
        if return_error:
            self.set_error(
                sqrt(
                    A_sigma**2 * abs(0.5 / A.data) ** 2
                    + B_sigma**2 * abs(0.5 / B.data) ** 2
                )
                / 2
                / pi
                / dt
            )
        else:
            self.set_error(None)
        return self

    def diff(self, thisaxis, backwards=False):
        if backwards is True:
            self.data = self[thisaxis, ::-1].data
        self.run_nopop(mydiff, thisaxis)
        if backwards is True:
            self.data = self[thisaxis, ::-1].data
        if len(self.axis_coords) > 0:
            t = self.getaxis(thisaxis)
            dt = t[1] - t[0]
            self.data /= dt
        return self

    def sum(self, axes):
        "calculate the sum along axes, also transforming error as needed"
        if isinstance(axes, str):
            axes = [axes]
        for j in range(0, len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except Exception:
                print("|-ERROR FINDING DIMENSION-----")
                print("| dimlabels is: ", self.dimlabels)
                print("| doesn't contain: ", axes[j])
                print("|-----------------------------")
                raise
            self.data = np.sum(self.data, axis=thisindex)
            if self.get_error() is not None:
                self.set_error(
                    np.sqrt(np.mean(self.get_error() ** 2, axis=thisindex))
                )
            self._pop_axis_info(thisindex)
        return self

    def sum_nopop(self, axes):
        if isinstance(axes, str):
            axes = [axes]
        for j in range(0, len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except Exception:
                print("error, dimlabels is: ", self.dimlabels)
                print("doesn't contain: ", axes[j])
                raise
            temp = list(self.data.shape)
            temp[thisindex] = 1
            self.data = np.sum(self.data, axis=thisindex)
            self.data = self.data.reshape(temp)
        return self

    # }}}
    # {{{ poly. fit
    def eval_poly(self, c, axis, inplace=False, npts=None):
        """Take `c` output (array of polynomial coefficents in ascending order)
        from :func:`~pyspecdata.nddata.polyfit`, and apply it along axis `axis`

        Parameters
        ----------
        c: ndarray
            polynomial coefficients in ascending polynomial order

        """
        if npts:
            temp = np.linspace(self[axis][0], self[axis][-1], npts)
            thisaxis = nddata(temp.copy(), [-1], [axis]).setaxis(axis, temp)
        else:
            thisaxis = self.fromaxis(axis)
        result = 0
        for j in range(len(c)):
            result += c[j] * thisaxis**j
        if inplace:
            self.data = result.data
            return self
        else:
            result.copy_props(self)
            result.set_units(axis, self.get_units(axis))
            return result

    def polyfit(self, axis, order=1, force_y_intercept=None):
        """polynomial fitting routine -- return the coefficients and the fit
        .. note:
            previously, this returned the fit data as a second argument called
            `formult`-- you very infrequently want it to be in the same size as
            the data, though;
            to duplicate the old behavior, just add the line
            ``formult = mydata.eval_poly(c,'axisname')``.

        .. seealso::
            :func:`~pyspecdata.nddata.eval_poly`

        Parameters
        ----------
        axis: str
            name of the axis that you want to fit along
            (not sure if this is currently tested for multi-dimensional data,
            but the idea should be that multiple fits would be returned.)
        order: int
            the order of the polynomial to be fit
        force_y_intercept: double or None
            force the y intercept to a particular value (e.g. 0)

        Returns
        -------
        c: np.ndarray
            a standard numpy np.array containing the coefficients (in ascending
            polynomial order)
        """
        x = self.getaxis(axis).copy().reshape(-1, 1)
        # {{{ make a copy of self with the relevant dimension second to last
        #     (i.e. rows)
        formult = self.copy()
        neworder = list(formult.dimlabels)
        neworder.pop(neworder.index(axis))
        if len(neworder) > 1:
            neworder = neworder[:-1] + [axis] + neworder[-1]
        else:
            neworder = [axis] + neworder
        formult.reorder(neworder)
        # }}}
        y = formult.data
        # {{{ now solve Lx = y, where x is appropriate for our polynomial
        startingpower = 0
        if force_y_intercept is not None:
            startingpower = 1
            L = [x**j for j in range(startingpower, order + 1)]
            L = np.concatenate(
                L, axis=1
            )  # note the totally AWESOME way in which this is done!
            y -= force_y_intercept
            c = np.dot(np.linalg.pinv(L), y)
            c = r_[force_y_intercept, c.ravel()]
        else:
            c = np.polyfit(
                x.ravel(), y, deg=order
            )  # better -- uses Hermite polys
            c = c[::-1]  # give in ascending order, as is sensible
        # }}}
        return c

    # }}}
    # {{{ max and mean
    def _wrapaxisfuncs(self, func):
        # {{{ for convenience, wrap the max and min functions
        if func == np.max:
            func = np.amax
        if func == np.min:
            func = np.amin
        if func == np.diff:
            func = mydiff
        return func
        # }}}

    def argmax(self, *args, **kwargs):
        r"""find the max along a particular axis, and get rid of that axis,
        replacing it with the index number of the max value

        Parameters
        ==========
        raw_index: bool
            return the raw (np.ndarray) numerical index, rather than the
            corresponding axis value Note that the result returned is still,
            however, an nddata (rather than numpy np.ndarray) object.
        """
        # {{{ process arguments
        axes = self._possibly_one_axis(*args)
        raw_index = False
        if "raw_index" in list(kwargs.keys()):
            raw_index = kwargs.pop("raw_index")
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:", repr(kwargs))
        if isinstance(axes, str):
            axes = [axes]
        # }}}
        for j in range(0, len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except Exception:
                print("error, dimlabels is: ", self.dimlabels)
                print("doesn't contain: ", axes[j])
                raise
            temp = self.data.copy()
            temp[~np.isfinite(temp)] = temp[np.isfinite(temp)].min()
            argmax_result = np.argmax(temp, axis=thisindex)
            if raw_index:
                self.data = argmax_result
            else:
                if self.axis_coords[thisindex] is None:
                    raise ValueError(
                        "It doesn't make sense to call argmax if you have"
                        " removed the axis coordinates! (getaxis yields None"
                        " for %s" % thisindex
                    )
                self.data = self.axis_coords[thisindex][argmax_result]
            self._pop_axis_info(thisindex)
        return self

    def argmin(self, *axes, **kwargs):
        r"""If `.argmin('axisname')` find the min along a particular axis, and
        get rid of that axis, replacing it with the index number of the max
        value.
        If `.argmin()`: return a dictionary giving the coordinates of the
        overall minimum point.

        Parameters
        ==========
        raw_index: bool
            Return the raw (np.ndarray) numerical index, rather than the
            corresponding axis value.
            Note that the result returned is still, however, an nddata (rather
            than numpy np.ndarray) object.
        """
        raw_index = process_kwargs([("raw_index", False)], kwargs)
        if len(axes) == 0:
            raw_indices = dict(
                zip(
                    self.dimlabels,
                    np.unravel_index(
                        self.data.ravel().argmin(), self.data.shape
                    ),
                )
            )
            if raw_index:
                return raw_indices
            else:
                return dict(
                    [(k, self.getaxis(k)[v]) for k, v in raw_indices.items()]
                )
        if type(axes) is str:
            axes = [axes]
        for j in range(0, len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except Exception:
                print("error, dimlabels is: ", self.dimlabels)
                print("doesn't contain: ", axes[j])
                raise
            temp = self.data.copy()
            temp[~np.isfinite(temp)] = temp[np.isfinite(temp)].max()
            argmin_result = np.argmin(temp, axis=thisindex)
            if raw_index:
                self.data = argmin_result
            else:
                self.data = self.axis_coords[thisindex][argmin_result]
            self._pop_axis_info(thisindex)
        return self

    def max(self):
        return self.data[np.isfinite(self.data)].max()

    def min(self):
        return self.data[np.isfinite(self.data)].min()

    def cdf(self, normalized=True, max_bins=500):
        """calculate the Cumulative Distribution Function for the data along
        `axis_name`

        only for 1D data right now

        Returns
        =======
        A new nddata object with an axis labeled `values`, and data
        corresponding to the CDF.
        """
        thisaxis = 0
        n_bins = self.data.shape[thisaxis]
        if n_bins > max_bins:
            n_bins = max_bins  # otherwise this takes a while
        bins, vals = np.histogram(self.data, bins=n_bins)
        retval = nddata(np.double(bins), [-1], ["values"]).labels(
            "values", vals[:-1] + (vals[1] - vals[0]) * 0.5
        )
        retval.run_nopop(np.cumsum, "values")
        if normalized:
            print("final value", retval["values", -1])
            retval /= retval["values", -1]
        return retval

    def mean_all_but(self, listofdims):
        "take the mean over all dimensions not in the list"
        for dimname in list(
            self.dimlabels
        ):  # I can't be popping from the list as I iterate over it
            if dimname not in listofdims:
                self.mean(dimname)
        return self

    def mean_weighted(self, axisname):
        r"""perform  the weighted mean along `axisname` (use :math:`\sigma`
        from :math:`\sigma = `self.get_error() do generate :math:`1/\sigma`
        weights) for now, it clears the error of `self`, though it would be
        easy to calculate the new error, since everything is linear

        unlike other functions, this creates working objects that are
        themselves nddata objects this strategy is easier than coding out the
        raw numpy math, but probably less efficient
        """
        # {{{ the weighted mean, pyspecdata style
        weight_matrix = self.copy().set_error(None)
        weight_matrix.data = 1.0 / self.get_error().copy()
        # {{{ find out where anything is nan, and set both error and weight to
        #     0
        nan_mask = np.isnan(self.data)
        nan_mask |= np.isnan(weight_matrix.data)
        weight_matrix.data[nan_mask] = 0
        self.data[nan_mask] = 0
        # }}}
        # {{{ make sure there are no infinite values, because I wouldn't be
        #     sure how to deal with this
        inf_mask = np.isinf(self.data)
        inf_mask |= np.isinf(weight_matrix.data)
        assert not np.any(inf_mask)
        # }}}
        normalization = weight_matrix.copy().run(np.sum, axisname)
        weight_matrix /= normalization
        self.data *= weight_matrix.data
        self.set_error(None)
        self.run(np.sum, axisname)
        # }}}
        return self

    def mean(self, *args, **kwargs):
        r"""Take the mean and (optionally) set the error to the standard
        deviation

        Parameters
        ----------
        std: bool
            whether or not to return the standard deviation as an error
        """
        logger.debug("entered the mean function")
        # {{{ process arguments
        if len(args) > 1:
            raise ValueError("you can't pass more than one argument!!")
        axes = self._possibly_one_axis(*args)
        if "return_error" in kwargs:
            raise ValueError(
                "return_error kwarg no longer used -- use std kwarg if you"
                " want to set the error to the std"
            )
        return_error = process_kwargs([("std", False)], kwargs)
        logger.debug(strm("return error is", return_error))
        if isinstance(axes, str):
            axes = [axes]
        # }}}
        for j in range(0, len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except Exception:
                logger.debug(strm("error, dimlabels is: ", self.dimlabels))
                logger.debug(strm("doesn't contain: ", axes[j]))
                raise
            if self.data_error is not None:
                this_axis_length = self.data.shape[thisindex]
                try:
                    self.data_error = sqrt(
                        np.sum(
                            (self.data * self.data_error) ** 2, axis=thisindex
                        )
                        / (this_axis_length**2)
                    )
                except Exception:
                    raise ValueError(
                        strm(
                            "shape of data",
                            np.shape(self.data),
                            "shape of data error",
                            np.shape(self.data_error),
                        )
                    )
            if return_error:  # since I think this is causing an error
                thiserror = np.std(self.data, axis=thisindex)
                if np.isscalar(thiserror):
                    thiserror = r_[thiserror]
            self.data = np.mean(self.data, axis=thisindex)
            if return_error:  # this needs to go after the data setting
                self.set_error(
                    thiserror
                )  # set the error to the standard deviation
            self._pop_axis_info(thisindex)
            logger.debug(strm("return error is", return_error))
        return self

    def mean_nopop(self, axis):
        self = self.run_nopop(np.mean, axis=axis)
        return self

    # }}}
    # {{{ running functions and popping dimensions
    def _pop_axis_info(self, thisindex):
        r"pop axis by index"
        self.dimlabels.pop(thisindex)
        if self.axis_coords != []:
            self.axis_coords.pop(thisindex)
            if (
                self.axis_coords_error is not None
                and len(self.axis_coords_error) > 0
            ):
                try:
                    self.axis_coords_error.pop(thisindex)
                except Exception:
                    raise RuntimeError(
                        strm(
                            "trying to pop",
                            thisindex,
                            "from",
                            self.axis_coords_error,
                        )
                    )
            if len(self.axis_coords_units) > 0:
                try:
                    self.axis_coords_units.pop(thisindex)
                except Exception:
                    raise IndexError(
                        strm(
                            "trying to pop",
                            thisindex,
                            "from",
                            self.axis_coords_units,
                        )
                    )
        return self

    def cov_mat(self, along_dim):
        """
        calculate covariance matrix for a 2D experiment

        Parameters
        ==========
        along_dim:  str
            the "observations" dimension of the data set (as opposed to
            the variable)
        """
        assert len(self.dimlabels) == 2, (
            "we are only calculating covariance matrices for datasets with one"
            " variable and on observation axis"
        )
        assert along_dim in self.dimlabels
        var_dim = list(set(self.dimlabels) - set([along_dim]))[0]
        var_dim_units = self.get_units(var_dim)
        if self.axn(along_dim) == 0:
            trans = False
        else:
            trans = True
        self.data = np.cov(self.data, rowvar=trans)
        self.setaxis(along_dim, self.getaxis(var_dim).copy())

        def add_subscript(start, sub):
            ismath = re.compile(r"\$(.*)\$")
            m = ismath.match(start)
            if m:
                (start,) = m.groups()
            # look for existing subscripts
            firsttry = re.compile("(.*)_{(.*)}")
            m = firsttry.match(start)
            if m:
                a, b = m.groups()
                return f"${a}_{{{b},{sub}}}$"
            secondtry = re.compile("(.*)_(.*)")
            m = secondtry.match(start)
            if m:
                a, b = m.groups()
                return f"${a}_{{{b},{sub}}}$"
            else:
                return f"${start}_{sub}$"

        firstdim = add_subscript(var_dim, "i")
        self.rename(along_dim, firstdim)
        self.rename(var_dim, add_subscript(var_dim, "j"))
        self.set_units(firstdim, var_dim_units)
        return self

    def popdim(self, dimname):
        thisindex = self.axn(dimname)
        thisshape = list(self.data.shape)
        if thisshape[thisindex] != 1:
            raise IndexError("trying to pop a dim that's not length 1")
        thisshape.pop(thisindex)
        self.data = self.data.reshape(thisshape)
        self._pop_axis_info(thisindex)
        return self

    def cropped_log(self, subplot_axes=None, magnitude=4):
        r"""For the purposes of plotting, this generates a copy where I take
        the log, spanning "magnitude" orders of magnitude.
        This is designed to be called as abs(instance).cropped_log(), so it
        doesn't make a copy
        """
        phaseinfo = None
        if self.data.dtype == np.complex128:
            absdata = abs(self)
            phaseinfo = self / absdata
            self.data = absdata.data
        self.run(np.log10)
        if subplot_axes is None:  # then do all
            self.data -= (
                self.data[np.isfinite(self.data)].flatten().max() - magnitude
            )  # span only 4 orders of magnitude
        else:
            print("smooshing along subplot_axes", subplot_axes)
            newdata = self.copy().smoosh(subplot_axes, dimname="subplot")
            print(ndshape(newdata))
            newdata.run(max, "subplot")
            print(newdata)
            newdata = self - newdata
            self.data = newdata.data + magnitude
        self.data[self.data < 0] = 0
        if phaseinfo is not None:
            self.data = self.data * phaseinfo.data
        return self

    def runcopy(self, *args):
        newdata = self.copy()
        func = args[0]
        func = self._wrapaxisfuncs(func)
        if len(args) > 1:
            axis = args[1]
            thisindex = newdata.dimlabels.index(axis)
            newdata.data = func(newdata.data, axis=thisindex)
            newdata._pop_axis_info(thisindex)
        else:
            newdata.data = func(newdata.data)
        return newdata

    def run(self, *args):
        """run a standard numpy function on the nddata:

        ``d.run(func,'axisname')`` will run function `func` (*e.g.* a
        lambda function) along axis named 'axisname'

        ``d.run(func)`` will run function `func` on the data

        **in general**: if the result of func reduces a dimension size to
        1, the 'axisname' dimension will be "popped" (it will not exist in
        the result) -- if this is not what you want, use ``run_nopop``
        """
        func = args[0]
        func = self._wrapaxisfuncs(func)
        if len(args) > 1:
            axis = args[1]
            try:
                thisindex = self.dimlabels.index(axis)
            except Exception:
                if not isinstance(axis, str):
                    raise ValueError(
                        'The format of run is run(func,"axisname"), but you'
                        " didn't give a string as the second argument --"
                        " maybe you fed the arguments backwards?"
                    )
                elif axis not in self.dimlabels:
                    raise ValueError(
                        "axis "
                        + axis
                        + " is not in dimlabels ("
                        + repr(self.dimlabels)
                        + ")"
                    )
                else:
                    raise Exception()
            logger.debug(
                strm(
                    "called function",
                    func,
                    "on axis",
                    axis,
                    "which has index",
                    thisindex,
                )
            )
            logger.debug(
                strm(
                    "data type before",
                    type(self.data),
                    "shape before",
                    self.data.shape,
                )
            )
            shape_before = self.data.shape
            self.data = func(self.data, axis=thisindex)
            shape_after = self.data.shape
            logger.debug(
                strm("data of type", type(self.data), "looks like", self.data)
            )
            if len(shape_after) == len(shape_before) - 1:
                self._pop_axis_info(thisindex)
                logger.debug("popping axis info")
            return self
        else:
            retval = func(self.data)
            if self.data.size == retval.size:
                self.data = retval
                return self
            else:
                raise ValueError(
                    "the function you've "
                    "chosen doesn't return data that's the "
                    "same size as what you started with"
                )

    def run_nopop(self, func, axis):
        func = self._wrapaxisfuncs(func)
        try:
            thisaxis = self.axn(axis)
        except Exception:
            raise IndexError(
                strm(
                    "I couldn't find the dimension",
                    axis,
                    "in the list of axes",
                    self.dimlabels,
                )
            )
        temp = list(self.data.shape)
        temp[thisaxis] = 1
        func_sig = inspect.signature(func)
        numnonoptargs = len([
            v.default
            for v in func_sig.parameters.values()
            if v.default == inspect.Parameter.empty
        ])
        kwargnames = [
            k
            for k, v in func_sig.parameters.items()
            if v.default != inspect.Parameter.empty
        ]
        if numnonoptargs == 2:
            if "axis" in kwargnames:
                self.data = func(self.getaxis(axis), self.data, axis=thisaxis)
            if "axes" in kwargnames:
                self.data = func(self.getaxis(axis), self.data, axes=thisaxis)
            else:
                raise ValueError(
                    "Your function doesn't have axis or axes as a keyword"
                    " argument!"
                )
        else:
            if numnonoptargs == 1:
                paramnames = [k for k in func_sig.parameters.keys()]
                if len(paramnames) == 1:
                    if (
                        func_sig.parameters[paramnames[0]].kind
                        == inspect.Parameter.VAR_POSITIONAL
                    ):
                        try:
                            self.data = func(self.data, axis=thisaxis)
                        except Exception:
                            self.data = func(self.data, axes=thisaxis)
                if "axis" in kwargnames:
                    self.data = func(self.data, axis=thisaxis)
                elif "axes" in kwargnames:
                    self.data = func(self.data, axes=thisaxis)
                else:
                    raise ValueError(
                        "Your function doesn't have axis or axes as a keyword"
                        " argument! The number of non-optional arguments are"
                        " %s. The keyword arguments are %s"
                        % (str(numnonoptargs), str(kwargnames))
                    )
            else:
                raise ValueError(
                    "you passed a function to run_nopop that doesn't"
                    "have either one or two arguments!"
                )
        # {{{ if the function doesn't rip out the dim, make sure we don't
        #     change the dims
        if len(self.data.shape) == len(temp):
            temp[thisaxis] = self.data.shape[thisaxis]
        # }}}
        self.data = self.data.reshape(temp)
        return self

    def item(self):
        r"like numpy item -- returns a number when zero-dimensional"
        try:
            return self.data.item()
        except Exception:
            raise ValueError(
                "your data has shape: "
                + str(ndshape(self))
                + " so you can't call item"
            )

    # }}}
    # {{{ ft-related functions
    def unitify_axis(self, axis_name, is_axis=True):
        "this just generates an axis label with appropriate units"
        if type(axis_name) in [int, np.int64]:
            axis_name = self.dimlabels[axis_name]
        if (
            self.get_prop("FT") is not None
            and axis_name in list(self.get_prop("FT").keys())
            and self.get_prop("FT")[axis_name]
        ):
            isft = True
        else:
            isft = False
        if is_axis:
            yunits = self.units_texsafe(axis_name)
            ph_re = re.compile("^ph([0-9]|_.*)")
            m = ph_re.match(axis_name)
            if m:
                first_group = m.groups()[0]
                if first_group.startswith("_"):
                    first_group = first_group[1:]
                if len(first_group) > 1:
                    first_group = "{%s}" % first_group
                if isft:
                    axis_name = "$\\Delta p_%s$" % first_group
                else:
                    axis_name = "$\\varphi_%s$" % first_group
            else:
                auto_underscore = re.compile("^([a-z]+)([0-9])")
                m = auto_underscore.match(axis_name)
                if m:
                    axis_name = "%s_%s" % (m.groups())
                j = axis_name.find("_")
                if j > -1:
                    prevword = axis_name[0:j]
                    if j + 1 < len(axis_name):
                        followword = axis_name[j + 1 :]
                    else:
                        followword = []
                    k = followword.find(" ")
                    if k > -1 and k < len(followword):
                        followword = followword[:k]
                    k = followword.find("_")
                    if len(followword) > 0:
                        if not (k > -1) and (
                            len(prevword) < 2 or len(followword) < 2
                        ):
                            if len(followword) > 1:
                                axis_name = (
                                    axis_name[: j + 1 + len(followword)]
                                    + "}$"
                                    + axis_name[j + 1 + len(followword) :]
                                )
                                axis_name = (
                                    axis_name[: j + 1]
                                    + "{"
                                    + axis_name[j + 1 :]
                                )
                            else:
                                axis_name = (
                                    axis_name[0 : j + 2]
                                    + "$"
                                    + axis_name[j + 2 :]
                                )
                            axis_name = "$" + axis_name
                if isft:
                    t_idx = axis_name.find("t")
                    if t_idx > -1:
                        if (
                            t_idx + 1 < len(axis_name)
                            and axis_name[t_idx + 1].isalpha()
                        ):
                            axis_name = r"F{" + axis_name + r"}"
                        else:
                            axis_name = axis_name.replace("t", "\\nu ")
                            if axis_name[0] != "$":
                                axis_name = "$" + axis_name + "$"
                            axis_name = axis_name.replace(" _", "_")
                    else:
                        axis_name = r"F{" + axis_name + r"}"
        else:
            yunits = self.units_texsafe()
        if yunits is not None:
            axis_name = axis_name + " / " + yunits
        return axis_name

    # {{{ the following are all in the desired format -- the repetition at the
    #     end is because each function is in its own file (module) of the same
    #     name
    _ft_conj = this_fourier._ft_conj
    ft = this_fourier.ft
    set_ft_prop = this_fourier.set_ft_prop
    get_ft_prop = this_fourier.get_ft_prop
    ft_state_to_str = this_fourier.ft_state_to_str
    ft_clear_startpoints = this_fourier.ft_clear_startpoints
    ift = this_fourier.ift
    _ft_shift = this_fourier._ft_shift
    ftshift = this_fourier.ftshift
    convolve = this_fourier.convolve
    extend_for_shear = this_fourier.extend_for_shear
    linear_shear = axis_manipulation.linear_shear
    inhomog_coords = axis_manipulation.inhomog_coords
    secsy_transform_manual = axis_manipulation.secsy_transform_manual
    secsy_transform = axis_manipulation.secsy_transform
    register_axis = axis_manipulation.register_axis
    fourier_shear = this_fourier.shear
    # }}}
    # }}}
    nnls = MM_nnls

    # {{{ interpolation and binning
    def run_avg(self, thisaxisname, decimation=20, centered=False):
        "a simple running average"
        temp = self.getaxis(thisaxisname).size % decimation
        decimation = int(decimation)
        if temp != 0:
            if centered:
                self = self[thisaxisname, temp // 2 : -int(temp / 2.0 + 0.5)]
            else:
                self = self[thisaxisname, 0:-temp]
        thisaxis = nddata(self.getaxis(thisaxisname), [-1], [thisaxisname])
        self.setaxis(thisaxisname, [])
        self.chunkoff(thisaxisname, ["avg"], [decimation])
        self.run(np.mean, "avg")
        thisaxis.chunkoff(thisaxisname, ["avg"], [decimation])
        thisaxis.run(np.mean, "avg")
        self.setaxis(thisaxisname, thisaxis.data)
        return self

    def spline_lambda(self, s_multiplier=None):
        """For 1D data, returns a lambda function to generate a Cubic Spline.

        Parameters
        ==========
        s_multiplier: float
            If this is specified, then use a
            smoothing BSpline, and set "s" in
            scipy to the
            `len(data)*s_multiplier`
        Returns
        =======
        nddata_lambda: lambda function
            Takes one argument, which is an
            array corresponding to the axis
            coordinates, and returns an
            nddata.
        """
        assert len(self.dimlabels) == 1, "currently only supports 1D data"
        if s_multiplier is not None:
            thefunc = lambda x, y, s=0: scipy.interpolate.BSpline(
                *scipy.interpolate.splrep(x, y, s=s)
            )
            kwargs = dict(s=len(self.dimlabels[0]) * s_multiplier)
        else:
            thefunc = scipy.interpolate.CubicSpline
            kwargs = {}
        myspline_re = thefunc(
            self.getaxis(self.dimlabels[0]), self.data.real, **kwargs
        )
        if np.iscomplexobj(self.data.dtype):
            myspline_im = thefunc(
                self.getaxis(self.dimlabels[0]), self.data.imag, **kwargs
            )
            nddata_lambda = (
                lambda x: nddata(
                    myspline_re(x) + 1j * myspline_im(x), self.dimlabels[0]
                )
                .setaxis(self.dimlabels[0], x)
                .set_units(
                    self.dimlabels[0], self.get_units(self.dimlabels[0])
                )
            )
        else:
            nddata_lambda = (
                lambda x: nddata(myspline_re(x), self.dimlabels[0])
                .setaxis(self.dimlabels[0], x)
                .set_units(
                    self.dimlabels[0], self.get_units(self.dimlabels[0])
                )
            )
        return nddata_lambda

    def interp(
        self, axis, axisvalues, past_bounds=None, return_func=False, **kwargs
    ):
        """interpolate data values given axis values

        Parameters
        ==========
        return_func : boolean
            defaults to False.  If True, it returns a function that accepts
            axis values and returns a data value.
        """
        oldaxis = self.getaxis(axis)
        if not return_func:
            if (isinstance(axisvalues, int)) or (
                isinstance(axisvalues, np.int32)
            ):
                axisvalues = np.linspace(oldaxis[0], oldaxis[-1], axisvalues)
            elif np.isscalar(axisvalues):
                axisvalues = r_[axisvalues]
            elif type(axisvalues) not in [np.ndarray, tuple]:
                raise ValueError(
                    "You passed a target axis of type"
                    + repr(type(axisvalues))
                    + "which I don't get"
                )
            if np.any(np.imag(axisvalues) > 1e-38):
                raise ValueError("I can't interpolate imaginary values")
            else:
                axisvalues = np.real(axisvalues)
            axisvalues_final = axisvalues
            if past_bounds is None:
                axisvalues = axisvalues.copy()
                axisvalues[axisvalues < oldaxis.min()] = oldaxis.min()
                axisvalues[axisvalues > oldaxis.max()] = oldaxis.max()
            elif not (past_bounds == "fail"):
                if isinstance(past_bounds, tuple):
                    if len(past_bounds) == 2:
                        axisvalues[axisvalues < oldaxis.min()] = past_bounds[0]
                        axisvalues[axisvalues > oldaxis.max()] = past_bounds[1]
                    else:
                        raise TypeError(
                            "If you pass axisvalues as a tuple, it must be of"
                            " length 2!"
                        )
                else:
                    axisvalues[axisvalues < oldaxis.min()] = past_bounds
                    axisvalues[axisvalues > oldaxis.max()] = past_bounds
        rdata = np.real(self.data)
        idata = np.imag(self.data)
        thiserror = self.get_error()
        if thiserror is not None:
            rerrvar = np.real(thiserror) ** 2
            if thiserror[0].dtype == "complex128":
                ierrvar = np.imag(thiserror) ** 2
        if "kind" in list(kwargs.keys()):
            thiskind = kwargs.pop("kind")
        else:
            thiskind = "cubic"
            if len(rdata) < 4:
                thiskind = "quadratic"
                if len(rdata) < 3:
                    thiskind = "linear"
        thisaxis = self.axn(axis)
        logger.debug(strm("Using %s interpolation" % thiskind))

        def local_interp_func(local_arg_data, kind=thiskind):
            interpfunc = interp1d(
                oldaxis, local_arg_data, kind=kind, axis=thisaxis
            )
            try:
                retval = interpfunc(axisvalues)
            except Exception:
                raise TypeError("dtype of axis is" + repr(axisvalues.dtype))
            return retval

        if return_func:
            rfunc = interp1d(
                oldaxis,
                rdata,
                kind=thiskind,
                axis=thisaxis,
                bounds_error=False,
                fill_value=tuple(rdata[r_[0, -1]].tolist()),
            )
            ifunc = interp1d(
                oldaxis,
                idata,
                kind=thiskind,
                axis=thisaxis,
                bounds_error=False,
                fill_value=tuple(idata[r_[0, -1]].tolist()),
            )
            return lambda x: rfunc(x) + 1j * ifunc(x)
        rdata = local_interp_func(rdata)
        idata = local_interp_func(idata)
        self.data = rdata + 1j * idata
        self.setaxis(axis, axisvalues_final)
        if thiserror is not None:
            rerrvar = local_interp_func(
                rerrvar, kind="linear"
            )  # calculate the error variance of the real part, use linear to
            #    avoid nan problems
            if thiserror[0].dtype == "complex128":
                ierrvar = local_interp_func(ierrvar, kind="linear")
                self.set_error(sqrt(rerrvar) + 1j * sqrt(ierrvar))
            else:
                self.set_error(sqrt(rerrvar))
            err_nanmask = np.isnan(self.get_error())
            self.data[err_nanmask] = nan
        return self

    def invinterp(self, axis, values, **kwargs):
        "interpolate axis values given data values"
        copy = process_kwargs(
            [
                ("copy", False),
            ],
            kwargs,
            pass_through=True,
        )

        if np.isscalar(values):
            values = r_[values]
        origdata = self.data.copy()
        origaxis = self.getaxis(axis).copy()
        if np.any(np.imag(values) > 1e-38):
            raise ValueError("I can't interpolate imaginary values")
        else:
            values = np.real(values)
        args = origdata.argsort()
        origdata = origdata[args]
        rdata = np.real(origaxis)
        idata = np.imag(origaxis)
        rdata = rdata[args]
        idata = idata[args]
        # {{{ determine type of interpolation
        if "kind" in list(kwargs.keys()):
            thiskind = kwargs.pop("kind")
        else:
            thiskind = "cubic"
            if len(rdata) < 4:
                thiskind = "quadratic"
                if len(rdata) < 3:
                    thiskind = "linear"
        # }}}
        interpfunc = interp1d(origdata, rdata, kind=thiskind, **kwargs)
        try:
            rdata = interpfunc(values)
        except Exception:
            raise ValueError(
                strm(
                    "You passed",
                    values,
                    "and the data spans from",
                    origdata.min(),
                    "to",
                    origdata.max(),
                )
            )
        interpfunc = interp1d(origdata, idata, kind=thiskind, **kwargs)
        idata = interpfunc(values)
        cdata = rdata + 1j * idata
        if copy:
            mynewaxis = "name data before interpolation"
            if self.name() is not None:
                mynewaxis = self.name()
            retval = nddata(cdata, [-1], [mynewaxis])
            retval.labels(mynewaxis, values)
            return retval
        else:
            self.data = values
            self.setaxis(axis, cdata)
            return self

    def contiguous(self, lambdafunc, axis=None):
        r"""Return contiguous blocks that satisfy the condition given by
        `lambdafunc`

        this function returns the start and stop positions along the
        axis for the contiguous blocks for which lambdafunc returns
        true
        **Currently only supported for 1D data**

        .. note::
            adapted from \
                    stackexchange post http://stackoverflow.com/questions/\
                    4494404/find-large-number-of-consecutive-values-fulfilling\
                    -condition-in-a-numpy-array

        Parameters
        ----------

        lambdafunc : types.FunctionType
            If only one argument (lambdafunc) is given,
            then lambdafunc is
            a function that accepts a copy of the current nddata object
            (`self`) as the argument.
            If two arguments are given,
            the second is `axis`, and lambdafunc has two arguments,
            `self` and the value of `axis`.
        axis : {None,str}
            the name of the axis along which you want to find contiguous
            blocks

        Returns
        -------
        retval : np.ndarray
            An :math:`N\times 2` matrix, where the :math:`N` rows correspond to
            pairs of axis label that give ranges over which `lambdafunc`
            evaluates to `True`.
            These are ordered according to descending range width.

        Examples
        --------

        .. code:: python

            sum_for_contiguous = abs(forplot).mean('t1')
            fl.next("test contiguous")
            forplot = sum_for_contiguous.copy().set_error(None)
            fl.plot(forplot,alpha = 0.25,linewidth = 3)
            print("this is what the max looks like",0.5*sum_for_contiguous.\
                    set_error(None).runcopy(max,'t2'))
            print(sum_for_contiguous > 0.5*sum_for_contiguous.\
                    runcopy(max,'t2'))
            retval = sum_for_contiguous.contiguous(quarter_of_max,'t2')
            print("contiguous range / 1e6:",retval/1e6)
            for j in range(retval.shape[0]):
                a,b = retval[j,:]
                fl.plot(forplot['t2':(a,b)])

        """
        logger.debug(
            strm("(contiguous) shape of self inside contiguous", ndshape(self))
        )
        if axis is None:
            if len(self.dimlabels) == 1:
                axis = self.dimlabels[0]
                mask = lambdafunc(self.copy()).data
            else:
                raise TypeError(
                    "If there is more than one dimension, `axis` must be set"
                    " to something other than ``None``"
                )
        else:
            mask = lambdafunc(self.copy(), axis).data
        logger.debug(strm("(contiguous) shape of mask", mask.shape))
        if axis is None:
            (idx,) = np.diff(
                mask
            ).nonzero()  # this gives a list of indices for the boundaries
            #              between true/false
        else:
            (idx,) = np.diff(
                mask, axis=self.axn(axis)
            ).nonzero()  # this gives a list of indices for the boundaries
            #              between true/false
        idx += 1  # because np.diff only starts on index #1 rather than #0
        if mask[
            0
        ]:  # because I want to indicate the boundaries of True, if I am
            # starting on True, I need to make 0 a boundary
            idx = np.r_[0, idx]
        if mask[
            -1
        ]:  # If the end of mask is True, then I need to add a boundary there
            # as well
            idx = np.r_[idx, mask.size - 1]  # Edit
        idx.shape = (-1, 2)  # idx is 2x2 np.array of start,stop
        logger.debug(strm("(contiguous) DEBUG idx is", idx))
        logger.debug(
            strm("(contiguous) diffs for blocks", np.diff(idx, axis=1))
        )
        block_order = np.diff(idx, axis=1).flatten().argsort()[::-1]
        logger.debug(
            strm(
                "(contiguous) in descending order, the blocks are therefore",
                idx[block_order, :],
            )
        )
        return self.getaxis(axis)[idx[block_order, :]]

    def to_ppm(self, axis="t2", freq_param="SFO1", offset_param="OFFSET"):
        """Function that converts from Hz to ppm using Bruker parameters

        Parameters
        ==========
        axis: str, 't2' default
            label of the dimension you want to convert from frequency to ppm
        freq_param: str
            name of the acquisition parameter that stores the carrier frequency
            for this dimension
        offset_param: str
            name of the processing parameter that stores the offset of the ppm
            reference (TMS, DSS, etc.)

        .. todo::

            Figure out what the units of PHC1 in Topspin are (degrees per
            *what??*), and apply those as well.

            make this part of an inherited bruker class
        """
        if self.get_units(axis) == "ppm":
            return
        offset = self.get_prop("proc")[offset_param]
        logger.debug(strm("offset", offset, "not used"))
        SF = self.get_prop("proc")["SF"]
        sfo1 = self.get_prop("acq")[freq_param]
        tms_hz = (SF - sfo1) * 1e6
        if not self.get_ft_prop(axis):
            self.ft(
                axis, shift=True
            )  # this fourier transforms along t2, overwriting the data that
            #    was in self
        if axis == "t2":
            self.setaxis(axis, lambda x: x - tms_hz)
            self.setaxis(axis, lambda x: x / SF)
            self.set_units(axis, "ppm")
            self.set_prop("x_inverted", True)
        elif axis == "t1":
            self.setaxis(axis, lambda x: x / sfo1)
            self.set_units(axis, "ppm")
            max_ppm = self.getaxis(axis).max()
            self.setaxis(axis, lambda x: (x - max_ppm + offset))
            self.set_prop("y_inverted", True)
        return self

    # }}}

    def repwlabels(self, axis):
        return None

    def add_noise(self, intensity):
        """Add Gaussian (box-muller) noise to the data.

        Parameters
        ----------
        intensity : double OR function
            If a double, gives the standard deviation of the noise.
            If a function, used to calculate the standard deviation of the
            noise from the data:
            *e.g.* ``lambda x: max(abs(x))/10.``

        """
        if isinstance(intensity, type(emptyfunction)):
            intensity = intensity(lambda x: self.data)
        return_complex = np.iscomplexobj(self.data)
        if return_complex:
            self.data += intensity * np.random.normal(
                size=self.data.shape
            ) + 1j * intensity * np.random.normal(size=self.data.shape)
        else:
            self.data += intensity * np.random.normal(size=self.data.shape)
        return self

    # {{{ functions to manipulate and return the axes
    def reorder(self, *axes, **kwargs):
        r"""Reorder the dimensions
        the first arguments are a list of dimensions

        Parameters
        ----------
        *axes : str
            Accept any number of arguments that gives the dimensions, in the
            order that you want thee.
        first : bool
            (default True)
            Put this list of dimensions first, while False puts them last
            (where they then come in the order given).
        """
        first = True
        if "first" in kwargs:
            first = kwargs.pop("first")
        if len(kwargs) > 0:
            raise ValueError("I don't understand your kwargs!")
        if len(axes) == 1:
            axes = axes[0]
        else:
            axes = axes
        if isinstance(axes, str):
            axes = [axes]
        if len(axes) < len(self.dimlabels):
            oldorder = list(self.dimlabels)
            for thisaxis in axes:
                oldorder.pop(oldorder.index(thisaxis))
            if first:
                axes = axes + oldorder
            else:
                axes = oldorder + axes
        try:
            neworder = list(map(self.dimlabels.index, axes))
        except ValueError:
            raise ValueError(strm("one of", axes, "not in", self.dimlabels))
        self.dimlabels = list(map(self.dimlabels.__getitem__, neworder))
        if len(self.axis_coords) > 0:
            try:
                self.axis_coords = list(
                    map(self.axis_coords.__getitem__, neworder)
                )
            except Exception:
                raise IndexError(
                    strm(
                        "problem mapping",
                        list(map(len, self.axis_coords)),
                        "onto",
                        neworder,
                    )
                )
            if len(self.axis_coords_units) > 0:
                try:
                    self.axis_coords_units = list(
                        map(self.axis_coords_units.__getitem__, neworder)
                    )
                except Exception:
                    raise IndexError(
                        strm(
                            "problem mapping",
                            list(map(len, self.axis_coords_units)),
                            "onto",
                            neworder,
                        )
                    )
        try:
            self.data = self.data.transpose(neworder)
        except ValueError:
            raise ValueError(
                strm("you can't reorder", self.dimlabels, "as", neworder)
            )
        if self.data_error is not None:
            self.data_error = self.data_error.transpose(neworder)
        return self

    def plot_labels(self, labels, fmt=None, **kwargs_passed):
        r"this only works for one axis now"
        axisname = self.dimlabels[0]
        if fmt is None:
            plot_label_points(
                self.getaxis(axisname), self.data, labels, **kwargs_passed
            )
        else:
            plot_label_points(
                self.getaxis(axisname),
                self.data,
                [fmt % j for j in labels],
                **kwargs_passed,
            )
        return

    def labels(self, *args):
        r"""label the dimensions, given in listofstrings with the axis labels
        given in listofaxes -- listofaxes must be a numpy np.array; you can
        pass either a dictionary or a axis name (string)/axis label (numpy
        np.array) pair
        """
        if len(args) == 2:
            listofstrings, listofaxes = args
        elif len(args) == 1 and isinstance(args[0], dict):
            listofstrings = list(args[0].keys())
            listofaxes = list(args[0].values())
        else:
            raise ValueError(
                strm("I can't figure out how to deal with the arguments", args)
            )
        for j in range(0, len(listofaxes)):
            if isinstance(listofaxes[j], list):
                listofaxes[j] = np.array(listofaxes[j])
        listofstrings = autostringconvert(listofstrings)
        if isinstance(listofstrings, str):
            listofstrings = [listofstrings]
            listofaxes = [listofaxes]
        if not isinstance(listofstrings, list):
            raise TypeError(
                "the arguments passed to the .labels() method must be a list"
                " of the axis names followed by the list of the axis arrays"
            )
        elif all(map((lambda x: isinstance(x, np.str_)), listofstrings)):
            listofstrings = list(map(str, listofstrings))
        elif not all(map((lambda x: isinstance(x, str)), listofstrings)):
            raise TypeError(
                "the arguments passed to the .labels() method must be a list"
                " of the axis names followed by the list of the axis arrays"
            )
        for j in range(0, len(listofstrings)):
            if listofaxes[j] is None:
                self.setaxis(listofstrings[j], None)
            else:
                # {{{ test that the axis is the right size
                if np.isscalar(listofaxes[j]):  # interpret as a timestep
                    listofaxes[j] = (
                        listofaxes[j] * r_[0 : ndshape(self)[listofstrings[j]]]
                    )
                if type(listofaxes[j]) not in [np.ndarray, list]:
                    raise TypeError(
                        "You passed an axis label of type "
                        + repr(type(listofaxes[j]))
                        + " for the axis "
                        + listofstrings[j]
                        + " to the labels method, which you can't do --> it"
                        " must be an nddata"
                    )
                if (
                    len(listofaxes[j]) != ndshape(self)[listofstrings[j]]
                ) and (len(listofaxes[j]) != 0):
                    raise IndexError(
                        "You're trying to attach an axis of len %d to the '%s'"
                        " dimension, which has %d data points (shape of self"
                        " is %s)"
                        % (
                            len(listofaxes[j]),
                            listofstrings[j],
                            ndshape(self)[listofstrings[j]],
                            repr(ndshape(self)),
                        )
                    )
                # }}}
                self.setaxis(listofstrings[j], listofaxes[j])
        return self

    def check_axis_coords_errors(self):
        if len(self.axis_coords_error) > len(self.dimlabels):
            raise ValueError(
                "this failed because there are more sets of axis errors than"
                " there are axes!\nlen(axis_coords_error) = %s\naxes = %s"
                % (repr(len(self.axis_coords_error)), repr(self.dimlabels))
            )

    def sort(self, axisname, reverse=False):
        whichaxis = self.dimlabels.index(axisname)
        if reverse:
            order = np.argsort(-1 * self.axis_coords[whichaxis])
        else:
            order = np.argsort(self.axis_coords[whichaxis])
        datacopy = self.copy()
        for j in range(
            0, len(order)
        ):  # do it this way, so that it deals with other dimensions correctly
            self.check_axis_coords_errors()
            self[axisname, j] = datacopy[axisname, order[j]]
        self.axis_coords[whichaxis] = self.axis_coords[whichaxis][order]
        return self

    def copyaxes(self, other):
        raise ValueError("use copy_axes")

    def copy_axes(self, other):
        # in the case that the dimensions match, and we want to copy the labels
        for thisdim in self.dimlabels:
            if thisdim in other.dimlabels:
                thisax = other.getaxis(thisdim)
                if thisax is not None:
                    thisax = thisax.copy()
                self.setaxis(thisdim, thisax)
                if other.get_error(thisdim) is not None:
                    self.set_error(thisdim, np.copy(other.get_error(thisdim)))
                if other.get_units(thisdim) is not None:
                    self.set_units(thisdim, other.get_units(thisdim))
        return self

    def axis(self, axisname):
        "returns a 1-D axis for further manipulation"
        return nddata(self.getaxis(axisname).copy(), [-1], [axisname]).labels(
            axisname, self.getaxis(axisname).copy()
        )

    def _axis_inshape(self, axisname):
        newshape = np.ones(len(self.data.shape), dtype="uint")
        thisaxis = self.axn(axisname)
        newshape[thisaxis] = self.data.shape[thisaxis]
        newshape = list(newshape)
        retval = self.getaxis(axisname)
        if retval is None:
            raise AttributeError(axisname + " does not have axis labels!")
        try:
            return retval.copy().reshape(newshape)
        except ValueError:
            raise ValueError(
                strm(
                    "Trying to reshape axis from",
                    retval.shape,
                    "to",
                    newshape,
                    "so I can manipulate it like data",
                )
            )

    def retaxis(self, axisname):
        thisaxis = self._axis_inshape(axisname)
        return nddata(thisaxis, thisaxis.shape, list(self.dimlabels)).labels(
            axisname, thisaxis.flatten()
        )

    def fromaxis(self, *args, **kwargs):
        """Generate an nddata object from one of the axis labels.

        Can be used in one of several ways:

        * ``self.fromaxis('axisname')``: Returns an nddata where `retval.data`
          consists of the given axis values.
        * ``self.fromaxis('axisname',inputfunc)``: use `axisname` as the input
          for `inputfunc`, and load the result into `retval.data`
        * ``self.fromaxis(inputsymbolic)``: Evaluate `inputsymbolic` and load
          the result into `retval.data`

        Parameters
        ==========
        axisname : str | list
            The axis (or list of axes) to that is used as the argument of
            `inputfunc` or the function represented by `inputsymbolic`.
            If this is the only argument, it cannot be a list.
        inputsymbolic : sympy.Expr
            A sympy expression whose only symbols are the names of axes.
            It is preferred, though not required, that this is passed
            without an `axisname` argument -- the axis names are then
            inferred from the symbolic expression.
        inputfunc : function
            A function (typically a lambda function) that taxes the values of
            the axis given by `axisname` as input.
        overwrite : bool
            Defaults to `False`. If set to `True`, it overwrites `self` with
            `retval`.
        as_array : bool
            Defaults to `False`. If set to `True`, `retval` is a properly
            dimensioned numpy ndarray rather than an nddata.

        Returns
        =======
        retval : nddata | ndarray
            An expression calculated from the axis(es) given by `axisname` or
            inferred from `inputsymbolic`.
        """
        overwrite, as_array = process_kwargs(
            [("overwrite", False), ("as_array", False)], kwargs
        )
        if len(args) == 1:
            if isinstance(args[0], str):
                axisname = args[0]
                # {{{ copied from old retaxis function, then added the
                #     overwrite capability
                if overwrite:
                    retval = self.retaxis(axisname)
                    thisaxis = self._axis_inshape(axisname)
                    self.data = thisaxis
                    return self
                else:
                    axis_data = self.getaxis(axisname).flatten()
                    # copy is needed here, or data and axis will be the same
                    # object
                    retval = nddata(
                        axis_data, axis_data.shape, [axisname]
                    ).setaxis(axisname, np.copy(axis_data))
                    retval.set_units(axisname, self.get_units(axisname))
                    retval.data_units = self.data_units
                    retval.name(self.name())
                    retval.copy_props(
                        self
                    )  # be sure to include info about ft startpoint
                    return retval
                # }}}
            else:
                if issympy(args[0]):
                    func = args[0]
                    symbols_in_func = func.atoms(sp.Symbol)
                    logger.debug(
                        strm(
                            "identified this as a sympy expression (",
                            func,
                            ") with symbols",
                            symbols_in_func,
                        )
                    )
                    symbols_not_in_dimlabels = set(
                        map(str, symbols_in_func)
                    ) - set(self.dimlabels)
                    if len(symbols_not_in_dimlabels) > 0:
                        raise ValueError(
                            "You passed a symbolic function, but the symbols"
                            + str(symbols_not_in_dimlabels)
                            + " are not axes"
                        )
                else:
                    raise ValueError(
                        "I don't know what to do with this type of argument!"
                    )
        elif len(args) == 2:
            axisnames = args[0]
            func = args[1]
            if not isinstance(axisnames, list):
                axisnames = [axisnames]
        else:
            raise ValueError(
                "Wrong number of arguments!! -- you passed "
                + repr(len(args))
                + " arguments!"
            )
        if issympy(func):
            logging.debug(
                strm(
                    "about to run sympy sp.utilities.lambdify,"
                    " symbols_in_func is",
                    symbols_in_func,
                )
            )
            try:
                lambdified_func = sp.utilities.lambdify(
                    list(symbols_in_func), func, modules=mat2array
                )
            except Exception:
                raise ValueError(
                    strm(
                        "Error parsing axis variables",
                        list(map(sp.core.var, axisnames)),
                        "that you passed and function",
                        func,
                        "that you passed",
                    )
                )
            func = lambdified_func
            axisnames = list(map(str, symbols_in_func))
        elif not hasattr(func, "__call__"):
            raise ValueError(
                "I can't interpret the second argument as a function! It is"
                " type "
                + str(type(func))
            )
        # I can't do the following for sympy, because the argument count is
        # always zero
        if not issympy(args[0]) and func.__code__.co_argcount != len(
            axisnames
        ):
            raise ValueError(
                strm(
                    "The axisnames you passed",
                    axisnames,
                    "and the argument count",
                    func.__code__.co_argcount,
                    "don't match",
                )
            )
        list_of_axes = [self._axis_inshape(x) for x in axisnames]
        retval = func(*list_of_axes)
        if issympy(retval):
            raise RuntimeError(
                "The sympy function that you passed doesn't match the"
                " automatically generated axis variables (obtained by mapping"
                " sympy.var onto the axis variables, without any kwargs). The"
                " atoms left over are:\n"
                + str(func.atoms)
            )
        logging.debug(strm("at this point, list of axes is:", list_of_axes))
        if len(list_of_axes) == 0:
            return nddata(float(func()))
        newshape = np.ones_like(list_of_axes[0].shape)
        for j in list_of_axes:
            newshape *= np.array(j.shape)
        if overwrite:
            self.data = retval.reshape(newshape)
            return self
        else:
            if as_array:
                return retval.reshape(newshape)
            else:
                retval = nddata(retval, newshape, list(self.dimlabels)).labels(
                    axisnames, [self.getaxis(x).copy() for x in axisnames]
                )
                retval.axis_coords_units = list(self.axis_coords_units)
                retval.data_units = self.data_units
                retval.name(self.name())
                retval.copy_props(
                    self
                )  # be sure to include info about ft startpoint
                return retval

    def getaxis(self, axisname):
        if self.axis_coords is None or len(self.axis_coords) == 0:
            return None
        else:
            retval = self.axis_coords[self.axn(axisname)]
        if retval is None:
            return None
        elif len(retval) > 0:
            return retval
        else:
            return None

    def extend(self, axis, extent, fill_with=0, tolerance=1e-5):
        r"""Extend the (domain of the) dataset and fill with a pre-set value.

        The coordinates associated with
        `axis` must be uniformly ascending with spacing :math:`dx`.
        The function will extend `self`
        by adding a point every :math:`dx` until the axis
        includes the point `extent`.  Fill the newly created datapoints with
        `fill_with`.

        Parameters
        ----------

        axis : str
            name of the axis to extend
        extent : double
            Extend the axis coordinates of `axis` out to this value.

            The value of `extent` must be less the smallest (most negative)
            axis coordinate or greater than the largest (most positive)
            axis coordinate.
        fill_with : double
            fill the new data points with this value (defaults to 0)
        tolerance : double
            when checking for ascending axis labels, etc.,
            values/differences must match to within tolerance
            (assumed to represent the actual precision, given
            various errors, etc.)
        """
        # check for uniformly ascending
        u = self.getaxis(axis)
        du = (u[-1] - u[0]) / (len(u) - 1.0)
        thismsg = (
            "In order to expand, the axis must be ascending (and equally"
            " spaced)"
        )
        assert all(abs(np.diff(u) - du) / du < tolerance), thismsg  # absolute
        # figure out how many points I need to add, and on which side of the
        # axis
        assert du > 0, thismsg  # ascending
        start_index = 0
        stop_index = len(u)  # this is the index at which the data
        #                     stops.  To start with, we assume the
        #                     data stops where the original
        #                     position stops, and if needed, we
        #                     added points with stop_index_addto
        logger.debug(
            strm(
                "attempting to extend axis that runs from",
                u[0],
                "to",
                u[-1],
                "out to",
                extent,
            )
        )
        if extent < u[0]:
            start_index = -int(
                (u[0] - extent) // du
            )  # the part after the negative is positive
            if (
                start_index * du + (u[0] - extent)
            ) / du < -tolerance:  # the first quantity here is negative
                start_index -= 1
        elif extent > u[-1]:
            stop_index_addto = int((extent - u[-1]) // du)
            if (
                (extent - u[-1]) - du * stop_index_addto
            ) / du > tolerance:  # the first quantity here is negative
                stop_index_addto += 1
            stop_index += stop_index_addto
        else:
            raise RuntimeError(
                "extent ({:g}) needs to be further than the bounds on '{:s}',"
                " which are {:g} and {:g}".format(extent, axis, u[0], u[-1])
            )
        # {{{ create a new np.array, and put self.data into it
        newdata = list(self.data.shape)
        newdata[self.axn(axis)] = stop_index - start_index
        if fill_with == 0:
            newdata = np.zeros(newdata, dtype=self.data.dtype)
        else:
            newdata = fill_with * np.ones(newdata, dtype=self.data.dtype)
        newdata_slice = [slice(None, None, None)] * len(newdata.shape)
        # since start_index is negative, -start_index points have been added to
        # the beginning of the data (and the original data is len(u) in length)
        newdata_slice[self.axn(axis)] = slice(
            -start_index, len(u) - start_index, None
        )
        newdata[tuple(newdata_slice)] = self.data
        self.data = newdata
        # }}}
        # construct the new axis
        new_u = u[0] + du * r_[start_index:stop_index]
        self.setaxis(axis, new_u)
        if start_index < 0:
            # if we are extending to negative values, we need to inform
            # the FT machinery!
            if self.get_ft_prop(axis):
                self.set_ft_prop(axis, ["start", "freq"], new_u[0])
            else:
                self.set_ft_prop(axis, ["start", "time"], new_u[0])
        return self

    def setaxis(self, *args):
        """set or alter the value of the coordinate axis

        Can be used in one of several ways:

        * ``self.setaxis('axisname', values)``: just sets the values
        * ``self.setaxis('axisname', '#')``: just
            number the axis in numerically increasing order,
            with integers,
            (e.g. if you have smooshed it from a couple
            other dimensions.)
        * ``self.fromaxis('axisname',inputfunc)``: take the existing function,
            apply inputfunc, and replace
        * ``self.fromaxis(inputsymbolic)``: Evaluate `inputsymbolic` and load
            the result into the axes, appropriately
        """
        if len(args) == 2:
            axis, value = args
            if np.isscalar(value) and value == "#":
                self.setaxis(axis, r_[0 : ndshape(self)[axis]])
                return self
        elif len(args) == 1 and issympy(args[0]):
            func = args[0]
            symbols_in_func = func.atoms(sp.Symbol)
            logger.debug(
                strm(
                    "identified this as a sympy expression (",
                    func,
                    ") with symbols",
                    symbols_in_func,
                )
            )
            symbols_not_in_dimlabels = set(map(str, symbols_in_func)) - set(
                self.dimlabels
            )
            if len(symbols_not_in_dimlabels) > 0:
                raise ValueError(
                    "You passed a symbolic function, but the symbols"
                    + str(symbols_not_in_dimlabels)
                    + " are not axes"
                )
            logging.debug(
                strm(
                    "about to run sympy sp.utilities.lambdify,"
                    " symbols_in_func is",
                    symbols_in_func,
                )
            )
            try:
                lambdified_func = sp.utilities.lambdify(
                    list(symbols_in_func), func, modules=mat2array
                )
            except Exception:
                raise ValueError(
                    strm(
                        "Error parsing axis variables",
                        list(symbols_in_func),
                        "that you passed and function",
                        func,
                        "that you passed",
                    )
                )
            value = lambdified_func
            axis = list(map(str, symbols_in_func))
            assert len(axis) == 1, (
                "currently only supported for 1 axis at a time -- if you want"
                " to do for more than one axis, please create a pull request"
                " with an example"
            )
            axis = axis[0]
        else:
            raise ValueError(
                "not a valid argument to setaxis -- look at the documentation!"
            )
        if axis == "INDEX":
            raise ValueError(
                "Axes that are called INDEX are special, and you are not"
                " allowed to label them!"
            )
        if isinstance(value, type(emptyfunction)):
            x = self.getaxis(axis)
            x[:] = value(x.copy())
            return self
        if type(value) in [float, int, np.double, np.float64]:
            value = np.linspace(0.0, value, self.axlen(axis))
        if isinstance(value, list):
            value = np.array(value)
        if self.axis_coords is None or len(self.axis_coords) == 0:
            self.axis_coords = [None] * len(self.dimlabels)
            self.axis_coords_error = [None] * len(self.dimlabels)
        if value is None:
            self.axis_coords[self.axn(axis)] = None
        else:
            a = len(value)
            b = self.data.shape[self.axn(axis)]
            assert a == b, (
                "Along the axis %s, the length of the axis you passed (%d)"
                " doesn't match the size of the data (%d)." % (axis, a, b)
            )
            self.axis_coords[self.axn(axis)] = value
        return self

    def shear(
        self,
        along_axis,
        propto_axis,
        shear_amnt,
        zero_fill=True,
        start_in_conj=False,
        method="linear",
    ):
        r"""Shear the data :math:`s`:

        :math:`s(x',y,z) = s(x+ay,y,z)`

        where :math:`x` is the `altered_axis` and :math:`y` is the
        `propto_axis`.  (Actually typically 2D, but :math:`z` included
        just to illustrate other dimensions that aren't involved)

        .. note: Unfortunately, currently, when the data is automatically
            extended, if both the start and endpoint of `along_axis` are on the
            same side of zero, some unnecessary padding will be added between
            the beginning of `along_axis` and zero.

        Parameters
        ----------

        method : {'fourier','linear'}

            fourier
                Use the Fourier shift theorem (*i.e.*, sinc interpolation).  A
                shear is equivalent to the following in the conjugate domain:

                ..math: `\tilde{s}(f_x,f'_y,z) = \tilde{s}(f_x,f_y-af_x,f_z)`

                Because of this, the algorithm **also**
                automatically `extend`s the data in `f_y` axis.
                Equivalently, it increases the resolution
                (decreases the interval between points) in the
                `propto_axis` dimension.  This prevents aliasing
                in the conjugate domain, which will corrupt the
                data *w.r.t.* successive transformations. It does
                this whether or not `zero_fill` is set
                (`zero_fill` only controls filling in the
                "current" dimension)

            linear
                Use simple linear interpolation.

        altered_axis : str

            The coordinate for which data is altered, *i.e.*
            ..math: `x` such that ..math: `f(x+ay,y)`.

        by_amount : double

            The amount of the shear (..math: `a` in the previous)

        propto_axis : str

            The shift along the `altered_axis` dimension is
            proportional to the shift along `propto_axis`.
            The position of data relative to the `propto_axis` is not
            changed.
            Note that by the shift theorem, in the frequency domain,
            an equivalent magnitude, opposite sign, shear is applied
            with the `propto_axis` and `altered_axis` dimensions
            flipped.

        start_in_conj : {False, True}, optional

            Defaults to False

            For efficiency, one can replace a double (I)FT call followed by a
            shear call with a single shear call where `start_in_conj` is set.

            `self` before the call is given in the conjugate domain  (*i.e.*,
            :math:`f` *vs.* :math:`t`) along both dimensions from the one
            that's desired.  This means: (1) `self` after the function call
            transformed into the conjugate domain from that before the call and
            (2) `by_amount`, `altered_axis`, and `propto_axis` all refer to the
            shear in the conjugate domain that the data is in at the end of the
            function call.
        """
        if not (
            self.get_ft_prop(along_axis) is None
            and self.get_ft_prop(propto_axis) is None
        ):
            if (
                self.get_ft_prop(along_axis) ^ self.get_ft_prop(propto_axis)
            ) ^ start_in_conj:
                if start_in_conj:
                    raise ValueError(
                        "if you pass start_in_conj, the two dimensions need to"
                        " be in conjugate domains, but you have: "
                        + self.ft_state_to_str(along_axis, propto_axis)
                    )
                else:
                    raise ValueError(
                        "(unless you intended to pass start_in_conj) the two"
                        " dimensions need to be in the same domain, but you"
                        " have: "
                        + self.ft_state_to_str(along_axis, propto_axis)
                    )
        if method == "fourier":
            return self.fourier_shear(
                along_axis, propto_axis, shear_amnt, zero_fill=zero_fill
            )
        elif method == "linear":
            return self.linear_shear(
                along_axis, propto_axis, shear_amnt, zero_fill=zero_fill
            )
        else:
            raise ValueError(
                "The shear method must be either linear or fourier"
            )

    def getaxisshape(self, axisname):
        thishape = np.ones(len(self.dimlabels))
        thisaxis = self.dimlabels.index(axisname)
        thishape[thisaxis] = self.data.shape[thisaxis]
        return thishape

    def circshift(self, axis, amount):
        if amount != 0:
            if abs(amount) > ndshape(self)[axis]:
                ValueError(
                    strm(
                        "Trying to circshift by ",
                        amount,
                        "which is bitter than the size of",
                        axis,
                    )
                )
            newdata = ndshape(self).alloc(dtype=self.data.dtype)
            newdata[axis, :-amount] = self[axis, amount:]
            newdata[axis, -amount:] = self[axis, :amount]
            self.data = newdata.data
        return self

    # }}}
    # {{{ breaking up and combining axes
    def smoosh(self, dimstocollapse, dimname=0, noaxis=False):
        r"""Collapse (smoosh) multiple dimensions into one dimension.

        Parameters
        ----------
        dimstocollapse : list of strings
            the dimensions you want to collapse to one result dimension
        dimname : None, string, integer (default 0)

            if dimname is:

            * None: create a new (direct product) name,
            * a number: an index to the ``dimstocollapse`` list.  The resulting
                smooshed dimension will be named ``dimstocollapse[dimname]``.
                Because the default is the number 0, the new dimname will be
                the first dimname given in the list.
            * a string: the name of the resulting smooshed dimension (can be
                part of the ``dimstocollapse`` list or not)

        noaxis : bool
            if set, then just skip calculating the axis for the new dimension,
            which otherwise is typically a complicated record array

        Returns
        -------
        self: nddata
            the dimensions `dimstocollapse` are smooshed into a single
            dimension, whose name is determined by `dimname`.
            The axis for the resulting, smooshed dimension is a structured
            np.array consisting of two fields that give the labels along the
            original axes.

        ..todo::
            when we transition to axes that are stored using a
            slice/linspace-like format,
            allow for smooshing to determine a new axes that is standard
            (not a structured np.array) and that increases linearly.
        """
        assert (type(dimstocollapse) in [list, tuple]) and len(
            dimstocollapse
        ) > 1, (
            "What?? You must try to collapse more than one dimension!! -- you"
            " claim you want to collapse '%s'"
            % str(dimstocollapse)
        )
        not_present = set(dimstocollapse) - set(self.dimlabels)
        if len(not_present) > 0:
            raise ValueError(
                strm(
                    not_present,
                    "was not found in the list of dimensions",
                    self.dimlabels,
                )
            )
        # {{{ first, put them all at the end, in order given here
        retained_dims = list(self.dimlabels)
        logger.debug(strm("old order", retained_dims))
        # {{{ if I'm using a dimension here, be sure to grab its current
        #     position
        if dimname is None:
            final_position = -1
            dimname = " $\\times$ ".join(dimstocollapse)
        else:
            if isinstance(dimname, int):
                dimname = dimstocollapse[dimname]
                final_position = self.axn(dimname)
            elif dimname in self.dimlabels:
                final_position = self.axn(dimname)
            else:
                final_position = -1
        # }}}
        # {{{ store the dictionaries for later use
        axis_coords_dict = self.mkd(self.axis_coords)
        axis_coords_error_dict = self.mkd(self.axis_coords_error)
        axis_coords_units_dict = self.mkd(self.axis_coords_units)
        # }}}
        old_units = []
        logger.debug(strm("dims to collapse", dimstocollapse))
        for this_name in dimstocollapse:
            this_idx = retained_dims.index(this_name)
            retained_dims.pop(this_idx)
            axis_coords_error_dict.pop(this_name)
            old_units.append(axis_coords_units_dict.pop(this_name))
            axis_coords_dict.pop(this_name)
            if this_idx < final_position:
                final_position -= 1
        logger.debug(strm("old units", old_units))
        new_units = list(set(old_units))
        if len(new_units) > 1:
            new_units = " ".join(map(str, new_units))
        elif new_units == 1:
            new_units = new_units[0]
        else:
            new_units = None
        # this might be sub-optimal, but put the dims to collapse at the end,
        # and move them back later if we want
        new_order = retained_dims + dimstocollapse
        self.reorder(new_order)
        logger.debug(strm("new order", new_order))
        # }}}
        # {{{ then, reshape the data (and error)
        logger.debug(strm("old shape", self.data.shape))
        new_shape = list(self.data.shape)[: -len(dimstocollapse)]
        logger.debug(strm("dimensions to keep", new_shape))
        dimstocollapse_shapes = np.array(
            self.data.shape[-len(dimstocollapse) :]
        )
        new_shape += [dimstocollapse_shapes.prod()]
        self.data = self.data.reshape(new_shape)
        if self.get_error() is not None:
            self.set_error(self.get_error().reshape(new_shape))
        logger.debug(strm("new shape", self.data.shape))
        # }}}
        # {{{ now for the tricky part -- deal with the axis labels
        # {{{ in order, make a list of the relevant axis names, dtypes, and
        #     sizes
        axes_with_labels = [
            j for j in dimstocollapse if self.getaxis(j) is not None
        ]  # specifically, I am only concerned with the np.ones I am collapsing
        #    that have labels
        if noaxis:
            logger.debug("noaxis was specified")
            axis_coords_dict[dimname] = None
            axis_coords_error_dict[dimname] = None
        else:
            logger.debug("starting construction of the smooshed axis")
            axes_with_labels_haserror = [
                self.get_error(j) is not None for j in axes_with_labels
            ]
            axes_with_labels_dtype = [
                (j, self.getaxis(j).dtype) for j in axes_with_labels
            ]  # an appropriate spec. for a structured np.array
            axes_with_labels_size = [
                self.getaxis(j).size for j in axes_with_labels
            ]
            # }}}
            logger.debug(
                strm("the dtype that I want is:", axes_with_labels_dtype)
            )
            logger.debug(
                strm("the axes that have labels are:", axes_with_labels)
            )
            logger.debug(
                strm(
                    "the axes that have labels have sizes:",
                    axes_with_labels_size,
                )
            )
            # {{{ we construct a multidimensional axis
            multidim_axis_error = None
            if len(axes_with_labels_dtype) > 0:
                # create a new axis of the appropriate shape and size
                multidim_axis_label = np.empty(
                    axes_with_labels_size, dtype=axes_with_labels_dtype
                )
                if np.any(axes_with_labels_haserror):
                    multidim_axis_error = np.empty(
                        axes_with_labels_size,
                        dtype=[
                            (
                                axes_with_labels[j],
                                self.getaxis(axes_with_labels[j]).dtype,
                            )
                            for j in range(len(axes_with_labels))
                            if axes_with_labels_haserror[j]
                        ],
                    )
                # one at a time index the relevant dimension, and load in the
                # information
                full_slice = [slice(None, None, None)] * len(
                    axes_with_labels_dtype
                )
                for this_index, thisdim in enumerate(axes_with_labels):
                    axis_for_thisdim = self.getaxis(thisdim)
                    if axes_with_labels_haserror[this_index]:
                        axis_error_for_thisdim = self.get_error(thisdim)
                    logger.debug(
                        strm("the axis for", thisdim, "is", axis_for_thisdim)
                    )
                    for j in range(axes_with_labels_size[this_index]):
                        this_slice = list(full_slice)
                        this_slice[this_index] = j  # set this element
                        multidim_axis_label[thisdim][tuple(this_slice)] = (
                            axis_for_thisdim[j]
                        )
                        if axes_with_labels_haserror[this_index]:
                            multidim_axis_error[thisdim][tuple(this_slice)] = (
                                axis_error_for_thisdim[j]
                            )
                logger.debug(
                    strm(
                        "shape of multidim_axis_label is now",
                        multidim_axis_label.shape,
                        "(",
                        axes_with_labels,
                        ")",
                    )
                )
                logger.debug(
                    strm(
                        "multidim_axis_label is:\n", repr(multidim_axis_label)
                    )
                )
                multidim_axis_label = (
                    multidim_axis_label.flatten()
                )  # then flatten the axis
                logger.debug(
                    strm(
                        "shape of multidim_axis_label is now",
                        multidim_axis_label.shape,
                    )
                )
                logger.debug(
                    strm(
                        "multidim_axis_label is:\n", repr(multidim_axis_label)
                    )
                )
            else:
                raise ValueError(
                    "You requested that smoosh generate an axis, but I don't"
                    " know what dtype to assign to it (what fields to use). "
                    " This is likely because you don't have axes assigned to"
                    " the dimensions you're trying to smoosh.  Consider"
                    " calling smoosh with noaxis=True, instead"
                )
            axis_coords_dict[dimname] = multidim_axis_label
            axis_coords_error_dict[dimname] = multidim_axis_error
            # }}}
        # {{{ update axis dictionary with the new info
        logger.debug(
            strm(
                "end up with axis_coords_dict (%d)" % len(axis_coords_dict),
                axis_coords_dict,
            )
        )
        logger.debug(
            strm(
                "end up with axis_coords_error_dict (%d)"
                % len(axis_coords_error_dict),
                axis_coords_error_dict,
            )
        )
        # }}}
        # }}}
        # {{{ make new dimlabels, and where relevant, project the new
        #     dictionary onto these dimlabels
        axis_coords_units_dict[dimname] = new_units
        self.dimlabels = retained_dims + [dimname]
        logger.debug(
            strm(
                "end up with dimlabels",
                self.dimlabels,
                "and shape",
                self.data.shape,
            )
        )
        self.axis_coords = self.fld(axis_coords_dict)
        self.axis_coords_error = self.fld(axis_coords_error_dict)
        self.axis_coords_units = self.fld(axis_coords_units_dict)
        logger.debug(
            strm(
                "new axis coords (%d)" % len(self.axis_coords),
                self.axis_coords,
            )
        )
        logger.debug(
            strm(
                "new axis coords errors (%d)" % len(self.axis_coords_error),
                self.axis_coords_error,
            )
        )
        logger.debug(
            strm(
                "new axis coords unitss (%d)" % len(self.axis_coords_units),
                self.axis_coords_units,
            )
        )
        # }}}
        # {{{ then deal with the units
        # }}}
        # {{{ finally, if I need to, reorder again to put the new dimension
        #     where I want it
        # }}}
        return self

    def chunk(self, axisin, *otherargs):
        r""" "Chunking" is defined here to be the opposite of taking a direct
        product, increasing the number of dimensions by the inverse of the
        process by which taking a direct product decreases the number of
        dimensions.  This function chunks axisin into multiple new axes
        arguments.:
            axesout -- gives the names of the output axes
            shapesout -- optional -- if not given, it assumes equal length --
            if given, one of the values can be -1, which is assumed length

        When there are axes, it assumes that the axes of the new dimensions
        are nested -- *e.g.*, it will chunk a dimension with axis:
        [1,2,3,4,5,6,7,8,9,10]
        into dimensions with axes:
        [0,1,2,3,4], [1,6]

        ..todo::
            Deal with this efficiently when we move to new-style axes
        """
        if len(otherargs) == 2:
            axesout, shapesout = otherargs
        elif len(otherargs) == 1:
            if isinstance(otherargs[0], list):
                axesout = otherargs[0]
                shapesout = ndshape(self)[axisin] ** (1.0 / len(axesout))
                if (
                    abs(shapesout - np.round(shapesout)) > 1e-15
                ):  # there is some kind of roundoff error here
                    raise ValueError(
                        "In order for chunk to be called with only a list of"
                        " axes, the shape of the dimension you are trying to"
                        " split (here %s) must be an Nth root of the original"
                        " dimension size (here: %d), where N (here %d) is the"
                        " number of dimensions you are trying to chunk into"
                        % (axisin, ndshape(self)[axisin], len(axesout))
                    )
                else:
                    shapesout = np.round(shapesout)
                shapesout = [shapesout] * len(axesout)
            elif isinstance(otherargs[0], dict):
                axesout, shapesout = list(otherargs[0].keys()), list(
                    otherargs[0].values()
                )
            else:
                raise ValueError("I don't know how to deal with this type!")
        else:
            raise ValueError("otherargs must be one or two arguments!")
        assert not np.any(
            [j in self.dimlabels for j in axesout if j != axisin]
        ), strm(
            "You are trying to create dimensions",
            [j for j in axesout if j != axisin],
            "one of which matches one of the existing labels",
            self.dimlabels,
        )
        if np.any([j == -1 for j in shapesout]):
            j = shapesout.index(-1)
            if j < len(shapesout) - 1:
                shapesout[j] = int(
                    np.round(
                        ndshape(self)[axisin]
                        / np.prod(r_[shapesout[0:j], shapesout[j + 1 :]])
                    )
                )
            else:
                shapesout[j] = int(
                    np.round(ndshape(self)[axisin] / np.prod(shapesout[0:j]))
                )
        if np.prod(shapesout) != ndshape(self)[axisin]:
            raise ValueError(
                "The size of the axis (%s) you're trying to split (%s) doesn't"
                " match the size of the axes you're trying to split it into"
                " (%s = %s)"
                % (
                    repr(axisin),
                    repr(ndshape(self)[axisin]),
                    repr(axesout),
                    repr(shapesout),
                )
            )
        thisaxis = self.axn(axisin)
        if self.getaxis(axisin) is not None and len(self.getaxis(axisin)) > 0:
            axes_tmp = self.getaxis(axisin).reshape(shapesout)
            new_axes = []
            for j in range(len(axes_tmp.shape)):
                this_slicer = [0] * len(axes_tmp.shape)
                this_slicer[j] = slice(None, None, None)
                new_axes.append(axes_tmp[tuple(this_slicer)])
        else:
            new_axes = None
        # {{{ if there is a list of axis coordinates, add in slots for the new
        #     axes
        if isinstance(self.axis_coords, list):
            if len(self.axis_coords) == 0:
                self.axis_coords = [None] * len(self.dimlabels)
            for j in range(len(axesout) - 1):
                self.axis_coords.insert(thisaxis, None)
        if isinstance(self.axis_coords_error, list):
            if len(self.axis_coords_error) == 0:
                self.axis_coords_error = [None] * len(self.dimlabels)
            for j in range(len(axesout) - 1):
                self.axis_coords_error.insert(thisaxis, None)
        if isinstance(self.axis_coords_units, list):
            if len(self.axis_coords_units) == 0:
                self.axis_coords_units = [None] * len(self.dimlabels)
            for j in range(len(axesout) - 1):
                self.axis_coords_units.insert(thisaxis, None)
        # }}}
        newshape = (
            list(self.data.shape[0:thisaxis])
            + shapesout
            + list(self.data.shape[thisaxis + 1 :])
        )
        newshape = list(map(int, newshape))
        newnames = (
            list(self.dimlabels[0:thisaxis])
            + axesout
            + list(self.dimlabels[thisaxis + 1 :])
        )
        self.data = self.data.reshape(newshape)
        orig_axis_units = self.get_units(axisin)
        self.dimlabels = newnames
        if new_axes is not None:
            for j in range(len(axesout)):
                self.setaxis(axesout[j], new_axes[j])
                self.set_units(axesout[j], orig_axis_units)
        return self

    def chunk_auto(self, axis_name, which_field=None, dimname=None):
        r"""assuming that axis "axis_name" is currently labeled with a
        structured np.array, choose one field ("which_field") of that
        structured np.array to generate a new dimension
        Note that for now, by definition, no error is allowed on the axes.
        However, once I upgrade to using structured arrays to handle axis and
        data errors, I will want to deal with that appropriately here.
        """

        def check_data(a):
            """we need this because other things expect dimlabels to be a list
            of strings"""
            if isinstance(a.dimlabels, np.recarray):
                a.dimlabels = [
                    str(j[0]) if len(j) == 1 else j
                    for j in a.dimlabels.tolist()
                ]
            return a

        if which_field is None:
            which_field = self[axis_name].dtype.names[0]
        axis_number = self.axn(axis_name)
        new_axis, indices = np.unique(
            self.getaxis(axis_name)[which_field], return_inverse=True
        )  # we are essentially creating a hash table for the axis.  According
        #    to numpy documentation, the hash indices that this returns should
        #    also be sorted sorted.
        logger.debug(strm("(chunk auto) indices look like this:", indices))
        # {{{ check that there are equal numbers of all the unique new_axis
        index_count = np.array(
            [np.count_nonzero(indices == j) for j in range(indices.max() + 1)]
        )
        if all(index_count == index_count[0]):
            logger.debug(
                strm(
                    "(chunk auto) Yes, there are equal numbers of all unique"
                    " new_axis! (Each element of the hash table has been"
                    " indexed the same number of times.)"
                )
            )
            # }}}
            # {{{ store the old shape and generate the new shape
            current_shape = list(self.data.shape)
            logger.debug(strm("(chunk auto) old shape -- ", current_shape))
            new_shape = np.insert(
                current_shape, axis_number + 1, len(new_axis)
            )
            new_shape[axis_number] /= len(
                new_axis
            )  # the indices of the hash table become the new dimension
            # }}}
            # {{{ actually reorder the data and error -- perhaps a view would
            #     be more efficient here
            old_data = self.data
            has_data_error = not (self.get_error() is None)
            self.data = np.empty(new_shape, dtype=self.data.dtype)
            if has_data_error:
                old_error = self.get_error()
                self.set_error(np.empty(new_shape, dtype=self.data.dtype))
            # }}}
            # {{{ adjust all the relevant axis information
            # {{{ generate an axis label along the axis I'm chunking that's
            #     stripped of the field that I'm creating a dimension from
            #     (i.e.  chunking off a new dimension based on) -- because I am
            #     independently manipulating the data, I don't use
            #     self.getaxis()
            x_strip_current_field = self.axis_coords[axis_number][[
                j
                for j in self.axis_coords[axis_number].dtype.names
                if j != which_field
            ]]
            # }}}
            # {{{ reshape the axis coordinate so that it becomes a 2D np.array
            #     with the new dimension chunked off
            self.axis_coords[axis_number] = np.empty(
                (len(x_strip_current_field) // len(new_axis), len(new_axis)),
                dtype=x_strip_current_field.dtype,
            )
            if not (self.get_error(axis_name) is None):
                raise ValueError(
                    "Until I do the structured np.array upgrade chunk_auto"
                    " will not be able to deal with an axis that has errors."
                )
            # }}}
            # {{{ everything is now ready to sort the data and residual axis
            #     into ordered slots
            # {{{ initialize the slices
            copy_to_slice = len(new_shape) * [
                slice(None, None, None)
            ]  # this is the memory address inside the new data (where stuff
            #    goes)
            copy_from_slice = len(current_shape) * [
                slice(None, None, None)
            ]  # this is the memory address inside the old data (where stuff
            #    comes from)
            # }}}
            if has_data_error:
                data_error_location = self.get_error()
            for j in range(len(new_axis)):  # j is the index in the hash table
                copy_to_slice[axis_number + 1] = j
                copy_from_slice[axis_number] = np.where(indices == j)[0]
                self.data[tuple(copy_to_slice)] = old_data[
                    tuple(copy_from_slice)
                ]
                if has_data_error:
                    data_error_location[tuple(copy_to_slice)] = old_error[
                        tuple(copy_from_slice)
                    ]
                logger.debug(
                    strm(
                        "(chunk auto) ",
                        j,
                        "matches at",
                        x_strip_current_field[copy_from_slice[axis_number]],
                    )
                )
                self.axis_coords[axis_number][:, j] = x_strip_current_field[
                    copy_from_slice[axis_number]
                ]
            # }}}
            logger.debug(
                strm(
                    "(chunk auto) new axis -- ", self.axis_coords[axis_number]
                )
            )
            logger.debug(strm("(chunk auto) new shape -- ", self.data.shape))
            # {{{ housekeeping for the various axes + data properties -- should
            #     perhaps be possible to do this first, then use .getaxis()
            self.dimlabels.insert(axis_number + 1, which_field)
            self.axis_coords.insert(axis_number + 1, new_axis)
            # {{{ by definition, axis can have neither errors nor units
            #     associated, for now.
            self.axis_coords_error.insert(axis_number + 1, None)
            self.axis_coords_units.insert(axis_number + 1, None)
            # }}}
            logger.debug(
                strm(
                    "(chunk auto) the dimensions of ",
                    self.dimlabels[axis_number],
                    "are (?? x ",
                    self.dimlabels[axis_number + 1],
                    ")=",
                    self.axis_coords[axis_number].shape,
                )
            )
            # }}}
            # }}}
            # {{{ deal appropriately with the "remainder axis" (axis_number)
            if dimname is None:
                remainder_axis_name = "_and_".join(
                    self.axis_coords[axis_number].dtype.names
                )
            else:
                remainder_axis_name = dimname
            # {{{ if everything is the same along the dimension that I've just
            # created (which is the second dimension), then get rid of the
            # duplicate labels
            test_axis = self.axis_coords[axis_number].T
            logger.debug(strm("(chunk auto) test axis -- ", test_axis))
            test_axis = (
                np.ascontiguousarray(test_axis)
                .flatten()
                .view([("", test_axis.dtype)] * test_axis.shape[1])
            )
            if all(test_axis == test_axis[0]):
                self.axis_coords[axis_number] = self.axis_coords[axis_number][
                    :, 0
                ].reshape(1, -1)
                logger.debug(
                    strm(
                        "(chunk auto) collapsed to",
                        self.axis_coords[axis_number],
                    )
                )
            # }}}
            if (
                self.axis_coords[axis_number].shape[0] == 1
            ):  # then this is a "valid" axis -- because, for each position of
                # the new axis, there is only one value of the remainder axis
                self.axis_coords[axis_number] = self.axis_coords[
                    axis_number
                ].reshape(-1)
                self.dimlabels[axis_number] = remainder_axis_name
                if (
                    len(self.axis_coords[axis_number].dtype) == 1
                ):  # only one field, which by the previous line will be named
                    # appropriately, so drop the structured np.array name
                    new_dtype = self.axis_coords[axis_number].dtype.descr[0][1]
                    self.axis_coords[axis_number] = np.array(
                        self.axis_coords[axis_number], dtype=new_dtype
                    )  # probably more efficiently done with a view, but leave
                    #    alone for now
                return check_data(self)
            else:
                # {{{ generate an index list to label the remainder axis, and
                #     generate a new nddata with the actual values (which are
                #     not copied across the new dimension that matches
                #     which_field)
                remainder_axis_index_list = r_[
                    0 : self.axis_coords[axis_number].shape[0]
                ]
                new_data = nddata(
                    self.axis_coords[axis_number],
                    self.axis_coords[axis_number].shape,
                    [remainder_axis_name, which_field],
                )
                self.axis_coords[axis_number] = remainder_axis_index_list
                new_data.labels(
                    [remainder_axis_name, which_field],
                    [
                        self.axis_coords[axis_number].copy(),
                        self.axis_coords[axis_number + 1].copy(),
                    ],
                )
                self.dimlabels[axis_number] = remainder_axis_name + "_INDEX"
                # }}}
                return check_data(self), check_data(new_data)
            # }}}
        else:
            raise ValueError(
                "Along the axis '"
                + axis_name
                + "', the field '"
                + which_field
                + "' does not represent an axis that is repeated one or more"
                " times!  The counts for how many times each element along"
                " the field is used is "
                + repr(index_count)
            )
            return

    def squeeze(self, return_dropped=False):
        r"""squeeze singleton dimensions

        Parameters
        ==========
        return_dropped: bool (default False)
           return a list of the dimensions that were dropped as a second
           argument
        Returns
        =======
        self

        return_dropped: list
            (optional, only if return_dropped is True)
        """
        mask = np.array(self.data.shape) > 1
        logger.debug(strm(list(zip(mask, self.dimlabels))))
        self.data = self.data.squeeze()
        retval = []
        if isinstance(self.axis_coords, list):
            for k, v in [
                (self.dimlabels[j], self.axis_coords[j])
                for j in range(len(self.dimlabels))
                if not mask[j]
            ]:
                retval.append(k)
                if v is not None:
                    self.set_prop(k, v[0])
        self.dimlabels = [v for j, v in enumerate(self.dimlabels) if mask[j]]
        if isinstance(self.axis_coords, list):
            self.axis_coords = [
                v for j, v in enumerate(self.axis_coords) if mask[j]
            ]
        if isinstance(self.axis_coords_error, list):
            self.axis_coords_error = [
                v for j, v in enumerate(self.axis_coords_error) if mask[j]
            ]
        if isinstance(self.axis_coords_units, list):
            self.axis_coords_units = [
                v for j, v in enumerate(self.axis_coords_units) if mask[j]
            ]
        if return_dropped:
            return self, retval
        else:
            return self

    # }}}
    # {{{ messing with data -- get, set, and copy
    def __getslice__(self, *args):
        raise ValueError(strm("getslice! ", args))

    def __setitem__(self, key, val):
        righterrors = None
        logger.debug(strm("key", key))
        if isinstance(key, nddata):
            logger.debug("initially, rightdata appears to be nddata")
            _, B = self.aligndata(key)
            key = B.data  # now the next part will handle this
        elif isinstance(key, np.ndarray):  # if selector is an np.ndarray
            logger.debug("initially, rightdata appears to be np.ndarray")
            if key.dtype is not np.dtype("bool"):
                raise ValueError(
                    "I don't know what to do with an np.ndarray subscript that"
                    " has dtype "
                    + repr(key.dtype)
                )
            if key.shape != self.data.shape:
                raise ValueError(
                    "The shape of your logical mask "
                    + repr(key.shape)
                    + " and the shape of your data "
                    + repr(self.data.shape)
                    + " are not compatible (matching or singleton) -- I really"
                    " don't think that you want to do this!"
                )
            self.data[key] = val
            return
        elif isinstance(key, str):
            logger.debug("setting the axis")
            self.setaxis(key, val)
            return self
        if isinstance(val, nddata):
            logger.debug(
                "rightdata appears to be nddata after initial treatment"
            )
            # {{{ reorder so the shapes match
            unshared_indices = list(set(val.dimlabels) ^ set(self.dimlabels))
            shared_indices = list(self.dimlabels)
            if "INDEX" in unshared_indices:
                unshared_indices.remove("INDEX")
            shared_indices = [
                j for j in shared_indices if j not in unshared_indices
            ]
            if len(val.dimlabels) != len(shared_indices) or (
                not all([
                    val.dimlabels[j] == shared_indices[j]
                    for j in range(0, len(shared_indices))
                ])
            ):
                val.reorder(shared_indices)
            # }}}
            rightdata = val.data
            righterrors = val.get_error()
        else:  # assume it's an np.ndarray
            logger.debug(
                "rightdata appears to be np.ndarray after initial treatment"
            )
            rightdata = val
            # {{{ if I just passed a function, assume that I'm applying some
            #     type of data-based mask
            if isinstance(key, type(emptyfunction)):
                thisfunc = key
                self.data[thisfunc(self.data)] = rightdata
                return
            # }}}
            if not isinstance(rightdata, np.ndarray):  # in case its a scalar
                rightdata = np.array([rightdata])
        slicedict, axesdict, errordict, unitsdict = self._parse_slices(
            key
        )  # pull left index list from parse slices
        leftindex = tuple(self.fld(slicedict))
        rightdata = rightdata.squeeze()
        logger.debug(
            strm("after squeeze, rightdata has shape", rightdata.shape)
        )
        if len(rightdata.shape) > 0:
            left_shape = np.shape(self.data[leftindex])
            try:
                self.data[leftindex] = rightdata.reshape(
                    left_shape
                )  # assign the data
            except Exception:
                raise IndexError(
                    strm(
                        "ERROR ASSIGNING NDDATA:\n",
                        "self.data.shape:",
                        self.data.shape,
                        "left index",
                        leftindex,
                        "\n",
                        "rightdata.shape:",
                        rightdata.shape,
                        "--> shape of left slice: ",
                        left_shape,
                    )
                )
        else:
            self.data[leftindex] = rightdata
        lefterror = self.get_error()
        if lefterror is not None:
            lefterror[leftindex] = righterrors.squeeze()
        return self

    # {{{ standard trig functions
    def __getattribute__(self, arg):
        fundict = {
            "exp": np.exp,
            "sin": np.sin,
            "cos": np.cos,
            "tan": np.tan,
            "sinh": np.sinh,
            "cosh": np.cosh,
            "tanh": np.tanh,
            "log": np.log,
            "log10": np.log10,
        }
        if arg in list(fundict.keys()):
            argf = fundict[arg]

            def retfun():
                retval = self.copy()
                retval.data = argf(retval.data)
                return retval

            return retfun
        elif arg == "shape":
            return ndshape(self)
        elif arg == "isfortran":
            raise ValueError(
                "you tried to call isfortran on an nddata object -- this"
                " probably means you're doing something wrong -- possibly that"
                " you are passing an nddata object when you should be passing"
                " a standard numpy ndarray"
            )
        else:
            return super().__getattribute__(arg)

    @property
    def C(self):
        """shortcut for copy

        btw, what we are doing is analogous to a ruby function with
        functioname!() modify result, and we can use the "out" keyword in
        numpy.

        ..todo::
            (new idea)
            This should just set a flag that says "Do not allow this data to be
            substituted in place,"
            so that if something goes to edit the data in place,
            it instead first makes a copy.

            also here, see
            `Definition of shallow and deep copy <https://docs.python.org/2/li\
                    brary/copy.html>`_

            (older idea)
            We should offer "N", which generates something like a copy,
            but which is sets the equivalent of "nopop".
            For example, currently, you need to do something like
            ``d.C.argmax('t2')``,
            which is very inefficient, since it copies the whole np.array.
            So, instead, we should do
            ``d.N.argmax('t2')``, which tells argmax and all other
            functions not to overwrite "self" but to return a new object.
            This would cause things like "run_nopop" to become obsolete.
        """
        return self.copy()

    @C.setter
    def C(self):
        raise ValueError(
            "You can't set the C property -- it's used to generate a copy"
        )

    @property
    def angle(self):
        """Return the angle component of the data.

        This has error, which is calculated even if there is no error in
        the original data -- in the latter case, a uniform error of 1 is
        assumed. (This is desirable since phase is a tricky beast!)
        """
        retval = self.copy(data=False)
        retval.data = np.angle(self.data)
        if np.isscalar(self.data):
            # if scalar, no error
            return retval
        else:
            # from ∂φ/∂A=-i n/2A
            # when ρe(iφ)=xAⁿ
            dangle_dA = np.empty_like(self.data)
            mask = self.data != 0
            dangle_dA[mask] = 1 / (2 * self.data[mask])
            dangle_dA[~mask] = np.nan
            A_sigma = self.get_error()
            A_sigma = 1 if A_sigma is None else A_sigma
            retval.set_error(abs(dangle_dA * A_sigma))
            return retval

    @angle.setter
    def angle(self):
        raise ValueError("Can't independently set the angle component yet")

    @property
    def imag(self):
        "Return the imag component of the data"
        retval = self.copy(data=False)
        # data=False excludes the error
        retval.data_error = self.data_error
        retval.data = self.data.imag
        return retval

    @imag.setter
    def imag(self):
        raise ValueError("Can't independently set the imag component yet")

    @property
    def real(self):
        "Return the real component of the data"
        retval = self.copy(data=False)
        # data=False excludes the error
        retval.data_error = self.data_error
        retval.data = self.data.real
        return retval

    @real.setter
    def real(self):
        raise ValueError("Can't independently set the real component yet")

    # }}}
    def copy(self, data=True):
        r"""Return a full copy of this instance.

        Because methods typically change the data in place, you might want to
        use this frequently.

        Parameters
        ----------
        data : boolean
            Default to True.
            False doesn't copy the data -- this is for internal use,
            *e.g.* when you want to copy all the metadata and perform a
            calculation on the data.

            The code for this also provides the definitive list of the
            nddata metadata.
        """
        if data:
            retval = deepcopy(self)
            retval.other_info = deepcopy(self.other_info)
            return retval
        else:
            retval = nddata(0)  # np.empty
            # {{{ data info
            retval.dimlabels = list(self.dimlabels)
            retval.data = None
            retval.data_error = None
            if hasattr(self, "data_units"):
                retval.data_units = deepcopy(self.data_units)
            if hasattr(self, "data_covariance"):
                retval.data_covariance = deepcopy(self.data_covariance)
            # }}}
            # {{{ axes
            retval.axis_coords = deepcopy(self.axis_coords)
            retval.axis_coords_error = deepcopy(self.axis_coords_error)
            retval.axis_coords_units = deepcopy(self.axis_coords_units)
            # }}}
            retval.other_info = deepcopy(self.other_info)
            return retval

    def set_to(self, otherinst):
        r"""Set data inside the current instance to that of the other instance.

        Goes through the list of attributes specified in copy,
        and assigns them to the element of the current instance.

        This is to be used:

        *   for constructing classes that inherit nddata with additional
            methods.
        *   for overwriting the current data with the result of a slicing
            operation
        """
        self.data = otherinst.data
        self.dimlabels = otherinst.dimlabels
        self.data_error = otherinst.data_error
        if hasattr(otherinst, "data_units"):
            self.data_units = otherinst.data_units
        if hasattr(otherinst, "data_covariance"):
            self.data_covariance = otherinst.data_covariance
        self.axis_coords = otherinst.axis_coords
        self.axis_coords_error = otherinst.axis_coords_error
        self.axis_coords_units = otherinst.axis_coords_units
        self.other_info = otherinst.other_info
        return self

    def like(self, value):
        r"""provide "zeros_like" and "ones_like" functionality

        Parameters
        ==========
        value: float
            1 is "ones_like" 0 is "zeros_like", etc.
        """
        retval = self.copy(data=False)
        retval.data = np.empty_like(self.data)
        retval.data[:] = value
        return retval

    def __getitem__(self, args):
        if isinstance(args, type(emptyfunction)):
            # {{{ just a lambda function operates on the data
            thisfunc = args
            newdata = self.copy()
            mask = thisfunc(newdata.data)
            newdata.data = newdata.data[mask]
            if len(newdata.dimlabels) == 1:
                x = newdata.getaxis(newdata.dimlabels[0])
                newdata.setaxis(newdata.dimlabels[0], x[mask])
            else:
                raise ValueError(
                    "I don't know how to do this for multidimensional data"
                    " yet!"
                )
            return newdata
            # }}}
        elif isinstance(args, nddata):
            # {{{ try the nddata mask
            A = args
            if isinstance(A, nddata) and A.data.dtype is np.dtype("bool"):
                thisshape = ndshape(A)
                nonsingleton = []
                for thisdim in A.dimlabels:
                    if thisshape[thisdim] != 1:
                        nonsingleton.append(thisdim)
                if len(nonsingleton) != 1:
                    raise ValueError(
                        "To index with an nddata, you must have only one"
                        " dimension"
                    )
                else:
                    self.setaxis(
                        nonsingleton[0],
                        self.getaxis(nonsingleton[0])[A.data.flatten()],
                    )
                _, B = self.aligndata(A)
                A = B.data  # now the next part will handle this
                if A.dtype is not np.dtype("bool"):
                    raise ValueError(
                        "I don't know what to do with an np.ndarray subscript"
                        " that has dtype "
                        + repr(A.dtype)
                    )
                if A.shape != self.data.shape:
                    temp = np.array(A.shape) == 1
                    if all(
                        np.array(A.shape)[temp]
                        == np.array(self.data.shape)[temp]
                    ):
                        pass
                    else:
                        raise ValueError(
                            "The shape of your logical mask "
                            + repr(A.shape)
                            + " and the shape of your data "
                            + repr(self.data.shape)
                            + " are not compatible (matching or singleton) --"
                            " I really don't think that you want to do this!"
                        )
                self.data = self.data[A]
                return self
            else:
                errmsg = "you passed a single argument of type " + repr(
                    type(A)
                )
                if isinstance(A, nddata):
                    errmsg += " with dtype " + repr(A.data.dtype)
                errmsg += " -- I don't know what to do with this"
                raise ValueError(errmsg)
            # }}}
        elif isinstance(args, str):
            return self.getaxis(args)
        else:
            if type(args) is not slice:
                if type(args) not in [tuple, list]:
                    raise ValueError(
                        "the first argument to your nddata slice/index is not"
                        " a string -- I don't understand that!  Are you trying"
                        " to pass nddata to a function that only accepts numpy"
                        " ndarrays?"
                    )
                elif type(args[0]) is not str:
                    raise ValueError(
                        "the first argument to your nddata slice/index is not"
                        " a string -- I don't understand that!  Are you trying"
                        " to pass nddata to a function that only accepts numpy"
                        " ndarrays?"
                    )
            slicedict, axesdict, errordict, unitsdict = self._parse_slices(
                args
            )
            if (
                not isinstance(args, slice)
                and isinstance(args[1], list)
                and isinstance(args[0], str)
                and len(args) == 2
            ):
                return concat([self[args[0], x] for x in args[1]], args[0])
            indexlist = tuple(self.fld(slicedict))
            newlabels = [
                x for x in self.dimlabels if not np.isscalar(slicedict[x])
            ]  # generate the new list of labels, in order, for all dimensions
            #    that are not indexed by a scalar
        # {{{ properly index the data error
        if self.data_error is not None:
            try:
                newerror = self.data_error[indexlist]
            except Exception:
                raise ValueError(
                    "Problem trying to index data_error"
                    + repr(self.data_error)
                    + " with",
                    repr(indexlist),
                )
        else:
            newerror = None
        # }}}
        if len(self.axis_coords) > 0:
            if errordict is not None:
                axis_coords_error = [errordict[x] for x in newlabels]
            else:
                axis_coords_error = None
            if unitsdict is not None:
                axis_coords_units = [unitsdict[x] for x in newlabels]
            else:
                axis_coords_units = None
            try:
                sliced_data = self.data[indexlist]
            except Exception:
                raise ValueError(
                    strm(
                        "the slice values that you've passed",
                        "don't seem to match the size of the data",
                        "the shape of the data is",
                        self.data.shape,
                        "and the index list (the slice indeces passed to the",
                        "underlying numpy data) I generate from this"
                        " command is",
                        indexlist,
                        "likely, one of the slice indeces is out of bounds for"
                        " the size of the data",
                    )
                )
            try:
                retval = nddata(
                    sliced_data,
                    sliced_data.shape,
                    newlabels,
                    axis_coords=[axesdict[x] for x in newlabels],
                    axis_coords_error=axis_coords_error,
                    data_error=newerror,
                    other_info=self.other_info,
                )
            except Exception:
                raise ValueError(
                    strm(
                        "likely some problem recasting the data when"
                        "trying to initialize a new nddata: shape of"
                        "self.data",
                        self.data.shape,
                        "indexlist",
                        indexlist,
                    )
                )
            retval.axis_coords_units = axis_coords_units
            retval.data_units = self.data_units
            return retval
        else:
            retval = nddata(
                self.data[indexlist],
                self.data[indexlist].shape,
                newlabels,
                other_info=self.other_info,
            )
            retval.axis_coords_units = self.axis_coords_units
            retval.data_units = self.data_units
            return retval

    def _possibly_one_axis(self, *args):
        if len(args) == 1:
            return args[0]
        if len(args) > 1:
            raise ValueError("you can't pass more than one argument!!")
        if len(args) == 0:
            if len(self.dimlabels) == 1:
                axes = self.dimlabels
            elif len(self.dimlabels) == 0:
                raise ValueError(
                    "You're trying to do something to data with no dimensions"
                )
            else:
                raise ValueError(
                    "If you have more than one dimension, you need to tell me"
                    " which one!!"
                )
        return axes

    def get_range(self, dimname, start, stop):
        """get raw indices that can be used to generate a slice for the start
        and (non-inclusive) stop

        Uses the same code as the standard slicing format (the 'range' option
        of parseslices)

        Parameters
        ==========
        dimname: str
            name of the dimension
        start: float
            the coordinate for the start of the range
        stop: float
            the coordinate for the stop of the range

        Return
        ======
        start: int
            the index corresponding to the start of the range
        stop: int
            the index corresponding to the stop of the range
        """
        axesdict = self.mkd(self.axis_coords)
        if len(axesdict) == 0:
            raise ValueError(
                f"possible that no axes are labeled? {ndshape(self)}"
            )
        if axesdict[dimname] is None:
            raise ValueError(
                "You passed a range-type slice"
                + " selection, but to do that, your axis coordinates need to"
                + f" be labeled! (The axis coordinates of {dimname} aren't"
                + " labeled)"
            )
        temp = np.diff(axesdict[dimname])
        if not all(temp * np.sign(temp[0]) > 0):
            raise ValueError(
                strm(
                    "you can only use the range format on data where the axis"
                    " is in consecutively increasing or decreasing order, and"
                    " the differences that I see are",
                    temp * np.sign(temp[0]),
                ),
                "if you like, you can still do this by first calling .sort( on"
                " the %s axis" % dimname,
            )
        if np.sign(temp[0]) == -1:
            thisaxis = axesdict[dimname][::-1]
        else:
            thisaxis = axesdict[dimname]
        if start is None:
            start = -inf
        if stop is None:
            stop = inf
        if start > stop:
            start, stop = stop, start
        # at this point, start is indeed the lower value, and stop indeed the
        # higher
        if start == inf:
            raise ValueError(
                strm(
                    "this is not going to work -- I interpret range",
                    start,
                    stop,
                    "I get to",
                    start,
                    ",",
                    stop,
                )
            )
        elif start == -inf:
            start = 0
        else:
            logger.debug(strm("looking for", start))
            start = np.searchsorted(thisaxis, start)
            if start >= len(thisaxis):
                raise ValueError(
                    'the lower value of your slice %s on the "%s" axis (which'
                    " runs from %g to %g) is higher than the highest value of"
                    " the axis coordinates!"
                    % (
                        (
                            str((start, stop)),
                            dimname,
                        )
                        + tuple(self.getaxis(dimname)[r_[0, -1]])
                    )
                )
        stop_float = stop
        if stop == inf:
            stop = len(
                thisaxis
            )  # not an exact match (inf doesn't match the index), so needs to
            #    be inclusive already
        elif stop == -inf:
            raise ValueError(
                strm(
                    "this is not going to work -- I I get to",
                    start,
                    ",",
                    stop,
                )
            )
        else:
            logger.debug(strm("looking for", stop))
            stop = np.searchsorted(thisaxis, stop)
        # at this point, the result is inclusive if stop is
        # not an exact match, but exclusive if it is
        if stop < len(thisaxis) and thisaxis[stop] == stop_float:
            stop += 1  # make it inclusive
        if np.sign(temp[0]) == -1:
            stop = len(thisaxis) - 1 - stop
            start = len(thisaxis) - 1 - start
            stop, start = start, stop
        del temp
        if start == stop:
            stop += 1
        return start, stop

    def _parse_slices(self, args):
        """This controls nddata slicing:
        it previously took
        \'axisname\',value
        pairs where value was an index or a lambda function.
        Now, it also takes
        \'axisname\':value
        and
        \'axisname\':(value1,value2)
        pairs, where the values give either a single value or an inclusive
        range on the axis, respectively
        """
        logger.debug(
            strm(
                "about to start parsing slices",
                args,
                "for data with axis_coords of length",
                len(self.axis_coords),
                "and dimlabels",
                self.dimlabels,
                "for ndshape of",
                ndshape(self),
            )
        )
        errordict = None  # in case it's not set
        if self.axis_coords_units is not None:
            unitsdict = self.mkd(self.axis_coords_units)
        axesdict = None  # in case it's not set
        if isinstance(args, slice):
            args = [args]
        # {{{ make a sensible list of tuples that's easier to understand
        sensible_list = []  # type, dimension, arguments
        testf = lambda x: x + 1
        j = 0
        while j < len(args):
            if isinstance(args[j], str):  # works for str and np.str_
                dimname = args[j]
                if isinstance(dimname, np.str_):
                    dimname = str(
                        dimname
                    )  # on upgrading + using on windows, this became
                    #    necessary, for some reason I don't understand
                if isinstance(args[j + 1], type(testf)):
                    sensible_list.append((hash("func"), dimname, args[j + 1]))
                else:
                    sensible_list.append((hash("np"), dimname, args[j + 1]))
                j += 2
            elif type(args[j]) is slice:
                dimname = args[j].start
                if isinstance(dimname, np.str_):
                    dimname = str(dimname)
                target = args[j].stop
                if np.isscalar(target):
                    sensible_list.append((hash("idx"), dimname, target))
                elif type(target) in [tuple, list]:
                    assert len(target) in [1, 2], strm(
                        "for",
                        args[j],
                        "I expected a 'dimname':(range_start,range_stop)",
                    )
                    if len(target) == 1:
                        sensible_list.append(
                            (hash("range"), dimname, target[0], None)
                        )
                    else:
                        sensible_list.append(
                            (hash("range"), dimname, target[0], target[1])
                        )
                elif isinstance(target, np.ndarray) and target.size == 2:
                    sensible_list.append(
                        (hash("range"), dimname, target[0], target[1])
                    )
                else:
                    raise ValueError(
                        "for part of your slice, you said for dimension",
                        dimname,
                        "you wanted",
                        target,
                        "but the second argument must be a tuple, list, or"
                        " array of length 2!",
                    )
                j += 1
            else:  # works for str and np.str_
                raise ValueError(
                    "I have read in slice argument",
                    args[:j],
                    "but then I get confused!",
                )

        def pprint(a):
            b = {hash(j): j for j in ["idx", "range", "np", "func"]}
            return (b[a[0]],) + a[1:]

        logger.debug(
            strm("Here is the sensible list:", [str(j) for j in sensible_list])
        )
        # }}}
        if type(args) in [float, np.int32, int, np.double]:
            raise ValueError(
                strm("You tried to pass just a nddata[", type(args), "]")
            )
        if isinstance(args[0], str) or isinstance(args[0], slice):
            # {{{ create a slicedict and errordict to store the slices
            slicedict = dict(
                list(
                    zip(
                        list(self.dimlabels),
                        [slice(None, None, None)] * len(self.dimlabels),
                    )
                )
            )  # initialize to all none
            if len(self.axis_coords) > 0:
                logger.debug(
                    strm(
                        "trying to make dictionaries from axis coords of len",
                        len(self.axis_coords),
                        "and axis_coords_error of len",
                        len(self.axis_coords_error),
                        "when dimlabels has len",
                        len(self.dimlabels),
                    )
                )
                axesdict = self.mkd(self.axis_coords)
                if len(self.axis_coords_error) > 0:
                    errordict = self.mkd(self.axis_coords_error)
            else:
                axesdict = self.mkd(self.axis_coords)
                logger.debug(
                    strm(
                        "length of axis_coords not greater than 0, generated"
                        " dictionary",
                        axesdict,
                    )
                )
            # }}}
            # {{{ map the slices onto the axis coordinates and errors
            for thistuple in sensible_list:
                thisop = thistuple[0]
                thisdim = thistuple[1]
                thisargs = thistuple[2:]
                # print "DEBUG, type of slice",x,"is",type(y)
                if thisop == hash("np"):
                    slicedict[thisdim] = thisargs[0]
                    if np.isscalar(thisargs[0]):
                        axesdict.pop(
                            thisdim
                        )  # pop the axes for all scalar dimensions
                    else:
                        if axesdict[thisdim] is not None:
                            axesdict[thisdim] = axesdict[thisdim][
                                slicedict[thisdim]
                            ]
                elif thisop == hash("func"):
                    mask = thisargs[0](axesdict[thisdim])
                    slicedict[thisdim] = mask
                    if axesdict[thisdim] is not None:
                        axesdict[thisdim] = axesdict[thisdim][mask]
                elif thisop == hash("range"):
                    if axesdict[thisdim] is None:
                        raise ValueError(
                            "You passed a range-type slice"
                            + " selection, but to do that, your axis"
                            " coordinates need to"
                            + " be labeled! (The axis coordinates of"
                            f" {thisdim} aren't"
                            + " labeled)"
                        )
                    temp = np.diff(axesdict[thisdim])
                    if not all(temp * np.sign(temp[0]) > 0):
                        raise ValueError(
                            strm(
                                "you can only use the range format on data"
                                " where the axis is in consecutively"
                                " increasing or decreasing order, and the"
                                " differences that I see are",
                                temp * np.sign(temp[0]),
                            ),
                            "if you like, you can still do this by first"
                            " calling .sort( on the %s axis" % thisdim,
                        )
                    if np.sign(temp[0]) == -1:
                        thisaxis = axesdict[thisdim][::-1]
                    else:
                        thisaxis = axesdict[thisdim]
                    if len(thisargs) > 2:
                        raise ValueError(
                            "range with more than two values not currently"
                            " supported"
                        )
                    elif len(thisargs) == 1:
                        temp_low = thisargs[0]
                        temp_high = inf
                    else:
                        temp_low = thisargs[0]
                        temp_high = thisargs[1]
                        if temp_low is None:
                            temp_low = -inf
                        if temp_high is None:
                            temp_high = inf
                        if temp_low > temp_high:
                            temp_low, temp_high = temp_high, temp_low
                    # at this point, temp_low is indeed the lower value, and
                    # temp_high indeed the higher
                    logger.debug(
                        strm(
                            "after initial processing, range is",
                            temp_low,
                            temp_high,
                        )
                    )
                    if temp_low == inf:
                        raise ValueError(
                            strm(
                                "this is not going to work -- I interpret"
                                " range",
                                thisargs,
                                "I get to",
                                temp_low,
                                ",",
                                temp_high,
                            )
                        )
                    elif temp_low == -inf:
                        temp_low = 0
                    else:
                        logger.debug(strm("looking for", temp_low))
                        temp_low = np.searchsorted(thisaxis, temp_low)
                        if temp_low >= len(thisaxis):
                            raise ValueError(
                                'the lower value of your slice %s on the "%s"'
                                " axis (which runs from %g to %g) is higher"
                                " than the highest value of the axis"
                                " coordinates!"
                                % (
                                    (
                                        str((thisargs[0], thisargs[1])),
                                        thisdim,
                                    )
                                    + tuple(self.getaxis(thisdim)[r_[0, -1]])
                                )
                            )
                        logger.debug(
                            strm(
                                "i found",
                                thisaxis[temp_low],
                                "for the low end of the slice",
                                thisargs,
                            )
                        )
                    temp_high_float = temp_high
                    if temp_high == inf:
                        temp_high = len(
                            thisaxis
                        )  # not an exact match (inf doesn't match the index),
                        #    so needs to be inclusive already
                    elif temp_high == -inf:
                        raise ValueError(
                            strm(
                                "this is not going to work -- I interpret"
                                " range",
                                thisargs,
                                "I get to",
                                temp_low,
                                ",",
                                temp_high,
                            )
                        )
                    else:
                        logger.debug(strm("looking for", temp_high))
                        temp_high = np.searchsorted(thisaxis, temp_high)
                    # at this point, the result is inclusive if temp_high is
                    # not an exact match, but exclusive if it is
                    if (
                        temp_high < len(thisaxis)
                        and thisaxis[temp_high] == temp_high_float
                    ):
                        temp_high += 1  # make it inclusive
                    logger.debug(
                        strm(
                            "before looking at direction of axis, I have",
                            temp_low,
                            temp_high,
                        )
                    )
                    if np.sign(temp[0]) == -1:
                        logger.debug("identified descending axis")
                        temp_high = len(thisaxis) - temp_high
                        temp_low = len(thisaxis) - 1 - temp_low
                        temp_high, temp_low = temp_low, temp_high
                    del temp
                    if temp_low == temp_high:
                        temp_high += 1
                    slicedict[thisdim] = slice(
                        temp_low, temp_high, None
                    )  # inclusive
                    axesdict[thisdim] = axesdict[thisdim][slicedict[thisdim]]
                elif thisop == hash("idx"):
                    if thisdim not in axesdict.keys():
                        raise ValueError(f"{thisdim} not in {axesdict.keys()}")
                    if axesdict[thisdim] is None:
                        raise ValueError(
                            "You passed a labeled index"
                            + " selection, but to do that, your axis"
                            " coordinates need to"
                            + " be labeled! (The axis coordinates of"
                            f" {thisdim} aren't"
                            + " labeled)"
                        )
                    temp = abs(axesdict[thisdim] - thisargs[0]).argmin()
                    slicedict[thisdim] = temp
                    axesdict.pop(thisdim)
            if errordict is not None and errordict != np.array(None):
                for x, y in slicedict.items():
                    if errordict[x] is not None:
                        if np.isscalar(y):
                            errordict.pop(x)
                        elif isinstance(y, type(emptyfunction)):
                            mask = y(axesdict[x])
                            errordict[x] = errordict[x][mask]
                        else:
                            try:
                                errordict[x] = errordict[x][y]  # default
                            except Exception:
                                raise IndexError(
                                    strm(
                                        "Trying to index",
                                        errordict,
                                        "-->",
                                        x,
                                        "=",
                                        errordict[x],
                                        "with",
                                        y,
                                        "error started as",
                                        self.axis_coords_error,
                                    )
                                )
            if unitsdict is not None and unitsdict != np.array(None):
                for x, y in slicedict.items():
                    if unitsdict[x] is not None:
                        if np.isscalar(y):
                            unitsdict.pop(x)
            logger.debug(
                strm(
                    "Here is the slice dict:",
                    slicedict,
                    "and the axes dict",
                    axesdict,
                )
            )
            return slicedict, axesdict, errordict, unitsdict
            # }}}
        else:
            raise ValueError(
                strm(
                    "label your freaking dimensions! (type of args[0] is ",
                    type(args[0]),
                    "and it should be str!)",
                )
            )

    # }}}
    # {{{ hdf5 write
    def hdf5_write(self, h5path, directory="."):
        r"""Write the nddata to an HDF5 file.

        `h5path` is the name of the file followed by the node path where
        you want to put it -- it does **not** include the directory where
        the file lives.
        The directory can be passed to the `directory` argument.

        You can use either :func:`~pyspecdata.find_file` or
        :func:`~pyspecdata.nddata_hdf5` to read the data, as shown below.
        When reading this, please note that HDF5 files store *multiple*
        datasets,
        and each is named (here, the name is `test_data`).

        .. code-block:: python

            from pyspecdata import *
            init_logging('debug')
            a = nddata(r_[0:5:10j], 'x')
            a.name('test_data')
            try:
                a.hdf5_write('example.h5',getDATADIR(exp_type='Sam'))
            except Exception:
                print("file already exists, not creating again -- delete the
                file or node if wanted")
            # read the file by the "raw method"
            b = nddata_hdf5('example.h5/test_data',
                    getDATADIR(exp_type='Sam'))
            print("found data:",b)
            # or use the find file method
            c = find_file('example.h5', exp_type='Sam',
                    expno='test_data')
            print("found data:",c)

        Parameters
        ----------
        h5path : str
            The name of the file followed by the node path where
            you want to put it -- it does **not** include the directory where
            the file lives.
            (Because HDF5 files contain an internal directory-like group
            structure.)
        directory : str
            the directory where the HDF5 file lives.
        """
        for thisax in self.dimlabels:
            if self.getaxis(thisax) is None or len(self.getaxis(thisax)) == 0:
                raise ValueError(
                    strm(
                        "The axis",
                        thisax,
                        "appears not to have a label!  I refuse to save data"
                        " to HDF5 if you do not label all your axes!!",
                    )
                )
        # {{{ add the final node based on the name stored in the nddata
        #     structure
        if h5path[-1] != "/":
            h5path += "/"  # make sure it ends in a slash first
        try:
            thisname = self.get_prop("name")
        except Exception:
            raise ValueError(
                strm(
                    "You're trying to save an nddata object which",
                    "does not yet have a name, and you can't do this! Run",
                    "yourobject.name('setname')",
                )
            )
        if isinstance(thisname, str):
            h5path += thisname
        else:
            raise ValueError(
                strm(
                    "problem trying to store HDF5 file; you need to",
                    "set the ``name'' property of the nddata object to a"
                    " string",
                    "first!",
                )
            )
        h5file, bottomnode = h5nodebypath(
            h5path, directory=directory
        )  # open the file and move to the right node
        try:
            # print 'DEBUG 1: bottomnode is',bottomnode
            # }}}
            # {{{ print out the attributes of the data
            myattrs = normal_attrs(self)
            # {{{ separate them into data and axes
            mydataattrs = list(filter((lambda x: x[0:4] == "data"), myattrs))
            myotherattrs = list(filter((lambda x: x[0:4] != "data"), myattrs))
            myotherattrs = [
                x
                for x in myotherattrs
                if x not in ["C", "sin", "cos", "exp", "log10"]
            ]
            myaxisattrs = list(
                filter((lambda x: x[0:4] == "axis"), myotherattrs)
            )
            myotherattrs = list(
                filter((lambda x: x[0:4] != "axis"), myotherattrs)
            )
            logger.debug(
                strm(
                    lsafe(
                        "data attributes:",
                        list(
                            zip(
                                mydataattrs,
                                [
                                    type(self.__getattribute__(x))
                                    for x in mydataattrs
                                ],
                            )
                        ),
                    ),
                    "\n\n",
                )
            )
            logger.debug(
                strm(
                    lsafe(
                        "axis attributes:",
                        list(
                            zip(
                                myaxisattrs,
                                [
                                    type(self.__getattribute__(x))
                                    for x in myaxisattrs
                                ],
                            )
                        ),
                    ),
                    "\n\n",
                )
            )
            logger.debug(
                strm(
                    lsafe(
                        "other attributes:",
                        list(
                            zip(
                                myotherattrs,
                                [
                                    type(self.__getattribute__(x))
                                    for x in myotherattrs
                                ],
                            )
                        ),
                    ),
                    "\n\n",
                )
            )
            # }}}
            # }}}
            # {{{ write the data table
            if "data" in mydataattrs:
                if (
                    "data_error" in mydataattrs
                    and self.get_error() is not None
                    and len(self.get_error()) > 0
                ):
                    thistable = np.core.rec.fromarrays(
                        [self.data, self.get_error()], names="data,error"
                    )
                    mydataattrs.remove("data_error")
                else:
                    thistable = np.core.rec.fromarrays(
                        [self.data], names="data"
                    )
                mydataattrs.remove("data")
                datatable = h5table(bottomnode, "data", thistable)
                # print 'DEBUG 2: bottomnode is',bottomnode
                # print 'DEBUG 2: datatable is',datatable
                logger.debug(strm("Writing remaining axis attributes\n\n"))
                if len(mydataattrs) > 0:
                    h5attachattributes(datatable, mydataattrs, self)
            else:
                raise ValueError(
                    "I can't find the data object when trying to save the HDF5"
                    " file!!"
                )
            # }}}
            # {{{ write the axes tables
            if "axis_coords" in myaxisattrs:
                if len(self.axis_coords) > 0:
                    # {{{ create an 'axes' node
                    axesnode = h5child(
                        bottomnode,  # current node
                        "axes",  # the child
                        create=True,
                    )
                    # }}}
                    for j, axisname in enumerate(
                        self.dimlabels
                    ):  # make a table for each different dimension
                        myaxisattrsforthisdim = dict([
                            (x, self.__getattribute__(x)[j])
                            for x in list(myaxisattrs)
                            if len(self.__getattribute__(x)) > 0
                        ])  # collect the attributes for this dimension and
                        #    their values
                        logger.debug(
                            strm(
                                lsafe(
                                    "for axis",
                                    axisname,
                                    "myaxisattrsforthisdim=",
                                    myaxisattrsforthisdim,
                                )
                            )
                        )
                        if (
                            "axis_coords" in list(myaxisattrsforthisdim.keys())
                            and myaxisattrsforthisdim["axis_coords"]
                            is not None
                        ):
                            if (
                                "axis_coords_error"
                                in list(myaxisattrsforthisdim.keys())
                                and myaxisattrsforthisdim["axis_coords_error"]
                                is not None
                                and len(
                                    myaxisattrsforthisdim["axis_coords_error"]
                                )
                                > 0
                            ):  # this is needed to avoid all errors, though I
                                # guess I could use try/except
                                thistable = np.core.rec.fromarrays(
                                    [
                                        myaxisattrsforthisdim["axis_coords"],
                                        myaxisattrsforthisdim[
                                            "axis_coords_error"
                                        ],
                                    ],
                                    names="data,error",
                                )
                                myaxisattrsforthisdim.pop("axis_coords_error")
                            else:
                                thistable = np.core.rec.fromarrays(
                                    [myaxisattrsforthisdim["axis_coords"]],
                                    names="data",
                                )
                            myaxisattrsforthisdim.pop("axis_coords")
                        datatable = h5table(axesnode, axisname, thistable)
                        # print 'DEBUG 3: axesnode is',axesnode
                        logger.debug(
                            strm(
                                "Writing remaining axis attributes for",
                                axisname,
                                "\n\n",
                            )
                        )
                        if len(myaxisattrsforthisdim) > 0:
                            h5attachattributes(
                                datatable,
                                list(myaxisattrsforthisdim.keys()),
                                list(myaxisattrsforthisdim.values()),
                            )
            # }}}
            # {{{ Check the remaining attributes.
            logger.debug(
                strm(
                    lsafe(
                        "other attributes:",
                        list(
                            zip(
                                myotherattrs,
                                [
                                    type(self.__getattribute__(x))
                                    for x in myotherattrs
                                ],
                            )
                        ),
                    ),
                    "\n\n",
                )
            )
            logger.debug(strm("Writing remaining other attributes\n\n"))
            if len(myotherattrs) > 0:
                # print 'DEBUG 4: bottomnode is',bottomnode
                test = repr(
                    bottomnode
                )  # somehow, this prevents it from claiming that the
                #    bottomnode is None --> some type of bug?
                logger.debug(strm("bottomnode", test))
                h5attachattributes(
                    bottomnode,
                    [
                        j
                        for j in myotherattrs
                        if not self._contains_symbolic(j)
                    ],
                    self,
                )
                warnlist = [
                    j
                    for j in myotherattrs
                    if (not self._contains_symbolic(j))
                    and isinstance(self.__getattribute__(j), dict)
                ]
                # {{{ to avoid pickling, test that none of the attributes I'm
                #     trying to write are dictionaries or lists
                if len(warnlist) > 0:
                    print(
                        "WARNING!!, attributes", warnlist, "are dictionaries!"
                    )
                warnlist = [
                    j
                    for j in myotherattrs
                    if (not self._contains_symbolic(j))
                    and isinstance(self.__getattribute__(j), list)
                ]
                if len(warnlist) > 0:
                    print("WARNING!!, attributes", warnlist, "are lists!")
                # }}}
                logger.debug(
                    strm(
                        lsafe(
                            "other attributes:",
                            list(
                                zip(
                                    myotherattrs,
                                    [
                                        type(self.__getattribute__(x))
                                        for x in myotherattrs
                                    ],
                                )
                            ),
                        ),
                        "\n\n",
                    )
                )
            # }}}
        finally:
            h5file.close()

    # }}}


class testclass:
    def __getitem__(self, *args, **kwargs):
        print("you called __getitem__ with args", args, "and kwargs", kwargs)
        return

    def __getattribute__(self, *args, **kwargs):
        print(
            "you called __getattribute__ with args", args, "and kwargs", kwargs
        )
        return


class nddata_hdf5(nddata):
    def __repr__(self):
        if hasattr(self, "_node_children"):
            return repr(self.datanode)
        else:
            return nddata.__repr__(self)
        atexit.register(self._cleanup)

    def _cleanup(self):
        if hasattr(self, "_node_children"):
            self.h5file.close()
            del self.h5file
            del self.datanode
        return

    def __init__(self, pathstring, directory="."):
        self.pathstring = pathstring
        # try:
        self.h5file, self.datanode = h5nodebypath(
            pathstring, check_only=True, directory=directory
        )
        logger.debug("about to call _init_datanode")
        self._init_datanode(self.datanode)
        atexit.register(self._cleanup)

    def _init_datanode(self, datanode, **kwargs):
        datadict = h5loaddict(datanode)
        # {{{ load the data, and pop it from datadict
        try:
            datarecordarray = datadict["data"][
                "data"
            ]  # the table is called data, and the data of the table is called
            #    data
            mydata = datarecordarray["data"]
        except Exception:
            raise ValueError("I can't find the nddata.data")
        try:
            kwargs.update({"data_error": datarecordarray["error"]})
        except Exception:
            logger.debug(strm("No error found\n\n"))
        datadict.pop("data")
        # }}}
        # {{{ be sure to load the dimlabels
        mydimlabels = [j.decode("utf-8") for j in datadict["dimlabels"]]
        if len(mydimlabels) == 1:
            if len(mydimlabels[0]) == 1:
                mydimlabels = list(
                    [mydimlabels[0][0]]
                )  # for some reason, think I need to do this for length 1
        # }}}
        # {{{ load the axes and pop them from datadict
        datadict.pop("dimlabels")
        if "axes" in list(datadict.keys()):
            myaxiscoords = [None] * len(mydimlabels)
            myaxiscoordserror = [None] * len(mydimlabels)
            myaxis_units = [None] * len(mydimlabels)
            logger.debug(
                strm(
                    "about to read out the various axes:",
                    list(datadict["axes"].keys()),
                )
            )
            for axisname in list(datadict["axes"].keys()):
                try:
                    axisnumber = mydimlabels.index(axisname)
                except AttributeError:
                    raise AttributeError(
                        strm(
                            "mydimlabels is not in the right format!\nit looks"
                            " like this:\n",
                            mydimlabels,
                            type(mydimlabels),
                        )
                    )
                except ValueError:
                    raise ValueError(
                        strm(
                            "mydimlabels is not in the right format!\nit looks"
                            " like this:\n",
                            mydimlabels,
                            type(mydimlabels),
                        )
                    )
                recordarrayofaxis = datadict["axes"][axisname]["data"]
                if "axis_coords_units" in datadict["axes"][axisname].keys():
                    myaxis_units[axisnumber] = datadict["axes"][axisname][
                        "axis_coords_units"
                    ]
                else:
                    if ("Scans" in axisname) or ("ph" in axisname):
                        pass
                    else:
                        print(
                            "You didn't set units for %s before saving the"
                            " data!!!" % axisname
                        )
                myaxiscoords[axisnumber] = recordarrayofaxis["data"]
                if "error" in recordarrayofaxis.dtype.names:
                    myaxiscoordserror[axisnumber] = recordarrayofaxis["error"]
                datadict["axes"][axisname].pop("data")
                for k in list(datadict["axes"][axisname].keys()):
                    logger.debug(
                        strm(
                            "Warning, attribute",
                            k,
                            "of axis table",
                            axisname,
                            "remains, but the code to load this is not yet"
                            " supported",
                        )
                    )
                datadict["axes"].pop(axisname)
            kwargs.update({"axis_coords": myaxiscoords})
            kwargs.update({"axis_coords_units": myaxis_units})
            kwargs.update({"axis_coords_error": myaxiscoordserror})
        elif len(mydimlabels) > 1:
            raise ValueError(
                "The current version uses the axis labels to"
                "figure out the shape of the data\nBecause you stored"
                "unlabeled data, I can't figure out the shape of the"
                "data!!"
            )
            # the reshaping this refers to is done below
        # }}}
        logger.debug(
            strm(
                "about to initialize data with shape",
                mydata.shape,
                "labels",
                mydimlabels,
                "and kwargs",
                kwargs,
            )
        )
        nddata.__init__(self, mydata, mydata.shape, mydimlabels, **kwargs)
        # {{{ reshape multidimensional data to match the axes
        if len(mydimlabels) > 1:
            det_shape = []
            for thisdimlabel in mydimlabels:
                try:
                    temp = self.getaxis(thisdimlabel)
                except Exception:
                    temp = -1  # no axis is given
                if isinstance(temp, np.ndarray):
                    temp = len(temp)
                det_shape.append(temp)
            try:
                self.data = self.data.reshape(det_shape)
            except Exception:
                raise RuntimeError(
                    strm(
                        "The data is of shape",
                        self.data.shape,
                        "and I try to reshape it into",
                        det_shape,
                        "corresponding to the dimensions",
                        mydimlabels,
                        "--> this fails!",
                    )
                )
        # }}}
        for remainingattribute in list(datadict.keys()):
            self.__setattr__(remainingattribute, datadict[remainingattribute])
        self.h5file.close()
        del self.h5file
        del self.datanode
        return


# }}}


class ndshape(ndshape_base):
    r"""A class for describing the shape and dimension names of nddata objects.

    A main goal of this class is to allow easy generation (allocation) of new
    arrays -- see :func:`alloc`.

    """

    def alloc(self, dtype="complex128", labels=False, format=0):
        r"""Use the shape object to allocate an empty nddata object.

        Parameters
        ----------
        labels :
            Needs documentation
        format : 0, 1, or None
            What goes in the allocated array.
            `None` uses numpy empty.

        Example
        -------

        If you want to create new empty array that's 10x3 with dimensions "x"
        and "y":

        >>> result = ndshape([10,3],['x','y']).alloc(format=None)

        You can also do things like creating a new array based on the size of
        an existing array (create a new array without dimension x, but with new
        dimension z)

        >>> myshape = ndshape(mydata)
        >>> myshape.pop('x')
        >>> myshape + (10,'z')
        >>> result = myshape.alloc(format=None)
        """
        try:
            if format == 0:
                try:
                    emptyar = np.zeros(tuple(self.shape), dtype=dtype)
                except TypeError:
                    raise TypeError(
                        "You passed a type of "
                        + repr(dtype)
                        + ", which was likely not understood (you also passed"
                        " a shape of "
                        + repr(tuple(self.shape))
                        + ")"
                    )
            elif format == 1:
                emptyar = np.ones(tuple(self.shape), dtype=dtype)
            elif format is None:
                emptyar = np.empty(tuple(self.shape), dtype=dtype)
            else:
                emptyar = format * np.ones(tuple(self.shape), dtype=dtype)
        except TypeError:
            raise TypeError(
                strm(
                    "Wrong type for self.shape",
                    list(map(type, self.shape)),
                    "this probably means that you swapped the size and name"
                    " arguments -- ",
                    self.shape,
                    "should be numbers, not names",
                )
            )
        retval = nddata(emptyar, self.shape, self.dimlabels)
        if labels:
            retval.labels(
                self.dimlabels, [np.double(r_[0:x]) for x in self.shape]
            )
        return retval


image = this_plotting.image.image
# }}}


# {{{ fitdata
class fitdata(nddata):
    r"""Inherits from an nddata and enables curve fitting through use of a
    sympy expression.

    The user creates a fitdata class object from an existing nddata
    class object, and on this fitdata object can define the
    :func:`functional_form` of the curve it would like to fit to the
    data of the original nddata.
    This functional form must be provided as a sympy expression, with
    one of its variables matching the name of the dimension that the
    user would like to fit to.
    The user provides fit coefficients using :func:`fit_coeff` and
    obtains output using :func:`fit` and :func:`eval`.

    If you haven't done this before,
    create a jupyter notebook (not checked in, just for your own playing
    around) with:
    ```
    import sympy as s
    s.init_printing()
    ```
    you can then use `s.symbols(` to create symbols/variables that
    allow you to build the mathematical expression for your fitting
    function
    """

    def __init__(self, *args, **kwargs):
        # {{{ manual kwargs
        fit_axis = None
        if "fit_axis" in list(kwargs.keys()):
            fit_axis = kwargs.pop("fit_axis")
        # }}}
        if isinstance(args[0], nddata):
            myattrs = normal_attrs(args[0])
            for j in range(0, len(myattrs)):
                self.__setattr__(
                    myattrs[j], args[0].__getattribute__(myattrs[j])
                )
            # nddata.__init__(self,
            #        args[0].data,
            #        args[0].data.shape,
            #        args[0].dimlabels,
            #        axis_coords = args[0].axis_coords,
            #        data_error = args[0].data_error,
            #        axis_coords_error = args[0].axis_coords_error,
            #        axis_coords_units = args[0].axis_coords_units,
            #        data_units = args[0].data_units,
            #        other_info = args[0].other_info,
            #        **kwargs)
        else:
            # self.__base_init(*args,**kwargs)
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
        # {{{ in the class, only store the forced values and indices they are
        #     set to
        self.set_to = None
        self.set_indices = None
        self.active_indices = None
        # }}}
        return

    def parameter_derivatives(self, xvals, set=None, set_to=None):
        r"""return a matrix containing derivatives of the parameters, can set
        dict set, or keys set, vals set_to"""
        logger.debug(strm("parameter derivatives is called!"))
        if np.iscomplex(self.data.flatten()[0]):
            print(lsafen("Warning, taking only real part of fitting data!"))
        if isinstance(set, dict):
            set_to = list(set.values())
            set = list(set.keys())
        solution_list = dict([
            (
                (self.symbolic_dict[k], set_to[j])
                if k in set
                else (self.symbolic_dict[k], self.output(k))
            )
            for j, k in enumerate(self.symbol_list)
        ])  # load into the solution list
        number_of_i = len(xvals)
        parameters = self._active_symbols()
        mydiff_sym = [[]] * len(self.symbolic_vars)
        x = self.symbolic_x
        fprime = np.zeros([len(parameters), number_of_i])
        for j in range(0, len(parameters)):
            thisvar = self.symbolic_dict[parameters[j]]
            mydiff_sym[j] = sp.core.diff(self.symbolic_func, thisvar)
            try:
                mydiff = mydiff_sym[j].subs(solution_list)
            except Exception:
                raise ValueError(
                    strm(
                        "error trying to substitute",
                        mydiff_sym[j],
                        "with",
                        solution_list,
                    )
                )
            try:
                fprime[j, :] = np.array([
                    complex(mydiff.subs(x, xvals[k]))
                    for k in range(0, len(xvals))
                ])
            except ValueError:
                raise ValueError(
                    strm(
                        "Trying to set index",
                        j,
                        "np.shape(fprime)",
                        np.shape(fprime),
                        "np.shape(xvals)",
                        np.shape(xvals),
                        "the thing I'm trying to",
                        "compute looks like this",
                        [
                            mydiff.subs(x, xvals[k])
                            for k in range(0, len(xvals))
                        ],
                    )
                )
            except Exception:
                raise ValueError(
                    strm(
                        "Trying to set index",
                        j,
                        "np.shape(fprime)",
                        np.shape(fprime),
                        "np.shape(xvals)",
                        np.shape(xvals),
                    )
                )
        return fprime

    @property
    def function_string(self):
        r"""A property of the fitdata class which stores a string
        output of the functional form of the desired fit expression
        provided in func:`functional_form` in LaTeX format"""
        retval = sp.printing.latex(self.symbolic_expr).replace("$", "")
        return (
            r"$f(%s)=" % (sp.printing.latex(sp.Symbol(self.fit_axis)))
            + retval
            + r"$"
        )

    @function_string.setter
    def function_string(self):
        raise ValueError(
            "You cannot set the string directly -- change the functional_form"
            " property instead!"
        )

    @property
    def functional_form(self):
        r"""A property of the fitdata class which is set by the user,
        takes as input a sympy expression of the desired fit
        expression"""
        print("Getting symbolic function")
        return self.symbolic_expr

    @functional_form.setter
    def functional_form(self, sym_expr):
        r"""The functional form, given as a sympy expression, to
        which you would like to fit the data."""
        assert issympy(
            sym_expr
        ), "for now, the functional form must be a sympy expression!"
        self.symbolic_expr = sym_expr
        # {{{ adapted from fromaxis, trying to adapt the variable
        symbols_in_expr = self.symbolic_expr.atoms(sp.Symbol)
        logger.debug(
            strm(
                "identified this as a sympy expression (",
                self.symbolic_expr,
                ") with symbols",
                symbols_in_expr,
            )
        )
        symbols_in_expr = set(map(str, symbols_in_expr))
        # the next are the parameters
        self.fit_axis = set(self.dimlabels) & symbols_in_expr
        if len(self.fit_axis) == 0:
            raise ValueError(
                "I can't find any variables that might correspond to a"
                " dimension you want to fit along.The variables are",
                symbols_in_expr,
                "and the dimensions are",
                self.dimlabels,
            )
        elif len(self.fit_axis) > 1:
            raise ValueError(
                "currently only 1D fitting is supported, though this should be"
                " easyto change -- I see potential fit axes %s"
                % str(self.fit_axis)
            )
        # the next line gives the parameters
        self.symbolic_vars = symbols_in_expr - self.fit_axis
        # this gets used later in p_ini
        self.number_of_parameters = len(self.symbolic_vars)
        # }}}
        self.fit_axis = list(self.fit_axis)[0]
        # redefine as real to avoid weird piecewise derivatives
        self.fit_axis_sym = sp.core.var(self.fit_axis, real=True)
        self.symbolic_vars = list(self.symbolic_vars)
        self.symbolic_vars.sort()  # so I get consistent behavior
        self.symbolic_vars = [
            sp.core.var(j, real=True) for j in self.symbolic_vars
        ]
        self.symbol_list = [str(j) for j in self.symbolic_vars]
        args = self.symbolic_vars + [self.fit_axis]
        self.fitfunc_multiarg = sp.utilities.lambdify(
            tuple(args), self.symbolic_expr, modules=mat2array
        )

        def fn(p, x):
            p = self.add_inactive_p(p)
            assert len(p) == len(self.symbolic_vars), (
                "length of parameter passed to fitfunc doesn't match number of"
                " symbolic parameters"
            )
            return self.fitfunc_multiarg(*tuple(list(p) + [x]))

        self.fitfunc = fn
        # leave the gradient for later

    def analytical_covariance(self):
        r"""Not up to date"""
        covarmatrix = np.zeros([len(self._active_symbols())] * 2)
        # {{{ try this ppt suggestion --> his V is my fprime, but
        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis))
        dirproductform = False
        if dirproductform:
            sigma = self.get_error()
            f1 = fprime.shape[0]
            f2 = fprime.shape[1]
            fprime1 = fprime.reshape(f1, 1, f2)  # j index
            fprime2 = fprime.reshape(1, f1, f2)  # k index
            fprime_prod = fprime1 * fprime2
            fprime_prod = fprime_prod.reshape(-1, f2).T  # direct product form
            try:
                covarmat = np.dot(
                    np.linalg.pinv(fprime_prod), (sigma**2).reshape(-1, 1)
                )
            except ValueError:
                raise ValueError(
                    strm(
                        "shape of fprime_prod",
                        np.shape(fprime_prod),
                        "shape of inverse",
                        np.shape(np.linalg.pinv(fprime_prod)),
                        "shape of sigma",
                        np.shape(sigma),
                    )
                )
            covarmatrix = covarmat.reshape(f1, f1)
            for l in range(0, f1):
                for m in range(0, f1):
                    if l != m:
                        covarmatrix[l, m] /= 2
        else:
            # Note that lots of commented out code was deleted when this
            # comment was added
            sigma = self.get_error()
            J = np.matrix(fprime.T)
            S = np.matrix(self.get_covariance())
            Omegainv = S**-1
            print("a")
            minimizer = scipy.linalg.inv(J.T * Omegainv * J) * J.T * Omegainv
            covarmatrix = minimizer * S * minimizer.T
        # }}}
        return covarmatrix

    def copy(
        self,
    ):  # for some reason, if I don't override this with the same thing, it
        # doesn't override
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
        logger.debug("*** *** *** *** ***")
        logger.debug(str(this_set))
        logger.debug("*** *** *** *** ***")
        set_indices = list(
            map(self.symbol_list.index, this_set)
        )  # calculate indices once for efficiency
        active_mask = np.ones(len(self.symbol_list), dtype=bool)
        active_mask[set_indices] = (
            False  # generate the mask of indices that are actively fit
        )
        return set_indices, set_to, active_mask

    def remove_inactive_p(self, p):
        return p[self.active_mask]

    def add_inactive_p(self, p):
        if self.set_indices is not None:
            # {{{ uncollapse the function
            temp = p.copy()
            p = np.zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            # }}}
            p[self.set_indices] = (
                self.set_to
            )  # then just set the forced values to their given values
        return p

    def fitfunc(self, p, x):
        r"""This is the function that does the actual fitting, and takes a
        properly ordered list of parameters as well as an np.ndarray for the x
        axis."""
        raise ValueError(
            "this should have been overwritten when you set functional_form!"
        )

    def residual(self, p, x, y, sigma):
        """just the residual (error) function"""
        fit = self.fitfunc(p, x)
        normalization = np.sum(1.0 / sigma)
        # print 'DEBUG: y=',y,'\nfit=',fit,'\nsigma=',sigma,'\n\n'
        sigma[sigma == 0.0] = 1
        try:
            # as noted here: https://stackoverflow.com/questions/6949370/sc\
            #                ipy-scipy.optimize.leastsq-dfun-usage
            # this needs to be fit - y, not vice versa
            retval = (fit - y) / sigma * normalization
        except ValueError:
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
            )
        return retval

    def pinv(self, *args, **kwargs):
        retval = self.linear(*args, **kwargs)
        y = retval.data
        yerr = retval.get_error()
        x_axis = retval.dimlabels[0]
        x = retval.getaxis(x_axis)
        nopowerindex = np.argmax(x)
        mask = np.logical_not(r_[0 : len(x)] == nopowerindex)
        y = y[mask]
        yerr = yerr[mask]
        x = x[mask]
        L = c_[x.reshape((-1, 1)), np.ones((len(x), 1))]
        retval = np.dot(np.linalg.pinv(L, rcond=1e-17), y)
        logger.debug(
            strm(
                r"\label{fig:pinv_figure_text}y=",
                y,
                "yerr=",
                yerr,
                "%s=" % x_axis,
                x,
                "L=",
                L,
            )
        )
        logger.debug("\n\n")
        logger.debug(strm("recalc y = ", np.dot(L, retval)))
        logger.debug(strm("recalc E = ", 1.0 - 1.0 / np.dot(L, retval)))
        logger.debug(strm("actual E = ", self.data))
        return retval

    def linear(self, *args, **kwargs):
        r"""return the linear-form function, either smoothly along the fit
        function, or on the raw data, depending on whether or not the taxis
        argument is given
        can take optional arguments and pass them on to eval"""
        # print "DEBUG called linear"
        if len(args) == 1:
            taxis = self._taxis(args[0])  # handle integer as well
            return self.linfunc(
                taxis, self.eval(taxis, **kwargs).data
            )  # if we pass an argument, return the function across the entire
            #    time axis passed
        else:
            return self.linfunc(
                self.getaxis(self.fit_axis),
                self.data,
                yerr=self.get_error(),
                xerr=self.get_error(self.fit_axis),
            )  # otherwise, return the raw data

    def output(self, *name):
        r"""give the fit value of a particular symbol, or a dictionary of all
        values.

        Parameters
        ----------
        name: str (optional)
            name of the symbol.
            If no name is passed, then output returns a dictionary of the
            resulting values.

        Returns
        -------
        retval: dict or float
            Either a dictionary of all the values, or the value itself.
        """
        if not hasattr(self, "fit_coeff") or self.fit_coeff is None:
            return None
        p = self.fit_coeff.copy()
        if self.set_indices is not None:
            # {{{ uncollapse the function
            temp = p.copy()
            p = np.zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            # }}}
            p[self.set_indices] = (
                self.set_to
            )  # then just set the forced values to their given values
        # this should also be generic
        if len(name) == 1:
            try:
                return p[self.symbol_list.index(name[0])]
            except Exception:
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

    def _pn(self, name):
        return self.symbol_list.index(name)

    def _active_symbols(self):
        if not hasattr(self, "active_symbols"):
            if self.set_indices is not None:
                self.active_symbols = [
                    x
                    for x in self.symbol_list
                    if self.active_mask[self._pn(x)]
                ]
            else:
                self.active_symbols = list(self.symbol_list)
        return self.active_symbols

    def _pn_active(self, name):
        return self._active_symbols().index(name)

    def covar(self, *names):
        r"""give the covariance for the different symbols"""
        if len(names) == 1:
            names = [names[0], names[0]]
        if self.covariance is not None:
            return self.covariance[
                self._pn_active(names[0]), self._pn_active(names[1])
            ].copy()
        else:
            return None

    def covarmat(self, *names):
        if (len(names) == 1) and (names[0] == "recarray"):
            if hasattr(self, "active_mask"):
                active_symbols = [
                    x
                    for x in self.symbol_list
                    if self.active_mask[self._pn(x)]
                ]
            else:
                active_symbols = list(self.symbol_list)
            if len(active_symbols) != self.covariance.shape[0]:
                raise ValueError(
                    strm(
                        "length of active symbols",
                        active_symbols,
                        "doesnt match covariance matrix size(",
                        self.covariance.shape[0],
                        ")!",
                    )
                )
            recnames = ["labels"] + active_symbols
            recdata = []
            for j in range(0, self.covariance.shape[0]):
                thisdata = [active_symbols[j]] + list(
                    np.double(self.covariance[j, :].copy())
                )  # the first index is the row
                recdata.append(make_rec(thisdata, recnames))
            return r_[tuple(recdata)]
        if len(names) > 0:
            indices = list(
                map(self._pn_active, names)
            )  # slice out only these rows and columns
            return self.covariance[r_[indices], :][:, r_[indices]].copy()
        else:
            try:
                return self.covariance.copy()
            except Exception:
                return np.zeros([len(self.fit_coeff)] * 2, dtype="double")

    def latex(self):
        r"""show the latex string for the function, with all the symbols
        substituted by their values"""
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
            symbol = sp.printing.latex(self.symbolic_vars[j]).replace("$", "")
            logger.debug(strm('DEBUG: replacing symbol "', symbol, '"'))
            location = retval.find(symbol)
            while location != -1:
                if retval[location - 1] == "-":
                    newstring = (
                        retval[: location - 1]
                        + dp(-1 * p[j])
                        + retval[location + len(symbol) :]
                    )  # replace the symbol in the written function with the
                    #    appropriate number
                else:
                    newstring = (
                        retval[:location]
                        + dp(p[j])
                        + retval[location + len(symbol) :]
                    )  # replace the symbol in the written function with the
                    #    appropriate number
                logger.debug(
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
        logger.debug(
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

    def settoguess(self):
        "a debugging function, to easily plot the initial guess"
        self.fit_coeff = np.real(self.guess())
        return self

    def _taxis(self, taxis):
        r"""You can enter None, to get the fit along the same range as the
        data, an integer to give the number of points, or a range of data,
        which will return with 300 points"""
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

    def eval(self, taxis, set_what=None, set_to=None):
        """calculate the fit function along the axis taxis.

        Parameters
        ----------
        taxis: ndarray, int
            :if ndarray: the new axis coordinates along which we want to
                calculate the fit.
            :if int: number of evenly spaced points along the t-axis along the
                fit
        set_what: 'str', optional
            forcibly sets a specific symbol
        set_to: double, optional
            the specific value (int) you are assigning the symbol you included

        Returns
        -------
        self: nddata
            the fit function evaluated along the axis coordinates that were
            passed
        """
        if isinstance(set_what, dict):
            set_to = list(set_what.values())
            set_what = list(set_what.keys())
        taxis = self._taxis(taxis)
        if hasattr(self, "fit_coeff") and self.fit_coeff is not None:
            p = self.fit_coeff.copy()
        else:
            p = np.array([np.nan] * len(self.symbol_list))
        # {{{ LOCALLY apply any forced values
        # changed line below from set to set_what, and now it works
        if set_what is not None:
            if self.set_indices is not None:
                raise ValueError(
                    "you're trying to set indices in an eval"
                    " function for a function that was fit constrained; this"
                    " is not currently supported"
                )
            set_indices, set_to, active_mask = self.gen_indices(
                set_what, set_to
            )
            p[set_indices] = set_to
        # }}}
        # {{{ make a new, blank np.array with the fit axis expanded to fit
        #     taxis
        newdata = ndshape(self)
        newdata[self.fit_axis] = np.size(taxis)
        newdata = newdata.alloc()
        newdata.set_plot_color(self.get_plot_color())
        # }}}
        # {{{ keep all axis labels the same, except the expanded one
        newdata.axis_coords = list(newdata.axis_coords)
        newdata.labels([self.fit_axis], list([taxis]))
        # }}}
        newdata.data[:] = self.fitfunc(p, taxis).flatten()
        return newdata

    def makereal(self):
        self.data = np.real(self.data)
        return

    def rename(self, previous, new):
        if previous == self.fit_axis:
            self.fit_axis = new
        nddata.rename(self, previous, new)
        return self

    def fit(self, set_what=None, set_to=None, force_analytical=False):
        r"""actually run the fit"""
        if isinstance(set_what, dict):
            set_to = list(set_what.values())
            set_what = list(set_what.keys())
        x = self.getaxis(self.fit_axis)
        if np.iscomplex(self.data.flatten()[0]):
            logger.debug(
                strm("Warning, taking only real part of fitting data!")
            )
        y = np.real(self.data)
        sigma = self.get_error()
        if sigma is None:
            print(
                "{\\bf Warning:} You have no error associated with your plot,"
                " and I want to flag this for now\n\n"
            )
            warnings.warn(
                "You have no error associated with your plot, and I want to"
                " flag this for now",
                Warning,
            )
            sigma = np.ones(np.shape(y))
        if set_what is None:
            p_ini = self.guess()
        if set_what is not None:
            self.set_indices, self.set_to, self.active_mask = self.gen_indices(
                set_what, set_to
            )
            p_ini = self.remove_inactive_p(p_ini)
        leastsq_args = (self.residual, p_ini)
        leastsq_kwargs = {
            "args": (x, y, sigma),
            "full_output": True,
        }  # 'maxfev':1000*(len(p_ini)+1)}
        p_out, this_cov, infodict, mesg, success = scipy.optimize.leastsq(
            *leastsq_args, **leastsq_kwargs
        )
        try:
            p_out, this_cov, infodict, mesg, success = scipy.optimize.leastsq(
                *leastsq_args, **leastsq_kwargs
            )
        # {{{ just give various explicit errors
        except TypeError:
            if not isinstance(x, np.ndarray) and not isinstance(y, np.ndarray):
                raise TypeError(
                    strm(
                        "scipy.optimize.leastsq failed because the two arrays",
                        "aren't of the right",
                        "type",
                        "type(x):",
                        type(x),
                        "type(y):",
                        type(y),
                    )
                )
            else:
                if np.any(np.shape(x) != np.shape(y)):
                    raise RuntimeError(
                        strm(
                            "scipy.optimize.leastsq failed because the two"
                            " arrays donot match in size size",
                            "np.shape(x):",
                            np.shape(x),
                            "np.shape(y):",
                            np.shape(y),
                        )
                    )
            raise TypeError(
                strm(
                    "scipy.optimize.leastsq failed because of a type error!",
                    "type(x):",
                    showtype(x),
                    "type(y):",
                    showtype(y),
                    "type(sigma)",
                    showtype(sigma),
                    "np.shape(x):",
                    np.shape(x),
                    "np.shape(y):",
                    np.shape(y),
                    "np.shape(sigma)",
                    np.shape(sigma),
                    "p_ini",
                    type(p_ini),
                    p_ini,
                )
            )
        except ValueError as err:
            raise ValueError(
                strm(
                    'scipy.optimize.leastsq failed with "',
                    err,
                    '", maybe there is something wrong with the input:',
                    self,
                )
            )
        except Exception:
            raise ValueError("scipy.optimize.leastsq failed; I don't know why")
        # }}}
        if success not in [1, 2, 3, 4]:
            # {{{ up maximum number of evals
            if mesg.find("maxfev"):
                leastsq_kwargs.update({"maxfev": 50000})
                (
                    p_out,
                    this_cov,
                    infodict,
                    mesg,
                    success,
                ) = scipy.optimize.leastsq(*leastsq_args, **leastsq_kwargs)
                if success != 1:
                    if mesg.find("two consecutive iterates"):
                        print(
                            r"{\Large\color{red}{\bf Warning data is not"
                            r" fit!!! output shown for debug purposes only!}}",
                            "\n\n",
                        )
                        print(
                            r"{\color{red}{\bf Original message:}",
                            lsafe(mesg),
                            "}",
                            "\n\n",
                        )
                        infodict_keys = list(infodict.keys())
                        infodict_vals = list(infodict.values())
                        if "nfev" in infodict_keys:
                            infodict_keys[infodict_keys.index("nfev")] = (
                                "nfev, number of function calls"
                            )
                        if "fvec" in infodict_keys:
                            infodict_keys[infodict_keys.index("fvec")] = (
                                "fvec, the function evaluated at the output"
                            )
                        if "fjac" in infodict_keys:
                            infodict_keys[infodict_keys.index("fjac")] = (
                                "fjac, A permutation of the R matrix of a QR"
                                " factorization of the final approximate"
                                " Jacobian matrix, stored column wise."
                                " Together with ipvt, the covariance of the"
                                " estimate can be approximated."
                            )
                        if "ipvt" in infodict_keys:
                            infodict_keys[infodict_keys.index("ipvt")] = (
                                "ipvt, an integer np.array of length N which"
                                " defines a permutation matrix, p, such that"
                                " fjac*p = q*r, where r is upper triangular"
                                " with diagonal elements of nonincreasing"
                                " magnitude.  Column j of p is column ipvt(j)"
                                " of the identity matrix"
                            )
                        if "qtf" in infodict_keys:
                            infodict_keys[infodict_keys.index("qtf")] = (
                                "qtf, the vector (transpose(q)*fvec)"
                            )
                        for k, v in zip(infodict_keys, infodict_vals):
                            print(r"{\color{red}{\bf %s:}%s}" % (k, v), "\n\n")
                        # self.fit_coeff = None
                        # self.settoguess()
                        # return
                    else:
                        raise RuntimeError(
                            strm(
                                "scipy.optimize.leastsq finished with an error"
                                " message:",
                                mesg,
                            )
                        )
                    # }}}
            else:
                raise RuntimeError(
                    strm(
                        "scipy.optimize.leastsq finished with an error"
                        " message:",
                        mesg,
                    )
                )
        else:
            logger.debug(
                "Fit finished successfully with a code of %d and a message"
                " ``%s''" % (success, mesg)
            )
        self.fit_coeff = p_out  # note that this is stored in HIDDEN form
        dof = len(x) - len(p_out)
        if hasattr(self, "symbolic_x") and force_analytical:
            self.covariance = self.analytical_covariance()
        else:
            if force_analytical:
                raise RuntimeError(
                    strm(
                        "I can't take the analytical",
                        "covariance!  This is problematic.",
                    )
                )
            if this_cov is None:
                print(
                    r"{\color{red}"
                    + lsafen(
                        "cov is none! why?!, x=",
                        x,
                        "y=",
                        y,
                        "sigma=",
                        sigma,
                        "p_out=",
                        p_out,
                        "success=",
                        success,
                        "output:",
                        p_out,
                        this_cov,
                        infodict,
                        mesg,
                        success,
                    ),
                    "}\n",
                )
            self.covariance = this_cov
        if self.covariance is not None:
            try:
                self.covariance *= (
                    np.sum(infodict["fvec"] ** 2) / dof
                )  # scale by chi_v "RMS of residuals"
            except TypeError:
                raise TypeError(
                    strm(
                        "type(self.covariance)",
                        type(self.covariance),
                        "type(infodict[fvec])",
                        type(infodict["fvec"]),
                        "type(dof)",
                        type(dof),
                    )
                )
        logger.debug(
            strm(
                "at end of fit covariance is shape",
                np.shape(self.covariance),
                "fit coeff shape",
                np.shape(self.fit_coeff),
            )
        )
        return

    def bootstrap(
        self, points, swap_out=np.exp(-1.0), minbounds={}, maxbounds={}
    ):
        print(r"\begin{verbatim}")
        fitparameters = list(self.symbol_list)
        recordlist = np.array(
            [tuple([0] * len(fitparameters))] * points,
            {
                "names": tuple(fitparameters),
                "formats": tuple(["double"] * len(fitparameters)),
            },
        )  # make an instance of the recordlist
        for runno in range(0, points):
            success = False  # because sometimes this doesn't work
            while success is False:
                thiscopy = self.copy()
                # {{{ discard datapoints
                origsizecheck = np.double(np.size(thiscopy.data))
                mask = thiscopy.random_mask(
                    thiscopy.fit_axis, threshold=swap_out
                )
                thiscopy.data = thiscopy.data[mask]
                derr = thiscopy.get_error()
                x = thiscopy.getaxis(thiscopy.fit_axis)
                x = x[mask]  # note that x is probably no longer a pointer
                derr = derr[mask]
                # }}}
                # {{{ now extend
                number_to_replace = origsizecheck - thiscopy.data.size
                # print 'DEBUG: number_to_replace',number_to_replace
                random_indices = np.int32(
                    (
                        np.random.rand(number_to_replace)
                        * (thiscopy.data.size - 1.0)
                    ).round()
                )
                thiscopy.data = r_[
                    thiscopy.data, thiscopy.data.copy()[random_indices]
                ]
                thiscopy.labels(
                    [thiscopy.fit_axis], [r_[x, x.copy()[random_indices]]]
                )
                thiscopy.set_error(r_[derr, derr.copy()[random_indices]])
                # }}}
                try:
                    thiscopy.fit()
                    success = True
                    if len(minbounds) > 0:
                        for k, v in minbounds.items():
                            if thiscopy.output(k) < v:
                                success = False
                    if len(maxbounds) > 0:
                        for k, v in maxbounds.items():
                            if thiscopy.output(k) > v:
                                success = False
                except Exception:
                    # print 'WARNING, didn\'t fit'
                    success = False
                # here, use the internal routines, in case there are
                # constraints, etc
                if success is True:
                    for (
                        name
                    ) in thiscopy.symbol_list:  # loop over all fit coeff
                        recordlist[runno][name] = thiscopy.output(name)
        print(r"\end{verbatim}")
        return recordlist  # collect into a single recordlist np.array

    def set_guess(self, dict_of_values):
        """sets parameters to guess/estimated value to compare fit.

        Parameters
        ----------
        dict_of_values: dict
            dictionary of values set to parameters in fit equation.
            Allows for the setting of multiple variables depending on
            what's defined in this dictionary. The keys of the dictionary
            must be sympy symbols

        Returns
        -------
        self: nddata
            The modified nddata
        """

        input_guesses = set(dict_of_values.keys())
        print(input_guesses)
        symbols_not_present = input_guesses - set(self.symbolic_vars)
        if len(symbols_not_present) > 0:
            raise ValueError(
                strm(
                    "You specified the symbol(s)",
                    symbols_not_present,
                    "but I can't find this in the symbols for the fitting"
                    " function, which are",
                    self.symbolic_vars,
                )
            )
        symbols_not_set = set(self.symbolic_vars) - input_guesses
        self.guess_dict = dict_of_values
        self.guess_dict.update({k: 1 for k in symbols_not_set})
        return self

    def guess(self, use_pseudoinverse=False):
        r"""old code that I am preserving here -- provide the guess for our
        parameters; by default, based on pseudoinverse"""
        if use_pseudoinverse:
            self.has_grad = False
            if np.iscomplex(self.data.flatten()[0]):
                print(
                    lsafen("Warning, taking only real part of fitting data!")
                )
            y = np.real(self.data)
            # I ended up doing the following, because as it turns out
            # T1 is a bad fit function, because it has a singularity!
            # this is probably why it freaks out if I set this to zero
            # on the other hand, setting a value of one seems to be
            # bad for very short T1 samples
            which_starting_guess = 0
            thisguess = self.starting_guesses[which_starting_guess]
            numguesssteps = 20
            # {{{ for some reason (not sure) adding a dimension to y
            new_y_shape = list(y.shape)
            new_y_shape.append(1)
            y = y.reshape(tuple(new_y_shape))
            # }}}
            # {{{ evaluate f, fprime and residuals
            guess_dict = dict(list(zip(self.symbol_list, list(thisguess))))
            fprime = self.parameter_derivatives(
                self.getaxis(self.fit_axis), set=guess_dict
            )
            f_at_guess = np.real(self.eval(None, set=guess_dict).data)
            try:
                f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
            except Exception:
                raise ValueError(
                    strm(
                        "trying to reshape f_at_ini_guess from",
                        f_at_guess.shape,
                        "to",
                        new_y_shape,
                    )
                )
            thisresidual = sqrt((y - f_at_guess) ** 2).sum()
            # }}}
            lastresidual = thisresidual
            for j in range(0, numguesssteps):
                logger.debug(
                    strm(
                        "\n\n.core.guess) " + r"\begin{verbatim} fprime = \n",
                        fprime,
                        "\nf_at_guess\n",
                        f_at_guess,
                        "y=\n",
                        y,
                        "\n",
                        r"\end{verbatim}",
                    )
                )
                logger.debug(
                    strm(
                        "\n\n.core.guess) shape of parameter derivatives",
                        np.shape(fprime),
                        "shape of output",
                        np.shape(y),
                        "\n\n",
                    )
                )
                regularization_bad = True
                alpha_max = 100.0
                alpha_mult = 2.0
                alpha = 0.1  # maybe I can rather estimate this based on the
                #              change in the residual, similar to in L-M?
                logger.debug(
                    strm(
                        "\n\n.core.guess) value of residual before"
                        " regularization %d:" % j,
                        thisresidual,
                    )
                )
                while regularization_bad:
                    newguess = np.real(
                        np.array(thisguess)
                        + np.dot(
                            pinvr(fprime.T, alpha), (y - f_at_guess)
                        ).flatten()
                    )
                    mask = newguess < self.guess_lb
                    newguess[mask] = self.guess_lb[mask]
                    mask = newguess > self.guess_ub
                    newguess[mask] = self.guess_ub[mask]
                    if np.any(np.isnan(newguess)):
                        logger.debug(
                            strm(
                                "\n\n.core.guess) Regularization blows up"
                                " $\\rightarrow$ increasing $\\alpha$ to"
                                " %0.1f\n\n" % alpha
                            )
                        )
                        alpha *= alpha_mult
                    else:
                        # {{{ evaluate f, fprime and residuals
                        guess_dict = dict(
                            list(zip(self.symbol_list, list(newguess)))
                        )
                        # only evaluate fprime once we know this is good, below
                        f_at_guess = np.real(
                            self.eval(None, set=guess_dict).data
                        )
                        try:
                            f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
                        except Exception:
                            raise IndexError(
                                strm(
                                    "trying to reshape f_at_ini_guess from",
                                    f_at_guess.shape,
                                    "to",
                                    new_y_shape,
                                )
                            )
                        thisresidual = sqrt((y - f_at_guess) ** 2).sum()
                        # }}}
                        if (thisresidual - lastresidual) / lastresidual > 0.10:
                            alpha *= alpha_mult
                            logger.debug(
                                strm(
                                    "\n\n.core.guess) Regularized Pinv gave a"
                                    " step uphill $\\rightarrow$ increasing"
                                    " $\\alpha$ to %0.1f\n\n" % alpha
                                )
                            )
                        else:  # accept the step
                            regularization_bad = False
                            thisguess = newguess
                            lastresidual = thisresidual
                            fprime = self.parameter_derivatives(
                                self.getaxis(self.fit_axis), set=guess_dict
                            )
                    if alpha > alpha_max:
                        print(
                            "\n\n.core.guess) I can't find a new guess without"
                            " increasing the alpha beyond %d\n\n" % alpha_max
                        )
                        if (
                            which_starting_guess
                            >= len(self.starting_guesses) - 1
                        ):
                            print(
                                "\n\n.core.guess) {\\color{red} Warning!!!}"
                                " ran out of guesses!!!%d\n\n" % alpha_max
                            )
                            return thisguess
                        else:
                            which_starting_guess += 1
                            thisguess = self.starting_guesses[
                                which_starting_guess
                            ]
                            print(
                                "\n\n.core.guess) try a new starting guess:",
                                lsafen(thisguess),
                            )
                            j = 0  # restart the loop
                            # {{{ evaluate f, fprime and residuals for the new
                            #     starting guess
                            guess_dict = dict(
                                list(zip(self.symbol_list, list(thisguess)))
                            )
                            fprime = self.parameter_derivatives(
                                self.getaxis(self.fit_axis), set=guess_dict
                            )
                            f_at_guess = np.real(
                                self.eval(None, set=guess_dict).data
                            )
                            try:
                                f_at_guess = f_at_guess.reshape(
                                    tuple(new_y_shape)
                                )
                            except Exception:
                                raise RuntimeError(
                                    strm(
                                        "trying to reshape f_at_ini_guess"
                                        " from",
                                        f_at_guess.shape,
                                        "to",
                                        new_y_shape,
                                    )
                                )
                            thisresidual = sqrt((y - f_at_guess) ** 2).sum()
                            # }}}
                            regularization_bad = False  # jump out of this loop
                logger.debug(
                    strm(
                        "\n\n.core.guess) new value of guess after"
                        " regularization:",
                        lsafen(newguess),
                    )
                )
                logger.debug(
                    strm(
                        "\n\n.core.guess) value of residual after"
                        " regularization:",
                        thisresidual,
                    )
                )
            return thisguess
        else:
            if hasattr(self, "guess_dict"):
                return [self.guess_dict[k] for k in self.symbolic_vars]
            else:
                return [1.0] * self.number_of_parameters


# }}}
def sqrt(arg):
    if isinstance(arg, nddata):
        return arg**0.5
    elif isinstance(arg, sp.Symbol):
        return sympy_sqrt(arg)
    else:
        return np.sqrt(arg)


# {{{ determine the figure style, and load the appropriate modules
# this must come at end to prevent circular imports
if _figure_mode_setting == "latex":
    from .fornotebook import figlistl, obs

    figlist_var = figlistl
    lsafe = orig_lsafe
elif _figure_mode_setting == "standard":
    from .figlist import figlist

    def obsn(*x):  # because this is used in fornotebook, and I want it defined
        print("".join(x), "\n")

    def obs(*x):  # because this is used in fornotebook, and I want it defined
        print("".join(map(repr, x)))

    def lrecordarray(*x, **kwargs):
        return repr(
            x
        )  # if I'm not using tex, it's easier to not use the formatting

    def lsafe(*string, **kwargs):
        "replacement for normal lsafe -- no escaping"
        if len(string) > 1:
            lsafewkargs = lambda x: lsafe(x, **kwargs)
            return " ".join(list(map(lsafewkargs, string)))
        else:
            string = string[0]
        # {{{ kwargs
        if "wrap" in list(kwargs.keys()):
            wrap = kwargs.pop("wrap")
        else:
            wrap = None
        # }}}
        if not isinstance(string, str):
            string = repr(string)
        if wrap is True:
            wrap = 60
        if wrap is not None:
            string = "\n".join(textwrap.wrap(string, wrap))
        return string

    figlist_var = figlist
else:
    raise ValueError(
        "I don't understand the figures mode " + _figure_mode_setting
    )
# }}}
