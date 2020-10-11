r'''Provides the core components of pyspecdata.
Currently, this is a very large file that we will slowly break down into separate modules or packages.

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

There are many helper and utility functions that need to be sorted an documented by JMF,
and can be ignored.
These are somewhat wide-ranging in nature.
For example, :func:`box_muller` is a helper function (based on numerical recipes) used by :func:`nddata.add_noise`,
while h5 functions are helper functions for using pytables in a fashion that
will hopefull be intuitive to those familiar with SQL, etc.
'''
from .datadir import _my_config
from sys import exc_info
from os import listdir,environ
from os.path import sep as path_sep

# {{{ determine the figure style, and load the appropriate modules
_figure_mode_setting = _my_config.get_setting('figures', section='mode', environ='pyspecdata_figures')
if _figure_mode_setting is None:
    print("Warning!  Figure mode is not set, so I'm going to set it to standard by default!!!")
    _figure_mode_setting = 'standard'
    _my_config.set_setting('mode','figures','standard')
if _figure_mode_setting == 'latex':
    environ['ETS_TOOLKIT'] = 'qt4'
    import matplotlib; matplotlib.use('Agg')
# }}} -- continued below
from .general_functions import inside_sphinx
if not inside_sphinx():
    from pylab import *
else:
    pi = 3.14
from types import FunctionType as function
import textwrap
import atexit
import matplotlib
import matplotlib.transforms as mtransforms
from distutils.version import LooseVersion
from numpy import sqrt as np_sqrt
from numpy.lib.recfunctions import rename_fields,drop_fields
from mpl_toolkits.mplot3d import axes3d
from matplotlib.collections import PolyCollection
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D
from scipy.interpolate import griddata as scipy_griddata
if not inside_sphinx():
    import tables
import warnings
import re
from inspect import ismethod,signature,Parameter
from numpy.core import rec
from matplotlib.pyplot import cm
from copy import deepcopy 
import traceback
import sympy# doesn't like to be imported from fornotebook as part of a *
from scipy.optimize import leastsq
from scipy.signal import fftconvolve
import scipy.sparse as sparse
import numpy.lib.recfunctions as recf
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from .datadir import getDATADIR,log_fname
from . import fourier as this_fourier
from . import axis_manipulation
from . import nnls as this_nnls
from . import plot_funcs as this_plotting
from .general_functions import *
from .ndshape import ndshape_base
#rc('image',aspect='auto',interpolation='bilinear') # don't use this, because it gives weird figures in the pdf
rc('image',aspect='auto',interpolation='nearest')
#rcParams['text.usetex'] = True
rc('font', family='Arial')# I need this to render unicode
rcParams['xtick.direction'] = 'out'
rcParams['xtick.major.size'] = 12
rcParams['xtick.minor.size'] = 6
rcParams['ytick.direction'] = 'out'
rcParams['ytick.major.size'] = 12
rcParams['ytick.minor.size'] = 6
#rcParams['lines.linewidth'] = 3.0
rcParams['legend.fontsize'] = 12
rcParams['axes.grid'] = False
rcParams['font.size'] = 18
rcParams['image.cmap'] = 'jet'
rcParams['figure.figsize']=(16,12)
if inside_sphinx():
    mat2array = []
else:
    mat2array = [{'ImmutableMatrix': array}, 'numpy']# for sympy returns arrays rather than the stupid matrix class
logger = logging.getLogger('pyspecdata.core')
#{{{ constants
k_B = 1.380648813e-23
mu_0 = 4e-7*pi
mu_B = 9.27400968e-24#Bohr magneton
epsilon_0 = 8.854187817e-12
hbar = 6.6260695729e-34/2./pi
N_A = 6.02214179e23
gammabar_H = 4.258e7
gammabar_e = 2.807e10 # this is for a nitroxide
#}}}
def apply_oom(average_oom,numbers,prev_label=''):
    """scale numbers by the order of magnitude average_oom and change the
    name of the units by adding the appropriate SI prefix

    Parameters
    ----------
    average_oom: int or float
        the average order of magnitude to use
    numbers: ndarray
        The numbers to be scaled by average_oom.
        The array is modified in-place.
    prev_label: str
        a string representing the units

    Returns
    -------
    new_label: str
        prev_label is prefixed by the appropriate SI prefix
    """
    oom_names =   ['T' , 'G' , 'M' , 'k' , '' , 'm' , '\\mu ' , 'n' , 'p']
    oom_values = r_[12 , 9   , 6   , 3   , 0  , -3  , -6     , -9  , -12]
    eq = oom_values == average_oom
    if not any(eq):
        if all(average_oom < oom_values):
            oom_index = len(oom_values)-1
        elif all(average_oom > oom_values):
            oom_index = 0
        else:
            raise ValueError(strm("you passed",average_oom,"which I can't find a prefix for"))
    else:
        oom_index = nonzero(eq)[0][0]
    if numbers.dtype in ['int32','int64']:
        # this is not necessary unless we have an integer type
        logger.warning("you are trying to determine the SI prefix of a"
                "set of numbers that are described by integers.  This is"
                "probably not a good idea!!")
        new_values = numbers / 10.**oom_values[oom_index]
        numbers[:] = new_values.astype(numbers.dtype)
    else:
        numbers[:] /= 10.**oom_values[oom_index]
    return oom_names[oom_index]+prev_label
def mybasicfunction(first_figure = None):
    r'''this gives the format for doing the image thing
note also
nextfigure(fl,'name')
and
nextfigure({'lplotproperty':value})
'''
    fl = figlistini_old(first_figure)
    return figlistret(first_figure,figurelist,other)

def issympy(x):
    'tests if something is sympy (based on the module name)'
    return isinstance(x,sympy.Expr)

#{{{ function trickery
def mydiff(data,axis = -1):
    '''this will replace diff with a version that has the same number of indices, with the last being the copy of the first'''
    newdata = zeros(shape(data),dtype = data.dtype)
    indices = [slice(None,None,None)]*len(data.shape)
    indices[axis] = slice(None,-1,None)
    newdata[indices] = diff(data,axis = axis)
    #setfrom = list(indices)
    #indices[axis] = -1
    #setfrom[axis] = 0
    #newdata[indices] = newdata[setfrom]
    return newdata
#}}}
def normal_attrs(obj):
    myattrs = [x for x in dir(obj) if not ismethod(obj.__getattribute__(x))]
    myattrs = [x for x in myattrs if not x[0:2] == '__']
    # next line filters out properties
    myattrs = [x for x in myattrs if x not in ['C','angle','imag','real']]
    return myattrs
def showtype(x):
    if isinstance(x, ndarray):
        return ndarray,x.dtype
    else:
        return type(x)
def emptyfunction():
    pass
#{{{ structured array helper functions
def make_bar_graph_indices(mystructarray,list_of_text_fields,
        recursion_depth = 0,
        spacing = 0.1):
    r"This is a recursive function that is used as part of textlabel_bargraph; it does NOT work without the sorting given at the beginning of that function"
    #{{{ if there are still text fields left, then break down the array further, otherwise, just return the indices for this subarray
    if len(list_of_text_fields) > 0:
        unique_values = unique(mystructarray[list_of_text_fields[0]])# the return_index argument doesn't do what it's supposed to all the time, so I have to manually find the start indices, as given in the following line
        start_indices = [nonzero(mystructarray[list_of_text_fields[0]] == val)[0][0] for val in unique_values]
        # find the structured array for the unique value
        index_values = []
        label_values = []
        start_indices = r_[start_indices,len(mystructarray)] # I add this so I can do the next step
        logger.debug(strm('recursion depth is',recursion_depth,'and I am analyzing',list_of_text_fields[0],': '))
        logger.debug(strm('I found these unique values:',unique_values,'at these start indices:',start_indices[:-1]))
        for k in range(0,len(start_indices)-1):
            logger.debug(strm('recursion depth is',recursion_depth,'and I am analyzing',list_of_text_fields[0],': '))
            logger.debug(strm('trying to extract unique value',unique_values[k],'using the range',start_indices[k],start_indices[k+1]))
            logger.debug(strm('which has this data'))
            indiv_struct_array = mystructarray[start_indices[k]:start_indices[k+1]]
            logger.debug(strm(lsafen(indiv_struct_array)))
            these_index_values,these_labels = make_bar_graph_indices(indiv_struct_array,list_of_text_fields[1:],recursion_depth = recursion_depth+1)
            index_values.append(these_index_values)
            label_values.append([str(unique_values[k])+','+j for j in these_labels])
        #{{{ scale the result of each call down to the equal size (regardless of number of elements), shift by the position in this array, and return
        logger.debug(strm('recursion depth is',recursion_depth,'and I just COMPLETED THE LOOP, which gives a list of index values like this',index_values))
        max_indices = max(array(list(map(len,index_values)),dtype='double'))# the maximum width of the array inside
        index_values = [x+(max_indices-len(x))/2.0 for x in index_values]# if the bar is less than max indices, shift it over, so it's still in the center
        logger.debug(strm('recursion depth is',recursion_depth,'and I centered each set like this',index_values))
        index_values = [x/max_indices*(1-spacing)+(1-spacing)/2 for x in index_values]# scale down, so the width from left edge of bar to right edge of largest bar runs 0--> 1
        logger.debug(strm('recursion depth is',recursion_depth,'and I scaled down so each runs zero to one*(1-spacing) (centered) like this',index_values))
        # this adds an index value, and also collapses down to a single dimension list
        retval_indices = [x+num for num,val in enumerate(index_values) for x in val]
        # now collapse labels down to a single dimension
        retval_labels = [k for j in label_values for k in j]
        logger.debug(strm('recursion depth is',recursion_depth,'and I am passing up indices',retval_indices,'and labels',retval_labels))
        return retval_indices,retval_labels
        #}}}
    else:
        logger.debug(strm('recursion depth is',recursion_depth))
        N = len(mystructarray)
        logger.debug(strm('hit innermost (no text labels left) and passing up a list of indices that looks like this:',r_[0:N]))
        return r_[0:N],['']*N
    #}}}
def textlabel_bargraph(mystructarray,othersort = None,spacing = 0.1,ax = None,tickfontsize = 8):
    if ax is None:
        thisfig = gcf()
        ax = thisfig.add_axes([0.2,0.5,0.8,0.5])
        try:
            ax.tick_params(axis = 'both',which = 'major',labelsize = tickfontsize)
            ax.tick_params(axis = 'both',which = 'minor',labelsize = tickfontsize)
        except:
            print('Warning, in this version I can\'t set the tick params method for the axis')
    #{{{ find the text fields, put them first, and sort by them
    mystructarray = mystructarray.copy()
    list_of_text_fields = [str(j[0]) for j in mystructarray.dtype.descr if j[1][0:2] == '|S']
    mystructarray = sorted(mystructarray[list_of_text_fields + [x[0]
        for x in mystructarray.dtype.descr
        if x[0] not in list_of_text_fields]])
    logger.debug(strm('test --> now, it has this form:',lsafen(mystructarray)))
    #}}}
    error_fields = [str(j) for j in mystructarray.dtype.names if j[-6:] == '_ERROR']
    if len(error_fields) > 0:
        mystructarray_errors = mystructarray[error_fields]
        logger.debug("found error fields:",mystructarray_errors)
    mystructarray = mystructarray[[str(j) for j in mystructarray.dtype.names if j not in error_fields]]
    if othersort is not None:
        list_of_text_fields.append(othersort)
    logger.debug(strm('list of text fields is',lsafen(list_of_text_fields)))
    indices,labels = make_bar_graph_indices(mystructarray,list_of_text_fields,spacing = spacing)
    temp = list(zip(indices,labels))
    logger.debug(strm('(indices,labels) (len %d):'%len(temp),lsafen(temp)))
    logger.debug(strm('I get these labels (len %d):'%len(labels),labels,'for the data (len %d)'%len(mystructarray),lsafen(mystructarray)))
    indices = array(indices)
    indiv_width = min(diff(indices))*(1-spacing)
    remaining_fields = [x for x in mystructarray.dtype.names if x not in list_of_text_fields] # so they are in the right order, since set does not preserve order
    logger.debug(strm('The list of remaining (i.e. non-text) fields is',lsafen(remaining_fields)))
    colors = ['b','g','r','c','m','k']
    rects = []
    for j,thisfield in enumerate(remaining_fields):
        field_bar_width = indiv_width/len(remaining_fields)
        thiserror = None
        if thisfield+'_ERROR' in error_fields:
            thiserror = mystructarray_errors[thisfield+'_ERROR']
        try:
            rects.append(ax.bar(indices+j*field_bar_width,
                    mystructarray[thisfield],
                    field_bar_width,color = colors[j],
                    yerr = thiserror,#just to test
                    ecolor = 'k',
                    label = '$%s$'%thisfield))
        except Exception as e:
            raise RuntimeError(strm('Problem with bar graph: there are %d indices, but %d pieces of data'%(len(indices),
                len(mystructarray[thisfield])),
                'indices:', indices, 'data',mystructarray[thisfield])+explain_error(e))
    ax.set_xticks(indices+indiv_width/2)
    ax.set_xticklabels(labels)
    ax.legend([j[0] for j in rects],
            ['$%s$'%j for j in remaining_fields],loc = 'best')
    return
def lookup_rec(A,B,indexpair):
    r'''look up information about A in table B (i.e. chemical by index, etc)
    indexpair is either the name of the index
    or -- if it's differently named -- the pair of indices
    given in (A,B) respectively
    
    This will just drop any fields in B that are also in A,
    and the output uses the first indexname
    
    note that it it seems like the join_rec function above may be more efficient!!'''
    raise RuntimeError('You should now use decorate_rec!!')
    if type(indexpair) not in [tuple,list]:
        indexpair = (indexpair,indexpair)
    Bini = copy(B)
    B = recf.drop_fields(B,( set(B.dtype.names) & set(A.dtype.names) ) - {indexpair[1]}) # indexpair for B gets dropped later anyways
    joined = []
    for j in A:
        matchedrows =  B[B[indexpair[1]] == j[indexpair[0]]]
        for matchedrow in matchedrows:
            joined.append((j,matchedrow))
    if len(joined) == 0:
        raise IndexError(strm('Unable to find any matches between',
            A[indexpair[0]], 'and', B[indexpair[1]], '!'))
    whichisindex = joined[0][1].dtype.names.index(indexpair[1])
    allbutindex = lambda x: list(x)[0:whichisindex]+list(x)[whichisindex+1:]
    joined = concatenate([array(tuple(list(j[0])+allbutindex(j[1])),
                    dtype = dtype(j[0].dtype.descr+allbutindex(j[1].dtype.descr))).reshape(1) for j in joined])
    return joined
def reorder_rec(myarray,listofnames,first = True):
    try:
        indices_to_move = [myarray.dtype.names.index(j) for j in listofnames]
    except Exception as e:
        stuff_not_found = [j for j in listofnames if j not in myarray.dtype.names]
        if len(stuff_not_found) > 0:
            raise IndexError(strm(stuff_not_found,'is/are in the list you passed,',
            'but not one of the fields, which are', myarray.dtype.names))
        else:
            raise RuntimeError('unknown problem' + explain_error(e))
    old_type = list(myarray.dtype.descr)
    new_type = [old_type[j] for j in indices_to_move] + [old_type[j] for j in range(0,len(old_type)) if j not in indices_to_move]
    new_list_of_data = [myarray[j[0]] for j in new_type]
    return rec.fromarrays(new_list_of_data,dtype = new_type)
def lambda_rec(myarray,myname,myfunction,*varargs):
    r'''make a new field "myname" which consists of "myfunction" evaluated with the fields given by "myargs" as arguments
    the new field is always placed after the last argument name
    if myname is in myargs, the original row is popped'''
    if len(varargs) == 1:
        myargs = varargs[0]
    elif len(varargs) == 0:
        myargs = [myname]
    else:
        raise IndexError("For the fourth argument, you must pass either a list"
                " with the names of the arguments, or nothing (to use the field"
                " itself as an argument)")
    myargs = autostringconvert(myargs)
    if isinstance(myargs, str):
        myargs = (myargs,)
    if not isinstance(myargs, tuple):
        myargs = tuple(myargs)
    argdata = list(map((lambda x: myarray[x]),myargs))
    try:
        newrow = myfunction(*tuple(argdata))
    except TypeError:
        newrow = array([myfunction(*tuple([x[rownumber] for x in argdata])) for rownumber in range(0,len(argdata[0]))])
    if isinstance(newrow, list) and isinstance(newrow[0], str):
        newrow = array(newrow,dtype = '|S100')
    try:
        new_field_type = list(newrow.dtype.descr[0])
    except AttributeError as e:
        raise IndexError(strm("evaluated function on", argdata, "and got back",
            newrow, "which appears not to be a numpy array")+explain_error(e))
    new_field_type[0] = myname
    starting_names = myarray.dtype.names
    #{{{ make the dtype
    new_dtype = list(myarray.dtype.descr)
    #{{{ determine if I need to pop one of the existing rows due to a name conflict
    eliminate = None
    if myname in myargs:
        eliminate = myname
        insert_before = starting_names.index(myname) # if we are replacing, we want it in the same place
        new_dtype = [j for j in new_dtype if j[0] != eliminate]
    #}}}
    # if we haven't already eliminated, determine where to put it
    if eliminate is None:
        insert_before = starting_names.index(myargs[-1])+1
    # insert the new field where it goes
    new_dtype.insert(insert_before,tuple(new_field_type))
    #}}}
    #{{{ separate starting_names and ending_names
    if eliminate is None:
        ending_names = starting_names[insert_before:]
        starting_names = starting_names[:insert_before]
    else: # if I'm eliminating, I don't want to include the eliminated one
        ending_names = starting_names[insert_before+1:]
        starting_names = starting_names[:insert_before]
    #}}}
    return rec.fromarrays([myarray[x] for x in starting_names if x != eliminate]+[newrow]+[myarray[x] for x in ending_names if x != eliminate],dtype = new_dtype)
def join_rec(xxx_todo_changeme, xxx_todo_changeme1):
    (A,a_ind) = xxx_todo_changeme
    (B,b_ind) = xxx_todo_changeme1
    raise RuntimeError('You should now use decorate_rec!!')
def decorate_rec(xxx_todo_changeme2, xxx_todo_changeme3,drop_rows = False):
    r'''Decorate the rows in A with information in B --> if names overlap,
    keep the ones in A
    b_ind and a_ind can be either a single key, or a list of keys;
    if more than one element in B matches that in A, include both options!!'''
    (A,a_ind) = xxx_todo_changeme2
    (B,b_ind) = xxx_todo_changeme3
    dropped_rows = None
    # first find the list of indices that give us the data we want
    #{{{ process the arguments
    if (isinstance(b_ind, str)) and (isinstance(a_ind, str)):
        b_ind = [b_ind]
        a_ind = [a_ind]
    if ((isinstance(b_ind, list)) and (isinstance(a_ind, list))) and (len(b_ind) == len(a_ind)):
        pass
    else:
        raise ValueError('If you call a list for b_ind and/or a_ind, they must match in length!!!')
    if any([x not in B.dtype.names for x in b_ind]):
        problem_index = [x for x in b_ind if x not in B.dtype.names]
        raise ValueError(repr(problem_index)+' not in second argument, which has fields'+repr(B.dtype.names)+'!!!')
    if any([x not in A.dtype.names for x in a_ind]):
        problem_index = [x for x in a_ind if x not in A.dtype.names]
        raise ValueError(repr(problem_index)+' not in first argument, which has fields'+repr(A.dtype.names)+'!!!')
    #}}}
    B_reduced = B[b_ind] # a version of B reduced to only include the keys
    B_reduced = reorder_rec(B_reduced,b_ind)# again, because it doesn't do this just based on the indexing
    A_reduced = A[a_ind] # same for A
    A_reduced = reorder_rec(A_reduced,a_ind)# again, because it doesn't do this just based on the indexing
    # now, I need to generate a mapping from the b_ind to a_ind
    field_mapping = dict(list(zip(b_ind,a_ind)))
    # now I change the names so they match and I can compare them
    B_reduced.dtype.names = tuple([field_mapping[x] for x in B_reduced.dtype.names])
    #{{{ now find the list of indices for B that match each value of A
    old_B_reduced_names,old_B_reduced_types = tuple(zip(*tuple(B_reduced.dtype.descr)))
    B_reduced.dtype = dtype(list(zip(A_reduced.dtype.names,old_B_reduced_types)))
    if A_reduced.dtype != B_reduced.dtype:
        B_reduced.dtype = dtype(list(zip(old_B_reduced_names,old_B_reduced_types)))
        raise TypeError(strm('The datatype of A_reduced=', A_reduced.dtype,
            'and B_reduced=', B_reduced.dtype,
            'are not the same,  which is going to create problems!'))
    try:
        list_of_matching = [nonzero(B_reduced == j)[0] for j in A_reduced]
    except Exception as e:
        raise RuntimeError(strm('When trying to decorate,  A_reduced=', A_reduced,
            'with B_reduced=', B_reduced,
            'one or more of the following is an empty tuple,  which is wrong!:',
            [nonzero(B_reduced == j) for j in A_reduced])+explain_error(e))
    logger.debug(strm("(decorate\\_rec):: original list of matching",list_of_matching))
    length_of_matching = array([len(j) for j in list_of_matching])
    logger.debug(strm("(decorate\\_rec):: length of matching is",length_of_matching))
    if any(length_of_matching == 0):
        if drop_rows:
            if drop_rows == 'return':
                dropped_rows = A[length_of_matching == 0].copy()
            else:
                dropped_rows = A_reduced[length_of_matching == 0]
                print(r'{\color{red}Warning! decorate\_rec dropped fields in the first argument',lsafen(repr(list(zip(A_reduced.dtype.names * len(dropped_rows),dropped_rows.tolist())))),r'}')
            #{{{ now, remove all trace of the dropped fields
            A = A[length_of_matching != 0]
            list_of_matching = [j for j in list_of_matching if len(j)>0]
            length_of_matching = [len(j) for j in list_of_matching]
            #}}}
        else:
            raise ValueError(strm('There is no data in the second argument that has',
                b_ind,'fields to match the',a_ind,
                'fields of the first argument for the following records:',
                A_reduced[length_of_matching == 0],
                "if this is correct, you can set the drop_rows = True",
                "keyword argument to drop these fields"))
    # now, do a neat trick of stackoverflow to collapse a nested list
    # this gives just the indices in B that match the values of A
    list_of_matching = [j for i in list_of_matching for j in i]
    #}}}
    logger.debug(strm("(decorate\\_rec):: list of matching is",list_of_matching))
    # now grab the data for these rows
    add_data = B[list_of_matching]
    #{{{ finally, smoosh the two sets of data together
    #{{{ Now, I need to replicate the rows that have multiple matchesjk
    if any(length_of_matching > 1):
        index_replication_vector = [k for j in range(0,len(length_of_matching))
                for k in [j]*length_of_matching[j]]
        retval = A[index_replication_vector]
    else:
        retval = A.copy()
    #}}}
    #{{{ add the new fields
    new_dtypes = [j for j in B.dtype.descr if j[0] not in A.dtype.names]
    logger.debug(strm("(decorate\\_rec):: new dtypes:",repr(new_dtypes)))
    try:
        retval = newcol_rec(retval,new_dtypes)
    except Exception as e:
        raise ValueError(strm("Problem trying to add new columns with the dtypes",
            new_dtypes)+explain_error(e))
    #}}}
    logger.debug(strm("(decorate\\_rec):: add data:",repr(add_data)))
    for name in dtype(new_dtypes).names:
        logger.debug(strm("(decorate\\_rec):: trying to add data for",name,':',add_data[name][:]))
        retval[name][:] = add_data[name][:]
    #}}}
    if drop_rows == 'return':
        return retval,dropped_rows
    else:
        return retval
def newcol_rec(A,new_dtypes):
    r'''add new, empty (i.e. random numbers) fields to A, as given by new_dtypes
    --> note that there are deeply nested numpy functions to do this, but the options are confusing, and I think the way these work is efficient'''
    if isinstance(new_dtypes, dtype):
        new_dtypes = new_dtypes.descr
    elif isinstance(new_dtypes, tuple):
        new_dtypes = [new_dtypes]
    elif isinstance(new_dtypes, list):
        if not isinstance(new_dtypes[0], tuple):
            new_dtypes = [tuple(new_dtypes)]
    retval = empty(A.shape,dtype = A.dtype.descr + new_dtypes)
    for name in A.dtype.names:
        retval[name][:] = A[name][:]
    return retval
def applyto_rec(myfunc,myarray,mylist):
    r'apply myfunc to myarray with the intention of collapsing it to a smaller number of values'
    if not isinstance(mylist, list) and isinstance(mylist, str):
        mylist = [mylist]
    combined = []
    j = 0
    #{{{ make the list "combined", which I later concatenate
    while len(myarray) > 0:
        thisitem = myarray[0] # always grab the first row of what's left
        #{{{ initialize the empty new row
        if j == 0:
            newrow = thisitem.reshape(1)
        newrow = newrow.copy()
        #}}}
        #{{{ make a mask for all items that are identified as the same data
        # and copy the identical data to newrow
        mask = myarray[mylist[0]] == thisitem[mylist[0]]
        newrow[mylist[0]] = thisitem[mylist[0]]
        for k in range(1,len(mylist)):
            mask &= myarray[mylist[k]] == thisitem[mylist[k]]
            newrow[mylist[k]] = thisitem[mylist[k]]
        #}}}
        logger.debug(strm(lsafen('(applyto rec): for row %d, I select these:'%j)))
        myarray_subset = myarray[mask]
        logger.debug(strm(lsafen('(applyto rec): ',repr(myarray_subset))))
        other_fields = set(mylist)^set(thisitem.dtype.names)
        logger.debug(strm(lsafen('(applyto rec): other fields are:',other_fields)))
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = myfunc(myarray_subset[thisfield])
            except Exception as e:
                raise ValueError(strm("error in applyto_rec:  You usually get this",
                    "when one of the fields that you have NOT passed in the",
                    "second argument is a string.  The fields and types",
                    "are:",repr(myarray_subset.dtype.descr)) + explain_error(e))
        logger.debug(strm(lsafen("(applyto rec): for row %d, I get this as a result:"%j,newrow)))
        combined.append(newrow) # add this row to the list
        myarray = myarray[~mask] # mask out everything I have used from the original matrix
        logger.debug(strm(lsafen("(applyto rec): the array is now",repr(myarray))))
        j += 1
    #}}}
    combined = concatenate(combined)
    logger.debug(strm(lsafen("(applyto rec): final result",repr(combined),"has length",len(combined))))
    return combined
def meanstd_rec(myarray,mylist,standard_error = False):
    r'this is something like applyto_rec, except that it applies the mean and creates new rows for the "error," where it puts the standard deviation'
    if not isinstance(mylist, list) and isinstance(mylist, str):
        mylist = [mylist]
    combined = []
    other_fields = set(mylist)^set(myarray.dtype.names)
    logger.debug(strm('(meanstd_rec): other fields are',lsafen(other_fields)))
    newrow_dtype = [[j,('%s_ERROR'%j[0],)+j[1:]] if j[0] in other_fields else [j] for j in myarray.dtype.descr]
    newrow_dtype = [k for j in newrow_dtype for k in j]
    logger.debug(strm(lsafen('(meanstd rec): other fields are:',other_fields)))
    #{{{ make the list "combined", which I later concatenate
    j = 0
    while len(myarray) > 0:
        thisitem = myarray[0] # always grab the first row of what's left
        #{{{ initialize the empty new row
        newrow = zeros(1,dtype = newrow_dtype)
        #}}}
        #{{{ make a mask for all items that are identified as the same data
        # and copy the identical data to newrow
        mask = myarray[mylist[0]] == thisitem[mylist[0]]
        newrow[mylist[0]] = thisitem[mylist[0]]
        for k in range(1,len(mylist)):
            mask &= myarray[mylist[k]] == thisitem[mylist[k]]
            newrow[mylist[k]] = thisitem[mylist[k]]
        #}}}
        logger.debug(strm(lsafen('(meanstd rec): for row %d, I select these:'%j)))
        myarray_subset = myarray[mask]
        logger.debug(strm(lsafen('(meanstd rec): ',repr(myarray_subset))))
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = mean(myarray_subset[thisfield])
                if standard_error:
                    newrow[thisfield+"_ERROR"] = std(myarray_subset[thisfield])/sqrt(len(myarray_subset[thisfield]))
                else:
                    newrow[thisfield+"_ERROR"] = std(myarray_subset[thisfield])
            except:
                raise RuntimeError("error in meanstd_rec:  You usually get this",
                        "when one of the fields that you have NOT passed in the",
                        "second argument is a string.  The fields and types",
                        "are:",repr(myarray_subset.dtype.descr))
            #print 'for field',lsafe(thisfield),'I find',lsafen(newrow[thisfield])
        logger.debug(strm(lsafen("(meanstd rec): for row %d, I get this as a result:"%j,newrow)))
        combined.append(newrow) # add this row to the list
        myarray = myarray[~mask] # mask out everything I have used from the original matrix
        logger.debug(strm(lsafen("(meanstd rec): the array is now",repr(myarray))))
        j += 1
    #}}}
    combined = concatenate(combined)
    logger.debug(strm(lsafen("(meanstd rec): final result",repr(combined),"has length",len(combined))))
    return combined
def make_rec(*args,**kwargs):
    r'input,names or a single argument, which is a dictionary\nstrlen = 100 gives length of the strings (which need to be specified in record arrays)\nyou can also specify (especially useful with the dictionary format) the list order = [str1,str2,...] which orders the output records with the field containing str1 first, then the field containing str2, then any remaining fields'
    strlen,order,zeros_like = process_kwargs([('strlen',100),
        ('order',None),
        ('zeros_like',False)],kwargs)
    if len(args) == 1 and (isinstance(args[0], dict)):
        names = list(args[0].keys())
        input = list(args[0].values())
    elif len(args) == 2:
        input = args[0]
        names = args[1]
    else:
        raise ValueError(strm("I don't understand the arguments you passed to",
                "make_rec!!!\nshould be (list of values, list of field names),",
                "or a dictionary"))
    #{{{ apply the order kwarg
    if order is not None:
        newindices = []
        for orderitem in order:
            newindices += [j for j,k in enumerate(names) if (k.find(orderitem)>-1 and j not in newindices)]
        newindices += [j for j,k in enumerate(names) if j not in newindices]
        names = [names[j] for j in newindices]
        input = [input[j] for j in newindices]
    #}}}
    if not (isinstance(input, list) and isinstance(names, list)):
        raise TypeError('you must enter a list for both')
    types = list(map(type,input))
    shapes = list(map(shape,input))
    if all([j == shapes[0] for j in shapes]):
        if shapes[0] == ():# if it's one dimensional
            equal_shapes = False
            shapes = [(1)]*len(shapes)
        else:
            equal_shapes = True
            shape_of_array = shapes[0]
            shapes = [()]*len(shapes)
    else:
        equal_shapes = False
    for j,k in enumerate(input):
        if isinstance(k, list) and equal_shapes:
            k = k[0]
        if isinstance(k, str):
            types[j] = '|S%d'%strlen
        if isinstance(k, ndarray):
            types[j] = k.dtype
    try:
        mydtype = dtype(list(zip(names,types,shapes)))
    except Exception as e:
        raise ValueError(strm('problem trying to make names',names,' types',
            types,'shapes',shapes)+explain_error(e))
    if zeros_like:
        retval = zeros(zeros_like,dtype = mydtype)
        return retval
    if equal_shapes:
        retval = empty(shape_of_array,dtype = mydtype)
        for j,thisname in enumerate(names):
            try:
                retval[thisname][:] = input[j][:]
            except Exception as e:
                raise RuntimeError("error trying to load input for '"+thisname
                        +"' of shape "+repr(shape(input[j]))+" into retval field of shape "
                        +repr(shape(retval[thisname])))
        return retval
    else:
        try:
            return array([tuple(input)],dtype = mydtype)
        except:
            raise ValueError(strm('problem trying to assign data of type',list(map(type,input)),
                '\nvalues',input,'\nonto',mydtype,'\ndtype made from tuple:',
                list(zip(names,types,shapes))))
#{{{ convert back and forth between lists, etc, and ndarray
def make_ndarray(array_to_conv,name_forprint = 'unknown'): 
    if type(array_to_conv) in [int,int32,double,float,complex,complex128,float,bool,bool_]: # if it's a scalar
        pass
    elif isinstance(array_to_conv, str):
        pass
    elif type(array_to_conv) in [list,ndarray] and len(array_to_conv) > 0:
        array_to_conv = rec.fromarrays([array_to_conv],names = 'LISTELEMENTS') #list(rec.fromarrays([b])['f0']) to convert back
    elif type(array_to_conv) in [list,ndarray] and len(array_to_conv) == 0:
        array_to_conv = None
    elif array_to_conv is  None:
        pass
    else:
        raise TypeError(strm('type of value (',type(array_to_conv),') for attribute name',
            name_forprint,'passed to make_ndarray is not currently supported'))
    return array_to_conv
def unmake_ndarray(array_to_conv,name_forprint = 'unknown'): 
    r'Convert this item to an ndarray'
    if (isinstance(array_to_conv, recarray)) or (isinstance(array_to_conv, ndarray) and array_to_conv.dtype.names is not None and len(array_to_conv.dtype.names)>0):
        #{{{ if it's a record/structured array, it should be either a list or dictionary
        if 'LISTELEMENTS' in array_to_conv.dtype.names:
            if array_to_conv.dtype.names == tuple(['LISTELEMENTS']):
                retval = list(array_to_conv['LISTELEMENTS'])
            else:
                raise ValueError(strm('Attribute',name_forprint,
                    'is a recordarray with a LISTELEMENTS field, but it',
                    'also has other dimensions:',array_to_conv.dtype.names,
                    'not',tuple(['LISTELEMENTS'])))
        elif len(array_to_conv)==1:
            thisval = dict(list(zip(a.dtype.names,a.tolist()[0])))
        else: raise ValueError('You passed a structured array, but it has more',
                "than one dimension, which is not yet supported\nLater, this",
                'should be supported by returning a dictionary of arrays')
        #}}}
    elif isinstance(array_to_conv, ndarray) and len(array_to_conv)==1:
        #{{{ if it's a length 1 ndarray, then return the element
        retval = array_to_conv.tolist()
        logger.debug(strm(name_forprint,"=",type(array_to_conv),"is a numpy array of length one"))
        #}}}
    elif type(array_to_conv) in [string_,int32,float64,bool_]:
        #{{{ map numpy strings onto normal strings
        retval = array_to_conv.tolist()
        logger.debug(strm("name_forprint","=",type(array_to_conv),"is a numpy scalar"))
        #}}}
    elif isinstance(array_to_conv, list):
        #{{{ deal with lists
        logger.debug(strm(name_forprint,"is a list"))
        typeofall = list(map(type,array_to_conv))
        if all([x is string_ for x in typeofall]):
            logger.debug(strm(name_forprint,"=",typeofall,"are all numpy strings"))
            retval = list(map(str,array_to_conv))
        else:
            logger.debug(strm(name_forprint,"=",typeofall,"are not all numpy string"))
            retval = array_to_conv
        #}}}
    else:
        logger.debug(strm(name_forprint,"=",type(array_to_conv),"is not a numpy string or record array"))
        retval = array_to_conv
    return retval
#}}}
#}}}
def emptytest(x): # test is it is one of various forms of empty
   if type(x) in [list,array]:
       if len(x) == 0:
           return True
       elif x is array(None):
           return True
       elif len(x) > 0:
           return False
       #don't want the following, because then I may need to pop, etc
       #if type(x) is list and all(map(lambda x: x is None,x)): return True
   if size(x) == 1 and x is None: return True
   if size(x) == 0: return True
   return False
def lsafen(*string,**kwargs):
    "see lsafe, but with an added double newline"
    string = list(string)
    string += ['\n\n']
    return lsafe(*tuple(string),**kwargs)
def lsafe(*string,**kwargs):
    "Output properly escaped for latex"
    if len(string) > 1:
        lsafewkargs = lambda x: lsafe(x,**kwargs)
        return ' '.join(list(map(lsafewkargs,string)))
    else:
        string = string[0]
    #{{{ kwargs
    spaces = False
    if 'spaces' in list(kwargs.keys()):
        spaces = kwargs.pop('spaces')
    if 'wrap' in list(kwargs.keys()):
        wrap = kwargs.pop('wrap')
    else:
        wrap = None
    #}}}
    if not isinstance(string, str):
        string = str(string)
    if wrap is True:
        wrap = 60
    if wrap is not None:
        string = '\n'.join(textwrap.wrap(string,wrap))
    string = string.replace('\\','\\textbackslash ')
    if spaces:
        string = string.replace(' ','\\ ')
    string = string.replace('\n\t','\n\n\\quad ')
    string = string.replace('\t','\\quad ')
    string = string.replace('_',r'\_')
    string = string.replace('{',r'\{')
    string = string.replace('}',r'\}')
    string = string.replace('$$',r'ACTUALDOUBLEDOLLAR')
    string = string.replace(']',r'$]$')
    string = string.replace('[',r'$[$')
    string = string.replace('<',r'$<$')
    string = string.replace('>',r'$>$')
    string = string.replace('$$',r'')
    string = string.replace('ACTUALDOUBLEDOLLAR',r'$$')
    string = string.replace('^',r'\^')
    string = string.replace('#',r'\#')
    string = string.replace('%',r'\%')
    string = string.replace('&',r'\&')
    string = string.replace('+/-',r'\ensuremath{\pm}')
    string = string.replace('|',r'$|$')
    return string
#{{{ errors
def explain_error(e):
    '''Allows you to wrap existing errors with more explanation

    For example:

    >    except BaseException  as e:
    >        raise IndexError("I can't find the node "+pathstring+explain_error(e))
    >                + '\n'.join(['>\t'+j for j in str(e).split('\n')]))# this indents
    '''
    exc_type,exc_obj,exc_tb = exc_info()
    #code_loc = strm(os.path.relpath(exc_tb.tb_frame.f_code.co_filename,os.getcwd()), 'line', exc_tb.tb_lineno)
    code_loc = strm(exc_tb.tb_frame.f_code.co_filename, 'line', exc_tb.tb_lineno)
    return ('\n> Original error (%s -- %s):\n'%(exc_type,code_loc)
            + '\n'.join(['>\t'+j
                for j in str(e).split('\n')]))# this indents
class CustomError(Exception):
    def __init__(self, *value, **kwargs):
        raise NotImplementedError("You should get rid of CustomError and use explain_error instead")
        return
def copy_maybe_none(input):
    if input is None:
        return None
    else:
        if isinstance(input, list):
            return list(map(copy,input))
        else:
            return input.copy()
def maprep(*mylist):
    mylist = list(mylist)
    for j in range(0,len(mylist)):
        if not isinstance(mylist[j], str):
            mylist[j] = mylist[j].__repr__()
    return ' '.join(mylist)
#}}}
#{{{ HDF5 functions
#{{{ helper function for HDF5 search
def gensearch(labelname,format = '%0.3f',value = None,precision = None):
    'obsolete -- use h5gensearch'
    if value is None:
        raise ValueError('You must pass a value to gensearch')
    if precision is None:
        precision = value*0.01 # the precision is 1% of the value, if we don't give an argument
    searchstring_high = '(%s < %s + (%s))'%tuple([labelname]+[format]*2)
    #print "\n\nDEBUG check format:\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_high = searchstring_high%(value,precision)
    #print "\n\nDEBUG after substitution with",value,precision,":\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_low = '(%s > %s - (%s))'%tuple([labelname]+[format]*2)
    searchstring_low = searchstring_low%(value,precision)
    return searchstring_low + ' & ' + searchstring_high
def h5searchstring(*args,**kwargs):
    '''generate robust search strings
    :parameter fieldname,value:
    search AROUND a certain value (overcomes some type conversion issues) optional arguments are the format specifier and the fractional precision:
    **OR**
    :parameter field_and_value_dictionary:
    generate a search string that matches one or more criteria'''
    format,precision = process_kwargs([('format','%g'),
        ('precision',0.01)],
        kwargs)
    if len(args) == 2:
        fieldname,value = args
    elif len(args) == 1 and isinstance(args[0], dict):
        dict_arg = args[0]
        condlist = []
        for k,v in dict_arg.items():
            condlist.append(h5searchstring(k,v,format = format,precision = precision))
        return ' & '.join(condlist)
    else:
        raise ValueError("pass either field,value pair or a dictionary!")
    if isinstance(value, str):
        raise ValueError("string matching in pytables is broken -- search by hand and then use the index")
    precision *= value
    searchstring_high = '(%s < %s + (%s))'%tuple([fieldname]+[format]*2)
    #print "\n\nDEBUG check format:\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_high = searchstring_high%(value,precision)
    #print "\n\nDEBUG after substitution with",value,precision,":\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_low = '(%s > %s - (%s))'%tuple([fieldname]+[format]*2)
    searchstring_low = searchstring_low%(value,precision)
    return '(' + searchstring_low + ' & ' + searchstring_high + ')'
#}}}
def h5loaddict(thisnode):
    #{{{ load all attributes of the node
    retval = dict([(x,thisnode._v_attrs.__getattribute__(x))
        for x in thisnode._v_attrs._f_list('user')])
    #}}}
    for k,v in retval.items():#{{{ search for record arrays that represent normal lists
        retval[k]  = unmake_ndarray(v,name_forprint = k)
    if isinstance(thisnode, tables.table.Table):#{{{ load any table data
        logger.debug(strm("It's a table\n\n"))
        if 'data' in list(retval.keys()):
            raise AttributeError('There\'s an attribute called data --> this should not happen!')
        retval.update({'data':thisnode.read()})
    elif isinstance(thisnode, tables.group.Group):
        #{{{ load any sub-nodes as dictionaries
        mychildren = thisnode._v_children
        for thischild in list(mychildren.keys()):
            if thischild in list(retval.keys()):
                raise AttributeError('There\'s an attribute called ',thischild,' and also a sub-node called the',thischild,'--> this should not happen!')
            retval.update({thischild:h5loaddict(mychildren[thischild])})
        #}}}
    else:
        raise AttributeError(strm("I don't know what to do with this node:",thisnode))
    #}}}
    return retval
def h5child(thisnode,childname,clear = False,create = None):
    r'''grab the child, optionally clearing it and/or (by default) creating it'''
    #{{{ I can't create and clear at the same time
    if create and clear:
        raise ValueError("You can't call clear and create at the same time!\nJust call h5child twice, once with clear, once with create")
    if create is None:
        if clear == True:
            create = False
        else:
            create = True
    #}}}
    h5file = thisnode._v_file
    try:
        childnode = h5file.get_node(thisnode,childname)
        logger.debug(strm('found',childname))
        if clear:
            childnode._f_remove(recursive = True)
            childnode = None
    except tables.NoSuchNodeError as e:
        if create is False and not clear:
            raise RuntimeError('Trying to grab a node that does not exist with create = False'+explain_error(e))
        elif clear:
            childnode = None
        else:
            childnode = h5file.create_group(thisnode,childname)
            logger.debug('created',childname)
    return childnode
def h5remrows(bottomnode,tablename,searchstring):
    if isinstance(searchstring, dict):
        searchstring = h5searchstring(searchstring)
    try:
        thistable = bottomnode.__getattr__(tablename)
        counter = 0
        try:
            data = thistable.read_where(searchstring).copy()
        except Exception as e:
            raise RuntimeError(strm('Problem trying to remove rows using search string',
                searchstring, 'in', thistable, explain_error(e)))
        for row in thistable.where(searchstring):
            if len(thistable) == 1:
                thistable.remove()
                counter += 1
            else:
                try:
                    thistable.remove_rows(row.nrow - counter,row.nrow - counter + 1) # counter accounts for rows I have already removed.
                except:
                    print("you passed searchstring",searchstring)
                    print("trying to remove row",row)
                    print("trying to remove row with number",row.nrow)
                    print(help(thistable.remove_rows))
                    raise RuntimeError("length of thistable is "+repr(len(thistable))+" calling remove_rows with "+repr(row.nrow-counter))
                counter += 1
        return counter,data
    except tables.NoSuchNodeError:
        return False,None
def h5addrow(bottomnode,tablename,*args,**kwargs):
    '''add a row to a table, creating it if necessary, but don\'t add if the data matches the search condition indicated by `match_row`
    `match_row` can be either text or a dictionary -- in the latter case it's passed to h5searchstring
    '''
    match_row,only_last = process_kwargs([('match_row',None),('only_last',True)],kwargs)
    try: # see if the table exists
        mytable = h5table(bottomnode,tablename,None)
        tableexists = True
    except RuntimeError: # if table doesn't exist, create it
        newindex = 1
        tableexists = False
    if tableexists:
        #{{{ auto-increment "index"
        newindex = mytable.read()['index'].max() + 1
        #}}}
        # here is where I would search for the existing data
        if match_row is not None:
            if isinstance(match_row, dict):
                match_row = h5searchstring(match_row)
            logger.debug("trying to match row according to",lsafen(match_row))
            mytable.flush()
            try:
                matches = mytable.read_where(match_row)
            except NameError as e:
                raise NameError(' '.join(map(str,
                    [e,'\nYou passed',match_row,'\nThe columns available are',mytable.colnames,"condvars are",condvars])))
            except ValueError as e:
                raise NameError(' '.join(map(str,
                    [e,'\nYou passed',match_row,'\nThe columns available are',mytable.colnames])))
            if len(matches) > 0:
                if only_last:
                    logger.debug(strm(r'\o{',lsafen(len(matches),"rows match your search criterion, returning the last row"),'}'))
                    return mytable,matches['index'][-1]
                else:
                    return mytable,matches['index'][:]
            else:
                logger.debug("I found no matches")
    if len(args) == 1 and (isinstance(args[0], dict)):
        listofnames,listofdata = list(map(list,list(zip(*tuple(args[0].items())))))
    elif len(args) == 2 and isinstance(args[0], list) and isinstance(args[1], list):
        listofdata = args[0]
        listofnames = args[1]
    else:
        raise TypeError('h5addrow takes either a dictionary for the third argument or a list for the third and fourth arguments')
    try:
        listofdata = [newindex] + listofdata
    except:
        raise TypeError('newindex is'+repr(newindex)+'listofdata is'+repr(listofdata))
    listofnames = ['index'] + listofnames
    myrowdata = make_rec(listofdata,listofnames)
    if tableexists:
        try:
            mytable.append(myrowdata)
        except ValueError as e:
            print(lsafen("I'm about to flag an error, but it looks like there was an issue appending",myrowdata))
            tabledforerr = mytable.read()
            raise AttributeError(strm('Compare names and values table data vs. the row you are trying to add\n',
                '\n'.join(map(repr,list(zip(list(mytable.read().dtype.fields.keys()),
                list(mytable.read().dtype.fields.values()),
                list(myrowdata.dtype.fields.keys()),
                list(myrowdata.dtype.fields.values())))))),explain_error(e))
        mytable.flush()
    else:
        recorddata = myrowdata
        try:
            mytable = h5table(bottomnode,
                    tablename,
                    recorddata)
        except Exception as e:
            raise RuntimeError(strm('Error trying to write record array:',
                repr(recorddata),'from listofdata',listofdata,'and names',listofnames,
                explain_error(e)))
        mytable.flush()
    return mytable,newindex
def h5table(bottomnode,tablename,tabledata):
    'create the table, or if tabledata is None, just check if it exists'
    #{{{ save but don't overwrite the table
    h5file = bottomnode._v_file
    if tablename not in list(bottomnode._v_children.keys()):
        if tabledata is not None:
            if isinstance(tabledata, dict):
                tabledata = make_rec(tabledata)
            datatable = h5file.create_table(bottomnode,tablename,tabledata) # actually write the data to the table
        else:
            raise RuntimeError(' '.join(map(str,['You passed no data, so I can\'t create table',tablename,'but it doesn\'t exist in',bottomnode,'which has children',list(bottomnode._v_children.keys())])))
    else:
        if tabledata is not None:
            raise ValueError(strm('You\'re passing data to create the table,',tablename,' but the table already exists!'))
        else:
            pass
    return bottomnode._v_children[tablename]
    #}}}
def h5nodebypath(h5path,force = False,only_lowest = False,check_only = False,directory='.'):
    r'''return the node based on an absolute path, including the filename'''
    logger.debug(strm("DEBUG: called h5nodebypath on",h5path))
    h5path = h5path.split('/')
    #{{{ open the file / check if it exists
    logger.debug(strm(lsafen('h5path=',h5path)))
    logger.debug(strm('the h5path is',h5path))
    if h5path[0] in listdir(directory):
        logger.debug(strm('DEBUG: file exists\n\n'))
        log_fname('data_files',h5path[0],directory)
    else:
        if check_only:
            errmsg = log_fname('missing_data_files',h5path[0],directory,err=True)
            raise AttributeError("You're checking for a node in a file (%s) that does not exist"%(h5path[0])
                    +'\n'
                    +errmsg)
        logger.debug(strm('DEBUG: file does not exist\n\n'))
    mode = 'a'
    #if check_only: mode = 'r'
    logger.debug(strm('so I look for the file',h5path[0],'in directory',directory))
    try:
        h5file = tables.open_file(os.path.join(directory,h5path[0]),mode = mode,title = 'test file')
    except IOError as e:
        raise IOError('I think the HDF5 file has not been created yet, and there is a bug pytables that makes it freak out, but you can just run again.'+explain_error(e))
    #}}}
    currentnode = h5file.get_node('/') # open the root node
    logger.debug(strm("I have grabbe node",currentnode,"of file",h5file,'ready to step down search path'))
    for pathlevel in range(1,len(h5path)):#{{{ step down the path
            clear = False
            create = True
            if only_lowest or check_only:
                create = False
            if pathlevel == len(h5path)-1: # the lowest level
                if only_lowest:
                    create = True
                if force:
                    clear = True
            safetoleaveopen = False
            try:
                currentnode = h5child(currentnode, # current node
                        h5path[pathlevel], # the child
                        create = create,
                        clear = clear)
                logger.debug(strm(lsafen("searching for node path: descended to node",currentnode)))
                logger.debug(strm("searching for node path: descended to node",currentnode))
            except BaseException as e:
                logger.debug(strm(lsafen("searching for node path: got caught searching for node",h5path[pathlevel])))
                logger.debug(strm("searching for node path: got caught searching for node",h5path[pathlevel]))
                h5file.close()
                #print lsafen("DEBUG: Yes, I closed the file")
                raise IndexError(strm('Problem trying to load node ',h5path,explain_error(e)))
            #}}}
    return h5file,currentnode
def h5attachattributes(node,listofattributes,myvalues):
    listofattributes = [j for j in listofattributes # need to exclude the properties
            if j not in ['angle','real','imag']]
    if node is None:
        raise IndexError('Problem!, node passed to h5attachattributes: ',node,'is None!')
    h5file = node._v_file
    if isinstance(myvalues,nddata):
        attributevalues = [myvalues.__getattribute__(x) for x in listofattributes]
    elif isinstance(myvalues, list):
        attributevalues = myvalues
    else:
        raise TypeError("I don't understand the type of myvalues, which much be a list or a nddata object, from which the attribute values are retrieved")
    listout = list(listofattributes)
    for j,thisattr in enumerate(listofattributes):
        thisval = attributevalues[j]
        if type(thisval) in [dict]:
            dictnode = h5child(node,
                    thisattr,
                    clear = True)
            dictnode = h5child(node,
                    thisattr,
                    create = True)
            h5attachattributes(dictnode,
                    list(thisval.keys()),
                    list(thisval.values()))
            thisval = None
            listout.remove(thisattr)
        else:
            # {{{ pytables hates <U24 which is created from unicode
            if type(thisval) in [list,tuple]:
                if any([isinstance(x,str) for x in thisval]):
                    logger.debug(strm("going to convert",thisval,"to strings"))
                    thisval = [str(x) if isinstance(x,str) else x for x in thisval]
                    logger.debug(strm("now it looks like this:",thisval))
            thisval = make_ndarray(thisval,name_forprint = thisattr)
            # }}}
        if thisval is not None:
            try:
                node._v_attrs.__setattr__(thisattr,thisval)
            except Exception as e:
                raise RuntimeError("PyTables freaks out when trying to attach attribute "+repr(thisattr)+" with value "+repr(thisval)+"\nOriginal error was:\n"+str(e))
            listout.remove(thisattr)
    listofattributes[:] = listout # pointer
def h5inlist(columnname,mylist):
    'returns rows where the column named columnname is in the value of mylist'
    if isinstance(mylist, slice):
        if mylist.start is not None and mylist.stop is not None:
            return "(%s >= %g) & (%s < %g)"%(columnname,mylist.start,columnname,mylist.stop)
        elif mylist.stop is not None:
            return "(%s < %g)"%(columnname,mylist.stop)
        elif mylist.start is not None:
            return "(%s > %g)"%(columnname,mylist.start)
        else:
            raise ValueError()
    if isinstance(mylist, ndarray):
        mylist = mylist.tolist()
    if not isinstance(mylist, list):
        raise TypeError("the second argument to h5inlist must be a list!!!")
    if any([type(x) in [double,float64] for x in mylist]):
        if all([type(x) in [double,float64,int,int32,int64] for x in mylist]):
            return '('+'|'.join(["(%s == %g)"%(columnname,x) for x in mylist])+')'
    elif all([type(x) in [int,int,int32,int64] for x in mylist]):
        return '('+'|'.join(["(%s == %g)"%(columnname,x) for x in mylist])+')'
    elif all([isinstance(x, str) for x in mylist]):
        return '('+'|'.join(["(%s == '%s')"%(columnname,x) for x in mylist])+')'
    else:
        raise TypeError("I can't figure out what to do with this list --> I know what to do with a list of numbers or a list of strings, but not a list of type"+repr(list(map(type,mylist))))
def h5join(firsttuple,secondtuple,
    additional_search = '',
    select_fields = None,
    pop_fields = None):
    #{{{ process the first argument as the hdf5 table and indices, and process the second one as the structured array to join onto
    if not ((isinstance(firsttuple, tuple)) and (isinstance(secondtuple, tuple))):
        raise ValueError('both the first and second arguments must be tuples!')
    if not ((len(firsttuple) == 2) and (len(secondtuple) == 2)):
        raise ValueError('The length of the first and second arguments must be two!')
    tablenode = firsttuple[0]
    tableindices = firsttuple[1]
    logger.debug(strm('h5join tableindices looks like this:',tableindices))
    if not isinstance(tableindices, list):
        tableindices = [tableindices]
    logger.debug(strm('h5join tableindices looks like this:',tableindices))
    mystructarray = secondtuple[0].copy()
    mystructarrayindices = secondtuple[1]
    if not isinstance(mystructarrayindices, list):
        mystructarrayindices = [mystructarrayindices]
    #}}}
    #{{{ generate a search string to match potentially more than one key
    search_string = []
    if len(tableindices) != len(mystructarrayindices):
        raise ValueError('You must pass either a string or a list for the second element of each tuple!\nIf you pass a list, they must be of the same length, since the field names need to line up!')
    # this can't use h5inlist, because the and needs to be on the inside
    #{{{ this loop creates a list of lists, where the inner lists are a set of conditions that need to be satisfied
    # this is actually not causing  any trouble right now, but needs to be fixed, because of the way that it's doing the type conversion
    for thistableindex,thisstructarrayindex in zip(tableindices,mystructarrayindices):
        if thisstructarrayindex not in mystructarray.dtype.names:
            raise ValueError(repr(thisstructarrayindex)+" is not in "+repr(mystructarray.dtype.names))
        if isinstance(mystructarray[thisstructarrayindex][0],str):
            search_string.append(["(%s == '%s')"%(thistableindex,x) for x in mystructarray[thisstructarrayindex]])
        elif type(mystructarray[thisstructarrayindex][0]) in [int,double,float,float64,float32,int32,int64]:
            search_string.append(["(%s == %s)"%(thistableindex,str(x)) for x in mystructarray[thisstructarrayindex]])
            #print 'a g mapping for',[x for x in mystructarray[thisstructarrayindex]],'gives',search_string[-1],'\n\n'
        else:
            raise TypeError("I don't know what to do with a structured array that has a row of type"+repr(type(mystructarray[thisstructarrayindex][0])))
    #}}}
    search_string = [' & '.join(x) for x in zip(*tuple(search_string))] # this "and"s together the inner lists, since all conditions must be matched
    search_string = '('+'|'.join(search_string)+')' # then, it "or"s the outer lists, since I want to collect data for all rows of the table
    #}}}
    if len(additional_search) > 0:
        additional_search = " & (%s)"%additional_search
        search_string = search_string + additional_search
    logger.debug(strm('\n\nh5join generated the search string:',lsafen(search_string)))
    retval = tablenode.read_where(search_string)
    #{{{ then join the data together
    # here I'm debugging the join function, again, and again, and agin
    try:
        retval = decorate_rec((retval,tableindices),(mystructarray,mystructarrayindices)) # this must be the problem, since the above looks fine
    except Exception as e:
        raise Exception(strm('Some problems trying to decorate the table',
            retval, 'of dtype', retval.dtype, 'with the structured array',
            mystructarray, 'of dtype', mystructarray.dtype, explain_error(e)))
    if pop_fields is not None:
        if select_fields is not None:
            raise ValueError("It doesn't make sense to specify pop_fields and select_fields at the same time!!")
        select_fields = list(set(retval.dtype.names) ^ set(pop_fields))
    if select_fields is not None:
        logger.debug(strm('\n\nh5join original indices',lsafen(retval.dtype.names)))
        try:
            retval = retval[select_fields]
        except ValueError as e:
            raise ValueError(strm('One of the fields', select_fields, 'is not in',
                retval.dtype.names, explain_error(e)))
    #}}}
    return retval
#}}}
#{{{ indices to slice
#}}}
#{{{ old grid and tick
def gridandtick(ax,rotation=(0,0),precision=(2,2),
        labelstring=('',''),gridcolor=r_[0,0,0],
        formatonly = False,fixed_y_locator = None,
        logarithmic = False,use_grid = True,
        spines = None,y = True):
    #{{{ taken from matplotlib examples
    def adjust_spines(ax,spines):
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        for loc, spine in list(ax.spines.items()):
            if loc in spines:
                spine.set_position(('outward',5)) # outward by 5 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine
        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([],minor = False)
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([],minor = False)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    #}}}
    if spines is not None:
        adjust_spines(ax,spines = spines)
    if not formatonly:
        #{{{x ticks
        # determine the size
        width = abs(diff(ax.get_xlim()))
        if width==0:
            raise ValueError('x axis width is zero')
        widthexp = floor(log(width)/log(10.))-1
        scalefactor = 10**widthexp
        width /= scalefactor
        majorLocator = MultipleLocator(5*scalefactor)
        #majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[0]+'f'+labelstring[0])# labelstring can be used, for instance, for pi
        #ax.xaxis.set_major_formatter(majorFormatter)
        minorLocator   = MultipleLocator(1*scalefactor)
        ax.xaxis.set_major_locator(majorLocator)
        #for the minor ticks, use no labels; default NullFormatter
        ax.xaxis.set_minor_locator(minorLocator)
        #}}}
        if y:
            #{{{ y ticks
            width = abs(diff(ax.get_ylim()))
            if width==0:
                raise ValueError('y axis width is zero')
            widthexp = floor(log(width)/log(10.))-1
            scalefactor = 10**widthexp
            width /= scalefactor
            if fixed_y_locator is None:
                if logarithmic:
                    majorLocator = LogLocator(10)
                else:
                    majorLocator   = MultipleLocator(5*scalefactor)
            else:
                majorLocator   = MultipleLocator(fixed_y_locator[4::5])
            #majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[1]+'f'+labelstring[1])# labelstring can be used, for instance, for pi
            #ax.yaxis.set_major_formatter(majorFormatter)
            if fixed_y_locator is None:
                if logarithmic:
                    minorLocator = LogLocator(10,subs=r_[0:11])
                else:
                    minorLocator   = MultipleLocator(1*scalefactor)
            else:
                minorLocator   = FixedLocator(fixed_y_locator)
            ax.yaxis.set_major_locator(majorLocator)
            #for the minor ticks, use no labels; default NullFormatter
            ax.yaxis.set_minor_locator(minorLocator)
            #}}}
    ax.yaxis.grid(use_grid,which='major',color=gridcolor,alpha=0.15,linestyle='-')
    ax.xaxis.grid(use_grid,which='major',color=gridcolor,alpha=0.15,linestyle='-')
    ax.yaxis.grid(use_grid,which='minor',color=gridcolor,alpha=0.075,linestyle='-')
    ax.xaxis.grid(use_grid,which='minor',color=gridcolor,alpha=0.075,linestyle='-')
    labels = ax.get_xticklabels()
    setp(labels,rotation=rotation[0],fontsize=10)
    if y:
        labels = ax.get_yticklabels()
        setp(labels,rotation=rotation[1],fontsize=10)
    fig = gcf()
    fig.autofmt_xdate()
    return
def gridon(gridcolor=r_[0,0,0]):
    grid(True,which='major',color=gridcolor,alpha=0.1,linestyle='-')
    grid(True,which='minor',color=gridcolor,alpha=0.05,linestyle='-')
#}}}
#{{{ a better version?
def othergridandtick(ax,rotation=(0,0),precision=(2,2),labelstring=('',''),gridcolor=r_[0,0,0],y = True,x = True,spines = None):
    #{{{ taken from matplotlib examples
    def adjust_spines(ax,spines):
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        for loc, spine in list(ax.spines.items()):
            if loc in spines:
                spine.set_position(('outward',5)) # outward by 5 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine
        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([],minor = False)
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([],minor = False)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    #}}}
    if spines is not None:
        adjust_spines(gca(),spines = spines)
    if x:
        #{{{x ticks
        # determine the size
        ax.xaxis.set_major_locator(MaxNLocator(10)) # could use multiplelocator if it keeps try to do multiples of 2
        ax.xaxis.set_minor_locator(MaxNLocator(50))
        #}}}
    if y:
        #{{{ y ticks
        ax.yaxis.set_major_locator(MaxNLocator(10))
        ax.yaxis.set_minor_locator(MaxNLocator(50))
        #}}}
    grid(True,which='major',color=gridcolor,alpha=0.2,linestyle='-')
    grid(True,which='minor',color=gridcolor,alpha=0.1,linestyle='-')
    if x:
        labels = ax.get_xticklabels()
        setp(labels,rotation=rotation[0],fontsize=10)
    if y:
        labels = ax.get_yticklabels()
        setp(labels,rotation=rotation[1],fontsize=10)
    return
#}}}
#{{{ plot wrapper
global OLDplot
OLDplot = plot
global myplotfunc
myplotfunc = OLDplot
def whereblocks(a):
    """returns contiguous chunks where the condition is true
    but, see the "contiguous" method, which is more OO"""
    parselist = where(a)[0]
    jumps_at = where(diff(parselist)>1)[0]+1
    retlist = []
    lastjump = 0
    for jump in jumps_at:
        retlist += [parselist[lastjump:jump]]
        lastjump = jump
    retlist += [parselist[lastjump:]]
    return retlist
def autolegend(*args,**kwargs):
    #lg = legend(legendstr,'best'),loc = 2, borderaxespad = 0.)
    match_colors = False
    if 'match_colors' in list(kwargs.keys()):
        match_colors = kwargs.pop('match_colors')
    alpha = 0.45
    if 'alpha' in list(kwargs.keys()):
        alpha = kwargs.pop('alpha')
    if 'ax' in list(kwargs.keys()):
        ax_list = [kwargs.pop('ax')]
    else:
        ax_list = [gca()]
    if 'ax2' in list(kwargs.keys()):
        ax_list.append(kwargs.pop('ax2'))
    for ax in ax_list:
        if len(args)==0:
            lg = ax.legend(**kwargs)
        elif len(args)==1:
            lg = ax.legend(args[0],**kwargs)
        else:
            lg = ax.legend(args[0],args[1],**kwargs)
        if lg is None:
            raise ValueError("Warning! you called autolegend, but you don't seem to have anything labeled!!")
        else:
            lg.get_frame().set_alpha(alpha)
    if match_colors:
        for line, txt in zip(lg.get_lines(), lg.get_texts()): # from http://stackoverflow.com/questions/13828246/matplotlib-text-color-code-in-the-legend-instead-of-a-line 
                    txt.set_color(line.get_color())  
                    txt.set_alpha(line.get_alpha())  
    return lg
def autopad_figure(pad = 0.2,centered = False,figname = 'unknown'):
    #{{{ solve the axis issue --> this does just the left
    fig = gcf()
    ax = gca()
    labelsets = [] 
    #labelsets.append(('left',ax.get_yticklabels()))
    #labelsets.append(('left',ax.get_yticklines()))
    #labelsets.append(('right',ax.get_yticklines()))
    labelsets.append(('left',[ylabel(ax.get_ylabel())]))
    #labelsets.append(('bottom',ax.get_xticklabels()))
    #labelsets.append(('bottom',ax.get_xticklines()))
    if len(ax.get_xlabel()) > 0:
        labelsets.append(('bottom',[xlabel(ax.get_xlabel())]))
    #labelsets.append(('top',ax.get_xticklines()))
    if len(ax.get_title()) > 0:
        pass #labelsets.append(('top',[title(ax.get_title())]))
    compto = {}
    def on_draw(event):
        # find the sum of the widths of all things labeled with a 'y'
        spkwargs = {}
        compto['bottom'] = fig.subplotpars.bottom
        compto['left'] = fig.subplotpars.left
        compto['right'] = fig.subplotpars.right
        compto['top'] = fig.subplotpars.top
        for axisn in ['left','bottom','top','right']:
            bboxes = []
            labellist = [x[1] for x in labelsets if x[0] is axisn]
            for labels in labellist:
                for label in labels:
                    if isinstance(label, Line2D):
                        pass # just rely on the pad
                        #if any(map(lambda x: x == label.get_transform(),[ax.transData,ax.transAxes,fig.transFigure,None])):
                        #    print 'found it'
                        #else:
                        #    print 'didn not find it'
                        #bbox = label.get_window_extent(fig.canvas).inverse_transformed(ax.transData).inverse_transformed(fig.transFigure)
                    else:
                        try:
                            bbox = label.get_window_extent()
                        except Exception as e:
                            warnings.warn("I wasn't able to run autopad on figure"+figname+"\nGetting window extent throws error"+str(e))
                    # the figure transform goes from relative coords->pixels and we
                    # want the inverse of that
                    bboxes.append(bbox)
                # this is the bbox that bounds all the bboxes, again in relative
                # figure coords
            l = 0 
            if len(labellist):
                bbox = mtransforms.Bbox.union(bboxes)
                bboxi = bbox.inverse_transformed(fig.transFigure)
                if axisn in ['left','right']:
                    l = bboxi.width
                if axisn in ['top','bottom']:
                    l = bboxi.height
            l += pad
            if axisn in ['top','right']:
                l = 1-l
                if compto[axisn] > l:
                    spkwargs.update({axisn:l})
            else:
                if compto[axisn] < l:
                    spkwargs.update({axisn:l})
        if len(spkwargs) > 0:
            if centered and 'left' in list(spkwargs.keys()) and 'right' in list(spkwargs.keys()):
                big = max(r_[spkwargs['left'],1-spkwargs['right']])
                spkwargs.update({'left':big,'right':1-big})
            try:
                fig.subplots_adjust(**spkwargs) # pad a little
            except:
                raise RuntimeError('failed to adjust subplots spwargs = ',spkwargs)
            #print "adjusted to",spkwargs
            fig.canvas.draw()# recurse
        return False
    fig.canvas.mpl_connect('draw_event', on_draw)
    fig.subplots_adjust(left = 0, right = 1, top = 1, bottom =0)
    fig.canvas.draw()# it needs this to generate the 'renderers'
    fig.canvas.mpl_connect('draw_event', on_draw)
    fig.canvas.draw()
    return
    #}}}
def expand_x(*args):
    r'''expand the axes.  If an argument is passed, then it refers to the position relative to the current coordinates.  Values can be:
        :0: set this side of the axis to 0
        :None: leave this side of the axis alone
        :a double: rescale the distance from the center of the axis to this side by this number'''
    # this is matplotlib code to expand the x axis
    ax = gca()
    xlims = array(ax.get_xlim())
    width = abs(diff(xlims))
    thismean = mean(xlims)
    if len(args) > 0:
        if len(args) == 1 and isinstance(args, tuple):
            args = args[0]
        for j in range(2):
            if args[j] is None:
                pass
            elif args[j] == 0:
                xlims[j] = 0
            else:
                xlims[j] = args[j]*(xlims[j]-thismean) + thismean
    else:
        xlims[0] -= width/10
        xlims[1] += width/10
    ax.set_xlim(xlims)
def expand_y(*args):
    r'''expand the axes.  If an argument is passed, then it refers to the position relative to the current coordinates.  Values can be:
        :0: set this side of the axis to 0
        :None: leave this side of the axis alone
        :a double: rescale the distance from the center of the axis to this side by this number'''
    # this is matplotlib code to expand the x axis
    ax = gca()
    ylims = array(ax.get_ylim())
    width = abs(diff(ylims))
    thismean = mean(ylims)
    if len(args) > 0:
        if len(args) == 1 and isinstance(args, tuple):
            args = args[0]
        for j in range(2):
            if args[j] is None:
                pass
            elif args[j] == 0:
                ylims[j] = 0
            else:
                ylims[j] = args[j]*(ylims[j]-thismean) + thismean
    else:
        ylims[0] -= width/10
        ylims[1] += width/10
    ax.set_ylim(ylims)
def plot_label_points(x,y,labels,**kwargs_passed):
    kwargs = {'alpha':0.5,'color':'g','ha':'left','va':'center','rotation':0,'size':14}
    kwargs.update(kwargs_passed)
    for j in range(0,len(labels)):
        text(x[j],y[j],labels[j],**kwargs)
def addlabels(labelstring,x,y,labels):
    r'obsolete -- use plot_label_points'
    for j in range(0,len(labels)):
        text(x[j],y[j],labelstring%labels[j],alpha=0.5,color='g',ha='left',va='top',rotation=0)
def plot_color_counter(*args,**kwargs):
    """Try not to use this function any more -- the version-to-version support for capturing and setting color cycles in matplotlib is very very bad.  (And, the cycler object in newer versions of matplolib is confusing.) So, just import `cycle` from `itertools`, and use it to build a cycle that you directly call to set your properties.

    .. note::
        previous description:

        if passed an argument: make it so that the next line will have the properties given by the argument

        if not passed an argument: just return the current plot properties,so that I can cycle back to it"""
    ax = process_kwargs([('ax',gca())],kwargs)
    if len(args)>0:
        if LooseVersion(matplotlib.__version__) >= LooseVersion("1.5"):
            # {{{ find the element before the one we want
            retval = args[0]
            penultimate = next(ax._get_lines.prop_cycler)
            j = next(ax._get_lines.prop_cycler)
            not_in_list_counter = 1000.
            while j != args[0]:
                penultimate = j
                j = next(ax._get_lines.prop_cycler)
                not_in_list_counter -= 1
                if not_in_list_counter == 0:
                    raise ValueError("the value isn't in the cycler!")
            # }}}
            # {{{ now, set to the element before
            not_in_list_counter = 1000.
            while j != penultimate:
                j = next(ax._get_lines.prop_cycler)
                not_in_list_counter -= 1
                if not_in_list_counter == 0:
                    raise ValueError("the value isn't in the cycler!")
            # }}}
        else:
            try:
                ax._get_lines.count = args[0] # set the value of the color counter
            except:
                ax._get_lines.color_cycle = args[0] # set the value of the color counter
            retval = args[0]
    else:
        if LooseVersion(matplotlib.__version__) >= LooseVersion("1.5"):
            # {{{ I want to return the current element of the cycle
            one_too_far = next(ax._get_lines.prop_cycler)
            j = next(ax._get_lines.prop_cycler)
            not_in_list_counter = 1000.
            while j != one_too_far:
                penultimate = j
                j = next(ax._get_lines.prop_cycler)
                not_in_list_counter -= 1
                if not_in_list_counter == 0:
                    raise ValueError("the value isn't in the cycler!")
            retval = penultimate
            # }}}
        else:
            try: # this is different depending on the version of.core
                retval = ax._get_lines.count
            except:
                retval = ax._get_lines.color_cycle
    return retval
def contour_plot(xvals,yvals,zvals,color = 'k',alpha = 1.0,npts = 300,**kwargs):
    if 'inline_spacing' in list(kwargs.keys()):
        inline_spacing = kwargs.pop('inline_spacing')
    else:
        inline_spacing = 20
    xi = linspace(xvals.min(),xvals.max(),npts)
    yi = linspace(yvals.min(),yvals.max(),npts)
    #{{{ show the diffusivity
    #plot(array(xvals),array(yvals),'k')# to show where everything is
    zi = scipy_griddata((xvals,yvals),
        zvals,
        (xi[None,:],yi[:,None]))
    zi_min = zi[isfinite(zi)].min()
    zi_max = zi[isfinite(zi)].max()
    levels = r_[zi_min:zi_max:40j]
    CS = contour(xi,yi,zi,levels,colors = color,
            alpha = 0.25*alpha)
    oldspacing = levels[1]-levels[0]
    levels = r_[zi_min:zi_max:oldspacing*5]
    try:
        CS = contour(xi,yi,zi,levels,colors = color,
            alpha = alpha,**kwargs)
    except Exception as e:
        raise Exception(strm("Is there something wrong with your levels?:",levels,"min z",zi_min,"max z",zi_max,explain_error(e)))
    clabel(CS,fontsize = 9,inline = 1,
        #fmt = r'$k_\sigma/k_{\sigma,bulk} = %0.2f$',
        fmt = r'%0.2f',
        use_clabeltext = True,
        inline_spacing = inline_spacing,
        alpha = alpha)
    #}}}
def plot_updown(data,axis,color1,color2,symbol = '',**kwargs):
    if symbol == '':
        symbol = 'o'
    change = r_[1,diff(data.getaxis(axis))]
    changemask = change > 0
    if 'force_color' in list(kwargs.keys()) and kwargs['force_color'] == True:
        if hasattr(data,'other_info'):
            if 'plot_color' in data.get_prop():
                data.other_info.pop('plot_color')
    plot(data[axis,changemask],color1+symbol,**kwargs)
    if len(kwargs) > 0 and 'label' in list(kwargs.keys()): kwargs.pop('label') # if I'm doing a legend, I want it on the first
    plot(data[axis,~changemask],color2+symbol,**kwargs)
    return
def nextfigure(figurelist,name):
    'obsolete -- now use class'
    if isinstance(figurelist,figlist_var):
        figurelist.next(name)
        return figurelist
    else:
        print('Boo! not a new style name!')
    logger.debug(strm(lsafe('DEBUG figurelist, called with',name)))
    if name in figurelist:
        fig = figure(figurelist.index(name)+1)
        logger.debug(strm(lsafen('in',figurelist,'at figure',figurelist.index(name)+1,'switched figures')))
    else:
        fig = figure(len(figurelist)+1)
        fig.add_subplot(111)
        logger.debug(strm(lsafen('added, figure',len(figurelist)+1,'because not in figurelist',figurelist)))
        figurelist.append(name)
    return figurelist
def figlistret(first_figure,figure_list,*args,**kwargs):
    if 'basename' in list(kwargs.keys()):
        basename = kwargs['basename']
    else:
        basename = thisjobname()
    if first_figure is None:
        figure_list.show(basename+'.pdf')
        return args
    else:
        args += (figure_list,)
        if len(args) == 1:
            return args[0]
        else:
            return args
def figlistini(first_figure):
    r"""processes a figure list argument:
    typically, you want to have a figure_list keyword argument for every function, which is by default set to None, then call this on the argument -- it always returns a figure list, creating a new one if required
    similarly, somewhere I have another guy that processes the output, so that if it's set to None, it will by default dump and show the figure list,
    and not return a figure list in the output"""
    if first_figure is None:
        return figlist_var() 
    else:
        return first_figure
def figlistini_old(first_figure):
    if isinstance(first_figure,figlist_var):
        return first_figure
    else:
        print("Boo, not a new style name! (initialize)")
    logger.debug(strm(lsafe('DEBUG: initialize figlist')))
    if first_figure is None:
        logger.debug(strm(lsafen('empty')))
        return []
    else:
        logger.debug(strm(lsafen(first_figure.figurelist)))
        return first_figure
class figlist(object):
    r"""
    Attributes
    ----------
    basename : str
        A basename that can be changed to generate different sets of figures with different basenames.
        For example, this is useful if you are looping over different sets of data,
        and generating the same set of figures for each set of data (which would correspond to a basename).
    figurelist : list
        A list of the figure names
    figdict : dict
        A dictionary containing the figurelist and the figure numbers or objects that they correspond to.
        Keys of this dictionary must be elements of `figurelist`.
    propdict : dict
        Maintains various properties for each element in figurelist.
        Keys of this dictionary must be elements of `figurelist`.
    """
    def __init__(self,*arg,**kwargs):
        r"""Initialize a figure list, which can be used to generate a series of
        figures from the command line or prompt.  Then the same code (if
        `figlist_var` is used) can be included inside a ``python`` environment
        in a latex document.

        Parameters
        ----------
        black : double
            A fractional number giving how "black" "black" is. Typically 1.0 is
            actually too dark and makes things hard to see.
        mlab : object
            If you want to use mayavi, this should be the mlab (module?)
        file_name : str
            This is the argument passed to :func:`self.show`, and used to
            construct the file names.
        """
        self.black, self.env, self.mlab, self.file_name, self.line_spacing = process_kwargs([
            ('black',0.9),
            ('env',''),
            ('mlab','BLANK'),
            ('file_name','BLANK'),
            ('line_spacing','BLANK'),
            ],
                kwargs, pass_through=True)
        if len(kwargs) > 0:
            self.lplot_kwargs = kwargs
        if self.mlab == 'BLANK': del self.mlab
        if self.file_name == 'BLANK': del self.file_name
        if self.line_spacing == 'BLANK': del self.line_spacing
        logger.debug('DEBUG: initialize figlist')
        if len(arg) == 0:
            self.figurelist = []
        else:
            self.figurelist = arg[0]
        if len(kwargs) > 0:
            self.figurelist.append(kwargs)
        self.units = {}
        self.autolegend_list = {}
        self.twinx_list = {}
        self.basename = None
        return
    def twinx(self,autopad = False,orig = False,color = None):
        #self.figurelist.insert(self.get_fig_number(self.current)-1,{'autopad':False}) #doesn't work because it changes the figure number; I can get the number with fig = gcf(); fig.number, but I can't set it; it would be best to switch to using a list that contains all the figure numbers to match all their names -- or alternatively, note that matplotlib allows you to give them names, though I don't know how that works
        if self.current in list(self.twinx_list.keys()):
            ax1,ax2 = self.twinx_list[self.current]
            if color is not None:
                if 'twinx_color' not in list(self.propdict[self.current].keys()):
                        ax2.tick_params(axis = 'y',colors = color)
                        ax2.yaxis.label.set_color(color)
                        ax2.spines['right'].set_color(color)
                        self.propdict[self.current]['twinx_color'] = color
                else:
                    if color != self.propdict[self.current]['twinx_color']:
                        raise ValueError("conflicting values for the twinx color have been given!!")
        else:
            if autopad: autopad_figure()
            ax1 = gca()
            twinx()
            ax2 = gca()
            self.twinx_list[self.current] = (ax1,ax2)
            if color is not None:
                ax2.tick_params(axis = 'y',colors = color)
                ax2.yaxis.label.set_color(color)
                ax2.spines['right'].set_color(color)
                self.propdict[self.current]['twinx_color'] = color
        if orig:
            sca(ax1)
            return ax1
        else:
            sca(ax2)
            return ax2
    def use_autolegend(self,value = None):
        'No argument sets to true if it\'s not already set'
        if value is None:
            if not self.current in list(self.autolegend_list.keys()):
                self.autolegend_list.update({self.current:True})
            else: #leave it alone
                return
        else: #passed an explicit value
            self.autolegend_list.update({self.current:value})
            return
    def push_marker(self):
        """save the current plot to a "stack" so we can return to it with "pop_marker" """
        if hasattr(self,'current'): # if not, this is the first plot
            if not hasattr(self,'pushlist'):
                self.pushlist = []
            if not hasattr(self,'pushbasenamelist'):
                self.pushbasenamelist = []
            self.pushlist.append(self.current)
            self.pushbasenamelist.append(self.basename)
        return
    def pop_marker(self):
        """use the plot on the top of the "stack" (see push_marker) as the current plot"""
        if hasattr(self,'pushlist') and len(self.pushlist) > 0: # otherwise, we called push with no current plot
            if hasattr(self,'pushbasenamelist'):
                self.basename = self.pushbasenamelist.pop()
            self.next(self.pushlist.pop())
        return
    def get_num_figures(self):
        cleanlist = [x for x in self.figurelist if isinstance(x, str)]
        return len(cleanlist)
    def get_fig_number(self,name):
        cleanlist = [x for x in self.figurelist if isinstance(x, str)]
        try:
            return cleanlist.index(name)+1
        except ValueError:
            raise ValueError(strm("You are looking for",name,
                "which isn't in the list of figures",cleanlist))
    def next(self,input_name, legend=False,
            boundaries=None, twinx=None, fig=None,
            **kwargs):
        r"""Switch to the figure given by input_name, which is used not only as
        a string-based name for the figure, but also as a default title and as
        a base name for resulting figure files.

        **In the future, we actually want this to track the appropriate axis object!**

        Parameters
        ----------
        legend : bool
            If this is set, a legend is created *outside* the figure.
        twinx : {0,1}
            :1: plots on an overlayed axis (the matplotlib twinx) whose y axis
                is labeled on the right when you set this for the first time, you
                can also set a `color` kwarg that controls the coloring of the
                right axis. 
            :0: used to switch back to the left (default) axis
        boundaries :
            **need to add description**
        kwargs : dict
            Any other keyword arguments are passed to the matplotlib (mayavi)
            figure() function that's used to switch (create) figures.
        """
        # {{{ basic setup
        if not hasattr(self,'figdict'):
            self.figdict = {} # the dictionary of the various figures
        if not hasattr(self,'propdict'):
            self.propdict = {} # the properties belonging to those same figures
        logger.debug(strm("for plot",input_name,"basename is",self.basename))
        if (self.basename is not None #basename for groups of figures
                # I need to check that the basename hasn't already been added
                and not input_name.startswith(self.basename)):
            name = self.basename + ' ' + input_name
        else:
            logger.debug(strm("not using a basename",self.basename is not None))
            name = input_name
        # }}}
        if name.find('/') > 0:
            raise ValueError("don't include slashes in the figure name, that's just too confusing")
        logger.debug(strm('called with',name))
        if name in self.figurelist:# figure already exists
            if hasattr(self,'mlab'):
                # with this commit, I removed the kwargs and bgcolor, not sure why
                fig = self.mlab.figure(self.get_fig_number(name))
                fig.scene.render_window.aa_frames = 20
                fig.scene.anti_aliasing_frames = 20
            else:
                logging.debug(strm("I'm changing to figure",self.get_fig_number(name),"for",name))
                fig = self.figdict[name]
                figure(self.figdict[name].number)
            self.current = name
            #logging.debug(strm('in',self.figurelist,'at figure',self.get_fig_number(name),'switched figures'))
            if boundaries is not None:
                if 'boundaries' not in list(self.propdict[self.current].keys()) or self.propdict[self.current]['boundaries'] != boundaries:
                    raise ValueError("You're giving conflicting values for boundaries")
            if legend:
                if 'legend' not in list(self.propdict[self.current].keys()) or self.propdict[self.current]['legend'] != legend:
                    raise ValueError("You're giving conflicting values for legend")
        else:# figure doesn't exist yet
            num_figs_before_add = self.get_num_figures()
            self.current = name
            if self.current not in list(self.propdict.keys()):
                self.propdict[self.current] = {}
            if boundaries == False:
                self.propdict[self.current]['boundaries'] = False
                self.setprops(boundaries = False)
            if legend:
                self.propdict[self.current]['legend'] = True
                if 'figsize' not in list(kwargs.keys()):
                    kwargs.update({'figsize':(12,6)})
                if hasattr(self,'mlab'):
                    fig = self.mlab.figure(num_figs_before_add+1,bgcolor = (1,1,1),**kwargs)
                    fig.scene.render_window.aa_frames = 20
                    fig.scene.anti_aliasing_frames = 20
                else:
                    fig = figure(num_figs_before_add+1,**kwargs)
                fig.add_axes([0.075,0.2,0.6,0.7]) # l b w h
                self.use_autolegend('outside')
            else:
                self.propdict[self.current]['legend'] = False
                if fig is None:
                    if hasattr(self,'mlab'):
                        fig = self.mlab.figure(num_figs_before_add+1,bgcolor = (1,1,1),**kwargs)
                        fig.scene.render_window.aa_frames = 20
                        fig.scene.anti_aliasing_frames = 20
                    else:
                        fig = figure(num_figs_before_add+1,**kwargs)
                if twinx is not None:
                    fig.add_subplot(111)
            logger.debug(strm('added figure',len(self.figurelist)+1,'because not in figurelist',self.figurelist))
            self.figurelist.append(name)
            self.figdict.update({self.current:fig})
            if boundaries == False:
                self.setprops(boundaries = True)# set this back
        if twinx is not None:
            self.propdict[self.current]['twinx'] = True
            if twinx == 0:
                self.twinx(orig = True)
                fig = gcf()
            elif twinx == 1:
                self.twinx()
                fig = gcf()
            else:
                raise ValueError('If you pass twinx, pass 0 for the original or 1 for the right side')
            self.figdict.update({self.current:fig})
        return fig
    def plot(self,*args,**kwargs):
        r"""
        Parameters
        ----------
        linestyle: {':','--','.','etc.'}
            the style of the line
        plottype: {'semilogy','semilogx','loglog'}
            Select a logarithmic plotting style.
        nosemilog: True
            Typically, if you supply a log-spaced axis,
            a semilogx plot will be automatically selected.
            This overrides that behavior.
            Defaults to False.
        """
        if 'label' in kwargs.keys() or 'label_format_string' in kwargs.keys():
            self.use_autolegend()
        human_units = True
        if 'human_units' in list(kwargs.keys()):
            human_units = kwargs.pop('human_units')
        if human_units:
            firstarg = self.check_units(args[0],0,1) # check units, and if need be convert to human units, where x is the first dimension and y is the last
        else:
            firstarg = args[0]
        if 'label' not in list(kwargs.keys()) and isinstance(args[0],nddata):
            thisname = args[0].name()
            if thisname is not None:
                kwargs['label'] = thisname
        retval = plot(*tuple((firstarg,)+args[1:]),**kwargs)#just a placeholder for now, will later keep units + such
        ax = gca()
        if ax.get_title() is None or len(ax.get_title()) == 0:
            try:
                title(self.current)
            except:
                title('untitled')
        return retval
    def phaseplot_finalize(self):
        ("Performs plot decorations that are typically desired for a manual phasing"
        " plot.  This assumes that the ``y``-axis is given in units of half-cycles"
        " ($\pi$ radians).")
        ax = gca()
        ylim(-1,1)
        gridandtick(ax)
        ylabel(r'$\phi / \pi$')
        # now show the pi/2 lines
        axhline(y = 0.5,color = 'r',alpha = 0.5,linewidth = 2)
        axhline(y = -0.5,color = 'r',alpha = 0.5,linewidth = 2)
        return
    def check_units(self, testdata, x_index, y_index):
        logger.debug(strm("-"*30))
        logger.debug(strm("called check_units for figure",self.current))
        if isinstance(testdata,nddata):
            logger.debug(strm("(check_units) it's nddata"))
            testdata = testdata.copy().human_units()
            if len(testdata.dimlabels) > 1:
                logger.debug(strm("(check_units) more than one dimension"))
                if not hasattr(self,'current'):
                    raise ValueError("give your plot a name (using .next()) first! (this is used for naming the PDF's etc)")
                if self.current in list(self.units.keys()):
                        theseunits = (testdata.get_units(testdata.dimlabels[x_index]),testdata.get_units(testdata.dimlabels[y_index]))
                        if theseunits != self.units[self.current] and theseunits[0] != self.units[self.current]:
                                raise ValueError("the units don't match (old units %s and new units %s)! Figure out a way to deal with this!"%(theseunits,self.units[self.current]))
                else:
                    if isinstance(testdata,nddata):
                        self.units[self.current] = (testdata.get_units(testdata.dimlabels[x_index]),testdata.get_units(testdata.dimlabels[y_index]))
            else:
                logger.debug(strm("(check_units) only one dimension"))
                if not hasattr(self,'current'):
                    self.next('default')
                if self.current in list(self.units.keys()):
                    theseunits = (testdata.get_units(testdata.dimlabels[x_index]))
                    testunits = self.units[self.current]
                    if theseunits != testunits:
                        if isinstance(testunits, tuple) and testunits[1] is None:
                            pass
                        else:
                            raise ValueError("the units don't match (old units %s and new units %s)! Figure out a way to deal with this!"%(self.units[self.current],theseunits))
                else:
                    self.units[self.current] = (testdata.get_units(testdata.dimlabels[x_index]))
        logger.debug(strm("-"*30))
        return testdata
    def adjust_spines(self,spines):
        ax = gca()
        #{{{ taken from matplotlib examples
        for loc, spine in list(ax.spines.items()):
            if loc in spines:
                spine.set_position(('outward',10)) # outward by 10 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine

        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([])

        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([])
        #}}}
    def grid(self):
        ax = gca()
        if self.black:
            gridandtick(ax,gridcolor = r_[0.5,0.5,0.5])
        else:
            gridandtick(ax,gridcolor = r_[0,0,0])
        return
    image = this_plotting.image.fl_image
    def marked_text(self,marker,input_text="",sep='\n'):
        """Creates a named `marker` where we can place text.   If `marker`
        has been used, goes back and places text there."""
        if not hasattr(self,'textdict'):
            self.textdict = {}
        if marker in list(self.textdict.keys()):
            idx = self.textdict[marker]
            self.figurelist[idx]['print_string'] = (
                    self.figurelist[idx]['print_string']
                    + sep + input_text )
        else:
            self.setprops(print_string=input_text)
            idx = len(self.figurelist)-1
            self.textdict[marker] = idx
    def text(self,mytext):
        self.setprops(print_string = mytext)
    def setprops(self,**kwargs):
        self.figurelist.append(kwargs)
    def show_prep(self):
        for k,v in list(self.autolegend_list.items()):
            kwargs = {}
            if v:
                if isinstance(v, str):
                    if v[0:7] == 'colored':
                        kwargs.update(dict(match_colors = True))
                        v = v[7:]
                        if v == '':
                            v = True
                    if v == 'outside':
                        kwargs.update(dict(bbox_to_anchor=(1.05,1),loc = 2,borderaxespad=0.))
                self.next(k)
                logger.debug(strm("I am about to assign a legend for ",k,". Is it in the figurelist?:",k in self.figurelist))
                logger.debug(strm("print out the legend object:",gca().legend()))
                try:
                    autolegend(**kwargs)
                except:
                    try:
                        self.twinx(orig = True)
                    except Exception as e:
                        raise Exception(strm('error while trying to run twinx to place legend for',k,'\n\tfiglist is',self.figurelist,explain_error(e)))
                    try:
                        autolegend(**kwargs)
                    except Exception as e:
                        raise Exception(strm('error while trying to run autolegend function for',k,'\n\tfiglist is',self.figurelist,explain_error(e)))
    def show(self,*args,**kwargs):
        self.basename = None # must be turned off, so it can cycle through lists, etc, on its own
        if 'line_spacing' in list(kwargs.keys()): kwargs.pop('line_spacing')# for latex only
        if len(kwargs) > 0:
            raise ValueError("didn't understand kwargs "+repr(kwargs))
        logger.debug(strm("before show_prep, figlist is",self.figurelist))
        logger.debug(strm("before show_prep, autolegend list is",self.autolegend_list))
        self.show_prep()
        #{{{ just copy from fornnotebook to get the print string functionality
        kwargs = {}
        for figname in self.figurelist:
            logger.debug(strm("showing figure"+lsafen(figname)))
            if isinstance(figname, dict):
                kwargs.update(figname)
                if 'print_string' in kwargs:
                    print('\n\n')
                    print(kwargs.pop('print_string'))
                    print('\n\n')
        #}}}
        if len(args) == 1:
            if (args[0][:-4] == '.pdf') or (args[0][:-4] == '.png') or (args[0][:-4] == '.jpg'):
                print("you passed me a filename, but I'm just burning it")
        if hasattr(self,'mlab'):
            print("running mlab show!")
            self.mlab.show()
        else:
            #print "not running mlab show!"
            show()
    def label_point(self, data, axis, value, thislabel,
            show_point=True, xscale=1, **new_kwargs):
        """only works for 1D data: assume you've passed a single-point nddata, and label it

        xscale gives the unit scaling

        ..todo::

            Improve the unit scaling, so that this would also work.

            Allow it to include a format string that would use the value.
        Parameters
        ----------

        show_point : bool

            Defaults to `True`. Actually generate a point (circle), *vs.*
            just the label.
        """
        kwargs = {'alpha':0.5,'color':'k','ha':'left','va':'bottom','rotation':45,'size':14}
        kwargs.update(new_kwargs)
        y = double(data[axis:value].data)
        x_ind = argmin(abs(data.getaxis(axis)-value))
        x = data.getaxis(axis)[x_ind]
        text(x/xscale, y, thislabel, **kwargs)
        if show_point:
            plot(x/xscale, y, 'o', color=kwargs["color"],
                    alpha=kwargs["alpha"])
        return
    def header(self,number_above,input_string):
        header_list = ['\\section','\\subsection','\\subsubsection','\\paragraph','\\subparagraph']
        self.text(header_list[number_above+1]+'{%s}'%input_string)
        return number_above + 1
    def mesh(self,plotdata,Z_normalization = None,equal_scale = True,
            lensoffset = 1e-3,
            show_contours = False,
            grey_surf = False,
            **kwargs):
        plotdata = self.check_units(plotdata,0,1)
        if hasattr(self,'mlab'):
            fig = self.figdict[self.current]
            fig.scene.disable_render = True
            X,Y,Z,x_axis,y_axis = plotdata.matrices_3d(also1d = True)# return the axes, and also alter "plotdata" so it's downsampled
            X_normalization = X.max()
            X /= X_normalization
            if equal_scale:
                Y_normalization = X_normalization
            else:
                Y_normalization = Y.max()
            Y /= Y_normalization
            if Z_normalization is None:
                Z_normalization = Z.flatten().max()
            Z /= Z_normalization
            surf_kwargs = {}
            if grey_surf:
                surf_kwargs.update(color = (0.5,0.5,0.5))# opacity and the contour lines don't play well, otherwise I would like to make this transluscent
            self.mlab.surf(X,Y,Z,**surf_kwargs)
            if show_contours:
                contour_kwargs = {'line_width':24}
                contour_kwargs.update(opacity = 0.5)
                if not grey_surf:
                    contour_kwargs.update(color = (1,1,1))
                self.mlab.contour_surf(X,Y,Z+lensoffset,contours = r_[-1:1:10j].tolist(),**contour_kwargs)
                contour_kwargs.update(opacity = 0.1)
                self.mlab.contour_surf(X,Y,Z+lensoffset,contours = r_[-1:1:46j].tolist(),**contour_kwargs)# for some reason, 46 gives alignment (I think 9+1 and 9*5+1)
            if equal_scale:
                self.generate_ticks(plotdata,(x_axis,y_axis),X_normalization,Z_normalization)
            else:
                self.generate_ticks(plotdata,(x_axis,y_axis),X_normalization,Z_normalization,y_rescale = Y_normalization/X_normalization)
            fig.scene.disable_render = False
        else:
            # this should be upgraded, or rather moved to here
            plotdata.meshplot(alpha=1.0, cmap=cm.jet, **kwargs)
        return Z_normalization
    def generate_ticks(self,plotdata,axes,rescale,z_norm = None,y_rescale = 1,text_scale = 0.05,follow_surface = False,
            lensoffset = 0.5e-2,
            line_width = 1e-3,
            tube_radius = 1e-3,
            fine_grid = False,
            ):
        'generate 3d ticks and grid for mayavi'
        if follow_surface and z_norm is None:
            raise ValueError("if you choose to generate the mesh -- i.e. follow the surface -- then you need to pass the z normalization")
        x_axis,y_axis = axes
        x_dim = plotdata.dimlabels[0]
        y_dim = plotdata.dimlabels[1]
        def gen_list(thisaxis,desired_ticks = 7.):
            #{{{ out of the following list, choose the one that gives as close as possible to the desired ticks
            axis_span = thisaxis.max() - thisaxis.min()
            possible_iterators = r_[0.1,0.5,1,5,10,20,30,50,100,200,500,1000]
            iterator = possible_iterators[argmin(abs(axis_span/desired_ticks -
                possible_iterators))]
            #}}}
            logger.debug(strm('iterator is',iterator))
            return iterator,r_[ceil(thisaxis.min()/iterator):
                floor(thisaxis.max()/iterator)+1]*iterator
        #{{{ now, I need to get the list of multiples that falls inside the axis span
        xiterator,xlist = gen_list(x_axis)
        yiterator,ylist = gen_list(y_axis)
        logger.debug(strm('range of x ',x_axis.min(),x_axis.max()))
        logger.debug(strm('xlist',xlist))
        logger.debug(strm(plotdata.unitify_axis(0)))
        logger.debug(strm('range of y ',y_axis.min(),y_axis.max()))
        logger.debug(strm('ylist',ylist))
        logger.debug(strm(plotdata.unitify_axis(1)))
        #}}}
        if xiterator < 1:
            x_ticklabels = ['{:0.1f}'.format(j) for j in xlist]
        else:
            x_ticklabels = ['{:0.0f}'.format(j) for j in xlist]
        if yiterator < 1:
            y_ticklabels = ['{:0.1f}'.format(j) for j in ylist]
        else:
            y_ticklabels = ['{:0.0f}'.format(j) for j in ylist]
        #{{{ rescale absolutely everything
        xlist /= rescale
        ylist /= (rescale*y_rescale)
        x_axis /= rescale
        y_axis /= (rescale*y_rescale)
        #}}}
        x_range = r_[x_axis.min(),x_axis.max()]
        y_range = r_[y_axis.min(),y_axis.max()]
        extension_factor = text_scale * 3
        #{{{ y ticks
        if follow_surface:
            if fine_grid:
                dy = ylist[1]-ylist[0]
                finer_ylist = r_[ylist[0]-dy:ylist[-1]+dy:1j*((len(ylist)+2-1)*5+1)]
                finer_ylist = finer_ylist[finer_ylist>=y_axis.min()]
                finer_ylist = finer_ylist[finer_ylist<=y_axis.max()]
            else:
                finer_ylist = ylist
            for j,y in enumerate(finer_ylist):
                x_linedata = plotdata.getaxis(x_dim)/rescale
                z_linedata = plotdata[y_dim:(y*rescale)].data.flatten()/z_norm
                self.mlab.plot3d(x_linedata,y*ones_like(x_linedata),
                        z_linedata+lensoffset,
                        color = (0,0,0), line_width = line_width,
                        tube_radius = tube_radius)
        for j,y in enumerate(ylist):
            self.mlab.plot3d(x_range+extension_factor*r_[-1,1],
                    y*ones(2),zeros(2),
                    color = (0,0,0), line_width = line_width,
                    tube_radius = tube_radius)
            self.mlab.text3d(x_range[0]-2*extension_factor, y, 0,
                    y_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
            self.mlab.text3d(x_range[1]+2*extension_factor, y, 0,
                    y_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
        self.mlab.text3d(x_range[1] + 3 * extension_factor,y_range.mean(), 0,
                plotdata.unitify_axis(1), color = (0,0,0),
                scale = text_scale,
                orient_to_camera = False,
                orientation = (0,0,90))# the last angle appears to be rotaiton about z
        #}}}
        #{{{ x ticks
        if follow_surface:
            if fine_grid:
                dx = xlist[1]-xlist[0]
                finer_xlist = r_[xlist[0]-dx:xlist[-1]+dx:1j*((len(xlist)+2-1)*5+1)]
                finer_xlist = finer_xlist[finer_xlist>=x_axis.min()]
                finer_xlist = finer_xlist[finer_xlist<=x_axis.max()]
            else:
                finer_xlist = xlist
            for j,x in enumerate(finer_xlist):
                y_linedata = plotdata.getaxis(y_dim)/(rescale*y_rescale)
                z_linedata = plotdata[x_dim:(x*rescale)].data.flatten()/z_norm
                self.mlab.plot3d(x*ones_like(y_linedata),y_linedata,
                        z_linedata+lensoffset,
                        color = (0,0,0), line_width = line_width,
                        tube_radius = tube_radius)
        for j,x in enumerate(xlist):
            self.mlab.plot3d(x*ones(2),y_range+extension_factor*r_[-1,1],
                    zeros(2),
                    color = (0,0,0), line_width = line_width,
                    tube_radius = tube_radius)
            self.mlab.text3d(x, y_range[0]-2*extension_factor, 0,
                    x_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
            self.mlab.text3d(x, y_range[1]+2*extension_factor, 0,
                    x_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
        self.mlab.text3d(x_range.mean(), y_range[1] + 3 * extension_factor,
                0,
                plotdata.unitify_axis(0), color = (0,0,0),
                scale = text_scale,
                orient_to_camera = False,
                orientation = (0,0,180))# the last angle appears to be rotaiton about z
        #}}}
        return
    def __enter__(self):
        return self
    def __exit__(self, exception_type, exception_value, traceback):
        r'''show the plots, unless there are errors.

        Because this is executed before raising any errors, we want to avoid showing any plots if there are errors.
        Otherwise, it gets very confusing.
        '''
        if exception_type is None:
            if hasattr(self,'file_name'):
                if hasattr(self,'line_spacing'):
                    self.show(self.file_name,line_spacing = self.line_spacing)
                else:
                    self.show(self.file_name)
            else:
                self.show()
            return
def text_on_plot(x,y,thistext,coord = 'axes',**kwargs):
    ax = gca()
    if coord == 'axes':
        newkwargs = {'transform':ax.transAxes,'size':'x-large',"horizontalalignment":'center'}
    elif coord == 'data':
        print("Yes, I am using data transform")
        newkwargs = {'transform':ax.transData,'size':'small',"horizontalalignment":'right'}
    color = None
    if 'match_data' in list(kwargs.keys()):
        if isinstance(kwargs['match_data'], list):
            color = kwargs['match_data'][-1].get_color() # get the color of the last line
        elif kwargs['match_data'].get_plot_color() is not None:
            color = kwargs['match_data'].get_plot_color() # don't know when this works, but apparently, it does!
        if color is not None:
            newkwargs.update({'color':color})
        else:
            raise ValueError('You passed match_data to text_on_plot, but I can\'t find a color in the object')
        kwargs.pop('match_data')
    newkwargs.update(kwargs)
    return text(x,y,thistext,**newkwargs)
def plot(*args,**kwargs):
    """The base plotting function that wraps around matplotlib to do a couple convenient things.

    Parameters
    ----------
    label_format_string: str
        If supplied, it formats the values of the other dimension to turn them into a label string.
    human_units: bool
    """
    global myplotfunc
    has_labels = False
    #{{{ deal with axes and some other kwargs
    ax,human_units,label_format_string,normalize,noerr,longest_is_x = process_kwargs([('ax',gca()),
        ('human_units',False),
        ('label_format_string',None),
        ('normalize',False),
        ('noerr',False),
        ('longest_is_x',True),
        ],kwargs,pass_through = True)
    #}}}
    myplotfunc = ax.plot # default
    #{{{ all possible properties
    myformat = None 
    myxlabel = None
    myylabel = None
    myx = None
    myy = None
    #}}}
    #{{{assign all the possible combinations
    if len(args)==1:
        myy = args[0]
    elif (len(args)==2) and (isinstance(args[1], str)):
        myy = args[0]
        myformat = args[1]
    else:
        myx = args[0]
        myy = args[1]
    if len(args)==3:
        myformat = args[2]
    if isscalar(myx):
        myx = array([myx])
    if isscalar(myy):
        myy = array([myy])
    #}}}
    x_inverted = False
    #{{{ parse nddata
    if isinstance(myy,nddata):
        myy = myy.copy()
        # {{{ automatically reduce any singleton dimensions
        if not len(myy.dimlabels) == 1:
            if any(array(myy.data.shape) == 1):
                for singleton_dim in [lb for j,lb in enumerate(myy.dimlabels) if myy.data.shape[j] == 1]:
                    myy = myy[singleton_dim,0]
        # }}}
        if len(myy.data.shape)>1 and longest_is_x:
            longest_dim = argmax(myy.data.shape)
            all_but_longest = set(range(len(myy.data.shape)))^set((longest_dim,))
            if len(all_but_longest) > 0:
                last_not_longest = max(all_but_longest)
            else:
                last_not_longest = -1
            all_but_longest = list(all_but_longest) # seems to be sorted by default
        else:
            longest_dim = 0 # treat first as x, like before
            last_not_longest = -1
            if len(myy.data.shape)>1:
                all_but_longest = set(range(len(myy.data.shape)))^set((longest_dim,))
                all_but_longest = list(all_but_longest)
            else:
                all_but_longest = []
        if human_units: myy = myy.human_units()
        if myy.get_plot_color() is not None\
            and 'color' not in list(kwargs.keys()):# allow override
            kwargs.update({'color':myy.get_plot_color()})
        if myy.name() is not None:
            myylabel = myy.name()
        else:
            myylabel = 'data'
        myylabel = myy.unitify_axis(myylabel,is_axis = False)
        if (len(myy.dimlabels)>0):
            myxlabel = myy.unitify_axis(longest_dim)
        if myx is None:
            try:
                myx = myy.getaxis(myy.dimlabels[longest_dim])
            except:
                if len(myy.data.shape) == 0:
                    raise ValueError("I can't plot zero-dimensional data (typically arises when you have a dataset with one point)")
                myx = r_[0:myy.data.shape[longest_dim]]
        if not noerr and isinstance(myy.data_error, ndarray) and len(myy.data_error)>0: #then this should be an errorbar plot
            def thiserrbarplot(*tebargs,**tebkwargs):
                if isinstance(tebargs[-1], str):
                    tebkwargs.update({'fmt':tebargs[-1]})
                    return ax.errorbar(*tebargs[:-1],**tebkwargs)
                else:
                    return ax.errorbar(*tebargs,**tebkwargs)
            myplotfunc = thiserrbarplot
            #{{{ pop any singleton dims
            myyerror = myy.get_error()
            myyerror = squeeze(myyerror)
            #}}}
            kwargs.update({'yerr':myyerror})
            valueforxerr = myy.get_error(myy.dimlabels[longest_dim])
            if valueforxerr is not None: # if we have x errorbars too
                #print "DEBUG decided to assign to xerr:",valueforxerr
                kwargs.update({'xerr':valueforxerr})
        #{{{ deal with axis labels along y
        try:
            yaxislabels = myy.getaxis(myy.dimlabels[last_not_longest])
        except:
            pass
        # at this point, if there is no axis label, it will break and go to pass
        if yaxislabels is not None:
            if len(yaxislabels) > 0:
                if isinstance(yaxislabels[0], string_):
                    has_labels = True
                elif label_format_string is not None:
                    yaxislabels = [label_format_string%j for j in yaxislabels]
                    has_labels = True
        #}}}
        # {{{ add label if name is present, and squeeze -- could do this instead of ylabel, above
        if myy.get_prop('x_inverted'):
            x_inverted=True
        #myy_name = myy.name()
        if len(myy.data.shape) == 1:
            myy = myy.data
        else:
            myy = squeeze(myy.data.transpose([longest_dim]+all_but_longest))
        #if len(myy.data) == 1 and 'label' not in kwargs.keys() and myy_name is not None:
        #    kwargs.update('label',myy_name)
        # }}}
    #}}}
    #{{{ semilog where appropriate
    if (myx is not None) and (len(myx)>1) and all(myx>0): # by doing this and making myplotfunc global, we preserve the plot style if we want to tack on one point
        try:
            b = diff(log10(myx))
        except Exception as e:
            raise Exception(strm('likely a problem with the type of the x label, which is',myx,explain_error(e)))
        if (size(b)>3) and all(abs((b-b[0])/b[0])<1e-4) and not ('nosemilog' in list(kwargs.keys())):
            if 'plottype' not in list(kwargs.keys()):
                myplotfunc = ax.semilogx
    if ('nosemilog' in list(kwargs.keys())):
        #print 'this should pop nosemilog'
        kwargs.pop('nosemilog')
    if 'plottype' in list(kwargs.keys()):
        if kwargs['plottype'] == 'semilogy':
            myplotfunc = ax.semilogy
        elif kwargs['plottype'] == 'semilogx':
            myplotfunc = ax.semilogx
        elif kwargs['plottype'] == 'loglog':
            myplotfunc = ax.loglog
        elif kwargs['plottype'] == 'linear':
            myplotfunc = ax.plot
        else:
            raise ValueError(strm("plot type",kwargs['plottype'],"not allowed!"))
        kwargs.pop('plottype')
    #}}}
    #{{{ take care of manual colors
    if myformat is not None:
        colorpos = myformat.find('#')
        if  colorpos > -1:
            kwargs.update({'color':myformat[colorpos:colorpos+7]})
            myformat = myformat[0:colorpos] + myformat[colorpos+7:]
        ##kwargs.update({'fmt':myformat})
        linematched = False
        for linestyle in ['-','--','-.',':','None','  ']:
            if myformat.find(linestyle) > -1:
                linematched = True
                myformat.replace(linestyle,'')
                kwargs.update({'linestyle':linestyle})
        for markerlabel in ['o','.','d']:
            if myformat.find(markerlabel) > -1:
                if not linematched: kwargs.update({'linestyle':''})
                myformat.replace(markerlabel,'')
                kwargs.update({'marker':markerlabel})
        if len(myformat) == 0:
            myformat = None
    #}}}
    if normalize is not None and normalize:
        myy /= myy.max()
    #{{{ hsv plots when we have multiple lines
    if len(shape(myy.squeeze()))>1 and sum(array(shape(myy))>1):
        #{{{ hsv plots
        retval = []
        for j in range(0,myy.shape[1]):
            #{{{ this is the way to assign plot arguments
            plotargs = [k for k in (myx,myy[:,j],myformat) if k is not None]
            #}}}
            #{{{ here, i update the kwargs to include the specific color for this line
            newkwargs = kwargs.copy() # kwargs is a dict
            newkwargs.update({'color':cm.hsv(double(j)/double(myy.shape[1]))})
            #}}}
            #{{{ here, I update to use the labels
            if has_labels:
                newkwargs.update({'label':yaxislabels[j]})
            #}}}
            if any(isinf(myy)):
                myy[isinf(myy)] = NaN # added this to prevent an overflow error
            try:
                retval += [myplotfunc(*tuple(plotargs),**newkwargs)]
            except Exception as e:
                raise RuntimeError(strm("Error trying to plot using function",
                    myplotfunc, '\nwith',len(plotargs), "arguments",
                    '\nwhich were\n',plotargs, "\nand had len\n",
                    list(map(len, plotargs)), "and", len(newkwargs),
                    "\noptions", newkwargs, "of len",
                    ', '.join([str(type(j)) + " " + str(j) if isscalar(j)
                        else str(len(j)) for j in list(newkwargs.values())]),
                    explain_error(e)))
            if x_inverted:
                these_xlims = ax.get_xlim()
                ax.set_xlim((max(these_xlims),min(these_xlims)))
        #}}}
        #}}}
    else:
        plotargs = [j for j in [myx,real(myy),myformat] if j is not None]
        try:
            #print 'DEBUG plotting with args',plotargs,'and kwargs',kwargs,'\n\n'
            retval = myplotfunc(*plotargs,**kwargs)
        except Exception as e:
            raise RuntimeError(strm('error trying to plot',type(myplotfunc),'with value',myplotfunc,
                    '\nlength of the ndarray arguments:',['shape:'+str(shape(j)) if isinstance(j, ndarray) else j for j in plotargs],
                    '\nsizes of ndarray kwargs',dict([(j,shape(kwargs[j])) if isinstance(kwargs[j], ndarray) else (j,kwargs[j]) for j in list(kwargs.keys())]),
                    '\narguments = ',plotargs,
                    '\nkwargs =',kwargs)+explain_error(e))
        if x_inverted:
            these_xlims = ax.get_xlim()
            ax.set_xlim((max(these_xlims),min(these_xlims)))
    #{{{ attach labels and such
    if (myxlabel!=None):
        ax.set_xlabel(myxlabel)
    if (myylabel!=None):
        ax.set_ylabel(myylabel)
    try:
        ax.axis('tight')
    except Exception as e:
        raise Exception(strm('error trying to set axis tight after plot',
            myplotfunc, 'with arguments', plotargs, 'and kwargs', kwargs,
            '\nsizes of arguments:', [shape(j) for j in plotargs],
            '\nsizes of ndarray kwargs:',
            dict([(j, shape(kwargs[j])) for j in
                list(kwargs.keys()) if isinstance(kwargs[j], ndarray)])))
    #grid(True)
    #}}}
    return retval
#}}}
#{{{general functions
def box_muller(length, return_complex=True):
    r'''algorithm to generate normally distributed noise'''
    s1 = rand(length)
    s2 = rand(length)
    n1 = sqrt(-2*log(s1))*cos(2*pi*s2)
    if return_complex:
        n2 = sqrt(-2*log(s1))*sin(2*pi*s2)
        return (n1 + 1j * n2)*0.5
    else:
        return (n1)*0.5
#}}}

#{{{nddata
def dp(number,decimalplaces=2,scientific=False,max_front=3):
    """format out to a certain decimal places, potentially in scientific notation

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
        logger.debug(strm("trying to convert",number,"to scientific notation"))
        tenlog = int(floor(log10(abs(number))))
        number /= 10**tenlog
        fstring = '%0.'+'%d'%decimalplaces+r'f\times 10^{%d}'%tenlog
    else:
        fstring = '%0.'+'%d'%decimalplaces+'f'
        if len(fstring%number) > 1+decimalplaces+max_front:
            return dp(number, decimalplaces=decimalplaces, scientific=True)
    return fstring%number
#}}}
#{{{ concatenate datalist along dimname
def concat(datalist,dimname,chop = False):
    #{{{ allocate a new datalist structure  
    newdimsize = 0
    #print 'DEBUG: type(datalist)',type(datalist)
    try:
        shapes = list(map(ndshape,datalist))
    except Exception as e:
        if not isinstance(datalist, list):
            raise TypeError(strm('You didn\'t pass a list, you passed a',type(datalist)))
        raise RuntimeError(strm('Problem with what you passed to concat, list of types,',
            list(map(type,datalist)))+explain_error(e))
    other_info_out = datalist[0].other_info
    for j in range(0,len(datalist)):
        #{{{ make list for the shape to check, which contains the dimensions we are NOT concatting along
        if dimname in shapes[j].dimlabels:
            newdimsize += shapes[j][dimname]
            shapetocheck = list(shapes[j].shape)
            shapetocheck.pop(shapes[j].axn(dimname))
        else:
            newdimsize += 1
            shapetocheck = list(shapes[j].shape)
        #}}}
        if j == 0:
            shapetocheckagainst = shapetocheck
        else:
            if any(~(array(shapetocheck) == array(shapetocheckagainst))):
                if chop:
                    logger.debug(repr(shapetocheck),lsafen(repr(shapetocheckagainst)))
                    raise ValueError(strm('For item ',j,'in concat, ',
                        shapetocheck,'!=',shapetocheckagainst,
                        'where all the shapes of the things',
                        'you\'re trying to concat are:',
                        shapes))
                else:
                    raise ValueError(strm('For item ',j,'in concat, ',
                        shapetocheck,'!=',shapetocheckagainst,
                        'where all the shapes of the things you\'re trying to concat are:',
                        shapes))
    newdatalist = ndshape(datalist[-1])
    if dimname in newdatalist.dimlabels:
        newdatalist[dimname] = newdimsize
    else:
        newdatalist += ([newdimsize],[dimname])
    #print "DEBUG newdatalist is shaped like",newdatalist
    try:
        newdatalist = newdatalist.alloc()
    except:
        raise ValueError(strm("trying to alloc the newdatalist",newdatalist,
            "created a problem") + explain_error(e))
    if datalist[0].get_error() is not None:
        newdatalist.set_error(zeros(shape(newdatalist.data)))
    #}}}
    #{{{ actually contract the datalist
    newdimsize = 0 # now use it to track to position
    for j in range(0,len(datalist)):
        if dimname in shapes[j].dimlabels:
            newdatalist[dimname,newdimsize:newdimsize+shapes[j][dimname]] = datalist[j]
            newdimsize += shapes[j][dimname]
        else:
            newdatalist[dimname,newdimsize:newdimsize+1] = datalist[j]
            newdimsize += 1
    #}}}
    #{{{ pull the axis labels from the last item in the list
    if len(datalist[-1].axis_coords)>0:
        dimlabels = list(datalist[-1].dimlabels)
        axis_coords = list(datalist[-1].axis_coords)
        #print "axis_coords are",axis_coords,"for",dimlabels
        if dimname in dimlabels:
            thisindex = dimlabels.index(dimname)
            dimlabels.pop(thisindex)
            axis_coords.pop(thisindex)
        dimlabels += [dimname]
        axis_coords += [r_[0:newdimsize]]
        try:
            newdatalist.labels(dimlabels,axis_coords)
        except Exception as e:
            raise ValueError(strm("trying to attach axes of lengths",
                list(map(len,axis_coords)),"to",dimlabels)+explain_error(e))
    #}}}
    newdatalist.other_info = other_info_out
    return newdatalist
#}}}
class nddata (object):
    """This is the detailed API reference.
    For an introduction on how to use ND-Data, see the :ref:`Main ND-Data Documentation <nddata-summary-label>`.
    """
    want_to_prospa_decim_correct = False
    # {{{ initialization
    def __init__(self, *args, **kwargs):
        """initialize nddata -- several options.
        Depending on the information available, one of several formats can be used.

        3 arguments:
            ``nddata(inputarray, shape, dimlabels)``

            :inputarray:
                ndarray storing the data -- note that the size is ignored
                and the data is reshaped as needed
            :shape:
                a list (or array, *etc.*) giving the size of each dimension, in order
            :dimlabels:
                a list giving the names of each dimension, in order
        2 arguments:
            ``nddata(inputarray, dimlabels)``

            :inputarray:
                ndarray storing the data -- the data is *not* reshaped
            :dimlabels:
                a list giving the names of each dimension, in order
        2 arguments:
            ``nddata(inputarray, single_dimlabel)``

            :inputarray:
                ndarray storing the data -- must be 1D  
                inputarray is *also* used to label the single axis
            :single_dimlabel:
                a list giving the name of the single axis
        1 argument:
            ``nddata(inputarray, shape, dimlabels)``

            :inputarray:
                ndarray storing the data -- reduced to 1D  
                A single dimension, called "INDEX" is set.
                This suppresses the printing of axis labels.  
                This is used to store numbers and arrays
                that might have error and units,
                but aren't gridded data.
        keyword args
            these can be used to set the labels, etc, and are passed to :func:`__my_init__`

        """
        logger.debug('called init')
        if len(args) > 1:
            logger.debug('more than one argument')
            if len(args) == 2:
                if len(args[0].shape) == 1 and isinstance(args[1], str):
                    logger.debug('constructing 1D array')
                    self.__my_init__(args[0],[len(args[0])],[args[1]])
                    self.labels(args[1],args[0].copy())# needs to be a copy, or when we write data, we will change the axis
                elif all([isinstance(j, str) for j in args[1]]):
                    logger.debug('passed only axis labels')
                    self.__my_init__(args[0],
                            list(args[0].shape),args[1])
                else:
                    raise ValueError('You can pass two arguments only if you pass a 1d ndarray and a name for the axis') 
            elif len(args) == 3:
                self.__my_init__(args[0],args[1],args[2],**kwargs)
            else:
                raise ValueError(strm("You passed",len(args),"to nddata.  I don't know what to do with this."))
        else:
            logger.debug('only one argument')
            self.__my_init__(args[0],[-1],['INDEX'],**kwargs)
        return
    def __my_init__(self, data, sizes, dimlabels, axis_coords=[],
            ft_start_time=None, data_error=None, axis_coords_error=None,
            axis_coords_units=None, data_units=None, other_info={}):
        if ft_start_time is not None:
            raise ValueError('ft_start_time is obsolete -- you will want to pass a float value to the shift keyword argument of either .ft() or .ift()')
        self.genftpairs = False
        if not (isinstance(data, ndarray)):
            #if (type(data) is float64) or (type(data) is complex128) or (type(data) is list):
            if isscalar(data) or (isinstance(data, list)) or (isinstance(data, tuple)):
                data = array(data)
            else:
                raise TypeError(strm('data is not an array, it\'s',type(data),'!'))
        if not (isinstance(dimlabels, list)):
            raise TypeError(strm('you provided a multi-dimensional ndarray but a set of dimension labels of type',type(dimlabels),"if you want a 1D nddata, give a 1D array, or if you want a ND nddata, give a list of dimensions"))
        try:
            self.data = reshape(data,sizes)
        except:
            try:
                error_string = strm("While initializing nddata, you are trying trying to reshape a",data.shape,"array (",data.size,"data elements) with list of sizes",list(zip(dimlabels,sizes)),"(implying that there are ",prod(sizes),"data elements)")
            except TypeError:
                error_string = strm("While initializing nddata, you are trying trying to reshape a",data.shape,"array (",data.size,"data elements) with list of sizes",sizes)
            raise ValueError(error_string)
        self.dimlabels = dimlabels
        self.axis_coords = axis_coords
        #if len(axis_coords) > 0:
        #    testshape = data.shape
        #    if not all([len(axis_coords[j])==testshape[j] if axis_coords[j] is not None else True for j in range(0,len(axis_coords))]):
        #        raise IndexError('The length of your axis labels (axis_coords) (shape %s) and your axis data (shape %s) does not match!!!'%(repr([len(thiscoord) for thiscoord in axis_coords]),repr(data.shape)))
        self.data_error = data_error
        self.data_units = data_units
        self.other_info = deepcopy(other_info)
        if axis_coords_error is None:
            self.axis_coords_error = [None]*len(axis_coords)
        else:
            self.axis_coords_error = axis_coords_error
        if axis_coords_units is None:
            self.axis_coords_units = [None]*len(axis_coords)
        else:
            self.axis_coords_units = axis_coords_units 
        return
    # }}}
    def _contains_symbolic(self,string):
        return string[:9] == 'symbolic_' and hasattr(self,string)
    #{{{ for printing
    def __repr_pretty__(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))
    def __repr__(self):
        return str(self)
    def __str__(self):
        def show_array(x,indent = ''):
            x = repr(x)
            if x.startswith('array('):
                x = x.split('\n')
                # need to remove the "array(" and aligning spaces
                return ('\n'+indent).join(j[6:-1] for j in x)
            else: return x
        retval = show_array(self.data) 
        retval += '\n\t\t+/-'
        retval += show_array(self.get_error())
        if len(self.dimlabels) > 1 or len(self.dimlabels) == 0 or self.dimlabels[0] != "INDEX":
            retval += '\n\tdimlabels='
            retval += repr(self.dimlabels)
            retval += '\n\taxes='
            def rep_this_dict(starting_indent,thisdict,errordict):
                dictrep = []
                for k,v in thisdict.items():
                    dictrep.append('`'+k+'\':'+show_array(v,starting_indent)+starting_indent+'\t\t+/-'+repr(errordict[k]))
                return '{'+(','+starting_indent+'\t').join(dictrep)+'}' # separate with an extra comma, the existing indent, and a tab
            retval += rep_this_dict('\n\t',self.mkd(self.axis_coords),self.mkd(self.axis_coords_error))
        #retval += '\n\t\t+/-'
        #retval += rep_this_dict('\n\t\t',self.mkd(self.axis_coords_error))
        retval += '\n'
        return retval
    #}}}
    #{{{ for plotting
    def gnuplot_save(self,filename):
        x = self.getaxis(self.dimlabels[0])[:5]
        y = self.getaxis(self.dimlabels[1])[:5]
        z = self.data[:5,:5]
        print("size of x",size(x),"size of y",size(y),"size of z",size(z))
        print("x",x,"y",y,"z",z)
        data = empty((z.shape[0]+1,z.shape[1]+1))
        data[1:,1:] = z[:]
        data[0,0] = z.shape[1]
        data[0,1:] = y.flatten()
        data[1:,0] = x.flatten()
        print("data",data)
        fp = open('auto_figures/'+filename+'.dat','w')
        fp.write(float32(data).tostring())
        fp.write('\n')
        fp.close()
        return
    #{{{ sort and shape the data for 3d plotting
    def sort_and_xy(self):
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if len(self.dimlabels) > 2:
            raise ValueError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        x_axis = self.retaxis(x_dim).data
        y_axis = self.retaxis(y_dim).data
        #}}}
        return x_axis,y_axis
    def matrices_3d(self,also1d = False,invert = False,max_dimsize = 1024,downsample_self = False):
        ''' returns X,Y,Z,x_axis,y_axis
        matrices X,Y,Z, are suitable for a variety of mesh plotting, etc, routines
        x_axis and y_axis are the x and y axes
        '''
        this_size = array(self.data.shape)
        sortedself = self.copy()
        if any(this_size > max_dimsize):
            print(lsafen("Warning! The data is big (%s), so I'm automatically downsampling"%(ndshape(self))))
            for j in where(this_size > max_dimsize):
                downsampling = ceil(double(this_size[j]) / max_dimsize)
                print('downsampling',self.dimlabels[j],'by',downsampling)
                sortedself = sortedself[self.dimlabels[j],0::downsampling]
            print(lsafen("I reduced to a max of max_dimsize = %d so the data is now %s"%(max_dimsize,ndshape(sortedself))))

        x_axis,y_axis = sortedself.sort_and_xy()
        if invert:
            print("trying to invert meshplot-like data")
        X = x_axis*ones(shape(y_axis))
        Y = ones(shape(x_axis))*y_axis
        Z = real(sortedself.data)
        if invert:
            X = X[:,::-1]
            Y = Y[:,::-1]
            Z = Z[:,::-1]
        if downsample_self:
            self.data = sortedself.data
            self.setaxis(self.dimlabels[0],x_axis)
            self.setaxis(self.dimlabels[1],y_axis)
        if also1d:
            if invert:
                return X,Y,Z,x_axis[::-1],y_axis[::-1]
            else:
                return X,Y,Z,x_axis,y_axis
        else:
            return X,Y,Z
    #}}}
    def mayavi_surf(self):
        """use the mayavi surf function, assuming that we've already loaded mlab
        during initialization"""
        X,Y,Z = self.matrices_3d()
        s = self.mlab.surf(X,Y,Z)
        return s
    #{{{ 3D mesh plot
    def meshplot(self,stride = None,alpha = 1.0,onlycolor = False,light = None,rotation = None,cmap = cm.gray,ax = None,invert = False,**kwargs):
        r'''takes both rotation and light as elevation, azimuth
        only use the light kwarg to generate a black and white shading display'''
        X,Y,Z = self.matrices_3d()
        if light == True:
            light = [0,0]# I think this is 45 degrees up shining down from the left of the y axis
        if not onlycolor:
            if ax is None: 
                ax = self._init_3d_axis(ax,rotation = rotation)
            else:
                if rotation is not None:
                    raise ValueError("you can only set the rotation once! (you tried"+repr(rotation)+")")
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
            ls = LightSource(azdeg = light[1],altdeg = light[0])
            if cmap is not None:
                rgb = ls.shade(Z,cmap)
        else:
            mask = isfinite(Z.flatten())
            for_rgb = Z-Z.flatten()[mask].min()
            for_rgb /= for_rgb.flatten()[mask].max()
            if cmap is not None:
                rgb = cmap(for_rgb)
        if onlycolor:
            imshow(rgb)
        else:
            if light is None:
                if cmap is not None:
                    kwargs.update({'cmap':cmap})
                ax.plot_surface(X,Y,Z,
                        rstride = rstride,
                        cstride = cstride,
                        shade = True,
                        **kwargs)
            else:
                newkwargs = {}
                newkwargs['linewidth'] = 0.0
                newkwargs.update(kwargs)
                if cmap is not None:
                    newkwargs['facecolors'] = rgb
                ax.plot_surface(X,Y,Z,
                        rstride = rstride,
                        cstride = cstride,
                        alpha = alpha,
                        shade = False,
                        **newkwargs)
            ax.set_xlabel(x_dim)
            ax.set_ylabel(y_dim)
            ax.set_zlabel(self.name())
        if onlycolor:
            return
        else:
            return ax
    def contour(self,labels = True,**kwargs):
        """Contour plot -- kwargs are passed to the matplotlib
        `contour` function.

        See docstring of `figlist_var.image()` for an example

        Attributes
        ----------
        labels : boolean
            Whether or not the levels should be labeled.
            Defaults to True
        """
        x_axis,y_axis = self.dimlabels
        x = self.getaxis(x_axis)[:,None]
        y = self.getaxis(y_axis)[None,:]
        if 'levels' not in list(kwargs.keys()):
            levels = r_[self.data.min():self.data.max():30j]
        cs = contour(x*ones_like(y),ones_like(x)*y,self.data,**kwargs)
        if labels:
            clabel(cs,inline = 1,fontsize = 10)
        xlabel(self.unitify_axis(x_axis))
        ylabel(self.unitify_axis(y_axis))
        return cs
    def waterfall(self,alpha = 0.3,ax = None,rotation = None,color = 'b',edgecolor = 'k'):
        if ax is None: 
            ax = self._init_3d_axis(ax,rotation = rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        if len(self.dimlabels) > 2:
            raise ValueError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        try:
            x_axis = self.retaxis(x_dim).data
        except Exception as e:
            raise ValueError(strm('trying to get the info on axis', x_dim, 'which is',
                self.getaxis(x_dim))
                +explain_error(e))
        y_axis = self.retaxis(y_dim).data
        #}}}
        ax.set_xlabel(self.unitify_axis(x_dim))
        ax.set_ylabel(self.unitify_axis(y_dim))
        ax.set_zlabel(self.unitify_axis(self.name(),is_axis = False))
        verts = []
        xs = x_axis.flatten()
        xs = r_[xs[0],xs,xs[-1]] # add points for the bottoms of the vertices
        ys = y_axis.flatten()
        for j in range(0,len(ys)):
            zs = self[y_dim,j].data.flatten()
            zs = r_[0,zs,0]
            verts.append(list(zip(xs,zs))) # one of the faces
        poly = PolyCollection(verts, facecolors = [color]*len(verts), edgecolors = edgecolor) # the individual facecolors would go here
        poly.set_alpha(alpha)
        fig = gcf()
        #ax = fig.add_subplot(111,projection = '3d')
        ax.add_collection3d(poly,zs = ys, zdir = 'y')
        ax.set_zlim3d(self.data.min(),self.data.max())
        ax.set_xlim3d(xs.min(),xs.max())
        ax.set_ylim3d(ys.min(),ys.max())
        return ax
    def _init_3d_axis(self,ax,rotation = None):
        # other things that should work don't work correctly, so use this to initialize the 3D axis
        #ax.view_init(elev = rotation[0],azim = rotation[1])
        if rotation is None:
            rotation = [0,0]
        if ax is None:
            fig = gcf()
            ax = axes3d.Axes3D(fig)
            print("I'm trying to rotate to",rotation)
            #ax.view_init(20,-120)
            #ax.view_init(elev = 20 + rotation[1],azim = -120 + rotation[0])
            ax.view_init(azim = rotation[0],elev = rotation[1])
        return ax
    def oldtimey(self,alpha = 0.5,ax = None,linewidth = None,sclinewidth = 20.,light = True,rotation = None,invert = False,**kwargs):
        sortedself = self.copy()
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if invert:
            print("trying to invert oldtimey")
        if linewidth is None:
            linewidth = sclinewidth/sortedself.data.shape[1]
            print("setting linewidth to %0.1f"%linewidth)
        if ax is None: 
            ax = sortedself._init_3d_axis(ax,rotation = rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        ax = sortedself.meshplot(linewidth = 0,light = light,ax = ax,invert = invert)
        #return
        if len(sortedself.dimlabels) > 2:
            raise ValueError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = sortedself.dimlabels[0]
        y_dim = sortedself.dimlabels[1]
        x_axis = sortedself.retaxis(x_dim).data
        y_axis = sortedself.retaxis(y_dim).data
        #}}}
        verts = []
        xs = x_axis.flatten()
        ys = y_axis.flatten() # this is the depth dimension
        if invert:
            ys = ys[::-1]
        for j in range(0,len(ys)):
            zs = sortedself[y_dim,j].data.flatten() # pulls the data (zs) for a specific y slice
            if invert:
                zs = zs[::-1]
            ax.plot(xs,ones(len(xs))*ys[j],zs,'k',linewidth = linewidth)
        fig = gcf()
        ax.set_zlim3d(sortedself.data.min(),sortedself.data.max())
        ax.set_xlim3d(xs.min(),xs.max())
        #if invert:
        #    ax.set_ylim3d(ys.max(),ys.min())
        #else:
        ax.set_ylim3d(ys.min(),ys.max())
        return ax
    #}}}
    #}}}
    #{{{ error-related functions
    def normalize(self,axis,first_figure = None):#,whichpoint = slice(0,1,None)):
        x = self.data
        n = len(x)
        S = sparse.lil_matrix((n,n))
        S.setdiag((self.get_error())**2)
        self.set_error(None)
        first_point = self[axis,0:1].copy() # this makes another instance that contains just the first point, for error propagation
        B = sparse.lil_matrix((n,n))
        B.setdiag(1./x)
        B[0,:] = -x/(x[0]**2) # Sparse seems to only support row assignment, so make the transpose to give it what it wants
        B[0,0] = 0.0
        B = B.T
        E = B * S * (B.T) # verified that this is matrix multiplication
        self /= first_point # this gives the experimentally measured E
        #{{{ now, chop out the first point, which is meaningless
        self = self[axis,1:]
        E = E[1:,1:]
        #}}}
        self.set_error(sqrt(E.diagonal()))
        E.setdiag(zeros(n-1))
        self.data_covariance = E
        return self
    def get_covariance(self):
        '''this returns the covariance matrix of the data'''
        if hasattr(self,'data_covariance'):
            E = self.data_covariance.copy()
        else:
            n = size(self.data)
            E = sparse.lil_matrix((n,n))
        try:
            E.setdiag(self.get_error()**2)
        except Exception as e:
            raise ValueError(strm('Problem getting covariance because error is',self.get_error())+explain_error(e))
        return E.toarray()
    #}}}
    #{{{ shortcuts for axes
    def axlen(self,axis):
        r"""return the size (length) of an axis, by name
        
        Parameters
        ----------

        axis: str
            name of the axis whos length you are interested in
        """
        return shape(self.data)[self.axn(axis)]
    def axn(self,axis):
        r'''Return the index number for the axis with the name "axis"

        This is used by many other methods.
        As a simple example,
        self.:func:`axlen`(axis) (the axis length) returns
        ``shape(self.data)[self.axn(axis)]``

        Parameters
        ----------

        axis: str
            name of the axis
        '''
        try:
            return self.dimlabels.index(axis)
        except:
            raise ValueError(' '.join(map(repr,['there is no axis named',axis,'all axes are named',self.dimlabels])))
    def indices(self,axis_name,values):
        r'Return a string of indeces that most closely match the axis labels corresponding to values. Filter them to make sure they are unique.'
        x = self.getaxis(axis_name)
        retval = []
        for j in values:
            retval.append(argmin(abs(x - j)))
        retval = array(retval)
        return unique(retval)
    #}}}
    #{{{ dictionary functions -- these convert between two formats:
    # dictionary -- stuff labeled according the dimension label.
    # list -- same information, but it's assumed they are listed in the order given by "dimlabels"
    def mkd(self,*arg,**kwargs):
        'make dictionary format'
        give_None = process_kwargs([("give_None",True)],kwargs)
        if len(arg) == 1:
            input_list = arg[0]
            if emptytest(input_list):
                return dict(zip(self.dimlabels,
                    [None]*len(self.dimlabels)))
            if len(input_list) != len(self.dimlabels):
                print(r"{\color{red}WARNING! mkd error (John will fix this later):}")
                print("When making a dictionary with mkd, you must pass a list that has one element for each dimension!  dimlabels is "+repr(self.dimlabels)+" and you passed "+repr(arg)+'\n\n')
                raise ValueError("When making a dictionary with mkd, you must pass a list that has one element for each dimension!  dimlabels is "+repr(self.dimlabels)+" and you passed "+repr(arg))
            for i,v in enumerate(input_list):
                if isinstance(v, ndarray):
                    if v.shape == () and v.size == 0:
                        input_list[i] = None
                    if v.dtype.type in [str_, bytes_]:
                        input_list[i] = str(v)
            if give_None:
                return dict(zip(self.dimlabels,input_list))
            else:
                #{{{ don't return values for the things that are None
                mykeys = [self.dimlabels[j] for j in range(0,len(self.dimlabels)) if input_list[j] is not None]
                myvals = [input_list[j] for j in range(0,len(self.dimlabels)) if input_list[j] is not None]
                return dict(zip(mykeys,myvals))
                #}}}
        elif len(arg) == 0:
            if not give_None:
                raise ValueError("You can't tell me not to give none and then not pass me anything!!")
            return dict(zip(self.dimlabels,
                [None]*len(self.dimlabels)))
        else:
            raise ValueError(strm('.mkd() doesn\'t know what to do with %d arguments',len(arg)))
    def fld(self,dict_in,noscalar = False):
        'flatten dictionary -- return list'
        return [dict_in[x] for x in self.dimlabels]
    #}}}
    #{{{ set + get the error + units
    #{{{ set units
    def set_units(self,*args):
        if len(args) == 2:
            unitval = args[1] # later, have some type of processing bojive
            if self.axis_coords_units is None or len(self.axis_coords_units) == 0:
                self.axis_coords_units = [None] * len(self.dimlabels)
            self.axis_coords_units[self.axn(args[0])] = unitval
        elif len(args) == 1:
            unitval = args[0] # later, have some type of processing bojive
            self.data_units = unitval
        else:
            raise TypeError(".set_units() takes data units or 'axis' and axis units")
        return self
    def human_units(self):
        prev_label = self.get_units()
        # -- rescaling for y axis seems screwed up, so
        # just skip it
        #if prev_label is not None and len(prev_label)>0:
        #    #{{{ find the average order of magnitude, rounded down to the nearest power of 3
        #    average_oom = log10(abs(self.data))/3.
        #    average_oom = average_oom[isfinite(average_oom)].mean()
        #    #}}}
        #    logger.debug(strm("(human units): for data the average oom is",average_oom*3))
        #    if round(average_oom) == 0.0:
        #        average_oom = 0
        #    else:
        #        average_oom = 3*floor(average_oom)
        #    logger.debug(strm("(human units): for data I round this to",average_oom))
        #    this_str = apply_oom(average_oom,self.data,prev_label=prev_label) 
        #    self.set_units(this_str)
        #else:
        #    logger.debug(strm('data does not have a unit label'))
        for thisaxis in self.dimlabels:
            prev_label = self.get_units(thisaxis)
            if prev_label is not None and len(prev_label)>0:
                data_to_test = self.getaxis(thisaxis)
                logger.debug(strm("the axis",thisaxis,"looks like this:",data_to_test))
                if data_to_test is not None:
                    try:
                        data_to_test = data_to_test[isfinite(data_to_test)]
                    except:
                        raise ValueError(strm('data_to_test is',data_to_test,'isfinite is',isfinite(data_to_test)))
                    if len(data_to_test) == 0:
                        raise ValueError(strm("Your",thisaxis,"axis doesn't seem to have any sensible values!"))
                    #{{{ find the average order of magnitude, rounded down to the nearest power of 3

                    average_oom = log10(abs(data_to_test))/3.
                    logger.debug(strm("for axis: dtype",data_to_test.dtype))
                    logger.debug(strm("for axis: dtype",data_to_test))
                    logger.debug(strm("for axis: oom:",average_oom))
                    average_oom = average_oom[isfinite(average_oom)].mean()
                    #}}}
                    logger.debug(strm("for axis",thisaxis,"the average oom is",average_oom*3))
                    average_oom = 3*floor(average_oom)
                    logger.debug(strm("for axis",thisaxis,"I round this to",average_oom))
                    x = self.getaxis(thisaxis)
                    result_label = apply_oom(average_oom,x,prev_label=prev_label)
                    self.set_units(thisaxis,result_label)
                else:
                    logger.debug(strm(thisaxis,'does not have an axis label'))
            else:
                logger.debug(strm(thisaxis,'does not have a unit label'))
        return self
    #}}}
    #{{{ get units
    def units_texsafe(self,*args):
        retval = self.get_units(*args)
        if retval is None:
            return None
        if retval.find('\\') > -1:
            retval = '$'+retval+'$'
        return retval
    def replicate_units(self,other):
        for thisaxis in self.dimlabels:
            if other.get_units(thisaxis) is not None:
                self.set_units(thisaxis,other.get_units(thisaxis))
        if other.get_units() is not None:
            self.set_units(other.get_units(thisaxis))
        return self
    def get_units(self,*args):
        if len(args) == 1:
            if self.axis_coords_units is None:
                return None
            if len(self.axis_coords_units) == 0:
                return None
            try:
                return self.axis_coords_units[self.axn(args[0])]
            except:
                raise RuntimeError(strm('problem getting units for',args[0],'dimension',self.dimlabels,self.axis_coords_units))
        elif len(args) == 0:
            return self.data_units
        else:
            raise ValueError(".set_units() takes axis or nothing")
    #}}}
    #{{{ set error
    def set_error(self,*args):
        r'''set the errors: either
        
        `set_error('axisname',error_for_axis)` or `set_error(error_for_data)`

        `error_for_data` can be a scalar, in which case, **all** the data errors are set to `error_for_data`

        .. todo::
                several options below -- enumerate them in the documentation
        '''
        if (len(args) == 1) and isscalar(args[0]):
            if args[0] == 0:
                args = (zeros_like(self.data),)
            else:
                args = (ones_like(self.data) * args[0],)
        if (len(args) == 1) and (isinstance(args[0], ndarray)):
            self.data_error = reshape(args[0],shape(self.data))
        elif (len(args) == 1) and (isinstance(args[0], list)):
            self.data_error = reshape(array(args[0]),shape(self.data))
        elif (len(args) == 2) and (isinstance(args[0], str)) and (isinstance(args[1], ndarray)):
            self.axis_coords_error[self.axn(args[0])] = args[1]
        elif (len(args) == 2) and (isinstance(args[0], str)) and (isscalar(args[1])):
            self.axis_coords_error[self.axn(args[0])] = args[1]*ones_like(self.getaxis(args[0]))
        elif (len(args) == 1) and args[0] is None:
            self.data_error = None
        else:
            raise TypeError(' '.join(map(repr,['Not a valid argument to set_error:',list(map(type,args))])))
        return self
    #}}}
    #{{{ random mask -- throw out points
    def random_mask(self,axisname,threshold = exp(-1.0),inversion = False):
        r'''generate a random mask with about 'threshold' of the points thrown out'''
        if inversion:
            threshold = threshold / (1.0 - threshold)
        myr = rand(self.data.shape[self.axn(axisname)]) # random array same length as the axis
        return myr > threshold
    #}}}
    #{{{ get error
    def get_error(self,*args):
        '''get a copy of the errors\neither set_error('axisname',error_for_axis) or set_error(error_for_data)'''
        if (len(args) == 0):
            if self.data_error is None:
                return None
            else:
                return real(self.data_error)
        elif (len(args) == 1):
            thearg = args[0]
            if isinstance(thearg, str_):
                thearg = str(thearg) # like in the other spot, this became necessary with some upgrade, though I'm not sure that I should maybe just change the error functions to treat the numpy string in the same way
            if (isinstance(thearg, str)):
                if len(self.axis_coords_error) == 0: self.axis_coords_error = [None] * len(self.dimlabels) # is we have an empty axis_coords_error, need to fill with None's
                try:
                    errorforthisaxis = self.axis_coords_error[self.axn(thearg)]
                except Exception as e:
                    raise RuntimeError(strm('Problem trying to load error',self.axn(thearg),'for axis',thearg,'out of',self.axis_coords_error)
                            +explain_error(e))
                if errorforthisaxis is None:
                    return None
                else:
                    x = self.axis_coords_error[self.axn(thearg)]
                    if isinstance(x, ndarray):
                        if x.shape == ():
                            return None
                        else:
                            return real(self.axis_coords_error[self.axn(thearg)])
                    else:
                        return real(self.axis_coords_error[self.axn(thearg)])
        else:
            raise ValueError(strm('Not a valid argument to get_error: *args=',args,'map(type,args)=',list(map(type,args))))
        #}}}
    #}}}
    #{{{ match dims --
    def matchdims(self,other):
        r'add any dimensions to self that are not present in other'
        #print 'diagnose: matching',ndshape(self),'to',ndshape(other)
        addeddims =  list(set(self.dimlabels)^set(other.dimlabels))
        newdims = addeddims + self.dimlabels
        newshape = [1]*len(addeddims) + list(self.data.shape)
        #print 'diagnose: newshape',newshape,'newdims',newdims
        #{{{ reshape to the new dimensions  
        new_axis_coords = [r_[1]]*len(addeddims) + self.axis_coords
        self.data = self.data.reshape(newshape)
        self.dimlabels = newdims
        if len(self.axis_coords)>0:
            self.axis_coords = new_axis_coords
        #}}}
        #{{{ if we are adding dimensions, we will need to reorder to match the order of the other   
        if len(addeddims)>0:
            self.reorder(other.dimlabels)
        #}}}
        return self
    #}}}
    #{{{ rename
    def rename(self,previous,new):
        self.dimlabels[self.dimlabels.index(previous)] = new
        return self
    #}}}
    #{{{ display and other properties
    #{{{ set and get prop
    def unset_prop(self,arg):
        "remove a 'property'"
        self.other_info.pop(arg)
        if len(self.other_info) == 0:
            del self.other_info
        return self
    def set_prop(self,*args):
        r"""set a 'property' of the nddata
        This is where you can put all unstructured information (e.g. experimental parameters, etc)"""
        if len(args) == 2:
            propname,val = args
            self.other_info.update({propname:val})
        elif len(args) == 1 and isinstance(args[0], dict):
            self.other_info.update(args[0])
        else:
            raise ValueError("I don't know what you're passing to set prop!!!")
        return self
    def copy_props(self,other):
        r"""Copy all properties (see :func:`get_prop`) from another nddata
        object -- note that these include properties pertaining the the FT
        status of various dimensions."""
        self.other_info.update(deepcopy(other.other_info))
        return self
    def get_prop(self,propname=None):
        r'''return arbitrary ND-data properties (typically acquisition parameters *etc.*) by name (`propname`)
        
        In order to allow ND-data to store acquisition parameters and other info that accompanies the data,
        but might not be structured in a gridded format, nddata instances
        always have a `other_info` dictionary attribute,
        which stores these properties by name.

        If the property doesn't exist, this returns `None`.
        
        Parameters
        ----------
        propname: str
            Name of the property that you're want returned.
            If this is left out or set to "None" (not given), the names of the available
            properties are returned.
            If no exact match is found, and propname contains a . or * or [, it's
            assumed to be a regular expression.
            If several such matches are found, the error message is informative.

            .. todo::
                have it recursively search dictionaries (e.g. bruker acq)

        Returns
        -------
        The value of the property (can by any type) or `None` if the property doesn't exist.
        '''
        if propname is None:
            return self.other_info.keys()
        if propname not in self.other_info.keys():
            if '.' in propname or '*' in propname or '[' in propname:
                propname_re = re.compile(propname)
                matches = [j for j in self.other_info.keys() if propname_re.match(j)]
                if len(matches) == 0:
                    return None
                assert len(matches) == 1, "I found %d matches for regexp %s in properties: %s"%(len(matches),
                        propname,
                        ' '.join(matches))
                return self.other_info[matches[0]]
            else:
                return None
        return self.other_info[propname]
    def name(self,*arg):
        r"""args:
           .name(newname) --> Name the object (for storage, etc)
           .name() --> Return the name"""
        if len(arg) == 1:
            self.set_prop('name',arg[0])
            return self
        elif len(arg) == 0:
            return self.get_prop('name')
        else:
            raise ValueError("invalid number of arguments")
    #}}}
    #{{{ set and get plot color
    def set_plot_color(self,thiscolor):
        if thiscolor is None:
            return
        if thiscolor is str:
            colordict = {'r':[1,0,0],
                    'g':[0,1,0],
                    'b':[0,0,1],
                    'k':[0,0,0],
                    'y':[0.5,0.5,0],
                    'o':[0.75,0.25,0],
                    'c':[0,0.5,0.5]}
            try:
                thiscolor = colordict[thiscolor]
            except:
                raise ValueError(strm('Color',thiscolor,'not in dictionary'))
        self.other_info.update({'plot_color':thiscolor})
        return
    def get_plot_color(self):
        if 'plot_color' in self.get_prop():
            return self.other_info['plot_color']
        else:
            return None
    #}}}
    #}}}
    #{{{ arithmetic
    def along(self,dimname):
        """Specifies the dimension for the next matrix
        multiplication (represents the rows/columns)."""
        self._matmul_along = dimname
        return self
    def along(self,dimname):
        """Specifies the dimension for the next matrix
        multiplication (represents the rows/columns)."""
        self._matmul_along = dimname
        return self
    #@profile
    def dot(self,arg):
        """Tensor dot of self with arg -- dot all matching dimension labels.  This can be used to do matrix multiplication, but note that the order of doesn't matter, since the dimensions that are contracted are determined by matching the dimension names, not the order of the dimension.

        >>> a = nddata(r_[0:9],[3,3],['a','b'])
        >>> b = nddata(r_[0:3],'b')
        >>> print a.C.dot(b)
        >>> print a.data.dot(b.data)
        >>> a = nddata(r_[0:27],[3,3,3],['a','b','c'])
        >>> b = nddata(r_[0:9],[3,3],['a','b'])
        >>> print a.C.dot(b)
        >>> print tensordot(a.data,b.data,axes=((0,1),(0,1)))

        >>> a = nddata(r_[0:27],[3,3,3],['a','b','c'])
        >>> b = nddata(r_[0:9],[3,3],['a','d'])
        >>> print a.C.dot(b)
        >>> print tensordot(a.data,b.data,axes=((0),(0)))
        """
        A,B = self.aligndata(arg)
        matching_dims = list(set(self.dimlabels) & set(arg.dimlabels))
        assert len(matching_dims) > 0, "no matching dimensions!"
        # {{{ store the dictionaries for later use
        axis_coords_dict = A.mkd(A.axis_coords)
        axis_units_dict = A.mkd(A.axis_coords_units)
        axis_coords_error_dict = A.mkd(A.axis_coords_error)
        # }}}
        # manipulate "self" directly
        self.dimlabels = [j for j in A.dimlabels if j not in matching_dims]
        match_idx = [A.axn(j) for j in matching_dims]
        if (self.get_error() is not None) or (arg.get_error() is not None):
            raise ValueError("we plan to include error propagation here, but not yet provided")
        self.data = tensordot(A.data,B.data,axes=(match_idx,match_idx))
        logger.debug(strm("shape of A is",ndshape(A)))
        logger.debug(strm("shape of B is",ndshape(B)))
        logger.debug(strm("matching_dims are",matching_dims))
        newsize = [(A.data.shape[j] if A.data.shape[j] != 1 else B.data.shape[j])
                for j in range(len(A.data.shape)) if A.dimlabels[j] not in matching_dims]
        self.data = self.data.reshape(newsize)
        # {{{ use the dictionaries to reconstruct the metadata
        self.axis_coords = self.fld(axis_coords_dict)
        self.axis_coords_units = self.fld(axis_units_dict)
        self.axis_coords_error = self.fld(axis_coords_error_dict)
        # }}}
        return self
    def __add__(self,arg):
        if isscalar(arg):
            A = self.copy()
            if isinstance(arg, complex) and self.data.dtype not in [complex128,complex64]:
                A.data = complex128(A.data)
            A.data += arg
            # error does not change
            return A
        #{{{ shape and add
        A,B = self.aligndata(arg)
        logger.debug(strm('after alignment, right data looks like:',ndshape(B)))
        retval = A.copy()
        retval.data = A.data + B.data
        #}}}
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0
        if Aerr is not None:
            Rerr += (Aerr)**2
        if Berr is not None:
            Rerr += (Berr)**2
        Rerr = sqrt(real(Rerr)) # convert back to stdev
        if Aerr is None and Berr is None:
            Rerr = None
        retval.set_error(Rerr)
        return retval
    def __sub__(self,arg):
        return self.__add__(-1*arg)
    def __lt__(self,arg):
        if isinstance(arg, ndarray):
            retval = self.copy()
            retval.data = retval.data < arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data < B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data < arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    def __gt__(self,arg):
        if isinstance(arg, ndarray):
            retval = self.copy()
            retval.data = retval.data > arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data > B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data > arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    def __le__(self,arg):
        if isinstance(arg, ndarray):
            retval = self.copy()
            retval.data = retval.data <= arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data <= B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data <= arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    def __ge__(self,arg):
        if isinstance(arg, ndarray):
            retval = self.copy()
            retval.data = retval.data >= arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data >= B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data >= arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    #@profile
    def __matmul__(self,arg):
        assert type(arg) is nddata, "currently matrix multiplication only allowed if both are nddata"
        return self.C.dot(arg)
    def __mul__(self,arg):
        #{{{ do scalar multiplication
        if isscalar(arg):
            #print "multiplying",self.data.dtype,"with scalar of type",type(arg)
            A = self.copy()
            if isinstance(arg, complex) and self.data.dtype not in [complex128,complex64]:
                A.data = complex128(A.data)
            A.data *= arg
            if A.get_error() is not None:
                error = A.get_error()
                error *= abs(arg)
            return A
        #}}}
        #{{{ shape and multiply
        try:
            A,B = self.aligndata(arg)
        except Exception as e:
            if arg.name() is not None and self.name() is not None:
                raise ValueError(strm("Error aligning right (arg)", arg.name(),
                    "with left (self)", self.name())+explain_error(e))
            else:
                raise ValueError("Error aligning"+explain_error(e))
        retval = A.copy()
        retval.data = A.data * B.data
        #}}}
        #{{{ if we have error for both the sets of data, I should propagate that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0 # we can have error on one or both, so we're going to need to add up the variances
        if Aerr is not None:
            Rerr += (Aerr * B.data)**2
        if Berr is not None:
            Rerr += (Berr * A.data)**2
        Rerr = sqrt(real(Rerr)) # convert back to stdev
        if Aerr is None and Berr is None:
            Rerr = None
        #}}}
        retval.set_error(Rerr)
        return retval
    def __rpow__(self,arg):
        result = self.copy()
        result.set_error(None)
        logger.debug("error propagation for right power not currently supported (do you need this, really?)")
        assert isscalar(arg) or isinstance(arg, ndarray), "currently right power only supported for ndarray and scalars -- do you really need something else??"
        result.data = arg**self.data
        return result
    def __pow__(self,arg):
        if arg == -1:
            x = self.get_error()
            result = self.copy()
            result.data = 1.0/result.data
            if x is not None:
                result.set_error(abs(x.copy()/(self.data**2)))
            return result
        elif arg == 2:
            return self * self
        else:
            if self.get_error() is not None:
                raise ValueError(strm("nothing but -1 and 2 supported yet! (you tried to raise to a power of "+repr(arg)+")"))
            else:
                result = self.copy()
                result.data = result.data**arg
                return result
    def __truediv__(self,arg):
        return self.__div__(arg)
    def __div__(self,arg):
        if isscalar(arg):
            A = self.copy()
            A.data /= arg
            if A.get_error() is not None:
                error = A.get_error()
                error /= abs(arg)
            return A
        A,B = self.aligndata(arg)
        retval = A.copy()
        retval.data = A.data / B.data
        #{{{ if we have error for both the sets of data, I should propagate that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0 # we can have error on one or both, so we're going to need to add up the variances
        dt128 = dtype('complex128')
        if Aerr is not None:
            if (A.data.dtype is dt128) or (B.data.dtype is dt128):# this should avoid the error that Ryan gets
                Rerr += (complex128(Aerr)/complex128(B.data))**2
            else:
                Rerr += (Aerr/B.data)**2
        if Berr is not None:
            if (A.data.dtype is dt128) or (Berr.dtype is dt128) or (B.data.dtype is dt128):# this should avoid the error that Ryan gets
                Rerr += (complex128(A.data)*complex128(Berr)/(complex128(B.data)**2))**2
            else:
                try:
                    Rerr += (A.data*Berr/(B.data**2))**2
                except:
                    raise ValueError(strm('self was',self,
                        'arg was',arg,
                        'dtype of A.data',A.data.dtype,
                        'dtype of Berr',Berr.dtype,
                        'dtype of B.data',Berr) + explain_error(e))
        try:
            Rerr = sqrt(real(Rerr)) # convert back to stdev --> note that this has problems with complex numbers, hence the "abs" above
        except AttributeError as e:
            raise AttributeError(strm("Rerr gave an attribute error when you passed",Rerr) + explain_error(e))
        #print "DEBUG: step 3",Rerr
        #print "Rerr dtype",Rerr.dtype
        if Aerr is None and Berr is None:
            Rerr = None
        #}}}
        retval.set_error(Rerr)
        return retval
    def __invert__(self):
        if self.data.dtype is dtype('bool'):
            self.data = ~self.data
            return self
        else:
            raise ValueError('invert only implemented for boolean now')
    def __abs__(self):
        return self.runcopy(abs)
    __radd__ = __add__
    __rmul__ = __mul__
    def __rsub__(self,arg):
        return -1*(self-arg)
    def __neg__(self):
        return -1*self
    def __rdiv__(self,arg):
        return arg * (self**(-1))
    #def real(self):
    #    self.data = real(self.data)
    #    return self
    #}}}
    #{{{ align data
    def aligndata(self,arg):
        r'''This is a fundamental method used by all of the arithmetic operations.
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
        arg : nddata or ndarray
            The nddata that you want to align to `self`.
            If arg is an ndarray, it will try to match dimensions to self based
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
        '''
        #{{{ if zero dimensional, fake a singleton dimension and recurse
        #{{{ unless both are zero dimensional, in which case, just leave alone
        logger.debug(strm("starting aligndata"))
        if isscalar(arg) or isinstance(arg, ndarray):
            arg = nddata(arg)
            index_dims = [j for j in r_[0:len(arg.dimlabels)]
                     if arg.dimlabels[j]=='INDEX']
            for j in index_dims:# find dimension of matching length
                match_dims = nonzero(arg.data.shape[j]==array(self.data.shape))[0]
                if len(match_dims) > 0:
                    arg.dimlabels[j] = self.dimlabels[match_dims[0]]
                if len(match_dims) != len(index_dims):
                    raise ValueError("you seem to by multiplying by something with an 'INDEX' data and something that doesn't have that -- is this really what you want?  (this is commonly produced by multiplying a mismatched ndarray by an nddata)")
        if ndshape(self).zero_dimensional and ndshape(arg).zero_dimensional:
            logger.debug(strm("(1) yes, I found something zero dimensional"))
            return self.copy(),arg.copy()
        #}}}
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
        #}}}
        selfout = self.copy() # copy self
        assert len(selfout.data.shape) != 0 and len(arg.data.shape) != 0, ("neither"
         " self nor arg should be zero dimensional at this point (previous code"
         " should have taken care of that")
        # {{{create newdims, consisting of dimlabels for self, followed by the
        # names of the dimensions in arg that are not also in self -- order for
        # both is important; then create a matching selfshape
        augmentdims = [x for x in arg.dimlabels if x in
                set(self.dimlabels)^set(arg.dimlabels)] # dims in arg
        #                   but not self, ordered as they were in arg
        newdims = self.dimlabels + augmentdims
        selfshape = list(selfout.data.shape)+list(
                ones(len(augmentdims),dtype=uint64)) # there is no need to
        #       transpose self, since its order is preserved
        # }}}
        argout = arg.copy()
        # {{{ now create argshape for the reshaped argument
        new_arg_labels = [x for x in newdims if x in
                arg.dimlabels] #  only the labels valid for arg, ordered
        #                         as they are in newdims
        argshape = list(ones(len(newdims), dtype=int64))# should be a better solution
        logger.debug(strm("DEBUG 2: shape of self",ndshape(self),"self data shape",self.data.shape,"shape of arg",ndshape(arg),"arg data shape",arg.data.shape))
        logger.debug(strm("DEBUG 3: shape of selfout",ndshape(selfout),"selfout data shape",selfout.data.shape,"shape of argout",ndshape(argout),"argout data shape",argout.data.shape))
        #{{{ wherever the dimension already exists in arg, pull the shape from arg
        for j,k in enumerate(newdims):
            if k in argout.dimlabels:
                try:
                    argshape[j] = argout.data.shape[argout.axn(k)]
                except:
                    raise ValueError("There seems to be a problem because the" +
                            "shape of argout is now len:%d"%len(argout.data.shape),
                            argout.data.shape,"while the dimlabels is len:%d"%len(
                                argout.dimlabels),argout.dimlabels)
        # }}}
        # }}}
        # {{{ transpose arg to match newshape
        argorder = list(map(argout.dimlabels.index,new_arg_labels)) # for
        #          each new dimension, determine the position of the
        #          original dimension
        selfout.data = selfout.data.reshape(int64(selfshape)) # and reshape
        #          to its new shape
        selfout.dimlabels = newdims
        try:
            argshape = int64(argshape)
            argout.data = argout.data.transpose(argorder
                    ).reshape(argshape) # and reshape the data
        except ValueError as Argument:
            raise ValueError('the shape of the data is ' +
                    repr(argout.data.shape) + ' the transpose ' +
                    repr(argorder) + ' and the new shape ' +
                    repr(argshape) + ' original arg: ' +
                    repr(Argument))
        argout.dimlabels = newdims
        # }}}
        # {{{ transpose the data errors appropriately
        if selfout.get_error() is not None:
            try:
                temp = selfout.get_error().copy().reshape(selfshape)
            except ValueError as Argument:
                raise ValueError("The instance (selfout) has a shape of "
                        + repr(selfout.data.shape) +
                        " but its error has a shape of" +
                        repr(selfout.get_error().shape) +
                        "!!!\n\n(original argument:\n" +
                        repr(Argument) + "\n)")
            selfout.set_error(temp)
        if argout.get_error() is not None:
            try:
                temp = argout.get_error().copy().transpose(argorder).reshape(argshape)
            except ValueError as Argument:
                raise ValueError("The argument (argout) has a shape of "
                        + repr(argout.data.shape)
                        + " but its error has a shape of" +
                        repr(argout.get_error().shape) + "(it's " +
                        repr(argout.get_error()) +
                        ")!!!\n\n(original argument:\n" +
                        repr(Argument) + "\n)")
            argout.set_error(temp)
        # }}}
        if (len(selfout.axis_coords)>0) or (len(argout.axis_coords)>0):
            #{{{ transfer the errors and the axis labels
            #{{{ make dictionaries for both, and update with info from both, giving preference to self
            axesdict = selfout.mkd()
            #print "DEBUG 4: original mkd",axesdict
            errordict = selfout.mkd()
            # {{{ define a function that allows me to only update non-zero axes
            def non_empty_axes(input_data,ret_err = False):
                if ret_err:
                    temp_dict = input_data.mkd(input_data.axis_coords_error)
                else:
                    temp_dict = input_data.mkd(input_data.axis_coords)
                temp_dict = {k:v for k,v in temp_dict.items()
                        if v is not None and len(v) > 0}
                return temp_dict
            # }}}
            #{{{ add the axes and errors for B
            if isinstance(arg.axis_coords, list):
                if len(arg.axis_coords) > 0:
                    axesdict.update(non_empty_axes(arg))
            if isinstance(arg.axis_coords_error, list):
                if len(arg.axis_coords_error) > 0 and not all([x is None for x in arg.axis_coords_error]):
                    errordict.update(arg.mkd(arg.axis_coords_error))
            #}}}
            #{{{ add the axes and errors for A
            if isinstance(self.axis_coords, list):
                if len(self.axis_coords) > 0:
                    axesdict.update(non_empty_axes(self))
            if isinstance(self.axis_coords_error, list):
                if len(self.axis_coords_error) > 0 and not all([x is None for x in self.axis_coords_error]):
                    errordict.update(self.mkd(self.axis_coords_error))
            #}}}
            #}}}
            selfout.axis_coords_error = selfout.fld(errordict)
            argout.axis_coords_error = selfout.fld(errordict)
            selfout.axis_coords = selfout.fld(axesdict)
            argout.axis_coords = selfout.fld(axesdict)
            #}}}
            selfout.axis_coords_units = [None]*len(newdims)
            argout.axis_coords_units = [None]*len(newdims)
            for thisdim in newdims:
                if thisdim in self.dimlabels:
                    selfout.set_units(thisdim,self.get_units(thisdim))
                    argout.set_units(thisdim,self.get_units(thisdim))
                elif thisdim in arg.dimlabels:
                    selfout.set_units(thisdim,arg.get_units(thisdim))
                    argout.set_units(thisdim,arg.get_units(thisdim))
        return selfout,argout
    #}}}
    #{{{ integrate, differentiate, and sum
    def integrate(self, thisaxis, backwards=False, cumulative=False):
        r'''Performs an integration -- which is similar to a sum, except that it takes the axis into account, *i.e.*, it performs:
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
            '''
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        t = None
        if len(self.axis_coords)>0:
            t = self.getaxis(thisaxis)
            dt = t[1]-t[0]
        if t is None:
            raise ValueError("You can't call integrate on an unlabeled axis")
        if cumulative:
            self.run_nopop(cumsum,thisaxis)
            if backwards is True:
                self.data = self[thisaxis,::-1].data
        else:
            self.run(sum,thisaxis)
        self.data *= dt
        return self
    def diff(self,thisaxis,backwards = False):
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        self.run_nopop(mydiff,thisaxis)
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        if len(self.axis_coords)>0:
            t = self.getaxis(thisaxis)
            dt = t[1]-t[0]
            self.data /= dt
        return self
    def sum(self,axes):
        if (isinstance(axes, str)):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print('|-ERROR FINDING DIMENSION-----')
                print('| dimlabels is: ',self.dimlabels)
                print("| doesn't contain: ",axes[j])
                print('|-----------------------------')
                raise
            self.data = sum(self.data,
                    axis=thisindex)
            self._pop_axis_info(thisindex)
        return self
    def sum_nopop(self,axes):
        if (isinstance(axes, str)):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print('error, dimlabels is: ',self.dimlabels)
                print("doesn't contain: ",axes[j])
                raise
            temp = list(self.data.shape)
            temp[thisindex] = 1
            self.data = sum(self.data,
                    axis=thisindex)
            self.data = self.data.reshape(temp)
        return self
    #}}}
    #{{{ poly. fit
    def polyfit(self,axis,order=1,force_y_intercept = None):
        '''polynomial fitting routine -- return the coefficients and the fit
        ..note:
            later, should probably branch this off as a new type of fit class

        ..warning:
            for some reason, this version doesn't use orthogonal polynomials,
            as the numpy routine does -- we had diagnosed and determined that
            that creates noticeably different results, so fix that here.

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
        c: ndarray
            a standard numpy array containing the coefficients (in ascending polynomial order)
        formult: nddata
            an nddata containing the result of the fit
        '''
        x = self.getaxis(axis).copy().reshape(-1,1)
        #{{{ make a copy of self with the relevant dimension second to last (i.e. rows)
        formult = self.copy()
        neworder = list(formult.dimlabels)
        neworder.pop(neworder.index(axis))
        if len(neworder) > 1:
            neworder = neworder[:-1] + [axis] + neworder[-1]
        else:
            neworder = [axis] + neworder
        formult.reorder(neworder)
        #}}}
        y = formult.data
        #{{{ now solve Lx = y, where x is appropriate for our polynomial
        startingpower = 0
        if force_y_intercept is not None:
            startingpower = 1
        L =  concatenate([x**j for j in range(startingpower,order+1)],axis=1) # note the totally AWESOME way in which this is done!
        #print 'fitting to matrix',L
        if force_y_intercept is not None:
            y -= force_y_intercept
        c = dot(pinv(L),y)
        fity = dot(L,c)
        if force_y_intercept is not None:
            #print "\n\nDEBUG: forcing from",fity[0],"to"
            fity += force_y_intercept
            #print "DEBUG: ",fity[0]
            c = c_[force_y_intercept,c]
        #}}}
        #{{{ rather than have to match up everything, just drop the fit data into formult, which should be the same size, shape, etc
        formult.data = fity
        formult.set_error(None)
        #}}}
        return c,formult
    #}}}
    #{{{ max and mean
    def _wrapaxisfuncs(self,func):
        #{{{ for convenience, wrap the max and min functions
        if func == max:
            func = amax
        if func == min:
            func = amin
        if func == diff:
            func = mydiff
        return func
        #}}}
    def argmax(self,*args,**kwargs):
        r"""find the max along a particular axis, and get rid of that axis, replacing it with the index number of the max value
        
        Parameters
        ==========
        raw_index: bool
            return the raw (ndarray) numerical index, rather than the corresponding axis value
            Note that the result returned is still, however, an nddata (rather than numpy ndarray) object.
        """
        #{{{ process arguments
        axes = self._possibly_one_axis(*args)
        raw_index = False
        if 'raw_index' in list(kwargs.keys()):
            raw_index = kwargs.pop('raw_index')
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:",repr(kwargs))
        if (isinstance(axes, str)):
            axes = [axes]
        #}}}
        for j in range(0,len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except:
                print('error, dimlabels is: ',self.dimlabels)
                print("doesn't contain: ",axes[j])
                raise
            if raw_index:
                self.data = argmax(self.data,
                        axis=thisindex)
            else:
                if self.axis_coords[thisindex] is None:
                    raise ValueError("It doesn't make sense to call argmax if you have removed the axis coordinates! (getaxis yields None for %s"%thisindex)
                self.data = self.axis_coords[thisindex][argmax(self.data,
                    axis=thisindex)]
            self._pop_axis_info(thisindex)
        return self
    def argmin(self,*axes,**kwargs):
        r"""If `argmin('axisname')` find the min along a particular axis, and get rid of that
        axis, replacing it with the index number of the max value.
        If `argmin()`: return a dictionary giving the coordinates of the overall minimum point.
        
        Parameters
        ==========
        raw_index: bool
            Return the raw (ndarray) numerical index, rather than the corresponding axis value.
            Note that the result returned is still, however, an nddata (rather than numpy ndarray) object.
        """
        raw_index = process_kwargs([("raw_index",False)],kwargs)
        if len(axes) == 0:
            raw_indices = dict(zip(self.dimlabels,
                unravel_index(self.data.ravel().argmin(),self.data.shape)))
            if raw_index:
                return raw_indices
            else:
                return dict([(k,self.getaxis(k)[v]) for k,v in raw_indices.items()])
        if (type(axes) is str):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except:
                print('error, dimlabels is: ',self.dimlabels)
                print("doesn't contain: ",axes[j])
                raise
            if raw_index:
                self.data = argmin(self.data,
                        axis=thisindex)
            else:
                self.data = self.axis_coords[thisindex][argmin(self.data,
                    axis=thisindex)]
            self._pop_axis_info(thisindex)
        return self
    def max(self):
        return self.data.max()
    def min(self):
        return self.data.min()
    def cdf(self,normalized = True,max_bins = 500):
        """calculate the Cumulative Distribution Function for the data along `axis_name`

        only for 1D data right now
        
        Returns
        =======
        A new nddata object with an axis labeled `values`, and data corresponding to the CDF.
        """
        thisaxis = 0
        n_bins = self.data.shape[thisaxis]
        if n_bins > max_bins: n_bins = max_bins # otherwise this takes a while
        bins, vals = histogram(self.data,
                bins = n_bins)
        retval = nddata(double(bins),[-1],['values']).labels('values',
                vals[:-1]+(vals[1]-vals[0])*0.5)
        retval.run_nopop(cumsum,'values')
        if normalized:
            print('final value',retval['values',-1])
            retval /= retval['values',-1]
        return retval
    def mean_all_but(self,listofdims):
        'take the mean over all dimensions not in the list'
        for dimname in list(self.dimlabels):# I can't be popping from the list as I iterate over it
            if not dimname in listofdims:
                self.mean(dimname)
        return self
    def mean_weighted(self,axisname):
        r"""perform  the weighted mean along `axisname` (use $\sigma$ from $\sigma = $self.get_error() do generate $1/\sigma$ weights)
        for now, it clears the error of `self`, though it would be easy to calculate the new error, since everything is linear

        unlike other functions, this creates working objects that are themselves nddata objects
        this strategy is easier than coding out the raw numpy math, but probably less efficient"""
        #{{{ the weighted mean, pyspecdata style
        weight_matrix = self.copy().set_error(None)
        weight_matrix.data = 1. / self.get_error().copy()
        #{{{ find out where anything is nan, and set both error and weight to 0
        nan_mask = isnan(self.data)
        nan_mask |= isnan(weight_matrix.data)
        weight_matrix.data[nan_mask] = 0
        self.data[nan_mask] = 0
        #}}}
        #{{{ make sure there are no infinite values, because I wouldn't be sure how to deal with this
        inf_mask = isinf(self.data)
        inf_mask |= isinf(weight_matrix.data)
        assert not any(inf_mask)
        #}}}
        normalization = weight_matrix.copy().run(sum,axisname)
        weight_matrix /= normalization
        self.data *= weight_matrix.data
        self.set_error(None)
        self.run(sum,axisname)
        #}}}
        return self
    def mean(self,*args,**kwargs):
        r'''Take the mean and (optionally) set the error to the standard deviation

        Parameters
        ----------
        std: bool
            whether or not to return the standard deviation as an error
        '''
        logger.debug("entered the mean function")
        #{{{ process arguments
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        axes = self._possibly_one_axis(*args)
        if 'return_error' in kwargs:
            raise ValueError('return_error kwarg no longer used -- use std kwarg if you want to set the error to the std')
        return_error = process_kwargs([('std',False)],kwargs)
        logger.debug(strm("return error is",return_error))
        if (isinstance(axes, str)):
            axes = [axes]
        #}}}
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                logger.debug(strm('error, dimlabels is: ',self.dimlabels))
                logger.debug(strm("doesn't contain: ",axes[j]))
                raise
            if self.data_error is not None:
                this_axis_length = self.data.shape[thisindex]
                try:
                    self.data_error = sqrt(sum((self.data*self.data_error)**2,
                            axis=thisindex)/(this_axis_length**2))
                except:
                    raise ValueError(strm('shape of data',shape(self.data),'shape of data error',shape(self.data_error)))
            if return_error: # since I think this is causing an error
                thiserror = std(self.data,
                        axis=thisindex)
                if isscalar(thiserror):
                    thiserror = r_[thiserror]
            self.data = mean(self.data,
                    axis=thisindex)
            if return_error: # this needs to go after the data setting
                self.set_error(thiserror) # set the error to the standard deviation
            self._pop_axis_info(thisindex)
            logger.debug(strm("return error is",return_error))
        return self
    def mean_nopop(self,axis):
        self = self.run_nopop(mean,axis=axis)
        return self
    #}}}
    #{{{ running functions and popping dimensions
    def _pop_axis_info(self,thisindex):
        r'pop axis by index'
        self.dimlabels.pop(thisindex)
        if self.axis_coords!=[]:
            self.axis_coords.pop(thisindex)
            if self.axis_coords_error is not None and len(self.axis_coords_error) > 0:
                try:
                    self.axis_coords_error.pop(thisindex)
                except Exception as e:
                    raise RuntimeError(strm('trying to pop',thisindex,'from',self.axis_coords_error) + explain_error(e))
            if len(self.axis_coords_units) > 0:
                try:
                    self.axis_coords_units.pop(thisindex)
                except:
                    raise IndexError(strm('trying to pop',
                        thisindex, 'from', self.axis_coords_units))
        return self
    def popdim(self,dimname):
        thisindex = self.axn(dimname)
        thisshape = list(self.data.shape)
        if thisshape[thisindex]!=1:
            raise IndexError("trying to pop a dim that's not length 1")
        thisshape.pop(thisindex)
        self.data = self.data.reshape(thisshape)
        self._pop_axis_info(thisindex)
        return self
    def cropped_log(self,subplot_axes = None,magnitude = 4):
        r'''For the purposes of plotting, this generates a copy where I take the log, spanning "magnitude" orders of magnitude
        This is designed to be called as abs(instance).cropped_log(), so it doesn't make a copy'''
        phaseinfo = None
        if self.data.dtype == complex128:
            absdata = abs(self)
            phaseinfo = self/absdata
            self.data = absdata.data
        self.run(log10)
        if subplot_axes is None:# then do all
            self.data -= self.data.flatten().max() - magnitude # span only 4 orders of magnitude
        else:
            print("smooshing along subplot_axes",subplot_axes)
            newdata = self.copy().smoosh(subplot_axes,dimname = 'subplot')
            print(ndshape(newdata))
            newdata.run(max,'subplot')
            print(newdata)
            newdata = self - newdata
            self.data = newdata.data + magnitude
        self.data[self.data < 0] = 0
        if phaseinfo is not None:
            self.data = self.data * phaseinfo.data
        return self
    def runcopy(self,*args):
        newdata = self.copy()
        func = args[0]
        func = self._wrapaxisfuncs(func)
        if len(args)>1:
            axis = args[1]
            thisindex = newdata.dimlabels.index(axis)
            newdata.data = func(newdata.data,axis=thisindex)
            newdata._pop_axis_info(thisindex)
        else:
            newdata.data = func(newdata.data)
        return newdata
    def run(self,*args):
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
        if len(args)>1:
            axis = args[1]
            try:
                thisindex = self.dimlabels.index(axis)
            except Exception as e:
                if not isinstance(axis, str):
                    raise ValueError('The format of run is run(func,"axisname"), but you didn\'t give a string as the second argument -- maybe you fed the arguments backwards?')
                elif axis not in self.dimlabels:
                    raise ValueError("axis "+axis+
                            " is not in dimlabels ("+
                            repr(self.dimlabels)+")")
                else:
                    raise e
            self.data = func(self.data,axis=thisindex)
            self._pop_axis_info(thisindex)
            return self
        else:
            retval = func(self.data)
            if self.data.size == retval.size:
                self.data = retval
                return self
            else:
                return retval
    def run_nopop(self,func,axis):
        func = self._wrapaxisfuncs(func)
        try:
            thisaxis = self.axn(axis)
        except Exception as e:
            raise IndexError(strm("I couldn't find the dimension",axis,
                "in the list of axes",self.dimlabels))
        temp = list(self.data.shape)
        temp[thisaxis] = 1
        func_sig = signature(func)
        numnonoptargs = len([v.default for v in func_sig.parameters.values() if v.default==Parameter.empty])
        kwargnames =  [k for k,v in func_sig.parameters.items() if v.default!=Parameter.empty]
        if numnonoptargs == 2:
            if 'axis' in kwargnames:
                self.data = func(self.getaxis(axis),self.data,axis=thisaxis)
            if 'axes' in kwargnames:
                self.data = func(self.getaxis(axis),self.data,axes=thisaxis)
            else:
                raise ValueError("Your function doesn't have axis or axes as a keyword argument!")
        else:
            if numnonoptargs == 1:
                paramnames = [k for k in func_sig.parameters.keys()]
                if len(paramnames) == 1:
                    if func_sig.parameters[paramnames[0]].kind == Parameter.VAR_POSITIONAL:
                        try:
                            self.data = func(self.data,axis=thisaxis)
                        except:
                            self.data = func(self.data,axes=thisaxis)
                if 'axis' in kwargnames:
                    self.data = func(self.data,axis=thisaxis)
                elif 'axes' in kwargnames:
                    self.data = func(self.data,axes=thisaxis)
                else:
                    raise ValueError("Your function doesn't have axis or axes as a keyword argument! The number of non-optional arguments are %s. The keyword arguments are %s"%(str(numnonoptargs),str(kwargnames)))
            else:
                raise ValueError('you passed a function to run_nopop that doesn\'t'
                        'have either one or two arguments!')
        #{{{ if the function doesn't rip out the dim, make sure we don't change the dims
        if len(self.data.shape)==len(temp):
            temp[thisaxis] = self.data.shape[thisaxis]
        #}}}
        self.data = self.data.reshape(temp)
        return self
    def item(self):
        r"like numpy item -- returns a number when zero-dimensional"
        return self.data.item()
    #}}}
    #{{{ ft-related functions
    def unitify_axis(self,axis_name,
            is_axis=True):
        'this just generates an axis label with appropriate units'
        if type(axis_name) in [int,int64]:
            axis_name = self.dimlabels[axis_name]
        if self.get_prop('FT') is not None and axis_name in list(self.get_prop('FT').keys()) and self.get_prop('FT')[axis_name]:
            isft = True
        else:
            isft = False
        if is_axis:
            yunits = self.units_texsafe(axis_name)
            j = axis_name.find('_')
            if j > -1:
                prevword = axis_name[0:j]
                if j+1< len(axis_name):
                    followword = axis_name[j+1:]
                else:
                    followword = []
                k = followword.find(' ')
                if k > -1 and k < len(followword):
                    followword = followword[:k]
                k = followword.find('_')
                if len(followword) > 0:
                    if not (k > -1) and (len(prevword) < 2 or len(followword) < 2):
                        if len(followword) > 1:
                            axis_name = axis_name[:j+1+len(followword)]  + '}$' + axis_name[j+1+len(followword):]
                            axis_name = axis_name[:j+1] + '{' + axis_name[j+1:]
                        else:
                            axis_name = axis_name[0:j+2] + '$' + axis_name[j+2:]
                        axis_name = '$'+axis_name
            if isft:
                t_idx = axis_name.find('t')
                if t_idx>-1:
                    if t_idx+1 < len(axis_name) and axis_name[t_idx+1].isalpha():
                        axis_name = r'F{'+axis_name+r'}'
                    else:
                        axis_name = axis_name.replace('t','\\nu ')
                        if axis_name[0] != '$':
                            axis_name = '$' + axis_name + '$'
                #elif axis_name[:2] == 'ph':
                #    if len(axis_name) > 2:
                #        axis_name = r'$\Delta c_{'+axis_name[2:]+'}$'
                #    else:
                #        axis_name = r'$\Delta c$'
                else:
                    axis_name = r'F{'+axis_name+r'}'
        else:
            yunits = self.units_texsafe()
        if yunits is not None:
            axis_name = axis_name + ' / ' + yunits
        return axis_name
    #{{{ the following are all in the desired format -- the repetition at the end is because each function is in its own file (module) of the same name
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
    #}}}
    #}}}
    def nnls(self, dimname, newaxis_dict, kernel_func, l=0):
        r"""Perform regularized non-negative least-squares "fit" on self.

        Capable of solving for solution in 1 or 2 dimensions.

        We seek to minimize
        :math:`Q = \| Ax - b \|_2 + \|\lambda x\|_2`
        in order to obtain solution vector :math:`x` subject to non-negativity constraint
        given input matrix :math:`A`, the kernel, and input vector :math:`b`, the data.

        The first term assesses agreement between the fit :math:`Ax` and the data :math:`b`,
        and the second term accounts for noise with the regularization parameter :math:`\lambda`
        according to Tikhonov regularization.

        To perform regularized minimization in 2 dimensions, set `l` to :str:`BRD` and provide a
        tuple of parameters :str:`dimname`, :nddata:`newaxis_dict`, and :function:`kernel_func`.
        Algorithm described in Venkataramanan et al. 2002 is performed which determines optimal :math:`\lambda`
        for the data (DOI:10.1109/78.995059).
        
        See: `Wikipedia page on NNLS <https://en.wikipedia.org/wiki/Non-negative_least_squares>`_,
        `Wikipedia page on Tikhonov regularization <https://en.wikipedia.org/wiki/Tikhonov_regularization>`_
         
        Parameters
        ==========
        dimname: str
            Name of the "data" dimension that is to be replaced by a
            distribution (the "fit" dimension);
            *e.g.* if you are regularizing a set of functions
            :math:`\exp(-\tau*R_1)`, then this is :math:`\tau`
        newaxis_dict: dict or nddata
            a dictionary whose key is the name of the "fit" dimension
            (:math:`R_1` in the example above)
            and whose value is an array with the new axis labels.
            OR
            this can be a 1D nddata
            -- if it has an axis, the axis will be used to create the
            fit axis; if it has no axis, the data will be used
        kernel_func: function
            a function giving the kernel for the regularization.
            The first argument is the "data" variable
            and the second argument is the "fit" variable
            (in the example above, this would be something like
            ``lambda x,y: exp(-x*y)``)
        l : double (default 0) or str
            the regularization parameter
            :math:`lambda` -- if this is set to 0, the algorithm reverts to
            standard nnls.
            If this is set to :str:`BRD`, then algorithm expects tuple of each parameter
            described above in order to perform a 2-dimensional fit.

        Returns
        =======
        self:
            The regularized result.
            For future use, both the kernel (as an nddata, in a property called
            "nnls_kernel") and the residual (as an nddata, in a property called
            "nnls_residual") are stored as properties of the nddata.
            The regularized dimension is always last
            (innermost).
            If :str:`BRD` is specified, then the individual, uncompressed kernels :math:`K_{1}` and :math:`K_{2}` are returned as properties of the nddata "K1" and "K2" respectively. The number of singular values used to compressed each kernel is returned in properties of the nddata called, respectively, "s1" and "s2". 
        """
        logger.debug(strm('on first calling nnls, shape of the data is',ndshape(self),'is it fortran ordered?',isfortran(self.data)))
        tuple_syntax = False
        if isinstance(dimname, tuple):
            tuple_syntax = True
            assert len(dimname) == 2, "tuple of two dimension names only"
            assert type(dimname[0]) and isinstance(dimname[1], str), "first argument is tuple of two dimension names"
        else:
            assert isinstance(dimname, str), "first argument is dimension name or tuple of two dimension names"
        if isinstance(newaxis_dict, tuple):
            assert len(newaxis_dict) == 2, "tuple of two nddatas only"
            if isinstance(newaxis_dict[0],nddata) and isinstance(newaxis_dict[1],nddata):
                assert len(newaxis_dict[0].dimlabels) and len(newaxis_dict[1].dimlabels) == 1, "currently only set up for 1D"
        elif isinstance(newaxis_dict, dict):
            assert len(newaxis_dict) == 1, "currently only set up for 1D"
        elif isinstance(newaxis_dict,nddata):
            assert len(newaxis_dict.dimlabels) == 1, "currently only set up for 1D"
        else:
            raise ValueError("second argument is dictionary or nddata with new axis, or tuple of nddatas with new axes")
        if isinstance(kernel_func, tuple):
            assert callable(kernel_func[0]) and callable(kernel_func[1]), "third argument is tuple of kernel functions"
        else:
            assert callable(kernel_func), "third argument is kernel function"
        # construct the kernel
        # the kernel transforms from (columns) the "fit" dimension to (rows)
        # the "data" dimension
        if tuple_syntax:
            if isinstance(newaxis_dict[0],nddata):
                assert len(newaxis_dict[0].dimlabels) and len(newaxis_dict[1].dimlabels) == 1, "must be 1 dimensional!!"
                fitdim_name1 = newaxis_dict[0].dimlabels[0]
                fitdim_name2 = newaxis_dict[1].dimlabels[0]
                fit_axis1 = newaxis_dict[0].getaxis(fitdim_name1)
                fit_axis2 = newaxis_dict[1].getaxis(fitdim_name2)
                if fit_axis1 is None:
                    fit_axis1 = newaxis_dict[0].data
                if fit_axis2 is None:
                    fit_axis2 = newaxis_dict[1].data
        elif isinstance(newaxis_dict,nddata):
            assert len(newaxis_dict.dimlabels) == 1, "must be 1 dimensional!!"
            fitdim_name = newaxis_dict.dimlabels[0]
            fit_axis = newaxis_dict.getaxis(fitdim_name)
            if fit_axis is None:
                fit_axis = newaxis_dict.data
        else:
            fitdim_name = list(newaxis_dict.keys())[0]
            logger.debug(strm('shape of fit dimension is',newaxis_dict[fitdim_name].shape))
            fit_axis = newaxis_dict[fitdim_name]
        if tuple_syntax:
            fit_axis1 = nddata(fit_axis1,fitdim_name1)
            fit_axis2 = nddata(fit_axis2,fitdim_name2)
            data_axis1 = self.fromaxis(dimname[0])
            data_axis2 = self.fromaxis(dimname[1])
            data_axis1.squeeze()
            data_axis2.squeeze()
            data_axis1,fit_axis1 = data_axis1.aligndata(fit_axis1)
            data_axis2,fit_axis2 = data_axis2.aligndata(fit_axis2)
            K1 = kernel_func[0](data_axis1,fit_axis1).squeeze()
            K1_ret = K1
            logger.debug(strm('K1 dimlabels',K1.dimlabels,'and raw shape',K1.data.shape))
            K2 = kernel_func[1](data_axis2,fit_axis2).squeeze()
            K2_ret = K2
            logger.debug(strm('K2 dimlabels',K2.dimlabels,'and raw shape',K2.data.shape))
            # SVD and truncation of kernels
            U1,S1,V1 = np.linalg.svd(K1.data,full_matrices=False)
            U2,S2,V2 = np.linalg.svd(K2.data,full_matrices=False)
            logger.debug(strm('Uncompressed SVD K1:',[x.shape for x in (U1,S1,V1)]))
            logger.debug(strm('Uncompressed SVD K2:',[x.shape for x in (U2,S2,V2)]))
            default_cut = 1e-2
            s1 = where(S1 > default_cut)[0][-1]
            s2 = where(S2 > default_cut)[0][-1]
            U1 = U1[:,0:s1]
            S1 = S1[0:s1]
            V1 = V1[0:s1,:]
            U2 = U2[:,0:s2]
            S2 = S2[0:s2]
            V2 = V2[0:s2,:]
            S1 = S1*eye(s1)
            S2 = S2*eye(s2)
            logger.debug(strm('Compressed SVD of K1:',[x.shape for x in (U1,S1,V1)]))
            logger.debug(strm('Compressed SVD K2:',[x.shape for x in (U2,S2,V2)]))
            # would generate projected data here
            # compress data here
            K1 = S1.dot(V1)
            K2 = S2.dot(V2)
            K = K1[:,newaxis,:,newaxis]*K2[newaxis,:,newaxis,:]
            K = K.reshape(K1.shape[0]*K2.shape[0],K1.shape[1]*K2.shape[1])
            logger.debug(strm('Compressed K0, K1, and K2:',[x.shape for x in (K,K1,K2)]))

            data_compressed = U1.T.dot(self.data.dot(U2))
            logger.debug(strm('Compressed data:',data_compressed.shape))
            # data_for_nnls = nddata(data_compressed,[dimname[0],dimname[1]])
            # data_for_nnls.smoosh([dimname[0],dimname[1]],dimname=dimname[0])
            data_fornnls = empty(s1*s2)
            for s1_index in range(s1):
                for s2_index in range(s2):
                    temp = data_compressed[s1_index][s2_index]
                    data_fornnls[s1_index*s2+s2_index] = temp
            logger.debug(strm('Lexicographically ordered data:',data_fornnls.shape))
            if len(data_fornnls.shape) > 2:
                logger.debug(strm('Reshpaing data..'))
                data_fornnls = data_fornnls.reshape((prod(data_fornnls.shape[:-1]),data_fornnls.shape[-1]))
            
            if l == 'BRD':
                def chi(x_vec,val):
                    return 0.5*dot(x_vec.T,dot(dd_chi(G(x_vec),val**2),x_vec)) - dot(x_vec.T,data_fornnls[:,newaxis])
                def d_chi(x_vec,val):
                    return dot(dd_chi(G(x_vec),val**2),x_vec) - data_fornnls[:,newaxis]
                def dd_chi(G,val):
                    return G + (val**2)*eye(shape(G)[0])
                def G(x_vec):
                    return dot(K,dot(square_heaviside(x_vec),K.T))
                def H(product):
                    if product <= 0:
                        return 0
                    if product > 0:
                        return 1
                def square_heaviside(x_vec):
                    diag_heavi = []
                    for q in range(shape(K.T)[0]):
                        pull_val = dot(K.T[q,:],x_vec)
                        temp = pull_val[0]
                        temp = H(temp)
                        diag_heavi.append(temp)
                    diag_heavi = array(diag_heavi)
                    square_heavi = diag_heavi*eye(shape(diag_heavi)[0])
                    return square_heavi
                def optimize_alpha(input_vec,val):
                    alpha_converged = False
                    factor = sqrt(s1*s2)
                    T = linalg.inv(dd_chi(G(input_vec),val**2))
                    dot_product = dot(input_vec.T,dot(T,input_vec))
                    ans = dot_product*factor
                    ans = ans/linalg.norm(input_vec)/dot_product
                    tol = 1e-3
                    if abs(ans-val**2) <= tol:
                        logger.debug(strm('ALPHA HAS CONVERGED.'))
                        alpha_converged = True
                        return ans,alpha_converged
                    return ans,alpha_converged
                def newton_min(input_vec,val):
                    fder = dd_chi(G(input_vec),val)
                    fval = d_chi(input_vec,val)
                    return (input_vec + dot(linalg.inv(fder),fval))
                def mod_BRD(guess,maxiter=20):
                    smoothing_param = guess
                    alpha_converged = False
                    for iter in range(maxiter):
                        logger.debug(strm('ITERATION NO.',iter))
                        logger.debug(strm('CURRENT LAMBDA',smoothing_param))
                        retval,residual = this_nnls.nnls_regularized(K,data_fornnls,l=smoothing_param)
                        f_vec = retval[:,newaxis]
                        alpha = smoothing_param**2
                        c_vec = dot(K,f_vec) - data_fornnls[:,newaxis]
                        c_vec /= -1*alpha
                        c_update = newton_min(c_vec,smoothing_param)
                        alpha_update,alpha_converged = optimize_alpha(c_update,smoothing_param)
                        lambda_update = sqrt(alpha_update[0,0])
                        if alpha_converged:
                            logger.debug(strm('*** OPTIMIZED LAMBDA',lambda_update,'***'))
                            break
                        if not alpha_converged:
                            logger.debug(strm('UPDATED LAMBDA',lambda_update))
                            smoothing_param = lambda_update
                        if iter == maxiter-1:
                            logger.debug(strm('DID NOT CONVERGE.'))
                    return lambda_update
                retval, residual = this_nnls.nnls_regularized(K,data_fornnls,l=mod_BRD(guess=1.0))
            else:
                retval, residual = this_nnls.nnls_regularized(K,data_fornnls,l=l)
            logger.debug(strm('coming back from fortran, residual type is',type(residual))+ strm(residual.dtype if isinstance(residual, ndarray) else ''))
            newshape = []
            if not isscalar(l):
                newshape.append(len(l))
            logger.debug(strm('test***',list(self.data.shape)[:-1]))
            #newshape += list(self.data.shape)[:-1] # this would return parametric axis
            newshape.append(ndshape(fit_axis1)[fitdim_name1])
            newshape.append(ndshape(fit_axis2)[fitdim_name2])
            logger.debug(strm('before mkd, shape of the data is',ndshape(self),'len of axis_coords_error',len(self.axis_coords_error)))
            # {{{ store the dictionaries for later use
            axis_coords_dict = self.mkd(self.axis_coords)
            axis_units_dict = self.mkd(self.axis_coords_units)
            axis_coords_error_dict = self.mkd(self.axis_coords_error)
            # }}}
            retval = retval.reshape(newshape)
            self.data = retval
            # {{{ clear axis info
            self.axis_coords = None
            self.axis_coords_units = None
            self.axis_coords_error_dict = None
            # }}}
            # change the dimnesion names and data
            self.rename(dimname[0], fitdim_name1)
            self.rename(dimname[1], fitdim_name2)
            axis_coords_dict[fitdim_name1] = fit_axis1.getaxis(fitdim_name1)
            axis_units_dict[fitdim_name1] = None
            axis_coords_error_dict[fitdim_name1] = None
            axis_coords_dict[fitdim_name2] = fit_axis2.getaxis(fitdim_name2)
            axis_units_dict[fitdim_name2] = None
            axis_coords_error_dict[fitdim_name2] = None
            if not isscalar(l):
                self.dimlabels = ['lambda'] + self.dimlabels
                axis_coords_dict['lambda'] = l
                axis_units_dict['lambda'] = None
                axis_coords_error_dict['lambda'] = None
            self.data = retval
            if not isscalar(residual):
                # make the residual nddata as well
                residual_nddata = ndshape(self).pop(fitdim_name2).pop(fitdim_name1).alloc(dtype=residual.dtype)
                residual_nddata.data[:] = residual[:]
            else:
                residual_nddata = residual
            # store the kernel and the residual as properties
            self.set_prop('nnls_kernel',K)
            self.set_prop('s1',s1)
            self.set_prop('s2',s2)
            self.set_prop('nnls_residual',residual_nddata)
            self.set_prop('K1',K1_ret)
            self.set_prop('K2',K2_ret)
            # {{{ use the info from the dictionaries
            self.axis_coords = self.fld(axis_coords_dict)
            self.axis_coords_units = self.fld(axis_units_dict)
            self.axis_coords_error = self.fld(axis_coords_error_dict)
            return self
        else:
            fit_axis = nddata(fit_axis, fitdim_name)
            data_axis = self.fromaxis(dimname)
            data_axis.squeeze()
            data_axis, fit_axis = data_axis.aligndata(fit_axis)
            K = kernel_func(data_axis, fit_axis).squeeze()
            logger.debug(strm('K dimlabels',K.dimlabels,'and raw shape',K.data.shape))
            self.reorder(dimname, first=False) # make the dimension we will be regularizing innermost
            logger.debug(strm('shape of the data is',ndshape(self),'is it fortran ordered?',isfortran(self.data)))
            data_fornnls = self.data
            if len(data_fornnls.shape) > 2:
                data_fornnls = data_fornnls.reshape((prod(
                    data_fornnls.shape[:-1]),data_fornnls.shape[-1]))
            logger.debug(strm('shape of the data is',ndshape(self),"len of axis_coords_error",len(self.axis_coords_error)))
            retval, residual = this_nnls.nnls_regularized(K.data, data_fornnls, l=l)
            logger.debug(strm("coming back from fortran, residual type is",type(residual))+ strm(residual.dtype if isinstance(residual, ndarray) else ''))
            newshape = []
            if not isscalar(l):
                newshape.append(len(l))
            newshape += list(self.data.shape)[:-1] # exclude data dimension
            newshape.append(ndshape(fit_axis)[fitdim_name])
            logger.debug(strm('before mkd, shape of the data is',ndshape(self),"len of axis_coords_error",len(self.axis_coords_error)))
            # {{{ store the dictionaries for later use
            axis_coords_dict = self.mkd(self.axis_coords)
            axis_units_dict = self.mkd(self.axis_coords_units)
            axis_coords_error_dict = self.mkd(self.axis_coords_error)
            # }}}
            retval = retval.reshape(newshape)
            self.data = retval
            # {{{ clear all the axis info
            self.axis_coords = None
            self.axis_coords_units = None
            self.axis_coords_error_dict = None
            # }}}
            # change the dimension names and data
            self.rename(dimname, fitdim_name)
            # {{{ manipulate the dictionaries, and call fld below
            axis_coords_dict[fitdim_name] = fit_axis.getaxis(fitdim_name)
            axis_units_dict[fitdim_name] = None
            axis_coords_error_dict[fitdim_name] = None
            if not isscalar(l):
                self.dimlabels = ['lambda'] + self.dimlabels
                axis_coords_dict['lambda'] = l
                axis_units_dict['lambda'] = None
                axis_coords_error_dict['lambda'] = None
            # }}}
            self.data = retval
            if not isscalar(residual):
                # make the residual nddata as well
                residual_nddata = ndshape(self).pop(fitdim_name).alloc(dtype=residual.dtype)
                residual_nddata.data[:] = residual[:]
            else:
                residual_nddata = residual
            # store the kernel and the residual in the properties
            self.set_prop('nnls_kernel',K)
            self.set_prop('nnls_residual',residual_nddata)
            # {{{ use the info from the dictionaries
            self.axis_coords = self.fld(axis_coords_dict)
            self.axis_coords_units = self.fld(axis_units_dict)
            self.axis_coords_error = self.fld(axis_coords_error_dict)
            # }}}
            return self
    #{{{ interpolation and binning
    def run_avg(self,thisaxisname,decimation = 20,centered = False):
        'a simple running average'
        temp = self.getaxis(thisaxisname).size % decimation
        decimation = int(decimation)
        if temp != 0:
            if centered:
                self = self[thisaxisname,temp//2:-int(temp/2. + 0.5)]
            else:
                self = self[thisaxisname,0:-temp]
        thisaxis = nddata(self.getaxis(thisaxisname),[-1],[thisaxisname])
        self.setaxis(thisaxisname,[])
        self.chunkoff(thisaxisname,['avg'],[decimation])
        self.run(mean,'avg')
        thisaxis.chunkoff(thisaxisname,['avg'],[decimation])
        thisaxis.run(mean,'avg')
        self.setaxis(thisaxisname,thisaxis.data)
        return self
    def histogram(*args,**kwargs):
        per_hit = 1e-3
        if 'per_hit' in list(kwargs.keys()):
            per_hit = kwargs.pop('per_hit')
        if len(kwargs) > 0:
            raise ValueError("I don't understand the arguments:"+repr(kwargs))
        if len(args) == 1:
            if isinstance(args[0], ndarray):
                if args[0].shape[2] != len(self.dimlabels):
                    raise ValueError("You must pass an N x M array, where M is the number of dimensions in this array!")
            elif isinstance(args[0], dict):
                raise ValueError("should eventually support dictionaries (via fld), but doesn't yet")
            else:
                raise ValueError("I don't know what funny business you're up to passing me a"+repr(type(args[0])))
        else:
            raise ValueError("should eventually support array, label pair, but doesn't yet")
    def interp(self,axis,axisvalues, past_bounds=None, return_func=False, **kwargs):
        '''interpolate data values given axis values
        
        Parameters
        ==========
        return_func : boolean
            defaults to False.  If True, it returns a function that accepts
            axis values and returns a data value.
        '''
        oldaxis = self.getaxis(axis)
        if not return_func:
            if (isinstance(axisvalues, int)) or (isinstance(axisvalues, int32)):
                axisvalues = linspace(oldaxis[0],oldaxis[-1],axisvalues)
            elif isscalar(axisvalues):
                axisvalues = r_[axisvalues]
            elif (type(axisvalues) not in [ndarray,tuple]):
                raise ValueError("You passed a target axis of type"+repr(type(axisvalues))+"which I don't get")
            if any(imag(axisvalues) > 1e-38):
                raise ValueError("I can't interpolate imaginary values")
            else:
                axisvalues = real(axisvalues)
            if past_bounds is None:
                axisvalues[axisvalues<oldaxis.min()] = oldaxis.min()
                axisvalues[axisvalues>oldaxis.max()] = oldaxis.max()
            elif not (past_bounds == 'fail'):
                if isinstance(past_bounds, tuple):
                    if len(past_bounds) == 2:
                        axisvalues[axisvalues<oldaxis.min()] = past_bounds[0]
                        axisvalues[axisvalues>oldaxis.max()] = past_bounds[1]
                    else:
                        raise TypeError('If you pass axisvalues as a tuple, it must be of length 2!')
                else:
                    axisvalues[axisvalues<oldaxis.min()] = past_bounds
                    axisvalues[axisvalues>oldaxis.max()] = past_bounds
        rdata = real(self.data)
        idata = imag(self.data)
        thiserror = self.get_error()
        if thiserror is not None:
            rerrvar = real(thiserror)**2
            if thiserror[0].dtype == 'complex128':
                ierrvar = imag(thiserror)**2
        if 'kind' in list(kwargs.keys()):
            thiskind = kwargs.pop('kind')
        else:
            thiskind = 'cubic'
            if len(rdata) < 4:
                thiskind = 'quadratic'
                if len(rdata) < 3:
                    thiskind = 'linear'
        thisaxis = self.axn(axis)
        logger.debug(strm('Using %s interpolation'%thiskind))
        def local_interp_func(local_arg_data,kind = thiskind):
            interpfunc =  interp1d(oldaxis,local_arg_data,kind = kind,axis = thisaxis)
            try:
                retval = interpfunc(axisvalues)
            except:
                raise TypeError("dtype of axis is"+repr(axisvalues.dtype))
            return retval
        if return_func:
            rfunc = interp1d(oldaxis, rdata, kind=thiskind, axis=thisaxis, bounds_error=False, fill_value=tuple(rdata[r_[0,-1]].tolist()))
            ifunc = interp1d(oldaxis, idata, kind=thiskind, axis=thisaxis, bounds_error=False, fill_value=tuple(idata[r_[0,-1]].tolist()))
            return lambda x: rfunc(x) + 1j*ifunc(x)
        rdata = local_interp_func(rdata)
        idata = local_interp_func(idata)
        self.data = rdata + 1j * idata
        self.setaxis(axis,axisvalues)
        if thiserror is not None:
            rerrvar = local_interp_func(rerrvar,kind = 'linear') # calculate the error variance of the real part, use linear to avoid nan problems
            if thiserror[0].dtype == 'complex128':
                ierrvar = local_interp_func(ierrvar,kind = 'linear')
                self.set_error(sqrt(rerrvar) + 1j * sqrt(ierrvar))
            else:
                self.set_error(sqrt(rerrvar))
            err_nanmask = isnan(self.get_error())
            self.data[err_nanmask] = nan
        return self
    def invinterp(self,axis,values,**kwargs):
        'interpolate axis values given data values'
        copy = process_kwargs([('copy',False),
            ], kwargs, pass_through=True)
        
        if isscalar(values):
            values = r_[values]
        origdata = self.data.copy()
        origaxis = self.getaxis(axis).copy()
        if any(imag(values) > 1e-38):
            raise ValueError("I can't interpolate imaginary values")
        else:
            values = real(values)
        args = origdata.argsort()
        origdata = origdata[args]
        rdata = real(origaxis)
        idata = imag(origaxis)
        rdata = rdata[args]
        idata = idata[args]
        #{{{ determine type of interpolation
        if 'kind' in list(kwargs.keys()):
            thiskind = kwargs.pop('kind')
        else:
            thiskind = 'cubic'
            if len(rdata) < 4:
                thiskind = 'quadratic'
                if len(rdata) < 3:
                    thiskind = 'linear'
        #}}}
        interpfunc =  interp1d(origdata,rdata, kind=thiskind, **kwargs)
        try:
            rdata = interpfunc(values)
        except Exception as e:
            raise ValueError(strm('You passed',values,'and the data spans from',
                    origdata.min(),'to',origdata.max())+explain_error(e))
        interpfunc =  interp1d(origdata,idata, kind=thiskind, **kwargs)
        idata = interpfunc(values)
        cdata = rdata + 1j * idata
        if copy:
            newaxis = 'name data before interpolation'
            if self.name() is not None:
                newaxis = self.name()
            retval = nddata(cdata,[-1],[newaxis])
            retval.labels(newaxis,values)
            return retval
        else:
            self.data = values
            self.setaxis(axis,cdata)
            return self
    def contiguous(self,lambdafunc,axis = None):
        r"""Return contiguous blocks that satisfy the condition given by `lambdafunc`
        
        this function returns the start and stop positions along the
        axis for the contiguous blocks for which lambdafunc returns
        true
        **Currently only supported for 1D data**
        
        .. note::
            adapted from stackexchange post http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array

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
        retval : ndarray
            An :math:`N\times 2` matrix, where the :math:`N` rows correspond to pairs of axis
            label that give ranges over which `lambdafunc` evaluates to `True`.
            These are ordered according to descending range width.

        Examples
        --------

        .. code:: python

            sum_for_contiguous = abs(forplot).mean('t1')
            fl.next("test contiguous")
            forplot = sum_for_contiguous.copy().set_error(None)
            fl.plot(forplot,alpha = 0.25,linewidth = 3)
            print("this is what the max looks like",0.5*sum_for_contiguous.set_error(None).runcopy(max,'t2'))
            print(sum_for_contiguous > 0.5*sum_for_contiguous.runcopy(max,'t2'))
            retval = sum_for_contiguous.contiguous(quarter_of_max,'t2')
            print("contiguous range / 1e6:",retval/1e6)
            for j in range(retval.shape[0]):
                a,b = retval[j,:]
                fl.plot(forplot['t2':(a,b)])

        """
        logger.debug(strm("(contiguous) shape of self inside contiguous",ndshape(self)))
        if axis is None:
            if len(self.dimlabels) == 1:
                axis = self.dimlabels[0]
                mask = lambdafunc(self.copy()).data
            else:
                raise TypeError("If there is more than one dimension, `axis` must be set to something other than ``None``")
        else:
            mask = lambdafunc(self.copy(),axis).data
        logger.debug(strm("(contiguous) shape of mask",mask.shape))
        if axis is None:
            idx, = np.diff(mask).nonzero() # this gives a list of indices for the boundaries between true/false
        else:
            idx, = np.diff(mask, axis=self.axn(axis)).nonzero() # this gives a list of indices for the boundaries between true/false
        idx += 1 # because diff only starts on index #1 rather than #0
        if mask[0]: # because I want to indicate the boundaries of True, if I am starting on True, I need to make 0 a boundary
            idx = np.r_[0, idx]
        if mask[-1]: # If the end of mask is True, then I need to add a boundary there as well
            idx = np.r_[idx, mask.size-1] # Edit
        idx.shape = (-1,2) # idx is 2x2 array of start,stop
        logger.debug(strm('(contiguous) DEBUG idx is',idx))
        logger.debug(strm("(contiguous) diffs for blocks",diff(idx,axis=1)))
        block_order = diff(idx, axis=1).flatten().argsort()[::-1]
        logger.debug(strm("(contiguous) in descending order, the blocks are therefore",idx[block_order,:]))
        return self.getaxis(axis)[idx[block_order,:]]
    def to_ppm(self):
        """Function that converts from Hz to ppm using Bruker parameters

        .. todo::

            Figure out what the units of PHC1 in Topspin are (degrees per *what??*), and apply those as well.

            make this part of an inherited bruker class
        """
        if self.get_units('t2') == 'ppm': return
        offset = self.get_prop('proc')['OFFSET']
        sfo1 = self.get_prop('acq')['SFO1']
        if not self.get_ft_prop('t2'):
            self.ft('t2', shift=True) # this fourier transforms along t2, overwriting the data that was in self
        self.setaxis('t2', lambda x:
                x/sfo1).set_units('t2','ppm')
        max_ppm = self.getaxis('t2').max()
        self.setaxis('t2', lambda x:
                (x-max_ppm+offset))
        self.set_prop('x_inverted',True)
        return self
    #}}}
    def multimin(self,minfunc,axisname,filterwidth,numberofmins):
        cost = self.copy().convolve(axisname,filterwidth).run_nopop(minfunc)
        for j in range(0,numberofmins):
            #{{{ find the x value at which the minimum occurs
            xvalues = cost.copy().argmin(axisname,raw_index = True)
            #}}}
    def repwlabels(self,axis):
        return None
    def add_noise(self,intensity):
        '''Add Gaussian (box-muller) noise to the data.

        Parameters
        ----------
        intensity : double OR function
            If a double, gives the standard deviation of the noise.
            If a function, used to calculate the standard deviation of the noise from the data:
            *e.g.* ``lambda x: max(abs(x))/10.``
            
        '''
        if isinstance(intensity, type(emptyfunction)):
            intensity = intensity(lambda x: self.data)
        return_complex = iscomplexobj(self.data)
        self.data += box_muller(self.data.size, return_complex=return_complex).reshape(self.data.shape) * intensity
        return self
    #{{{ functions to manipulate and return the axes
    def reorder(self,*axes,**kwargs):
        r'''Reorder the dimensions
        the first arguments are a list of dimensions

        Parameters
        ----------
        *axes : str
            Accept any number of arguments that gives the dimensions, in the
            order that you want thee.
        first : bool
            (default True)
            Put this list of dimensions first, while False puts them last (where they then come in the order given).
        '''
        first = True
        if 'first' in kwargs:
            first = kwargs.pop('first')
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
            neworder = list(map(self.dimlabels.index,axes))
        except ValueError as e:
            raise ValueError(strm('one of',axes,'not in',self.dimlabels))
        self.dimlabels = list(map(self.dimlabels.__getitem__,neworder))
        if len(self.axis_coords)>0:
            try:
                self.axis_coords = list(map(self.axis_coords.__getitem__,neworder))
            except Exception as e:
                raise IndexError(strm('problem mapping',list(map(len,self.axis_coords)),'onto',neworder)
                        + explain_error(e))
            if len(self.axis_coords_units)>0:
                try:
                    self.axis_coords_units = list(map(self.axis_coords_units.__getitem__,neworder))
                except:
                    raise IndexError(strm('problem mapping',list(map(len,self.axis_coords_units)),
                        'onto',neworder)+explain_error(e))
        try:
            self.data = self.data.transpose(neworder)
        except ValueError as e:
            raise ValueError(strm('you can\'t reorder',self.dimlabels,'as',neworder)
                    +explain_error(e))
        if self.data_error is not None:
            self.data_error = self.data_error.transpose(neworder)
        return self
    def plot_labels(self,labels,fmt = None,**kwargs_passed):
        r'this only works for one axis now'
        axisname = self.dimlabels[0]
        if fmt is None:
            plot_label_points(self.getaxis(axisname),self.data,labels,**kwargs_passed)
        else:
            plot_label_points(self.getaxis(axisname),self.data,[fmt%j for j in labels],**kwargs_passed)
        return
    def labels(self,*args):
        r'''label the dimensions, given in listofstrings with the axis labels given in listofaxes -- listofaxes must be a numpy array;
        you can pass either a dictionary or a axis name (string)/axis label (numpy array) pair'''
        if len(args) == 2:
            listofstrings,listofaxes = args
        elif len(args) == 1 and isinstance(args[0], dict):
            listofstrings = list(args[0].keys())
            listofaxes = list(args[0].values())
        else:
            raise ValueError(strm("I can't figure out how to deal with the arguments",args))
        for j in range(0,len(listofaxes)):
            if isinstance(listofaxes[j], list):
                listofaxes[j] = array(listofaxes[j])
        listofstrings = autostringconvert(listofstrings)
        if isinstance(listofstrings, str):
            listofstrings = [listofstrings]
            listofaxes = [listofaxes]
        if not isinstance(listofstrings, list):
            raise TypeError("the arguments passed to the .labels() method must be a list of the axis names followed by the list of the axis arrays")
        elif all(map(( lambda x: isinstance(x, str_) ),listofstrings)):
            listofstrings = list(map(str,listofstrings))
        elif not all(map(( lambda x: isinstance(x,str) ),listofstrings)):
            raise TypeError("the arguments passed to the .labels() method must be a list of the axis names followed by the list of the axis arrays")
        for j in range(0,len(listofstrings)):
            if listofaxes[j] is None:
                self.setaxis(listofstrings[j],None)
            else:
                #{{{ test that the axis is the right size
                if isscalar(listofaxes[j]):#interpret as a timestep
                    listofaxes[j] = listofaxes[j] * r_[0:ndshape(self)[listofstrings[j]]]
                if type(listofaxes[j]) not in [ndarray,list]:
                    raise TypeError('You passed an axis label of type '+repr(type(listofaxes[j]))+' for the axis '+listofstrings[j]+' to the labels method, which you can\'t do --> it must be an nddata')
                if (len(listofaxes[j]) != ndshape(self)[listofstrings[j]]) and (len(listofaxes[j])!=0):
                    raise IndexError("You're trying to attach an axis of len %d to the '%s' dimension, which has %d data points (shape of self is %s)"%(len(listofaxes[j]),listofstrings[j],ndshape(self)[listofstrings[j]],repr(ndshape(self))))
                #}}}
                self.setaxis(listofstrings[j],listofaxes[j])
        return self
    def check_axis_coords_errors(self):
        if len(self.axis_coords_error) > len(self.dimlabels):
            raise ValueError("this failed because there are more sets of axis errors than there are axes!\nlen(axis_coords_error) = %s\naxes = %s"%(repr(len(self.axis_coords_error)),repr(self.dimlabels)))
    def sort(self,axisname,reverse = False):
        whichaxis = self.dimlabels.index(axisname)
        if reverse:
            order = argsort(-1*self.axis_coords[whichaxis])
        else:
            order = argsort(self.axis_coords[whichaxis])
        datacopy = self.copy()
        for j in range(0,len(order)): # do it this way, so that it deals with other dimensions correctly
            self.check_axis_coords_errors()
            self[axisname,j] = datacopy[axisname,order[j]]
        self.axis_coords[whichaxis] = self.axis_coords[whichaxis][order]
        return self
    def copyaxes(self,other):
        raise ValueError('use copy_axes')
    def copy_axes(self,other):
        # in the case that the dimensions match, and we want to copy the labels
        for thisdim in self.dimlabels:
            if thisdim in other.dimlabels:
                thisax = other.getaxis(thisdim)
                if thisax is not None:
                    thisax = thisax.copy()
                self.setaxis(thisdim,thisax)
                if other.get_error(thisdim) is not None:
                    self.set_error(thisdim, copy(other.get_error(thisdim)))
                if other.get_units(thisdim) is not None:
                    self.set_units(thisdim, other.get_units(thisdim))
        return self
    def axis(self,axisname):
        'returns a 1-D axis for further manipulation'
        thisaxis = self.axn(axisname)
        return nddata(self.getaxis(axisname).copy(),[-1],[axisname]).labels(axisname,self.getaxis(axisname).copy())
    def _axis_inshape(self,axisname):
        newshape = ones(len(self.data.shape),dtype = 'uint')
        thisaxis = self.axn(axisname)
        newshape[thisaxis] = self.data.shape[thisaxis]
        newshape = list(newshape)
        retval = self.getaxis(axisname)
        if retval is None:
            raise AttributeError(axisname+" does not have axis labels!")
        try:
            return retval.copy().reshape(newshape)
        except ValueError as e:
            raise ValueError(strm("Trying to reshape axis from",retval.shape,"to",newshape,"so I can manipulate it like data"))
    def retaxis(self,axisname):
        thisaxis = self._axis_inshape(axisname)
        return nddata(thisaxis,thisaxis.shape,list(self.dimlabels)).labels(axisname,thisaxis.flatten())
    def fromaxis(self,*args,**kwargs):
        '''Generate an nddata object from one of the axis labels.

        Can be used in one of several ways:

        * ``self.fromaxis('axisname')``: Returns an nddata where `retval.data` consists of the given axis values.
        * ``self.fromaxis('axisname',inputfunc)``: use `axisname` as the input for `inputfunc`, and load the result into `retval.data`
        * ``self.fromaxis(inputsymbolic)``: Evaluate `inputsymbolic` and load the result into `retval.data`

        Parameters
        ==========
        axisname : str | list
            The axis (or list of axes) to that is used as the argument of `inputfunc` or the function represented by `inputsymbolic`.
            If this is the only argument, it cannot be a list.
        inputsymbolic : sympy.Expr
            A sympy expression whose only symbols are the names of axes.
            It is preferred, though not required, that this is passed
            without an `axisname` argument -- the axis names are then
            inferred from the symbolic expression.
        inputfunc : function
            A function (typically a lambda function) that taxes the values of the axis given by `axisname` as input.
        overwrite : bool
            Defaults to `False`. If set to `True`, it overwrites `self` with `retval`.
        as_array : bool
            Defaults to `False`. If set to `True`, `retval` is a properly dimensioned numpy ndarray rather than an nddata.

        Returns
        =======
        retval : nddata | ndarray
            An expression calculated from the axis(es) given by `axisname` or inferred from `inputsymbolic`.
        '''
        overwrite,as_array = process_kwargs([('overwrite',False),
            ('as_array',False)],kwargs)
        if len(args) == 1:
            if isinstance(args[0], str):
                axisname = args[0]
                retval = self.retaxis(axisname)
                #{{{ copied from old retaxis function, then added the overwrite capability
                thisaxis = self._axis_inshape(axisname)
                if overwrite:
                    self.data = thisaxis
                    return self
                else:
                    retval = nddata(thisaxis,thisaxis.shape,list(self.dimlabels)).labels(axisname,thisaxis.flatten())
                    retval.axis_coords_units = list(self.axis_coords_units)
                    retval.data_units = self.data_units
                    retval.name(self.name())
                    return retval
                #}}}
            else:
                if issympy(args[0]):
                    func = args[0]
                    symbols_in_func = func.atoms(sympy.Symbol)
                    logger.debug(strm('identified this as a sympy expression (',func,') with symbols',symbols_in_func))
                    symbols_not_in_dimlabels = set(map(str,symbols_in_func))-set(self.dimlabels)
                    if len(symbols_not_in_dimlabels)>0:
                        raise ValueError("You passed a symbolic function, but the symbols"+str(symbols_not_in_dimlabels)+" are not axes")
                else:
                    raise ValueError("I don't know what to do with this type of argument!")
        elif len(args) == 2:
            axisnames = args[0]
            func = args[1]
            if not isinstance(axisnames, list):
                axisnames = [axisnames]
        else:
            raise ValueError('Wrong number of arguments!! -- you passed '+repr(len(args))+' arguments!')
        if issympy(func):
            logging.debug(strm("about to run sympy lambdify, symbols_in_func is",symbols_in_func))
            try:
                lambdified_func = sympy.lambdify(list(symbols_in_func), func,
                        modules=mat2array)
            except Exception as e:
                raise ValueError(strm('Error parsing axis variables',list(map(sympy.var,axisnames)),
                    'that you passed and function',func,'that you passed') +
                    explain_error(e))
            func = lambdified_func
            axisnames = list(map(str,symbols_in_func))
        elif not hasattr(func,'__call__'):
            raise ValueError("I can't interpret the second argument as a function! It is type "+str(type(func)))
        # I can't do the following for sympy, because the argument count is always zero
        if not issympy(args[0]) and func.__code__.co_argcount != len(axisnames):
            raise ValueError(strm("The axisnames you passed",axisnames,
                "and the argument count",func.__code__.co_argcount,"don't match"))
        list_of_axes = [self._axis_inshape(x) for x in axisnames]
        retval = func(*list_of_axes)
        if issympy(retval):
            raise RuntimeError("The sympy function that you passed doesn't match the automatically generated axis variables (obtained by mapping sympy.var onto the axis variables, without any kwargs). The atoms left over are:\n"+str(func.atoms))
        logging.debug(strm("at this point, list of axes is:",list_of_axes))
        if len(list_of_axes) == 0:
            return nddata(float(func()))
            #raise ValueError(strm("Doesn't seem like there are axes -- the axis names that you passed are",axisnames))
        newshape = ones_like(list_of_axes[0].shape)
        for j in list_of_axes:
            newshape *= array(j.shape)
        if overwrite:
            self.data = retval.reshape(newshape)
            return self
        else:
            if as_array:
                return retval.reshape(newshape)
            else:
                retval =  nddata(retval,
                        newshape,
                        list(self.dimlabels)).labels(axisnames,
                                [self.getaxis(x).copy() for x in axisnames])
                retval.axis_coords_units = list(self.axis_coords_units)
                retval.data_units = self.data_units
                retval.name(self.name())
                return retval
    def getaxis(self,axisname):
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
    def extend(self,axis,extent, fill_with=0, tolerance=1e-5):
        r"""If `axis` is uniformly ascending with spacing :math:`dx`,
        then extend by adding a point every :math:`dx` until the axis
        includes the point `extent`.  Fill the newly created datapoints with `fill_with`.

        Parameters
        ----------

        axis : str
            name of the axis to extend
        extent : double
            extend the axis `axis` out to this point
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
        thismsg = "In order to expand, the axis must be equally spaced (and ascending)"
        du = (u[-1] - u[0])/(len(u)-1.)
        assert all(abs(diff(u) - du)/du < tolerance), thismsg# absolute
        # figure out how many points I need to add, and on which side of the axis
        thismsg = "In order to expand, the axis must be ascending (and equally spaced)"
        assert du > 0, thismsg# ascending
        start_index = 0
        stop_index = len(u)
        if extent < u[0]:
            start_index = int(-(u[0] - extent) // du) # the part after the negative is positive
            if (start_index * du + (u[0] - extent))/du < -tolerance:# the first quantity here is negative
                start_index -= 1
        elif extent > u[-1]:
            stop_index_addto = int((extent - u[-1]) // du)
            if ((extent - u[-1]) - du * stop_index_addto)/du > tolerance:# the first quantity here is negative
                stop_index_addto += 1
            stop_index += stop_index_addto
        else:
            raise RuntimeError("extent ({:g}) needs to be further than the bounds on '{:s}', which are {:g} and {:g}".format(extent,axis,u[0],u[-1]))
        #{{{ create a new array, and put self.data into it
        newdata = list(self.data.shape)
        newdata[self.axn(axis)] = stop_index - start_index
        if fill_with == 0:
            newdata = zeros(newdata,dtype = self.data.dtype)
        else:
            newdata = fill_with * ones(newdata,dtype = self.data.dtype)
        newdata_slice = [slice(None,None,None)] * len(newdata.shape)
        newdata_slice[self.axn(axis)] = slice(-start_index,len(u)-start_index,None)
        logger.debug(strm("-------------------------"))
        logger.debug(strm("shape of newdata",newdata.shape))
        logger.debug(strm("shape of self.data",self.data.shape))
        logger.debug(strm("len of u",len(u)))
        logger.debug(strm("start index",start_index))
        logger.debug(strm("shape of slice",newdata[newdata_slice].shape))
        logger.debug(strm("-------------------------"))
        newdata[newdata_slice] = self.data
        self.data = newdata
        #}}}
        # construct the new axis
        new_u = u[0] + du * r_[start_index:stop_index]
        self.setaxis(axis,new_u)
        return self
    def setaxis(self,*args):
        """set or alter the value of the coordinate axis
        
        Can be used in one of several ways:

        * ``self.setaxis('axisname', values)``: just sets the values
        * ``self.setaxis('axisname', '#')``: just
            number the axis in numerically increasing order
            (e.g. if you have smooshed it from a couple
            other dimensions.)
        * ``self.fromaxis('axisname',inputfunc)``: take the existing function, apply inputfunc, and replace
        * ``self.fromaxis(inputsymbolic)``: Evaluate `inputsymbolic` and load the result into the axes, appropriately
        """
        if len(args) == 2:
            axis, value = args
            if isscalar(value) and value=='#':
                self.setaxis(axis,r_[0:ndshape(self)[axis]])
                return self
        elif len(args) == 1 and issympy(args[0]):
            func = args[0]
            symbols_in_func = func.atoms(sympy.Symbol)
            logger.debug(strm('identified this as a sympy expression (',func,') with symbols',symbols_in_func))
            symbols_not_in_dimlabels = set(map(str,symbols_in_func))-set(self.dimlabels)
            if len(symbols_not_in_dimlabels)>0:
                raise ValueError("You passed a symbolic function, but the symbols"+str(symbols_not_in_dimlabels)+" are not axes")
            logging.debug(strm("about to run sympy lambdify, symbols_in_func is",symbols_in_func))
            try:
                lambdified_func = sympy.lambdify(list(symbols_in_func), func,
                        modules=mat2array)
            except Exception as e:
                raise ValueError(strm('Error parsing axis variables',list(map(sympy.var,axisnames)),
                    'that you passed and function',func,'that you passed') +
                    explain_error(e))
            value = lambdified_func
            axis = list(map(str,symbols_in_func))
            assert len(axis)==1, "currently only supported for 1 axis at a time -- if you want to do for more than one axis, please create a pull request with an example"
            axis = axis[0]
        else:
            raise ValueError("not a valid argument to setaxis -- look at the documentation!")
        if axis == 'INDEX':
            raise ValueError("Axes that are called INDEX are special, and you are not allowed to label them!")
        if isinstance(value, type(emptyfunction)):
            x = self.getaxis(axis)
            x[:] = value(x.copy())
            return self
        if type(value) in [float,int,double,float64]:
           value = linspace(0.,value,self.axlen(axis))
        if isinstance(value, list):
            value = array(value)
        if self.axis_coords is None or len(self.axis_coords) == 0:
            self.axis_coords = [None]*len(self.dimlabels)
            self.axis_coords_error = [None]*len(self.dimlabels)
        if value is None:
            self.axis_coords[self.axn(axis)] = None
        else:
            a = len(value)
            b = self.data.shape[self.axn(axis)]
            assert  a == b, "Along the axis %s, the length of the axis you passed (%d) doesn't match the size of the data (%d)."%(axis,a,b)
            self.axis_coords[self.axn(axis)] = value
        return self
    def shear(self, along_axis, propto_axis, shear_amnt,
            zero_fill=True, start_in_conj=False, method='linear'):
        r'''Shear the data :math:`s`:

        :math:`s(x',y,z) = s(x+ay,y,z)`

        where :math:`x` is the `altered_axis` and :math:`y` is the
        `propto_axis`.  (Actually typically 2D, but :math:`z` included
        just to illustrate other dimensions that aren't involved)

        .. note: Unfortunately, currently, when the data is automatically extended,
            if both the start and endpoint of `along_axis` are on the same side
            of zero, some unnecessary padding will be added between the
            beginning of `along_axis` and zero.

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
            :math:`f` *vs.* :math:`t`) along both dimensions from the one that's
            desired.  This means: (1) `self` after the function call transformed
            into the conjugate domain from that before the call and (2)
            `by_amount`, `altered_axis`, and `propto_axis` all refer to the shear
            in the conjugate domain that the data is in at the end of the
            function call.
        '''
        if not (
                self.get_ft_prop(along_axis) is None
                and
                self.get_ft_prop(propto_axis) is None):
            if (self.get_ft_prop(along_axis) ^ self.get_ft_prop(propto_axis)) ^ start_in_conj:
                if start_in_conj:
                    raise ValueError("if you pass start_in_conj, the two dimensions need to be in conjugate domains, but you have: "+self.ft_state_to_str(along_axis,propto_axis))
                else:
                    raise ValueError("(unless you intended to pass start_in_conj) the two dimensions need to be in the same domain, but you have: "+self.ft_state_to_str(along_axis,propto_axis))
        if method == 'fourier':
            return self.fourier_shear(along_axis, propto_axis,
                    shear_amnt, zero_fill=zero_fill)
        elif method == 'linear':
            return self.linear_shear(along_axis, propto_axis,
                    shear_amnt, zero_fill=zero_fill)
        else:
            raise ValueError("The shear method must be either linear or fourier")
    def getaxisshape(self,axisname):
        thishape = ones(len(self.dimlabels))
        thisaxis = self.dimlabels.index(axisname) 
        thishape[thisaxis] = self.data.shape[thisaxis]
        return thishape
    def circshift(self,axis,amount):
        if amount!=0:
            if abs(amount) > ndshape(self)[axis]:
                ValueError(strm("Trying to circshift by ",amount,
                    "which is bitter than the size of",axis))
            newdata = ndshape(self).alloc(dtype=self.data.dtype)
            newdata[axis,:-amount] = self[axis,amount:]
            newdata[axis,-amount:] = self[axis,:amount]
            self.data = newdata.data
        return self
    #}}}
    #{{{ breaking up and combining axes
    def smoosh(self,dimstocollapse, dimname=0, noaxis=False):
        r'''Collapse (smoosh) multiple dimensions into one dimension.

        Parameters
        ----------
        dimstocollapse : list of strings
            the dimensions you want to collapse to one result dimension
        dimname : None, string, integer (default 0)

            if dimname is:

            * None: create a new (direct product) name,
            * a number: an index to the ``dimstocollapse`` list.  The resulting smooshed dimension will be named ``dimstocollapse[dimname]``. Because the default is the number 0, the new dimname will be the first dimname given in the list.
            * a string: the name of the resulting smooshed dimension (can be part of the ``dimstocollapse`` list or not)

        noaxis : bool
            if set, then just skip calculating the axis for the new dimension,
            which otherwise is typically a complicated record array

        Returns
        -------
        self: nddata
            the dimensions `dimstocollapse` are smooshed into a single dimension,
            whose name is determined by `dimname`.
            The axis for the resulting, smooshed dimension is a structured
            array consisting of two fields that give the labels along the
            original axes.

        ..todo::
            when we transition to axes that are stored using a
            slice/linspace-like format, 
            allow for smooshing to determine a new axes that is standard
            (not a structured array) and that increases linearly.
        '''
        assert (type(dimstocollapse) in [list,tuple]) and len(dimstocollapse)>1, "What?? You must try to collapse more than one dimension!! -- you claim you want to collapse '%s'"%str(dimstocollapse)
        not_present = set(dimstocollapse) - set(self.dimlabels)
        if len(not_present) > 0: raise ValueError(strm(not_present,"was not found in the list of dimensions",self.dimlabels))
        #{{{ first, put them all at the end, in order given here
        retained_dims = list(self.dimlabels)
        logger.debug(strm("old order",retained_dims))
        #{{{ if I'm using a dimension here, be sure to grab its current position
        if dimname is None:
            final_position = -1
            dimname = ' $\\times$ '.join(dimstocollapse)
        else:
            if isinstance(dimname, int):
                dimname = dimstocollapse[dimname]
                final_position = self.axn(dimname)
            elif dimname in self.dimlabels:
                final_position = self.axn(dimname)
            else:
                final_position = -1
        #}}}
        # {{{ store the dictionaries for later use
        axis_coords_dict = self.mkd(self.axis_coords)
        axis_coords_error_dict = self.mkd(self.axis_coords_error)
        axis_coords_units_dict = self.mkd(self.axis_coords_units)
        # }}}
        old_units = []
        logger.debug(strm("dims to collapse",dimstocollapse))
        for this_name in dimstocollapse:
            this_idx = retained_dims.index(this_name)
            retained_dims.pop(this_idx)
            axis_coords_error_dict.pop(this_name)
            old_units.append(axis_coords_units_dict.pop(this_name))
            axis_coords_dict.pop(this_name)
            if this_idx < final_position:
                final_position -= 1
        logger.debug(strm("old units",old_units))
        new_units = list(set(old_units))
        if len(new_units) > 1:
            new_units = ' '.join(map(str,new_units))
        elif new_units == 1:
            new_units = new_units[0]
        else:
            new_units = None
        # this might be sub-optimal, but put the dims to collapse at the end, and move them back later if we want
        new_order = retained_dims + dimstocollapse
        self.reorder(new_order)
        logger.debug(strm("new order",new_order))
        #}}}
        #{{{ then, reshape the data (and error)
        logger.debug(strm("old shape",self.data.shape))
        new_shape = list(self.data.shape)[:-len(dimstocollapse)]
        logger.debug(strm("dimensions to keep",new_shape))
        dimstocollapse_shapes = array(self.data.shape[-len(dimstocollapse):])
        new_shape += [dimstocollapse_shapes.prod()]
        self.data = self.data.reshape(new_shape)
        if self.get_error() is not None:
            self.set_error(self.get_error().reshape(new_shape))
        logger.debug(strm("new shape",self.data.shape))
        #}}}
        #{{{ now for the tricky part -- deal with the axis labels
        #{{{ in order, make a list of the relevant axis names, dtypes, and sizes
        axes_with_labels = [j for j in dimstocollapse if self.getaxis(j) is not None] # specifically, I am only concerned with the ones I am collapsing that have labels
        if noaxis:
            logger.debug('noaxis was specified')
            axis_coords_dict[dimname] = None
            axis_coords_error_dict[dimname] = None
        else:
            logger.debug('starting construction of the smooshed axis')
            axes_with_labels_haserror = [self.get_error(j) is not None for j in axes_with_labels]
            axes_with_labels_dtype = [(j,self.getaxis(j).dtype) for
                    j in axes_with_labels]# an appropriate spec. for a structured array
            axes_with_labels_size = [self.getaxis(j).size for j in axes_with_labels]
            #}}}
            logger.debug(strm("the dtype that I want is:",axes_with_labels_dtype))
            logger.debug(strm("the axes that have labels are:",axes_with_labels))
            logger.debug(strm("the axes that have labels have sizes:",axes_with_labels_size))
            # {{{ we construct a multidimensional axis
            multidim_axis_error = None
            if len(axes_with_labels_dtype) > 0:
                # create a new axis of the appropriate shape and size
                multidim_axis_label = empty(axes_with_labels_size,
                        dtype=axes_with_labels_dtype)
                if any(axes_with_labels_haserror):
                    multidim_axis_error = empty(axes_with_labels_size,
                            dtype = [(axes_with_labels[j],
                                self.getaxis(axes_with_labels[j]).dtype)
                                for j in range(len(axes_with_labels))
                                if axes_with_labels_haserror[j]])
                # one at a time index the relevant dimension, and load in the information
                full_slice = [slice(None,None,None)]*len(axes_with_labels_dtype)
                for this_index,thisdim in enumerate(axes_with_labels):
                    axis_for_thisdim = self.getaxis(thisdim)
                    if axes_with_labels_haserror[this_index]:
                        axis_error_for_thisdim = self.get_error(thisdim)
                    logger.debug(strm("the axis for",thisdim,"is",axis_for_thisdim))
                    for j in range(axes_with_labels_size[this_index]):
                        this_slice = list(full_slice)
                        this_slice[this_index] = j # set this element
                        multidim_axis_label[thisdim][tuple(this_slice)] = axis_for_thisdim[j]
                        if axes_with_labels_haserror[this_index]:
                            multidim_axis_error[thisdim][tuple(this_slice)] = axis_error_for_thisdim[j]
                logger.debug(strm("shape of multidim_axis_label is now",multidim_axis_label.shape,"(",axes_with_labels,")"))
                logger.debug(strm("multidim_axis_label is:\n",repr(multidim_axis_label)))
                multidim_axis_label = multidim_axis_label.flatten() # then flatten the axis
                logger.debug(strm("shape of multidim_axis_label is now",multidim_axis_label.shape))
                logger.debug(strm("multidim_axis_label is:\n",repr(multidim_axis_label)))
            else:
                raise ValueError("You requested that smoosh generate an axis, but I don't know what dtype to assign to it (what fields to use).  This is likely because you don't have axes assigned to the dimensions you're trying to smoosh.  Consider calling smoosh with noaxis=True, instead")
            axis_coords_dict[dimname] = multidim_axis_label
            axis_coords_error_dict[dimname] = multidim_axis_error
            # }}}
        #{{{ update axis dictionary with the new info
        logger.debug(strm("end up with axis_coords_dict (%d)"%len(axis_coords_dict),axis_coords_dict))
        logger.debug(strm("end up with axis_coords_error_dict (%d)"%len(axis_coords_error_dict),axis_coords_error_dict))
        #}}}
        #}}}
        #{{{ make new dimlabels, and where relevant, project the new dictionary onto these dimlabels
        axis_coords_units_dict[dimname] = new_units
        self.dimlabels = retained_dims + [dimname]
        logger.debug(strm("end up with dimlabels",self.dimlabels,"and shape",self.data.shape))
        self.axis_coords = self.fld(axis_coords_dict)
        self.axis_coords_error = self.fld(axis_coords_error_dict)
        self.axis_coords_units = self.fld(axis_coords_units_dict)
        logger.debug(strm("new axis coords (%d)"%len(self.axis_coords),self.axis_coords))
        logger.debug(strm("new axis coords errors (%d)"%len(self.axis_coords_error),self.axis_coords_error))
        logger.debug(strm("new axis coords unitss (%d)"%len(self.axis_coords_units),self.axis_coords_units))
        #}}}
        #{{{ then deal with the units
        #}}}
        #{{{ finally, if I need to, reorder again to put the new dimension where I want it
        #}}}
        return self
    def chunk(self,axisin,*otherargs):
        r'''"Chunking" is defined here to be the opposite of taking a direct product, increasing the number of dimensions by the inverse of the process by which taking a direct product decreases the number of dimensions.  This function chunks axisin into multiple new axes arguments.:
            axesout -- gives the names of the output axes
            shapesout -- optional -- if not given, it assumes equal length -- if given, one of the values can be -1, which is assumed length

        When there are axes, it assumes that the axes of the new dimensions
        are nested -- *e.g.*, it will chunk a dimension with axis: 
        [1,2,3,4,5,6,7,8,9,10]
        into dimensions with axes:
        [0,1,2,3,4], [1,6]

        ..todo::
            Deal with this efficiently when we move to new-style axes
        '''
        if len(otherargs) == 2:
            axesout,shapesout = otherargs
        elif len(otherargs) == 1:
            if isinstance(otherargs[0], list):
                axesout = otherargs[0]
                shapesout = ndshape(self)[axisin]**(1./len(axesout))
                if abs(shapesout-round(shapesout)) > 1e-15: # there is some kind of roundoff error here
                    raise ValueError('''In order for chunk to be called with
                            only a list of axes, the shape of the dimension you are
                            trying to split (here %s) must be an Nth root of
                            the original dimension size (here: %d), where N (here
                            %d) is the number of dimensions you are trying to chunk into'''%(axisin,ndshape(self)[axisin],len(axesout)))
                else:
                    shapesout = round(shapesout)
                shapesout = [shapesout] * len(axesout)
            elif isinstance(otherargs[0], dict):
                axesout,shapesout = list(otherargs[0].keys()),list(otherargs[0].values())
            else:
                raise ValueError("I don't know how to deal with this type!")
        else:
            raise ValueError("otherargs must be one or two arguments!")
        assert not any([j in self.dimlabels for j in axesout if j != axisin]), strm(
                "You are trying to create dimensions",[j for j in axesout if j!=axisin],
                "one of which matches one of the existing labels",self.dimlabels)
        if any([j == -1 for j in shapesout]):
            j = shapesout.index(-1)
            if j < len(shapesout)-1:
                shapesout[j] = int(round(ndshape(self)[axisin]/prod(r_[shapesout[0:j],shapesout[j+1:]])))
            else:
                shapesout[j] = int(round(ndshape(self)[axisin]/prod(shapesout[0:j])))
        if prod(shapesout) != ndshape(self)[axisin]:
            raise ValueError("The size of the axis (%s) you're trying to split (%s) doesn't match the size of the axes you're trying to split it into (%s = %s)"%(repr(axisin),repr(ndshape(self)[axisin]),repr(axesout),repr(shapesout)))
        thisaxis = self.axn(axisin)
        if self.getaxis(axisin) is not None and len(self.getaxis(axisin)) > 0:
            axes_tmp = self.getaxis(axisin).reshape(shapesout)
            new_axes = []
            for j in range(len(axes_tmp.shape)):
                this_slicer = [0]*len(axes_tmp.shape)
                this_slicer[j] = slice(None,None,None)
                new_axes.append(axes_tmp[this_slicer])
        else:
            new_axes = None
        #{{{ if there is a list of axis coordinates, add in slots for the new axes
        if isinstance(self.axis_coords, list):
            if len(self.axis_coords) == 0:
                self.axis_coords = [None] * len(self.dimlabels)
            for j in range(len(axesout)-1):
                self.axis_coords.insert(thisaxis,None)
        if isinstance(self.axis_coords_error, list):
            if len(self.axis_coords_error) == 0:
                self.axis_coords_error = [None] * len(self.dimlabels)
            for j in range(len(axesout)-1):
                self.axis_coords_error.insert(thisaxis,None)
        if isinstance(self.axis_coords_units, list):
            if len(self.axis_coords_units) == 0:
                self.axis_coords_units = [None] * len(self.dimlabels)
            for j in range(len(axesout)-1):
                self.axis_coords_units.insert(thisaxis,None)
        #}}}
        newshape = list(self.data.shape[0:thisaxis]) + shapesout + list(self.data.shape[thisaxis+1:])
        newshape = list(map(int,newshape))
        newnames = list(self.dimlabels[0:thisaxis]) + axesout + list(self.dimlabels[thisaxis+1:])
        self.data = self.data.reshape(newshape)
        orig_axis_units = self.get_units(axisin)
        self.dimlabels = newnames
        if new_axes is not None:
            for j in range(len(axesout)):
                self.setaxis(axesout[j],new_axes[j])
                self.set_units(axesout[j],orig_axis_units)
        return self
    def chunk_auto(self,axis_name,which_field,dimname = None):
        r'''assuming that axis "axis_name" is currently labeled with a structured array, choose one field ("which_field") of that structured array to generate a new dimension
        Note that for now, by definition, no error is allowed on the axes.
        However, once I upgrade to using structured arrays to handle axis and data errors, I will want to deal with that appropriately here.'''
        def check_data(a):
            "we need this because other things expect dimlabels to be a list of strings"
            if isinstance(a.dimlabels,recarray):
                a.dimlabels = [str(j[0]) if len(j) == 1 else j for j in a.dimlabels.tolist()]
            return a
        axis_number = self.axn(axis_name)
        new_axis,indices = unique(self.getaxis(axis_name)[which_field],
                return_inverse = True) # we are essentially creating a hash table for the axis.  According to numpy documentation, the hash indices that this returns should also be sorted sorted.
        logger.debug(strm("(chunk auto) indices look like this:",indices))
        #{{{ check that there are equal numbers of all the unique new_axis
        index_count = array([count_nonzero(indices == j) for j in range(indices.max()+1)])
        if all(index_count == index_count[0]):
            logger.debug(strm("(chunk auto) Yes, there are equal numbers of all unique new_axis! (Each element of the hash table has been indexed the same number of times.)"))
            #}}}
            #{{{ store the old shape and generate the new shape
            current_shape = list(self.data.shape)
            logger.debug(strm("(chunk auto) old shape -- ",current_shape))
            new_shape = insert(current_shape,axis_number + 1,len(new_axis))
            new_shape[axis_number] /= len(new_axis) # the indices of the hash table become the new dimension
            #}}}
            #{{{ actually reorder the data and error -- perhaps a view would be more efficient here
            old_data = self.data
            has_data_error = not (self.get_error() is None)
            self.data = empty(new_shape,dtype = self.data.dtype)
            if has_data_error:
                old_error = self.get_error()
                self.set_error(empty(new_shape,dtype = self.data.dtype))
            #}}}
            #{{{ adjust all the relevant axis information
            #{{{ generate an axis label along the axis I'm chunking that's stripped of the field that I'm creating a dimension from (i.e. chunking off a new dimension based on) -- because I am independently manipulating the data, I don't use self.getaxis()
            x_strip_current_field = self.axis_coords[axis_number][
                    [j for j in
                        self.axis_coords[axis_number].dtype.names
                        if j != which_field]]
            #}}}
            #{{{ reshape the axis coordinate so that it becomes a 2D array with the new dimension chunked off
            self.axis_coords[axis_number] = empty((len(x_strip_current_field)//len(new_axis),len(new_axis)),
                    dtype = x_strip_current_field.dtype)
            if not (self.get_error(axis_name) is None):
                raise ValueError("Until I do the structured array upgrade chunk_auto will not be able to deal with an axis that has errors.")
            #}}}
            #{{{ everything is now ready to sort the data and residual axis into ordered slots
            #{{{ initialize the slices
            copy_to_slice   = len(new_shape)*[slice(None,None,None)] # this is the memory address inside the new data (where stuff goes)
            copy_from_slice = len(current_shape)*[slice(None,None,None)] # this is the memory address inside the old data (where stuff comes from)
            #}}}
            if has_data_error:
                data_error_location = self.get_error()
            for j in range(len(new_axis)): # j is the index in the hash table
                copy_to_slice[axis_number + 1]     = j
                copy_from_slice[axis_number]       = where(indices == j)[0]
                self.data[copy_to_slice]           = old_data[copy_from_slice]
                if has_data_error:
                    data_error_location[copy_to_slice] = old_error[copy_from_slice]
                logger.debug(strm("(chunk auto) ",j,'matches at',x_strip_current_field[copy_from_slice[axis_number]]))
                self.axis_coords[axis_number][:,j] = x_strip_current_field[copy_from_slice[axis_number]]
            #}}}
            logger.debug(strm("(chunk auto) new axis -- ",self.axis_coords[axis_number]))
            logger.debug(strm("(chunk auto) new shape -- ",self.data.shape))
            #{{{ housekeeping for the various axes + data properties -- should perhaps be possible to do this first, then use .getaxis()
            self.dimlabels.insert(axis_number + 1,which_field)
            self.axis_coords.insert(axis_number + 1,new_axis)
            #{{{ by definition, axis can have neither errors nor units associated, for now.
            self.axis_coords_error.insert(axis_number + 1,None)
            self.axis_coords_units.insert(axis_number + 1,None)
            #}}}
            logger.debug(strm('(chunk auto) the dimensions of ',self.dimlabels[axis_number],'are (?? x ',self.dimlabels[axis_number+1],')=',self.axis_coords[axis_number].shape))
            #}}}
            #}}}
            #{{{ deal appropriately with the "remainder axis" (axis_number)
            if dimname is None:
                remainder_axis_name = '_and_'.join(self.axis_coords[axis_number].dtype.names)
            else:
                remainder_axis_name = dimname
            #{{{ if everything is the same along the dimension that I've just
            # created (which is the second dimension), then get rid of the
            # duplicate labels
            test_axis = self.axis_coords[axis_number].T
            logger.debug(strm("(chunk auto) test axis -- ",test_axis))
            test_axis = ascontiguousarray(test_axis).flatten().view([('',test_axis.dtype)]*test_axis.shape[1])
            if all(test_axis == test_axis[0]):
                self.axis_coords[axis_number] = self.axis_coords[axis_number][:,0].reshape(1,-1)
                logger.debug(strm("(chunk auto) collapsed to", self.axis_coords[axis_number]))
            #}}}
            if self.axis_coords[axis_number].shape[0] == 1:# then this is a "valid" axis -- because, for each position of the new axis, there is only one value of the remainder axis
                self.axis_coords[axis_number] = self.axis_coords[axis_number].reshape(-1)
                self.dimlabels[axis_number] = remainder_axis_name
                if len(self.axis_coords[axis_number].dtype) == 1: # only one field, which by the previous line will be named appropriately, so drop the structured array name
                    new_dtype = self.axis_coords[axis_number].dtype.descr[0][1]
                    self.axis_coords[axis_number] = array(self.axis_coords[axis_number],dtype = new_dtype) # probably more efficiently done with a view, but leave alone for now
                return check_data(self)
            else:
                #{{{ generate an index list to label the remainder axis, and generate a new nddata with the actual values (which are not copied across the new dimension that matches which_field)
                remainder_axis_index_list = r_[0:self.axis_coords[axis_number].shape[0]]
                new_data = nddata(self.axis_coords[axis_number],
                        self.axis_coords[axis_number].shape,
                        [remainder_axis_name,which_field])
                self.axis_coords[axis_number] = remainder_axis_index_list
                new_data.labels([remainder_axis_name,which_field],
                        [self.axis_coords[axis_number].copy(),self.axis_coords[axis_number + 1].copy()])
                self.dimlabels[axis_number] = remainder_axis_name + '_INDEX'
                #}}}
                return check_data(self),check_data(new_data)
            #}}}
        else:
            raise ValueError("Along the axis '"+axis_name+"', the field '"+which_field+"' does not represent an axis that is repeated one or more times!  The counts for how many times each element along the field is used is "+repr(index_count))
            return
    def squeeze(self,return_dropped=False):
        r'''squeeze singleton dimensions
        
        Parameters
        ==========
        return_dropped: bool (default False)
           return a list of the dimensions that were dropped as a second argument 
        Returns
        =======
        self

        return_dropped: list
            (optional, only if return_dropped is True)
        '''
        mask = array(self.data.shape) > 1
        logger.debug(strm(list(zip(mask,self.dimlabels))))
        self.data = self.data.squeeze()
        retval = []
        if isinstance(self.axis_coords, list):
            for k,v in [(self.dimlabels[j],self.axis_coords[j]) for j in range(len(self.dimlabels)) if not mask[j]]:
                retval.append(k)
                if v is not None:
                    self.set_prop(k,v[0])
        self.dimlabels = [v for j,v in enumerate(self.dimlabels) if mask[j]]
        if isinstance(self.axis_coords, list):
            self.axis_coords = [v for j,v in enumerate(self.axis_coords) if mask[j]]
        if isinstance(self.axis_coords_error, list):
            self.axis_coords_error = [v for j,v in enumerate(self.axis_coords_error) if mask[j]]
        if isinstance(self.axis_coords_units, list):
            self.axis_coords_units = [v for j,v in enumerate(self.axis_coords_units) if mask[j]]
        if return_dropped:
            return self, retval
        else:
            return self
    #}}}
    #{{{ messing with data -- get, set, and copy
    def __getslice__(self,*args):
        raise ValueError(strm('getslice! ',args))
    def __setitem__(self,key,val):
        righterrors = None
        logger.debug(strm('key',key))
        if isinstance(key, nddata):
            logger.debug("initially, rightdata appears to be nddata")
            _,B = self.aligndata(key)
            key = B.data # now the next part will handle this
        if isinstance(key, ndarray):# if selector is an ndarray
            logger.debug("initially, rightdata appears to be ndarray")
            if key.dtype is not dtype('bool'):
                raise ValueError("I don't know what to do with an ndarray subscript that has dtype "+repr(key.dtype))
            if key.shape != self.data.shape:
                raise ValueError("The shape of your logical mask "
                        +repr(key.shape)
                        +" and the shape of your data "
                        +repr(self.data.shape)
                        +" are not compatible (matching or singleton) -- I really don't think that you want to do this!")
            self.data[key] = val
            return
        if isinstance(val,nddata):
            logger.debug("rightdata appears to be nddata after initial treatment")
            #{{{ reorder so the shapes match
            unshared_indices = list(set(val.dimlabels) ^ set(self.dimlabels))
            shared_indices = list(self.dimlabels)
            if 'INDEX' in unshared_indices:
                unshared_indices.remove('INDEX')
            shared_indices = [j for j in shared_indices
                    if j not in unshared_indices]
            if len(val.dimlabels) != len(shared_indices) or (not all([val.dimlabels[j] == shared_indices[j] for j in range(0,len(shared_indices))])):
                val.reorder(shared_indices)
            #}}}
            rightdata = val.data
            righterrors = val.get_error()
        else: # assume it's an ndarray
            logger.debug("rightdata appears to be ndarray after initial treatment")
            rightdata = val
            #{{{ if I just passed a function, assume that I'm applying some type of data-based mask
            if isinstance(key, type(emptyfunction)):
                thisfunc = key
                self.data[thisfunc(self.data)] = rightdata
                return
            #}}}
            if (not isinstance(rightdata, ndarray)): # in case its a scalar
                rightdata = array([rightdata])
        slicedict,axesdict,errordict,unitsdict = self._parse_slices(key) # pull left index list from parse slices
        leftindex = tuple(self.fld(slicedict))
        rightdata = rightdata.squeeze()
        logger.debug(strm("after squeeze, rightdata has shape",rightdata.shape))
        if len(rightdata.shape) > 0:
            left_shape = shape(self.data[leftindex])
            try:
                self.data[leftindex] = rightdata.reshape(left_shape) # assign the data
            except:
                raise IndexError(strm('ERROR ASSIGNING NDDATA:\n',
                    'self.data.shape:',self.data.shape,
                    'left index',leftindex,'\n',
                    'rightdata.shape:',rightdata.shape,
                    '--> shape of left slice: ',left_shape))
        else:
            self.data[leftindex] = rightdata
        lefterror = self.get_error()
        if lefterror is not None:
            lefterror[leftindex] = righterrors.squeeze()
        return self
    # {{{ standard trig functions
    def __getattribute__(self,arg):
        fundict = {'exp':exp,
                'sin':sin,
                'cos':cos,
                'tan':tan,
                'sinh':sinh,
                'cosh':cosh,
                'tanh':tanh,
                'log':log,
                'log10':log10,
                }
        if arg in list(fundict.keys()):
            argf = fundict[arg]
            def retfun():
                retval = self.copy()
                retval.data = argf(retval.data)
                return retval
            return retfun
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
            This should just set a flag that says "Do not allow this data to be substituted in place,"
            so that if something goes to edit the data in place,
            it instead first makes a copy.

            also here, see `Definition of shallow and deep copy <https://docs.python.org/2/library/copy.html>`_

            (older idea)
            We should offer "N", which generates something like a copy,
            but which is sets the equivalent of "nopop".
            For example, currently, you need to do something like
            ``d.C.argmax('t2')``,
            which is very inefficient, since it copies the whole array.
            So, instead, we should do
            ``d.N.argmax('t2')``, which tells argmax and all other
            functions not to overwrite "self" but to return a new object.
            This would cause things like "run_nopop" to become obsolete.
        """
        return self.copy()
    @C.setter
    def C(self):
        raise ValueError("You can't set the C property -- it's used to generate a copy")
    @property
    def angle(self):
        "Return the angle component of the data"
        retval = self.copy(data=False)
        retval.data = angle(self.data)
        return retval
    @angle.setter
    def angle(self):
        raise ValueError("Can't independently set the angle component yet")
    @property
    def imag(self):
        "Return the imag component of the data"
        retval = self.copy(data=False)
        retval.data = self.data.imag
        return retval
    @imag.setter
    def imag(self):
        raise ValueError("Can't independently set the imag component yet")
    @property
    def real(self):
        "Return the real component of the data"
        retval = self.copy(data=False)
        retval.data = self.data.real
        return retval
    @real.setter
    def real(self):
        raise ValueError("Can't independently set the real component yet")
    # }}}
    def copy(self,data=True):
        r'''Return a full copy of this instance.
        
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
        '''
        if data:
            return deepcopy(self)
        else:
            retval =  nddata(0) # empty
            # {{{ data info
            retval.dimlabels = list(self.dimlabels)
            retval.data = None
            retval.data_error = None
            if hasattr(self,'data_units'):
                retval.data_units = deepcopy(self.data_units)
            if hasattr(self,'data_covariance'):
                retval.data_covariance = deepcopy(self.data_covariance)
            # }}}
            # {{{ axes
            retval.axis_coords = deepcopy(self.axis_coords)
            retval.axis_coords_error = deepcopy(self.axis_coords_error)
            retval.axis_coords_units = deepcopy(self.axis_coords_units)
            # }}}
            retval.other_info = deepcopy(self.other_info)
            return retval
    def set_to(self,otherinst):
        r'''Set data inside the current instance to that of the other instance.

        Goes through the list of attributes specified in copy,
        and assigns them to the element of the current instance.

        This is to be used:

        *   for constructing classes that inherit nddata with additional methods.
        *   for overwriting the current data with the result of a slicing operation
        '''
        self.data = otherinst.data
        self.dimlabels = otherinst.dimlabels
        self.data_error = otherinst.data_error
        if hasattr(otherinst,'data_units'):
            self.data_units = otherinst.data_units
        if hasattr(otherinst,'data_covariance'):
            self.data_covariance = otherinst.data_covariance
        self.axis_coords = otherinst.axis_coords
        self.axis_coords_error = otherinst.axis_coords_error
        self.axis_coords_units = otherinst.axis_coords_units
        self.other_info = otherinst.other_info
        return self
    def like(self,value):
        r'''provide "zeros_like" and "ones_like" functionality

        Parameters
        ==========
        value: float
            1 is "ones_like" 0 is "zeros_like", etc.
        '''
        retval = self.copy(data=False)
        retval.data = empty_like(self.data)
        retval.data[:] = value
        return retval
    def __getitem__(self,args):
        if isinstance(args, type(emptyfunction)):
            #{{{ just a lambda function operates on the data
            thisfunc = args
            newdata = self.copy()
            mask = thisfunc(newdata.data)
            newdata.data = newdata.data[mask]
            if len(newdata.dimlabels) == 1:
                x = newdata.getaxis(newdata.dimlabels[0])
                newdata.setaxis(newdata.dimlabels[0],x[mask])
            else:
                raise ValueError("I don't know how to do this for multidimensional data yet!")
            return newdata
            #}}}
        elif isinstance(args, nddata):
            #{{{ try the nddata mask
            A = args
            if isinstance(A,nddata) and A.data.dtype is dtype("bool"):
                thisshape = ndshape(A)
                nonsingleton = []
                for thisdim in A.dimlabels:
                    if thisshape[thisdim] != 1:
                        nonsingleton.append(thisdim)
                if len(nonsingleton) != 1:
                    raise ValueError("To index with an nddata, you must have only one dimension")
                else:
                    self.setaxis(nonsingleton[0],self.getaxis(nonsingleton[0])[A.data.flatten()])
                _,B = self.aligndata(A)
                A = B.data # now the next part will handle this
                if A.dtype is not dtype('bool'):
                    raise ValueError("I don't know what to do with an ndarray subscript that has dtype "+repr(A.dtype))
                if A.shape != self.data.shape:
                    temp = array(A.shape) == 1
                    if all( array(A.shape)[temp] == array(self.data.shape)[temp]):
                        pass
                    else:
                        raise ValueError("The shape of your logical mask "+repr(A.shape)+" and the shape of your data "+repr(self.data.shape)+" are not compatible (matching or singleton) -- I really don't think that you want to do this!")
                self.data = self.data[A]
                return self
            else:
                errmsg = "you passed a single argument of type "+repr(type(A))
                if isinstance(A, nddata):
                    errmsg += " with dtype "+repr(A.data.dtype)
                errmsg += " -- I don't know what to do with this"
                raise ValueError(errmsg)
            #}}}
        else:
            slicedict,axesdict,errordict,unitsdict = self._parse_slices(args)
            if not isinstance(args, slice) and isinstance(args[1], list) and isinstance(args[0], str) and len(args) == 2:
                return concat([self[args[0],x] for x in args[1]],args[0])
            indexlist = tuple(self.fld(slicedict))
            newlabels = [x for x in self.dimlabels if not isscalar(slicedict[x])] # generate the new list of labels, in order, for all dimensions that are not indexed by a scalar
        #{{{ properly index the data error
        if self.data_error is not None:
            try:
                newerror = self.data_error[indexlist]
            except:
                raise ValueError('Problem trying to index data_error'+repr(self.data_error)+' with',repr(indexlist))
        else:
            newerror = None
        #}}}
        if len(self.axis_coords)>0:
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
            except Exception as e:
                raise ValueError(strm("the slice values that you've passed",
                    "don't seem to match the size of the data",
                    "the shape of the data is",self.data.shape,
                    "and the index list (the slice indeces passed to the",
                    "underlying numpy data) I generate from this command is",
                    indexlist,
                    "likely, one of the slice indeces is out of bounds for the size of the data"))
            try:
                retval =  nddata(sliced_data,
                        sliced_data.shape,
                        newlabels,
                        axis_coords = [axesdict[x] for x in newlabels],
                        axis_coords_error = axis_coords_error,
                        data_error = newerror,
                        other_info = self.other_info)
            except Exception as e:
                raise ValueError(strm("likely some problem recasting the data when"
                    "trying to initialize a new nddata: shape of"
                    "self.data",self.data.shape,"indexlist",indexlist))
            retval.axis_coords_units = axis_coords_units
            retval.data_units = self.data_units
            return retval
        else:
            retval = nddata(self.data[indexlist],
                    self.data[indexlist].shape,
                    newlabels,
                    other_info = self.other_info)
            retval.axis_coords_units = self.axis_coords_units
            retval.data_units = self.data_units
            return retval
    def _possibly_one_axis(self,*args):
        if len(args) == 1:
            return args[0]
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        if len(args) == 0:
            if len(self.dimlabels) == 1:
                axes = self.dimlabels
            elif len(self.dimlabels) == 0:
                raise ValueError("You're trying to do something to data with no dimensions")
            else:
                raise ValueError("If you have more than one dimension, you need to tell me which one!!")
        return axes
    def _parse_slices(self,args):
        """This controls nddata slicing:
            it previously took
            \'axisname\',value
            pairs where value was an index or a lambda function.
            Now, it also takes
            \'axisname\':value
            and
            \'axisname\':(value1,value2)
            pairs, where the values give either a single value or an inclusive range on the axis, respectively"""
        logger.debug(strm("about to start parsing slices",args,"for data with axis_coords of length",len(self.axis_coords),"and dimlabels",self.dimlabels,"for ndshape of",ndshape(self)))
        errordict = None # in case it's not set
        if self.axis_coords_units is not None:
            unitsdict = self.mkd(self.axis_coords_units)
        axesdict = None # in case it's not set
        if isinstance(args, slice):
            args = [args]
        # {{{ make a sensible list of tuples that's easier to understand
        sensible_list = [] # type, dimension, arguments
        testf = lambda x: x+1
        j = 0
        while j < len(args):
            if isinstance(args[j],str):# works for str and str_
                dimname = args[j]
                if isinstance(dimname, str_):
                    dimname = str(dimname) # on upgrading + using on windows, this became necessary, for some reason I don't understand
                elif isinstance(args[j+1],type(testf)):
                    sensible_list.append((hash('func'),dimname,args[j+1]))
                else:
                    sensible_list.append((hash('np'),dimname,args[j+1]))
                j += 2
            elif type(args[j]) is slice:
                dimname = args[j].start
                if isinstance(dimname, str_):
                    dimname = str(dimname)
                target = args[j].stop
                if isscalar(target):
                    sensible_list.append((hash('idx'),dimname,target))
                elif type(target) in [tuple,list]:
                    assert len(target) in [1,2], strm("for",args[j],"I expected a 'dimname':(range_start,range_stop)")
                    if len(target) == 1:
                        sensible_list.append((hash('range'),dimname,target[0],None))
                    else:
                        sensible_list.append((hash('range'),dimname,target[0],target[1]))
                j += 1
            else:# works for str and str_
                raise ValueError("I have read in slice argument",args[:j],"but then I get confused!")
        def pprint(a):
            b = {hash(j):j for j in ['idx','range','np']}
            return (b[a[0]],)+a[1:]
        logger.debug(strm('Here is the sensible list:',[pprint(j) for j in sensible_list]))
        # }}}
        if type(args) in [float,int32,int,double]:
            raise ValueError(strm('You tried to pass just a nddata[',type(args),']'))
        if isinstance(args[0], str) or isinstance(args[0], slice):
            #{{{ create a slicedict and errordict to store the slices
            slicedict = dict(list(zip(list(self.dimlabels),[slice(None,None,None)]*len(self.dimlabels)))) #initialize to all none
            if len(self.axis_coords)>0:
                logger.debug(strm("trying to make dictionaries from axis coords of len",len(self.axis_coords),"and axis_coords_error of len",len(self.axis_coords_error),"when dimlabels has len",len(self.dimlabels)))
                axesdict = self.mkd(self.axis_coords)
                if len(self.axis_coords_error)>0:
                    errordict = self.mkd(self.axis_coords_error)
            else:
                axesdict = self.mkd(self.axis_coords)
                logger.debug(strm("length of axis_coords not greater than 0, generated dictionary",axesdict))
            #}}}
            #{{{ map the slices onto the axis coordinates and errors
            for thistuple in sensible_list:
                thisop = thistuple[0]
                thisdim = thistuple[1]
                thisargs = thistuple[2:]
                #print "DEBUG, type of slice",x,"is",type(y)
                if thisop == hash('np'):
                    slicedict[thisdim] = thisargs[0]
                    if isscalar(thisargs[0]):
                        axesdict.pop(thisdim) # pop the axes for all scalar dimensions
                    else:
                        if axesdict[thisdim] is not None:
                            axesdict[thisdim] = axesdict[thisdim][slicedict[thisdim]]
                elif thisop == hash('func'):
                    mask = thisargs[0](axesdict[thisdim])
                    slicedict[thisdim] = mask
                    if axesdict[thisdim] is not None:
                        axesdict[thisdim] = axesdict[thisdim][mask]
                elif thisop == hash('range'):
                    if axesdict[thisdim] is None:
                        raise ValueError("You passed a range-type slice"
                        +" selection, but to do that, your axis coordinates need to"
                        +f" be labeled! (The axis coordinates of {thisdim} aren't"
                        +" labeled)")
                    temp = diff(axesdict[thisdim]) 
                    if not all(temp*sign(temp[0])>0):
                        raise ValueError(strm("you can only use the range format on data where the axis is in consecutively increasing or decreasing order, and the differences that I see are",temp*sign(temp[0])),
                                "if you like, you can still do this by first calling .sort( on the %s axis"%thisdim)
                    if sign(temp[0]) == -1:
                        thisaxis = axesdict[thisdim][::-1]
                    else:
                        thisaxis = axesdict[thisdim]
                    if len(thisargs) > 2:
                        raise ValueError("range with more than two values not currently supported")
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
                            temp_low,temp_high = temp_high,temp_low
                    # at this point, temp_low is indeed the lower value, and temp_high indeed the higher
                    if temp_low == inf:
                        raise ValueError(strm("this is not going to work -- I interpret range",thisargs,"I get to",temp_low,",",temp_high))
                    elif temp_low == -inf:
                        temp_low = 0
                    else:
                        logger.debug(strm("looking for",temp_low))
                        temp_low = searchsorted(thisaxis,temp_low)
                        if temp_low >= len(thisaxis):
                            raise ValueError("the lower value of your slice on the %s axis is higher than the highest value of the axis coordinates!"%thisdim)
                        logger.debug(strm("i found",thisaxis[temp_low],"for the low end of the slice",
                            thisargs))
                    temp_high_float = temp_high
                    if temp_high == inf:
                        temp_high = len(thisaxis)-1
                    elif temp_high == -inf:
                        raise ValueError(strm("this is not going to work -- I interpret range",thisargs,"I get to",temp_low,",",temp_high))
                    else:
                        logger.debug(strm("looking for",temp_high))
                        temp_high = searchsorted(thisaxis,temp_high)
                    # at this point, the result is inclusive if temp_high is
                    # not an exact match, but exclusive if it is
                    if temp_high<len(thisaxis) and thisaxis[temp_high] == temp_high_float:
                        temp_high += 1 # make it inclusive
                    if sign(temp[0]) == -1:
                        temp_high = len(thisaxis) -1 -temp_high
                        temp_low = len(thisaxis) -1 -temp_low
                        temp_high, temp_low = temp_low, temp_high
                    del temp
                    if temp_low == temp_high:
                        temp_high += 1
                    slicedict[thisdim] = slice(temp_low,temp_high,None) # inclusive
                    axesdict[thisdim] = axesdict[thisdim][slicedict[thisdim]]
                elif thisop == hash('idx'):
                    if axesdict[thisdim] is None:
                        raise ValueError("You passed a labeled index"
                        +" selection, but to do that, your axis coordinates need to"
                        +f" be labeled! (The axis coordinates of {thisdim} aren't"
                        +" labeled)")
                    temp = abs(axesdict[thisdim] - thisargs[0]).argmin()
                    slicedict[thisdim] = temp
                    axesdict.pop(thisdim)
            if errordict is not None and errordict != array(None):
                for x,y in slicedict.items():
                    if errordict[x] is not None:
                        if isscalar(y):
                            errordict.pop(x)
                        elif isinstance(y, type(emptyfunction)):
                            mask = y(axesdict[x])
                            errordict[x] = errordict[x][mask]
                        else:
                            try:
                                errordict[x] = errordict[x][y] # default
                            except:
                                raise IndexError(strm('Trying to index',
                                        errordict,'-->',x,'=',errordict[x],'with',y,
                                        'error started as',self.axis_coords_error))
            if unitsdict is not None and unitsdict != array(None):
                for x,y in slicedict.items():
                    if unitsdict[x] is not None:
                        if isscalar(y):
                            unitsdict.pop(x)
            logger.debug(strm('Here is the slice dict:',slicedict,"and the axes dict",axesdict))
            return slicedict,axesdict,errordict,unitsdict
            #}}}
        else:
            raise ValueError(strm('label your freaking dimensions! (type of args[0] is ',
                type(args[0]),'and it should be str!)'))
    #}}}
    #{{{ hdf5 write
    def hdf5_write(self, h5path, directory='.'):
        r"""Write the nddata to an HDF5 file.

        `h5path` is the name of the file followed by the node path where
        you want to put it -- it does **not** include the directory where
        the file lives.
        The directory can be passed to the `directory` argument.

        You can use either :func:`~pyspecdata.find_file` or
        :func:`~pyspecdata.nddata_hdf5` to read the data, as shown below.
        When reading this, please note that HDF5 files store *multiple* datasets,
        and each is named (here, the name is `test_data`).

        .. code-block:: python

            from pyspecdata import *
            init_logging('debug')
            a = nddata(r_[0:5:10j], 'x')
            a.name('test_data')
            try:
                a.hdf5_write('example.h5',getDATADIR(exp_type='Sam'))
            except:
                print("file already exists, not creating again -- delete the file or node if wanted")
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
                raise ValueError(strm("The axis",thisax,"appears not to have a label!  I refuse to save data to HDF5 if you do not label all your axes!!"))
        #{{{ add the final node based on the name stored in the nddata structure
        if h5path[-1] != '/': h5path += '/' # make sure it ends in a slash first
        try:
            thisname = self.get_prop('name')
        except:
            raise ValueError(strm("You're trying to save an nddata object which",
                    "does not yet have a name, and you can't do this! Run",
                    "yourobject.name('setname')"))
        if isinstance(thisname, str):
            h5path += thisname
        else:
            raise ValueError(strm("problem trying to store HDF5 file; you need to",
                "set the ``name'' property of the nddata object to a string",
                "first!"))
        h5file,bottomnode = h5nodebypath(h5path, directory=directory) # open the file and move to the right node
        try:
            #print 'DEBUG 1: bottomnode is',bottomnode
            #}}}
            #{{{ print out the attributes of the data
            myattrs = normal_attrs(self)
            #{{{ separate them into data and axes
            mydataattrs = list(filter((lambda x: x[0:4] == 'data'),myattrs))
            myotherattrs = list(filter((lambda x: x[0:4] != 'data'),myattrs))
            myotherattrs = [x for x in myotherattrs if x not in ['C','sin','cos','exp','log10']]
            myaxisattrs = list(filter((lambda x: x[0:4] == 'axis'),myotherattrs))
            myotherattrs = list(filter((lambda x: x[0:4] != 'axis'),myotherattrs))
            logger.debug(strm(lsafe('data attributes:',list(zip(mydataattrs,[type(self.__getattribute__(x)) for x in mydataattrs]))),'\n\n'))
            logger.debug(strm(lsafe('axis attributes:',list(zip(myaxisattrs,[type(self.__getattribute__(x)) for x in myaxisattrs]))),'\n\n'))
            logger.debug(strm(lsafe('other attributes:',list(zip(myotherattrs,[type(self.__getattribute__(x)) for x in myotherattrs]))),'\n\n'))
            #}}}
            #}}}
            #{{{ write the data table
            if 'data' in mydataattrs:
                if 'data_error' in mydataattrs and self.get_error() is not None and len(self.get_error()) > 0:
                    thistable = rec.fromarrays([self.data,self.get_error()],names='data,error')
                    mydataattrs.remove('data_error')
                else:
                    thistable = rec.fromarrays([self.data],names='data')
                mydataattrs.remove('data')
                datatable = h5table(bottomnode,'data',thistable)
                #print 'DEBUG 2: bottomnode is',bottomnode
                #print 'DEBUG 2: datatable is',datatable
                logger.debug(strm("Writing remaining axis attributes\n\n"))
                if len(mydataattrs) > 0:
                    h5attachattributes(datatable,mydataattrs,self)
            else:
                raise ValueError("I can't find the data object when trying to save the HDF5 file!!")
            #}}}
            #{{{ write the axes tables
            if 'axis_coords' in myaxisattrs:
                if len(self.axis_coords) > 0:
                    #{{{ create an 'axes' node
                    axesnode = h5child(bottomnode, # current node
                            'axes', # the child
                            create = True)
                    #}}}
                    for j,axisname in enumerate(self.dimlabels): # make a table for each different dimension
                        myaxisattrsforthisdim = dict([(x,self.__getattribute__(x)[j])
                            for x in list(myaxisattrs) if len(self.__getattribute__(x)) > 0]) # collect the attributes for this dimension and their values
                        logger.debug(strm(lsafe('for axis',axisname,'myaxisattrsforthisdim=',myaxisattrsforthisdim)))
                        if 'axis_coords' in list(myaxisattrsforthisdim.keys()) and myaxisattrsforthisdim['axis_coords'] is not None:
                            if 'axis_coords_error' in list(myaxisattrsforthisdim.keys()) and myaxisattrsforthisdim['axis_coords_error'] is not None and len(myaxisattrsforthisdim['axis_coords_error']) > 0: # this is needed to avoid all errors, though I guess I could use try/except
                                thistable = rec.fromarrays([myaxisattrsforthisdim['axis_coords'],myaxisattrsforthisdim['axis_coords_error']],names='data,error')
                                myaxisattrsforthisdim.pop('axis_coords_error')
                            else:
                                thistable = rec.fromarrays([myaxisattrsforthisdim['axis_coords']],names='data')
                            myaxisattrsforthisdim.pop('axis_coords')
                        datatable = h5table(axesnode,axisname,thistable)
                        #print 'DEBUG 3: axesnode is',axesnode
                        logger.debug(strm("Writing remaining axis attributes for",axisname,"\n\n"))
                        if len(myaxisattrsforthisdim) > 0:
                            h5attachattributes(datatable,list(myaxisattrsforthisdim.keys()),list(myaxisattrsforthisdim.values()))
            #}}}
            #{{{ Check the remaining attributes.
            logger.debug(strm(lsafe('other attributes:',list(zip(myotherattrs,[type(self.__getattribute__(x)) for x in myotherattrs]))),'\n\n'))
            logger.debug(strm("Writing remaining other attributes\n\n"))
            if len(myotherattrs) > 0:
                #print 'DEBUG 4: bottomnode is',bottomnode
                test = repr(bottomnode) # somehow, this prevents it from claiming that the bottomnode is None --> some type of bug?
                h5attachattributes(bottomnode,
                    [j for j in myotherattrs if not self._contains_symbolic(j)],
                    self)
                warnlist = [j for j in myotherattrs if (not self._contains_symbolic(j)) and isinstance(self.__getattribute__(j), dict)]
                #{{{ to avoid pickling, test that none of the attributes I'm trying to write are dictionaries or lists
                if len(warnlist) > 0:
                    print("WARNING!!, attributes",warnlist,"are dictionaries!")
                warnlist = [j for j in myotherattrs if (not self._contains_symbolic(j)) and isinstance(self.__getattribute__(j), list)]
                if len(warnlist) > 0:
                    print("WARNING!!, attributes",warnlist,"are lists!")
                #}}}
                logger.debug(strm(lsafe('other attributes:',list(zip(myotherattrs,[type(self.__getattribute__(x)) for x in myotherattrs]))),'\n\n'))
            #}}}
        finally:
            h5file.close()
    #}}}
class testclass:
    def __getitem__(self,*args,**kwargs):
        print("you called __getitem__ with args",args,"and kwargs",kwargs)
        return
    def __getattribute__(self,*args,**kwargs):
        print("you called __getattribute__ with args",args,"and kwargs",kwargs)
        return
class nddata_hdf5 (nddata):
    def __repr__(self):
        if hasattr(self,'_node_children'):
            return repr(self.datanode)
        else:
            return nddata.__repr__(self)
        atexit.register(self._cleanup)
    def _cleanup(self):
        if hasattr(self,'_node_children'):
            self.h5file.close()
            del self.h5file
            del self.datanode
        return
    def __init__(self,pathstring,directory='.'):
        self.pathstring = pathstring
        #try:
        self.h5file,self.datanode = h5nodebypath(pathstring,
                check_only=True, directory=directory)
        #except BaseException  as e:
        #    raise IndexError("I can't find the node "+pathstring+explain_error(e))
        logger.debug("about to call _init_datanode")
        self._init_datanode(self.datanode)
        atexit.register(self._cleanup)
    def _init_datanode(self,datanode,**kwargs):
        datadict = h5loaddict(datanode)
        #{{{ load the data, and pop it from datadict
        try:
            datarecordarray = datadict['data']['data'] # the table is called data, and the data of the table is called data
            mydata = datarecordarray['data']
        except:
            raise ValueError("I can't find the nddata.data")
        try:
            kwargs.update({'data_error':datarecordarray['error']})
        except:
            logger.debug(strm("No error found\n\n"))
        datadict.pop('data')
        #}}}
        #{{{ be sure to load the dimlabels
        mydimlabels = [j.decode('utf-8') for j in datadict['dimlabels']]
        if len(mydimlabels) == 1:
            if len(mydimlabels[0]) == 1:
                mydimlabels = list([mydimlabels[0][0]]) # for some reason, think I need to do this for length 1
        #}}}
        #{{{ load the axes and pop them from datadict
        datadict.pop('dimlabels')
        if 'axes' in list(datadict.keys()):
            myaxiscoords = [None]*len(mydimlabels)
            myaxiscoordserror = [None]*len(mydimlabels)
            logger.debug(strm("about to read out the various axes:",list(datadict['axes'].keys())))
            for axisname in list(datadict['axes'].keys()):
                try:
                    axisnumber = mydimlabels.index(axisname)
                except AttributeError as e:
                    raise AttributeError(strm('mydimlabels is not in the right format!\nit looks like this:\n',
                        mydimlabels,type(mydimlabels))+explain_error(e))
                except ValueError as e:
                    raise ValueError(strm('mydimlabels is not in the right format!\nit looks like this:\n',
                        mydimlabels,type(mydimlabels))+explain_error(e))
                recordarrayofaxis = datadict['axes'][axisname]['data']
                myaxiscoords[axisnumber] = recordarrayofaxis['data']
                if 'error' in recordarrayofaxis.dtype.names:
                    myaxiscoordserror[axisnumber] = recordarrayofaxis['error']
                datadict['axes'][axisname].pop('data')
                for k in list(datadict['axes'][axisname].keys()):
                    logger.debug(strm("Warning, attribute",k,"of axis table",axisname,"remains, but the code to load this is not yet supported"))
                datadict['axes'].pop(axisname)
            kwargs.update({"axis_coords":myaxiscoords})
            kwargs.update({"axis_coords_error":myaxiscoordserror})
        elif len(mydimlabels)>1:
            raise ValueError("The current version uses the axis labels to"
                    "figure out the shape of the data\nBecause you stored"
                    "unlabeled data, I can\'t figure out the shape of the"
                    "data!!")
            # the reshaping this refers to is done below
        #}}}
        logger.debug(strm("about to initialize data with shape",mydata.shape,"labels",mydimlabels,"and kwargs",kwargs))
        nddata.__init__(self,
                mydata,
                mydata.shape,
                mydimlabels,
                **kwargs)
        #{{{ reshape multidimensional data to match the axes
        if len(mydimlabels)>1:
            det_shape = []
            for thisdimlabel in mydimlabels:
                try:
                    temp = self.getaxis(thisdimlabel)
                except:
                    temp = -1 # no axis is given
                if isinstance(temp, ndarray):
                    temp = len(temp)
                det_shape.append(temp)
            try:
                self.data = self.data.reshape(tuple([len(self.getaxis(x)) for x in mydimlabels]))
            except:
                raise RuntimeError(strm("The data is of shape", self.data.shape,
                    "and I try to reshape it into", tuple([len(self.getaxis(x))
                        for x in mydimlabels]), "corresponding to the dimensions",mydimlabels,"--> this fails!"))
        #}}}
        for remainingattribute in list(datadict.keys()):
            self.__setattr__(remainingattribute,datadict[remainingattribute])
        self.h5file.close()
        del self.h5file
        del self.datanode
        return
#}}}

class ndshape (ndshape_base):
    r'''The ndshape class, including the allocation method''' 
    def alloc(self,dtype='complex128',labels = False,format = 0):
        r'''Use the shape object to allocate an empty nddata object.

        Parameters
        ----------
        labels : 
            Needs documentation
        format : 0, 1, or None
            What goes in the allocated array.
            `None` uses numpy empty.
        '''
        try:
            if format == 0:
                try:
                    emptyar = zeros(tuple(self.shape),dtype=dtype)
                except TypeError:
                    raise TypeError("You passed a type of "+repr(dtype)+", which was likely not understood (you also passed a shape of "+repr(tuple(self.shape))+")")
            elif format == 1:
                emptyar = ones(tuple(self.shape),dtype=dtype)
            elif format is None:
                emptyar = empty(tuple(self.shape),dtype=dtype)
            else:
                emptyar = format*ones(tuple(self.shape),dtype=dtype)
        except TypeError as e:
            raise TypeError(strm('Wrong type for self.shape',list(map(type,self.shape)),'this probably means that you swapped the size and name arguments -- ',self.shape,'should be numbers, not names'))
        retval = nddata(emptyar,self.shape,self.dimlabels)
        if labels:
            retval.labels(self.dimlabels,[double(r_[0:x]) for x in self.shape])
        return retval

#{{{subplot_dim
class subplot_dim():
    def __init__(self,firstdim,seconddim):
        self.num = r_[firstdim,seconddim,0]
    def set(self,args,x='',g=True,y='',t='',a=''):
        if isinstance(args, int):
            number = args
            ax = subplot(*tuple(self.num+r_[0,0,number]))
            xlabel(x)
            ylabel(y)
            title(t)
            grid(g)
        elif (isinstance(args, tuple)) and (len(args) == 3):
            # the second value passed is 
            whichsmall = args[2]
            break_into = args[1]
            number = args[0]
            mydims = self.num*r_[1,break_into,1]+r_[
                    0,0,break_into*(number-1)+whichsmall]
            try:
                ax = subplot(*tuple(mydims))
            except:
                print('failed trying subplots: ', mydims)
                raise
            xlabel(x)
            ylabel(y)
            title(t)
            grid(g)
        else:
            print("problem, need to pass either 1 or 3 arguments to set")
            print('type of args: ',type(args))
        return ax
#}}}
def fa(input,dtype='complex128'):# make a fortran array
    return array(input,order='F',dtype=dtype) # will need transpose reverses the dimensions, since the bracketing still works in C order (inner is last index), but F tells it to store it appropriately in memory
def ndgrid(*input):
    thissize = list([1])
    thissize = thissize * len(input)
    output = list()
    for j in range(0,len(input)):
        tempsize = copy(thissize)
        tempsize[j] = input[j].size
        output.append(input[j].reshape(tempsize))
    return output
def pinvr(C,alpha):
    U,S,V = svd(C,full_matrices=0)
    #print 'U S V shapes:'
    #print U.shape
    #print S.shape
    #print V.shape
    if any(~isfinite(U)):
        raise ValueError('pinvr error, U is not finite')
    if any(~isfinite(V)):
        raise ValueError('pinvr error, V is not finite')
    if any(~isfinite(S)):
        raise ValueError('pinvr error, S is not finite')
    S = diag(S / (S**2 + alpha**2))
    if any(~isfinite(S)):
        raise ValueError('pinvr error, problem with S/(S^2+alpha^2) --> set your regularization higher')
    return dot(conj(transpose(V)),
            dot(S,conj(transpose(U))))
def sech(x):
    return 1./cosh(x)
def spectrogram(waveform,f_start,f_stop,npoints_fdom=40,tdom_div=2):
    #npoints_tdom = int(round(double(waveform.len)/double(npoints_fdom)))*npoints_tdom_mult
    npoints_tdom = waveform.len/tdom_div # this seems to be more legible than above 
    resolution = diff(waveform.x[0:2])

    sigma = abs(f_start-f_stop)/double(npoints_fdom)
    #print "sigma = %f resolution = %f"%(sigma,resolution)
    if sigma<4*resolution:
        sigma = 4*resolution

    waveform.def_filter(sigma,npoints_tdom)# define the filter and number of points for the spectrogram windowing (define the filter such that the points are spaced sigma apart)

    # go through and apply the filter for some range of points

    f_axis = linspace(f_start,f_stop,npoints_fdom)

    specgram = zeros((npoints_fdom,npoints_tdom),dtype="complex128")

    for j in range(0,npoints_fdom):

        t_axis, specgram[j,:] = waveform.do_filter(f_axis[j])
        #plot(t_axis,abs(specgram[j,:])) # leave this in for testing what it does in the fdom
    #image(specgram,y=f_axis/1e6,x=t_axis*1e6) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
    imshow(abs(specgram),extent=(t_axis[0]*1e6,t_axis[-1]*1e6,f_axis[-1]/1e6,f_axis[0]/1e6)) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
    return gca()
image = this_plotting.image.image
def colormap(points,colors,n=256):
    r = interp(linspace(0,1,n),points,colors[:,0].flatten())
    g = interp(linspace(0,1,n),points,colors[:,1].flatten())
    b = interp(linspace(0,1,n),points,colors[:,2].flatten())
    return reshape(r_[r,g,b],(3,n)).T
def myfilter(x,center = 250e3,sigma = 100e3):
    x = (x-center)**2
    x /= sigma**2
    return exp(-x)
#}}}

#{{{ fitdata
class fitdata(nddata):
    r'''Inherits from an nddata and enables curve fitting through use of a sympy expression.

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
    create a jupyter notebook (not checked in, just for your own playing around) with:
    ```
    import sympy as s
    s.init_printing()
    ```
    you can then use `s.symbols(` to create symbols/variables that
    allow you to build the mathematical expression for your fitting
    function
    '''
    def __init__(self,*args,**kwargs):
        #{{{ manual kwargs
        fit_axis = None
        if 'fit_axis' in list(kwargs.keys()):
            fit_axis = kwargs.pop('fit_axis')
        #}}}
        if isinstance(args[0],nddata):
            myattrs = normal_attrs(args[0])
            for j in range(0,len(myattrs)):
                self.__setattr__(myattrs[j],args[0].__getattribute__(myattrs[j]))
            #nddata.__init__(self,
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
            #self.__base_init(*args,**kwargs)
            nddata.__init__(self,*args,**kwargs)
        if fit_axis is None:
            if len(self.dimlabels) == 1:
                fit_axis = self.dimlabels[0]
            else:
                raise IndexError("Right now, we can only auto-determine the fit axis if there is a single axis")
        self.fit_axis = fit_axis
        #{{{ in the class, only store the forced values and indices they are set to
        self.set_to = None
        self.set_indices = None
        self.active_indices = None
        #}}}
        return
    def parameter_derivatives(self,xvals,set = None,set_to = None):
        r'return a matrix containing derivatives of the parameters, can set dict set, or keys set, vals set_to'
        logger.debug(strm('parameter derivatives is called!'))
        if iscomplex(self.data.flatten()[0]):
            print(lsafen('Warning, taking only real part of fitting data!'))
        if isinstance(set, dict):
            set_to = list(set.values())
            set = list(set.keys())
        solution_list = dict([(self.symbolic_dict[k],set_to[j])
            if k in set
            else (self.symbolic_dict[k],self.output(k))
            for j,k in enumerate(self.symbol_list)]) # load into the solution list
        number_of_i = len(xvals)
        parameters = self._active_symbols()
        mydiff_sym = [[]] * len(self.symbolic_vars)
        x = self.symbolic_x
        fprime = zeros([len(parameters),number_of_i])
        for j in range(0,len(parameters)):
            thisvar = self.symbolic_dict[parameters[j]]
            mydiff_sym[j] = sympy.diff(self.symbolic_func,thisvar)
            #print r'$\frac{\partial %s}{\partial %s}=%s$'%(self.function_name,repr(thisvar),sympy.latex(mydiff).replace('$','')),'\n\n'
            try:
                mydiff = mydiff_sym[j].subs(solution_list)
            except Exception as e:
                raise ValueError(strm('error trying to substitute', mydiff_sym[j],
                    'with', solution_list) + explain_error(e))
            try:
                fprime[j,:] = array([complex(mydiff.subs(x,xvals[k])) for k in range(0,len(xvals))])
            except ValueError as e:
                raise ValueError(strm('Trying to set index',j,
                    'shape(fprime)',shape(fprime),
                    'shape(xvals)',shape(xvals),'the thing I\'m trying to',
                    'compute looks like this',
                    [mydiff.subs(x,xvals[k]) for k in range(0,len(xvals))]))
            except Exception as e:
                raise ValueError(strm('Trying to set index',j,
                    'shape(fprime)',shape(fprime),
                    'shape(xvals)',shape(xvals))+explain_error(e))
        return fprime
    @property
    def function_string(self):
        r'''A property of the fitdata class which stores a string
        output of the functional form of the desired fit expression
        provided in func:`functional_form` in LaTeX format'''
        retval = sympy.latex(self.symbolic_expr).replace('$','')
        return r'$f(%s)='%(sympy.latex(self.fit_axis)) + retval + r'$'
    @function_string.setter
    def function_string(self):
        raise ValueError("You cannot set the string directly -- change the functional_form property instead!")
    @property
    def functional_form(self):
        r'''A property of the fitdata class which is set by the user,
        takes as input a sympy expression of the desired fit
        expression'''
        print("Getting symbolic function")
        return self.symbolic_expr
    @functional_form.setter
    def functional_form(self,sym_expr):
        r''' The functional form, given as a sympy expression, to
        which you would like to fit the data.'''
        assert issympy(sym_expr), "for now, the functional form must be a sympy expression!"
        self.symbolic_expr = sym_expr
        #{{{ adapted from fromaxis, trying to adapt the variable
        symbols_in_expr = self.symbolic_expr.atoms(sympy.Symbol)
        #logger.debug(strm('identified this as a sympy expression (',self.symbolic_expr,') with symbols',symbols_in_expr))
        print('identified this as a sympy expression (',self.symbolic_expr,') with symbols',symbols_in_expr)
        symbols_in_expr = set(map(str,symbols_in_expr))
        # the next are the parameters
        self.fit_axis = set(self.dimlabels) & symbols_in_expr
        if len(self.fit_axis) == 0:
            raise ValueError("I can't find any variables that might correspond to a dimension you want to fit along."
                    "The variables are", symbols_in_expr,
                    "and the dimensions are", self.dimlabels)
        elif len(self.fit_axis) > 1:
            raise ValueError("currently only 1D fitting is supported, though this should be easy"
                    "to change -- I see potential fit axes %s"%str(self.fit_axis))
        # the next line gives the parameters
        self.symbolic_vars = symbols_in_expr-self.fit_axis
        # this gets used later in p_ini
        self.number_of_parameters = len(self.symbolic_vars)
        #}}}
        self.fit_axis = list(self.fit_axis)[0]
        # redefine as real to avoid weird piecewise derivatives
        self.fit_axis_sym = sympy.var(self.fit_axis,real=True) 
        self.symbolic_vars = list(self.symbolic_vars)
        self.symbolic_vars.sort() # so I get consistent behavior
        self.symbolic_vars = [sympy.var(j,real=True) for j in self.symbolic_vars]
        self.symbol_list = [str(j) for j in self.symbolic_vars]
        args = self.symbolic_vars + [self.fit_axis]
        self.fitfunc_multiarg = sympy.lambdify(tuple(args), self.symbolic_expr, modules=mat2array)
        def raw_fn(p,x):
            assert len(p)==len(self.symbolic_vars), "length of parameter passed to fitfunc_raw doesn't match number of symbolic parameters"
            return self.fitfunc_multiarg(
                    *tuple([p[j] for j in range(len(self.symbolic_vars))] + [x]))
        self.fitfunc_raw = raw_fn
        # leave the gradient for later
    def analytical_covariance(self):
        r'''Not up to date'''
        covarmatrix = zeros([len(self._active_symbols())]*2)
        #{{{ try this ppt suggestion --> his V is my fprime, but 
        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis))
        dirproductform = False
        if dirproductform:
            sigma = self.get_error()
            f1 = fprime.shape[0]
            f2 = fprime.shape[1]
            fprime1 = fprime.reshape(f1,1,f2) # j index
            fprime2 = fprime.reshape(1,f1,f2) # k index
            fprime_prod = fprime1 * fprime2
            fprime_prod = fprime_prod.reshape(-1,f2).T # direct product form
            try:
                covarmat = dot(pinv(fprime_prod),(sigma**2).reshape(-1,1))
            except ValueError as e:
                raise ValueError(strm('shape of fprime_prod', shape(fprime_prod),
                    'shape of inverse', shape(pinv(fprime_prod)),
                    'shape of sigma', shape(sigma))+explain_error(e))
            covarmatrix = covarmat.reshape(f1,f1)
            for l in range(0,f1): 
                for m in range(0,f1): 
                    if l != m:
                        covarmatrix[l,m] /= 2
        else:
            sigma = self.get_error()
            #covarmatrix = dot(pinv(f),
            #        dot(diag(sigma**2),pinv(f.T)))
            J = matrix(fprime.T)
            #W = matrix(diag(1./sigma**2))
            #S = matrix(diag(sigma**2))
            #if hasattr(self,'data_covariance'):
            #    print "covariance data is present"
            S = matrix(self.get_covariance())
            Omegainv = S**-1
            #S = matrix(diag(sigma**2))
            #G = matrix(diag(1./sigma))
            #G = S**(-1/2) # analog of the above
            #covarmatrix = ((J.T * W * J)**-1) * J.T * W
            print('a')
            minimizer = inv(J.T * Omegainv * J) * J.T * Omegainv
            covarmatrix = minimizer * S * minimizer.T
            #covarmatrix = array(covarmatrix * S * covarmatrix.T)
            #covarmatrix = array((J.T * G.T * G * J)**-1 * J.T * G.T * G * S * G.T * G * J * (J.T * G.T * G * J)**-1)
            #try:
            #    betapremult = (J.T * Omegainv * J)**-1 * J.T * Omegainv
            #except:
            #    print 'from sigma','\n\n',diag(sigma**2),'\n\n','from covarmatrix','\n\n',S,'\n\n'
            #    raise RuntimeError('problem generating estimator (word?)')
            #covarmatrix = array( betapremult * S * betapremult.T)
        #print "shape of fprime",shape(fprime),"shape of fprime_prod",shape(fprime_prod),'sigma = ',sigma,'covarmat=',covarmatrix,'\n'
        #}}}
        # note for this code, that it depends on above code I later moved to  parameter_derivatives
        #for j in range(0,shape(covarmatrix)[0]):
        #    for k in range(0,shape(covarmatrix)[0]):
        #        #mydiff_second = sympy.diff(mydiff_sym[j],self.symbolic_vars[k]).subs(solution_list)
        #        #fdprime = array([mydiff_second.subs(x,xvals[l])/sigma[l] for l in range(0,len(xvals))]) # only divide by sigma once, since there is only one f
        #        #try:
        #        temp = 1.0/(fprime[j,:] * fprime[k,:])
        #        mask = isinf(temp)
        #        covarmatrix[j,k] = mean(sigma[~mask]**2 * temp[~mask])# + 2. * mean(fminusE * fdprime)
        #        #except:
        #        #    raise RuntimeError(strm('Problem multiplying covarmatrix', 'shape(fprime[j,:])',shape(fprime[j,:]), 'shape(fminusE)',shape(fminusE), 'shape(fdprime)',shape(fdprime)))
        #        #if j != k:
        #        #    covarmatrix[j,k] *= 2
        return covarmatrix
    def copy(self): # for some reason, if I don't override this with the same thing, it doesn't override
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
            raise ValueError(strm('length of this_set=', this_set,
                'and set_to', set_to, 'are not the same!'))
        logger.debug("*** *** *** *** ***")
        logger.debug(str(this_set))
        logger.debug("*** *** *** *** ***")
        set_indices = list(map(self.symbol_list.index,this_set)) # calculate indices once for efficiency
        active_mask = ones(len(self.symbol_list),dtype = bool)
        active_mask[set_indices] = False # generate the mask of indices that are actively fit
        return set_indices,set_to,active_mask
    def remove_inactive_p(self,p):
        return p[self.active_mask]
    def add_inactive_p(self,p):
        if self.set_indices is not None:
            #{{{ uncollapse the function
            temp = p.copy()
            p = zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            #}}}
            p[self.set_indices] = self.set_to # then just set the forced values to their given values
        return p
    def fitfunc(self,p,x):
        r"this wraps fitfunc_raw (which gives the actual form of the fit function) to take care of forced variables"
        p = self.add_inactive_p(p)
        return self.fitfunc_raw(p,x)
    def residual(self,p,x,y,sigma):
        '''just the error function'''
        fit = self.fitfunc(p,x)
        #normalization = sum(1.0/sigma)
        #print 'DEBUG: y=',y,'\nfit=',fit,'\nsigma=',sigma,'\n\n'
        sigma[sigma == 0.0] = 1
        try:
            # as noted here: https://stackoverflow.com/questions/6949370/scipy-leastsq-dfun-usage
            # this needs to be fit - y, not vice versa
            retval = (fit-y)/sigma #* normalization
        except ValueError as e:
            raise ValueError(strm('your error (',shape(sigma),
                    ') probably doesn\'t match y (',
                    shape(y),') and fit (',shape(fit),')')
                    + explain_error(e))
        return retval
    def pinv(self,*args,**kwargs):
        retval = self.linear(*args,**kwargs)
        y = retval.data
        yerr = retval.get_error()
        x_axis = retval.dimlabels[0]
        x = retval.getaxis(x_axis)
        nopowerindex = argmax(x)
        mask = logical_not(r_[0:len(x)] == nopowerindex)
        y = y[mask]
        yerr = yerr[mask]
        x = x[mask]
        L = c_[x.reshape((-1,1)),ones((len(x),1))]
        retval = dot(pinv(L,rcond = 1e-17),y)
        logger.debug(r'\label{fig:pinv_figure_text}y=',y,'yerr=',yerr,'%s='%x_axis,x,'L=',L)
        logger.debug('\n\n')
        logger.debug('recalc y = ',dot(L,retval))
        logger.debug('recalc E = ',1.0-1.0/dot(L,retval))
        logger.debug('actual E = ',self.data)
        return retval
    def linear(self,*args,**kwargs):
        r'''return the linear-form function, either smoothly along the fit function, or on the raw data, depending on whether or not the taxis argument is given
        can take optional arguments and pass them on to eval'''
        #print "DEBUG called linear"
        if len(args) == 1:
            taxis = self._taxis(args[0]) # handle integer as well
            return self.linfunc(taxis,self.eval(taxis,**kwargs).data) # if we pass an argument, return the function across the entire time axis passed
        else:
            return self.linfunc(self.getaxis(self.fit_axis),self.data,yerr = self.get_error(),xerr = self.get_error(self.fit_axis)) # otherwise, return the raw data
    def output(self,*name):
        r'''give the fit value of a particular symbol, or a dictionary of all values.

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
        '''
        if not hasattr(self,'fit_coeff') or self.fit_coeff is None:
            return None
        p = self.fit_coeff.copy()
        if self.set_indices is not None:
            #{{{ uncollapse the function
            temp = p.copy()
            p = zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            #}}}
            p[self.set_indices] = self.set_to # then just set the forced values to their given values
            #print "DEBUG trying to uncollapse in fitfunc w/ ",self.symbol_list,"; from",temp,"to",p
        # this should also be generic
        if len(name) == 1:
            try:
                return p[self.symbol_list.index(name[0])]
            except:
                raise ValueError(strm("While running output: couldn't find",
                    name,"in",self.symbol_list))
        elif len(name) == 0:
            return {self.symbol_list[j]:p[j] for j in range(len(p))}
        else:
            raise ValueError(strm("You can't pass",len(name),"arguments to .output()"))
    def _pn(self,name):
        return self.symbol_list.index(name)
    def _active_symbols(self):
        if not hasattr(self,'active_symbols'):
            if self.set_indices is not None:
                self.active_symbols = [x for x in self.symbol_list if self.active_mask[self._pn(x)]]
            else:
                self.active_symbols = list(self.symbol_list)
        return self.active_symbols
    def _pn_active(self,name):
        return self._active_symbols().index(name)
    def covar(self,*names):
        r'''give the covariance for the different symbols'''
        if len(names) == 1:
            names = [names[0],names[0]]
        if self.covariance is not None:
            return self.covariance[self._pn_active(names[0]),
                    self._pn_active(names[1])].copy()
        else:
            return None
    def covarmat(self,*names):
        if (len(names) == 1) and (names[0] == 'recarray'):
            if hasattr(self,'active_mask'):
                active_symbols = [x for x in self.symbol_list if self.active_mask[self._pn(x)]]
            else:
                active_symbols = list(self.symbol_list)
            if len(active_symbols) != self.covariance.shape[0]:
                raise ValueError(strm('length of active symbols',active_symbols,
                    'doesnt match covariance matrix size(',
                    self.covariance.shape[0],')!'))
            recnames = ['labels'] + active_symbols
            recdata = []
            for j in range(0,self.covariance.shape[0]): 
                thisdata = [active_symbols[j]] + list(double(self.covariance[j,:].copy())) # the first index is the row
                recdata.append(make_rec(thisdata,recnames))
            return r_[tuple(recdata)]
        if len(names) > 0:
            indices = list(map(self._pn_active,names)) # slice out only these rows and columns
            return self.covariance[r_[indices],:][:,r_[indices]].copy()
        else:
            try:
                return self.covariance.copy()
            except:
                return zeros([len(self.fit_coeff)]*2,dtype = 'double')
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
            symbol = sympy.latex(self.symbolic_vars[j]).replace('$','')
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
            '\n',retval,'\n',[allsymb[x] for x in argsort(locations)],
            '\n',printfargs))
        return retval
    def settoguess(self):
        'a debugging function, to easily plot the initial guess'
        self.fit_coeff = real(self.guess())
        return self
    def _taxis(self,taxis):
        r'You can enter None, to get the fit along the same range as the data, an integer to give the number of points, or a range of data, which will return with 300 points'
        if taxis is None:
            taxis = self.getaxis(self.fit_axis).copy()
        elif isinstance(taxis, int):
            taxis = linspace(self.getaxis(self.fit_axis).min(),
                    self.getaxis(self.fit_axis).max(),
                    taxis)
        elif not isscalar(taxis) and len(taxis) == 2:
            taxis = linspace(taxis[0],taxis[1],300)
        return taxis
    def eval(self,taxis,set_what = None,set_to = None):
        r'''after we have fit, evaluate the fit function along the axis taxis
        set_what and set_to allow you to forcibly set_what a specific symbol to a
        specific value --> however, this does not affect the class, but only
        the return value'''
        if isinstance(set_what, dict):
            set_to = list(set_what.values())
            set_what = list(set_what.keys())
        taxis = self._taxis(taxis)
        if hasattr(self,'fit_coeff') and self.fit_coeff is not None:
            p = self.fit_coeff.copy()
        else:
            p = array([NaN]*len(self.symbol_list))
        #{{{ LOCALLY apply any forced values
        # changed line below from set to set_what, and now it works
        if set_what is not None:
            if self.set_indices is not None:
                raise ValueError("you're trying to set indices in an eval"
                        " function for a function that was fit constrained; this"
                        " is not currently supported")
            set_indices,set_to,active_mask = self.gen_indices(set_what,set_to)
            p[set_indices] = set_to
        #}}}
        #{{{ make a new, blank array with the fit axis expanded to fit taxis
        newdata = ndshape(self)
        newdata[self.fit_axis] = size(taxis)
        newdata = newdata.alloc()
        newdata.set_plot_color(self.get_plot_color())
        #}}}
        #{{{ keep all axis labels the same, except the expanded one
        newdata.axis_coords = list(newdata.axis_coords)
        newdata.labels([self.fit_axis],list([taxis]))
        #}}}
        newdata.data[:] = self.fitfunc(p,taxis).flatten()
        return newdata
    def makereal(self):
        self.data = real(self.data)
        return
    def rename(self,previous,new):
        if previous == self.fit_axis:
            self.fit_axis = new
        nddata.rename(self,previous,new)
        return self
    def fit(self,set_what = None, set_to = None, force_analytical = False):
        r'''actually run the fit'''
        if isinstance(set_what, dict):
            set_to = list(set_what.values())
            set_what = list(set_what.keys())
        x = self.getaxis(self.fit_axis)
        if iscomplex(self.data.flatten()[0]):
            logger.debug(strm('Warning, taking only real part of fitting data!'))
        y = real(self.data)
        sigma = self.get_error()
        if sigma is None:
            print('{\\bf Warning:} You have no error associated with your plot, and I want to flag this for now\n\n')
            warnings.warn('You have no error associated with your plot, and I want to flag this for now',Warning)
            sigma = ones(shape(y))
        if set_what is None:
            p_ini = self.guess()
        if set_what is not None:
            self.set_indices,self.set_to,self.active_mask = self.gen_indices(set_what,set_to)
            p_ini = self.remove_inactive_p(p_ini)
        leastsq_args = (self.residual, p_ini)
        leastsq_kwargs = {'args':(x,y,sigma),
                    'full_output':True}# 'maxfev':1000*(len(p_ini)+1)}
        p_out,cov,infodict,mesg,success = leastsq(*leastsq_args,**leastsq_kwargs)
        try:
           p_out,cov,infodict,mesg,success = leastsq(*leastsq_args,**leastsq_kwargs)
        #{{{ just give various explicit errors
        except TypeError as err:
            if not isinstance(x, ndarray) and not isinstance(y, ndarray):
                raise TypeError(strm('leastsq failed because the two arrays',
                    "aren\'t of the right",
                    'type','type(x):',type(x),'type(y):',type(y)))
            else:
                if any(shape(x) != shape(y)):
                    raise RuntimeError(strm('leastsq failed because the two arrays do'
                            "not match in size size",
                            'shape(x):',shape(x),
                            'shape(y):',shape(y)))
            raise TypeError(strm('leastsq failed because of a type error!',
                'type(x):',showtype(x),'type(y):',showtype(y),
                'type(sigma)',showtype(sigma),'shape(x):',shape(x),
                'shape(y):',shape(y),'shape(sigma)',shape(sigma),
                'p_ini',type(p_ini),p_ini))
        except ValueError as err:
            raise ValueError(strm('leastsq failed with "',err,
                '", maybe there is something wrong with the input:',
                self))
        except Exception as e:
            raise ValueError('leastsq failed; I don\'t know why')
        #}}}
        if success not in [1,2,3,4]:
            #{{{ up maximum number of evals
            if mesg.find('maxfev'):
                leastsq_kwargs.update({ 'maxfev':50000 })
                p_out,cov,infodict,mesg,success = leastsq(*leastsq_args,**leastsq_kwargs)
                if success != 1:
                    if mesg.find('two consecutive iterates'):
                        print(r'{\Large\color{red}{\bf Warning data is not fit!!! output shown for debug purposes only!}}','\n\n')
                        print(r'{\color{red}{\bf Original message:}',lsafe(mesg),'}','\n\n')
                        infodict_keys = list(infodict.keys())
                        infodict_vals = list(infodict.values())
                        if 'nfev' in infodict_keys:
                            infodict_keys[infodict_keys.index('nfev')] = 'nfev, number of function calls'
                        if 'fvec' in infodict_keys:
                            infodict_keys[infodict_keys.index('fvec')] = 'fvec, the function evaluated at the output'
                        if 'fjac' in infodict_keys:
                            infodict_keys[infodict_keys.index('fjac')] = 'fjac, A permutation of the R matrix of a QR factorization of the final approximate Jacobian matrix, stored column wise. Together with ipvt, the covariance of the estimate can be approximated.'
                        if 'ipvt' in infodict_keys:
                            infodict_keys[infodict_keys.index('ipvt')] = 'ipvt, an integer array of length N which defines a permutation matrix, p, such that fjac*p = q*r, where r is upper triangular with diagonal elements of nonincreasing magnitude.  Column j of p is column ipvt(j) of the identity matrix'
                        if 'qtf' in infodict_keys:
                            infodict_keys[infodict_keys.index('qtf')] = 'qtf, the vector (transpose(q)*fvec)'
                        for k,v in zip(infodict_keys,infodict_vals):
                            print(r'{\color{red}{\bf %s:}%s}'%(k,v),'\n\n')
                        #self.fit_coeff = None
                        #self.settoguess()
                        #return
                    else:
                        raise RuntimeError(strm('leastsq finished with an error message:',mesg))
                    #}}}
            else:
                raise RuntimeError(strm('leastsq finished with an error message:',mesg))
        else:
            logger.debug("Fit finished successfully with a code of %d and a message ``%s''"%(success,mesg))
        self.fit_coeff = p_out # note that this is stored in HIDDEN form
        dof = len(x) - len(p_out)
        if hasattr(self,'symbolic_x') and force_analytical:
            self.covariance = self.analytical_covariance()
        else:
            if force_analytical: raise RuntimeError(strm("I can't take the analytical",
                "covariance!  This is problematic."))
            if cov is None:
                print(r'{\color{red}'+lsafen('cov is none! why?!, x=',x,'y=',y,'sigma=',sigma,'p_out=',p_out,'success=',success,'output:',p_out,cov,infodict,mesg,success),'}\n')
            self.covariance = cov
        if self.covariance is not None:
            try:
                self.covariance *= sum(infodict["fvec"]**2)/dof # scale by chi_v "RMS of residuals"
            except TypeError as e:
                raise TypeError(strm("type(self.covariance)",type(self.covariance),
                    "type(infodict[fvec])",type(infodict["fvec"]),
                    "type(dof)",type(dof)))
        logger.debug(strm("at end of fit covariance is shape",shape(self.covariance),"fit coeff shape",shape(self.fit_coeff)))
        return
    def bootstrap(self,points,swap_out = exp(-1.0),seedval = 10347,minbounds = {},maxbounds = {}):
        print(r'\begin{verbatim}')
        seed(seedval)
        fitparameters = list(self.symbol_list)
        recordlist = array([tuple([0]*len(fitparameters))]*points,
                {'names':tuple(fitparameters),'formats':tuple(['double']*len(fitparameters))}) # make an instance of the recordlist
        for runno in range(0,points):
            success = False # because sometimes this doesn't work
            while success is False:
                thiscopy = self.copy()
                #{{{ discard datapoints
                origsizecheck = double(size(thiscopy.data))
                mask = thiscopy.random_mask(thiscopy.fit_axis,threshold = swap_out)
                thiscopy.data = thiscopy.data[mask]
                derr = thiscopy.get_error()
                x = thiscopy.getaxis(thiscopy.fit_axis)
                x = x[mask] # note that x is probably no longer a pointer
                derr = derr[mask]
                #print 'DEBUG: size of data after cut',double(size(thiscopy.data))/origsizecheck,' (expected ',1.-swap_out,')'
                #}}}
                #{{{ now extend
                number_to_replace = origsizecheck - thiscopy.data.size
                #print 'DEBUG: number_to_replace',number_to_replace
                random_indices = int32((rand(number_to_replace)*(thiscopy.data.size-1.0)).round())
                thiscopy.data = r_[thiscopy.data,thiscopy.data.copy()[random_indices]]
                thiscopy.labels([thiscopy.fit_axis],[r_[x,x.copy()[random_indices]]])
                thiscopy.set_error(r_[derr,derr.copy()[random_indices]])
                #print 'DEBUG: size of data after extension',double(size(thiscopy.data))/origsizecheck
                #}}}
                try:
                    thiscopy.fit()
                    success = True
                    if len(minbounds) > 0:
                        for k,v in minbounds.items():
                            if thiscopy.output(k) < v:
                                success = False
                    if len(maxbounds) > 0:
                        for k,v in maxbounds.items():
                            if thiscopy.output(k) > v:
                                success = False
                except:
                    #print 'WARNING, didn\'t fit'
                    success = False
                # here, use the internal routines, in case there are constraints, etc
                if success is True:
                    for name in thiscopy.symbol_list: # loop over all fit coeff
                        recordlist[runno][name] = thiscopy.output(name)
        print(r'\end{verbatim}')
        return recordlist # collect into a single recordlist array
    def guess(self,use_pseudoinverse=False):
        r'''old code that I am preserving here -- provide the guess for our parameters; by default, based on pseudoinverse'''
        if use_pseudoinverse:
            self.has_grad = False
            if iscomplex(self.data.flatten()[0]):
                print(lsafen('Warning, taking only real part of fitting data!'))
            y = real(self.data)
            # I ended up doing the following, because as it turns out
            # T1 is a bad fit function, because it has a singularity!
            # this is probably why it freaks out if I set this to zero
            # on the other hand, setting a value of one seems to be
            # bad for very short T1 samples
            which_starting_guess = 0
            thisguess = self.starting_guesses[which_starting_guess]
            numguesssteps = 20
            #{{{ for some reason (not sure) adding a dimension to y
            new_y_shape = list(y.shape)
            new_y_shape.append(1)
            y = y.reshape(tuple(new_y_shape))
            #}}}
            #{{{ evaluate f, fprime and residuals
            guess_dict = dict(list(zip(self.symbol_list,list(thisguess))))
            fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
            f_at_guess = real(self.eval(None,set = guess_dict).data)
            try:
                f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
            except:
                raise ValueError(strm('trying to reshape f_at_ini_guess from',f_at_guess.shape,
                    'to',new_y_shape))
            thisresidual = sqrt((y-f_at_guess)**2).sum()
            #}}}
            lastresidual = thisresidual
            for j in range(0,numguesssteps):
                logger.debug('\n\n.core.guess) '+r'\begin{verbatim} fprime = \n',fprime,'\nf_at_guess\n',f_at_guess,'y=\n',y,'\n',r'\end{verbatim}')
                logger.debug('\n\n.core.guess) shape of parameter derivatives',shape(fprime),'shape of output',shape(y),'\n\n')
                regularization_bad = True
                alpha_max = 100.
                alpha_mult = 2.
                alpha = 0.1 # maybe I can rather estimate this based on the change in the residual, similar to in L-M?
                logger.debug(strm('\n\n.core.guess) value of residual before regularization %d:'%j,thisresidual))
                while regularization_bad:
                    newguess = real(array(thisguess) + dot(pinvr(fprime.T,alpha),(y-f_at_guess)).flatten())
                    mask = newguess < self.guess_lb
                    newguess[mask] = self.guess_lb[mask]
                    mask = newguess > self.guess_ub
                    newguess[mask] = self.guess_ub[mask]
                    if any(isnan(newguess)):
                        logger.debug(strm('\n\n.core.guess) Regularization blows up $\\rightarrow$ increasing $\\alpha$ to %0.1f\n\n'%alpha))
                        alpha *= alpha_mult
                    else:
                        #{{{ evaluate f, fprime and residuals
                        guess_dict = dict(list(zip(self.symbol_list,list(newguess))))
                        # only evaluate fprime once we know this is good, below
                        f_at_guess = real(self.eval(None,set = guess_dict).data)
                        try:
                            f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
                        except:
                            raise IndexError(strm('trying to reshape f_at_ini_guess from',
                                f_at_guess.shape,'to',new_y_shape))
                        thisresidual = sqrt((y-f_at_guess)**2).sum()
                        #}}}
                        if (thisresidual-lastresidual)/lastresidual > 0.10:
                            alpha *= alpha_mult
                            logger.debug(strm('\n\n.core.guess) Regularized Pinv gave a step uphill $\\rightarrow$ increasing $\\alpha$ to %0.1f\n\n'%alpha))
                        else: # accept the step
                            regularization_bad = False
                            thisguess = newguess
                            lastresidual = thisresidual
                            fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
                    if alpha > alpha_max:
                        print("\n\n.core.guess) I can't find a new guess without increasing the alpha beyond %d\n\n"%alpha_max)
                        if which_starting_guess >= len(self.starting_guesses)-1:
                            print("\n\n.core.guess) {\\color{red} Warning!!!} ran out of guesses!!!%d\n\n"%alpha_max)
                            return thisguess
                        else:
                            which_starting_guess += 1
                            thisguess = self.starting_guesses[which_starting_guess]
                            print("\n\n.core.guess) try a new starting guess:",lsafen(thisguess))
                            j = 0 # restart the loop
                            #{{{ evaluate f, fprime and residuals for the new starting guess
                            guess_dict = dict(list(zip(self.symbol_list,list(thisguess))))
                            fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
                            f_at_guess = real(self.eval(None,set = guess_dict).data)
                            try:
                                f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
                            except:
                                raise RuntimeError(strm('trying to reshape f_at_ini_guess from',
                                    f_at_guess.shape,'to',new_y_shape))
                            thisresidual = sqrt((y-f_at_guess)**2).sum()
                            #}}}
                            regularization_bad = False # jump out of this loop
                logger.debug(strm('\n\n.core.guess) new value of guess after regularization:',lsafen(newguess)))
                logger.debug(strm('\n\n.core.guess) value of residual after regularization:',thisresidual))
            return thisguess
        else:
            return [1.0]*self.number_of_parameters
#}}}
def sqrt(arg):
    if isinstance(arg,nddata):
        return arg**0.5
    elif isinstance(arg,sympy.symbol.Symbol):
        return sympy.sqrt(arg)
    else:
        return np_sqrt(arg)

# {{{ determine the figure style, and load the appropriate modules
if _figure_mode_setting == 'latex':
    from .fornotebook import *
    figlist_var = figlistl
elif _figure_mode_setting == 'standard':
    def obsn(*x): #because this is used in fornotebook, and I want it defined
        print(''.join(x),'\n')
    def obs(*x): #because this is used in fornotebook, and I want it defined
        print(''.join(map(repr,x)))
    def lrecordarray(*x,**kwargs):
        return repr(x) # if I'm not using tex, it's easier to not use the formatting
    def lsafe(*string,**kwargs):
        "replacement for normal lsafe -- no escaping"
        if len(string) > 1:
            lsafewkargs = lambda x: lsafe(x,**kwargs)
            return ' '.join(list(map(lsafewkargs,string)))
        else:
            string = string[0]
        #{{{ kwargs
        spaces = False
        if 'spaces' in list(kwargs.keys()):
            spaces = kwargs.pop('spaces')
        if 'wrap' in list(kwargs.keys()):
            wrap = kwargs.pop('wrap')
        else:
            wrap = None
        #}}}
        if not isinstance(string, str):
            string = repr(string)
        if wrap is True:
            wrap = 60
        if wrap is not None:
            string = '\n'.join(textwrap.wrap(string,wrap))
        return string
    figlist_var = figlist
else:
    raise ValueError("I don't understand the figures mode "+_figure_mode_setting)
# }}}
