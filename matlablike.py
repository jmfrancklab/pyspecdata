from pylab import *
import textwrap
import matplotlib.transforms as mtransforms
from numpy import sqrt as np_sqrt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.collections import PolyCollection
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D
import tables
import warnings
from inspect import ismethod
from numpy.core import rec
from matplotlib.pyplot import cm
import tables
from os import listdir,environ
from copy import deepcopy 
import traceback
import sympy
from scipy.optimize import leastsq
from scipy.signal import fftconvolve
import scipy.sparse as sparse
import numpy.lib.recfunctions as recf
from inspect import getargspec
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from datadir import getDATADIR
#rc('image',aspect='auto',interpolation='bilinear')
rc('image',aspect='auto',interpolation='nearest')
#rc('text',usetex=True) # this creates all sorts of other problems
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
DATADIR = getDATADIR() 
#{{{ constants
k_B = 1.380648813e-23
hbar = 6.6260695729e-34/2./pi
#}}}

def mybasicfunction(first_figure = None):
    r'''this gives the format for doing the image thing
note also
nextfigure(fl,'name')
and
nextfigure({'lplotproperty':value})
'''
    fl = figlistini_old(first_figure)
    return figlistret(first_figure,figurelist,other)

#{{{ function trickery
def mydiff(data,axis = -1):
    '''this will replace diff with a version that has the same number of indeces, with the last being the copy of the first'''
    newdata = zeros(shape(data),dtype = data.dtype)
    indeces = [slice(None,None,None)]*len(data.shape)
    indeces[axis] = slice(None,-1,None)
    newdata[indeces] = diff(data,axis = axis)
    #setfrom = list(indeces)
    #indeces[axis] = -1
    #setfrom[axis] = 0
    #newdata[indeces] = newdata[setfrom]
    return newdata
#}}}
def normal_attrs(obj):
    myattrs = filter(lambda x: not ismethod(obj.__getattribute__(x)),dir(obj))
    myattrs = filter(lambda x: not x[0:2] == '__',myattrs)
    return myattrs
def showtype(x):
    if type(x) is ndarray:
        return ndarray,x.dtype
    else:
        return type(x)
def emptyfunction():
    pass
#{{{ structured array helper functions
def make_bar_graph_indeces(mystructarray,list_of_text_fields,
        recursion_depth = 0,
        verbose = False,
        spacing = 0.1):
    r"This is a recursive function that is used as part of textlabel_bargraph; it does NOT work without the sorting given at the beginning of that function"
    #{{{ if there are still text fields left, then break down the array further, otherwise, just return the indeces for this subarray
    if len(list_of_text_fields) > 0:
        unique_values = unique(mystructarray[list_of_text_fields[0]])# the return_index argument doesn't do what it's supposed to all the time, so I have to manually find the start indeces, as given in the following line
        start_indeces = [nonzero(mystructarray[list_of_text_fields[0]] == val)[0][0] for val in unique_values]
        # find the structured array for the unique value
        index_values = []
        label_values = []
        start_indeces = r_[start_indeces,len(mystructarray)] # I add this so I can do the next step
        if verbose: print 'recursion depth is',recursion_depth,'and I am analyzing',list_of_text_fields[0],': ',
        if verbose: print 'I found these unique values:',unique_values,'at these start indeces:',start_indeces[:-1]
        for k in range(0,len(start_indeces)-1):
            if verbose: print 'recursion depth is',recursion_depth,'and I am analyzing',list_of_text_fields[0],': ',
            if verbose: print 'trying to extract unique value',unique_values[k],'using the range',start_indeces[k],start_indeces[k+1]
            if verbose: print 'which has this data'
            indiv_struct_array = mystructarray[start_indeces[k]:start_indeces[k+1]]
            if verbose: print lsafen(indiv_struct_array)
            these_index_values,these_labels = make_bar_graph_indeces(indiv_struct_array,list_of_text_fields[1:],recursion_depth = recursion_depth+1,verbose = verbose)
            index_values.append(these_index_values)
            label_values.append([str(unique_values[k])+','+j for j in these_labels])
        #{{{ scale the result of each call down to the equal size (regardless of number of elements), shift by the position in this array, and return
        if verbose: print 'recursion depth is',recursion_depth,'and I just COMPLETED THE LOOP, which gives a list of index values like this',index_values
        max_indeces = max(array(map(len,index_values),dtype='double'))# the maximum width of the array inside
        index_values = map(lambda x: x+(max_indeces-len(x))/2.0,index_values)# if the bar is less than max indeces, shift it over, so it's still in the center
        if verbose: print 'recursion depth is',recursion_depth,'and I centered each set like this',index_values
        index_values = map(lambda x: x/max_indeces*(1-spacing)+(1-spacing)/2,index_values)# scale down, so the width from left edge of bar to right edge of largest bar runs 0--> 1
        if verbose: print 'recursion depth is',recursion_depth,'and I scaled down so each runs zero to one*(1-spacing) (centered) like this',index_values
        # this adds an index value, and also collapses down to a single dimension list
        retval_indeces = [x+num for num,val in enumerate(index_values) for x in val]
        # now collapse labels down to a single dimension
        retval_labels = [k for j in label_values for k in j]
        if verbose: print 'recursion depth is',recursion_depth,'and I am passing up indeces',retval_indeces,'and labels',retval_labels
        return retval_indeces,retval_labels
        #}}}
    else:
        if verbose: print 'recursion depth is',recursion_depth,
        N = len(mystructarray)
        if verbose: print 'hit innermost (no text labels left) and passing up a list of indeces that looks like this:',r_[0:N]
        return r_[0:N],['']*N
    #}}}
def textlabel_bargraph(mystructarray,othersort = None,spacing = 0.1,verbose = False,ax = None,tickfontsize = 8):
    if ax is None:
        thisfig = gcf()
        ax = thisfig.add_axes([0.2,0.5,0.8,0.5])
        try:
            ax.tick_params(axis = 'both',which = 'major',labelsize = tickfontsize)
            ax.tick_params(axis = 'both',which = 'minor',labelsize = tickfontsize)
        except:
            print 'Warning, in this version I can\'t set the tick params method for the axis'
    #{{{ perform the necessary sorting!
    mystructarray = mystructarray.copy()
    list_of_text_fields = [str(j[0]) for j in mystructarray.dtype.descr if j[1][0:2] == '|S']
    mystructarray = mystructarray[list_of_text_fields + [x[0]
        for x in mystructarray.dtype.descr
        if x[0] not in list_of_text_fields]]
    mystructarray.sort()
    if verbose: print 'test --> now, it has this form:',lsafen(mystructarray)
    #}}}
    if othersort is not None:
        list_of_text_fields.append(othersort)
    #list_of_text_fields = ['chemical','run_number']
    if verbose: print 'list of text fields is',lsafen(list_of_text_fields)
    indeces,labels = make_bar_graph_indeces(mystructarray,list_of_text_fields,verbose = verbose,spacing = spacing)
    temp = zip(indeces,labels)
    if verbose: print '(indeces,labels) (len %d):'%len(temp),lsafen(temp)
    if verbose: print 'I get these labels (len %d):'%len(labels),labels,'for the data (len %d)'%len(mystructarray),lsafen(mystructarray)
    indeces = array(indeces)
    indiv_width = min(diff(indeces))*(1-spacing)
    remaining_fields = list(set(mystructarray.dtype.names)^set(list_of_text_fields))
    if verbose: print 'The list of remaining (i.e. non-text) fields is',lsafen(remaining_fields)
    colors = ['b','r','g']
    rects = []
    for j,thisfield in enumerate(remaining_fields):
        field_bar_width = indiv_width/len(remaining_fields)
        try:
            rects.append(ax.bar(indeces+j*field_bar_width,
                    mystructarray[thisfield],
                    field_bar_width,color = colors[j],
                    label = '$%s$'%thisfield))
        except:
            raise CustomError('Problem with bar graph: there are %d indeces, but %d pieces of data'%(len(indeces),len(mystructarray[thisfield])),'indeces:',indeces,'data',mystructarray[thisfield])
    ax.set_xticks(indeces+indiv_width/2)
    ax.set_xticklabels(labels)
    return
def lookup_rec(A,B,indexpair):
    r'''look up information about A in table B (i.e. chemical by index, etc)
    indexpair is either the name of the index
    or -- if it's differently named -- the pair of indeces
    given in (A,B) respectively
    
    This will just drop any fields in B that are also in A,
    and the output uses the first indexname
    
    note that it it seems like the join_rec function above may be more efficient!!'''
    raise RuntimeError('You should now use decorate_rec!!')
    if type(indexpair) not in [tuple,list]:
        indexpair = (indexpair,indexpair)
    Bini = copy(B)
    B = recf.drop_fields(B,( set(B.dtype.names) & set(A.dtype.names) ) - set([indexpair[1]])) # indexpair for B gets dropped later anyways
    joined = []
    for j in A:
        matchedrows =  B[B[indexpair[1]] == j[indexpair[0]]]
        for matchedrow in matchedrows:
            joined.append((j,matchedrow))
    if len(joined) == 0:
        raise CustomError('Unable to find any matches between',A[indexpair[0]],'and',B[indexpair[1]],'!')
    whichisindex = joined[0][1].dtype.names.index(indexpair[1])
    allbutindex = lambda x: list(x)[0:whichisindex]+list(x)[whichisindex+1:]
    joined = concatenate([array(tuple(list(j[0])+allbutindex(j[1])),
                    dtype = dtype(j[0].dtype.descr+allbutindex(j[1].dtype.descr))).reshape(1) for j in joined])
    return joined
def reorder_rec(myarray,listofnames,first = True):
    try:
        indeces_to_move = [myarray.dtype.names.index(j) for j in listofnames]
    except:
        stuff_not_found = [j for j in listofnames if j not in myarray.dtype.names]
        if len(stuff_not_found) > 0:
            raise CustomError(stuff_not_found,'is/are in the list you passed, but not one of the fields, which are',myarray.dtype.names)
        else:
            raise CustomError('unknown problem')
    old_type = list(myarray.dtype.descr)
    new_type = [old_type[j] for j in indeces_to_move] + [old_type[j] for j in range(0,len(old_type)) if j not in indeces_to_move]
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
        raise CustomError("For the fourth argument, you must pass either a list with the names of the arguments, or nothing (to use the field itself as an argument)")
    if type(myargs) is str:
        myargs = (myargs,)
    if type(myargs) is not tuple:
        myargs = tuple(myargs)
    argdata = map((lambda x: myarray[x]),myargs)
    try:
        newrow = myfunction(*tuple(argdata))
    except TypeError:
        newrow = array([myfunction(*tuple([x[rownumber] for x in argdata])) for rownumber in range(0,len(argdata[0]))])
    if type(newrow) is list and type(newrow[0]) is str:
        newrow = array(newrow,dtype = '|S100')
    try:
        new_field_type = list(newrow.dtype.descr[0])
    except AttributeError:
        raise CustomError("evaluated function on",argdata,"and got back",newrow,"which appears not to be a numpy array")
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
def join_rec((A,a_ind),(B,b_ind)):
    raise RuntimeError('You should now use decorate_rec!!')
def decorate_rec((A,a_ind),(B,b_ind),drop_rows = False,verbose = False):
    r'''Decorate the rows in A with information in B --> if names overlap,
    keep the ones in A
    b_ind and a_ind can be either a single key, or a list of keys;
    if more than one element in B matches that in A, include both options!!'''
    #A = A.copy() # because I play with it later
    dropped_rows = None
    # first find the list of indeces that give us the data we want
    #{{{ process the arguments
    if (type(b_ind) is str) and (type(a_ind) is str):
        b_ind = [b_ind]
        a_ind = [a_ind]
    if ((type(b_ind) is list) and (type(a_ind) is list)) and (len(b_ind) == len(a_ind)):
        pass
    else:
        raise ValueError('If you call a list for b_ind and/or a_ind, they must match in length!!!')
    if any([x not in B.dtype.names for x in b_ind]):
        raise ValueError(repr(b_ind)+' not in '+repr(B.dtype.names)+'!!!')
    if any([x not in A.dtype.names for x in a_ind]):
        raise ValueError(repr(a_ind)+' not in '+repr(B.dtype.names)+'!!!')
    #}}}
    B_reduced = B[b_ind] # a version of B reduced to only include the keys
    B_reduced = reorder_rec(B_reduced,b_ind)# again, because it doesn't do this just based on the indexing
    A_reduced = A[a_ind] # same for A
    A_reduced = reorder_rec(A_reduced,a_ind)# again, because it doesn't do this just based on the indexing
    # now, I need to generate a mapping from the b_ind to a_ind
    field_mapping = dict(zip(b_ind,a_ind))
    # now I change the names so they match and I can compare them
    B_reduced.dtype.names = tuple([field_mapping[x] for x in B_reduced.dtype.names])
    #{{{ now find the list of indeces for B that match each value of A
    old_B_reduced_names,old_B_reduced_types = tuple(zip(*tuple(B_reduced.dtype.descr)))
    B_reduced.dtype = dtype(zip(A_reduced.dtype.names,old_B_reduced_types))
    if A_reduced.dtype != B_reduced.dtype:
        B_reduced.dtype = dtype(zip(old_B_reduced_names,old_B_reduced_types))
        raise CustomError('The datatype of A_reduced=',A_reduced.dtype,'and B_reduced=',B_reduced.dtype,'are not the same, which is going to create problems!')
    try:
        list_of_matching = [nonzero(B_reduced == j)[0] for j in A_reduced]
    except:
        raise CustomError('When trying to decorate, A_reduced=',A_reduced,'with B_reduced=',B_reduced,'one or more of the following is an empty tuple, which is wrong!:',[nonzero(B_reduced == j) for j in A_reduced])
    if verbose: print "(decorate\\_rec):: original list of matching",list_of_matching
    length_of_matching = array([len(j) for j in list_of_matching])
    if verbose: print "(decorate\\_rec):: length of matching is",length_of_matching
    if any(length_of_matching == 0):
        if drop_rows:
            if drop_rows == 'return':
                dropped_rows = A[length_of_matching == 0].copy()
            else:
                dropped_rows = A_reduced[length_of_matching == 0]
                print r'{\color{red}Warning! decorate\_rec dropped fields in the first argument',lsafen(repr(zip(A_reduced.dtype.names * len(dropped_rows),dropped_rows.tolist()))),r'}'
            #{{{ now, remove all trace of the dropped fields
            A = A[length_of_matching != 0]
            list_of_matching = [j for j in list_of_matching if len(j)>0]
            length_of_matching = [len(j) for j in list_of_matching]
            #}}}
        else:
            raise CustomError('There is no data in the second argument that has',b_ind,'fields to match the',a_ind,'fields of the first argument for the following records:',A_reduced[length_of_matching == 0],'if this is correct, you can set the drop_rows = True keyword argument to drop these fields')
    # now, do a neat trick of stackoverflow to collapse a nested list
    # this gives just the indeces in B that match the values of A
    list_of_matching = [j for i in list_of_matching for j in i]
    #}}}
    if verbose: print "(decorate\\_rec):: list of matching is",list_of_matching
    # now grab the data for these rows
    add_data = B[list_of_matching]
    #{{{ finally, smash the two sets of data together
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
    if verbose: print "(decorate\\_rec):: new dtypes:",repr(new_dtypes)
    try:
        retval = newcol_rec(retval,new_dtypes)
    except:
        raise CustomError("Problem trying to add new columns with the dtypes",new_dtypes)
    #}}}
    if verbose: print "(decorate\\_rec):: add data:",repr(add_data)
    for name in dtype(new_dtypes).names:
        if verbose: print "(decorate\\_rec):: trying to add data for",name,':',add_data[name][:]
        retval[name][:] = add_data[name][:]
    #}}}
    if drop_rows == 'return':
        return retval,dropped_rows
    else:
        return retval
def newcol_rec(A,new_dtypes):
    r'''add new, empty (i.e. random numbers) fields to A, as given by new_dtypes
    --> note that there are deeply nested numpy functions to do this, but the options are confusing, and I think the way these work is efficient'''
    if type(new_dtypes) is dtype:
        new_dtypes = new_dtypes.descr
    elif type(new_dtypes) is tuple:
        new_dtypes = [new_dtypes]
    elif type(new_dtypes) is list:
        if type(new_dtypes[0]) is not tuple:
            new_dtypes = [tuple(new_dtypes)]
    retval = empty(A.shape,dtype = A.dtype.descr + new_dtypes)
    for name in A.dtype.names:
        retval[name][:] = A[name][:]
    return retval
def applyto_rec(myfunc,myarray,mylist,verbose = False):
    r'apply myfunc to myarray with the intention of collapsing it to a smaller number of values'
    if type(mylist) is not list and type(mylist) is str:
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
        if verbose: print lsafen('(applyto rec): for row %d, I select these:'%j)
        myarray_subset = myarray[mask]
        if verbose: print lsafen('(applyto rec): ',repr(myarray_subset))
        other_fields = set(mylist)^set(thisitem.dtype.names)
        if verbose: print lsafen('(applyto rec): other fields are:',other_fields)
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = myfunc(myarray_subset[thisfield])
            except:
                raise CustomError("error in applyto_rec:  You usually get this when one of the fields that you have NOT passed in the second argument is a string.  The fields and types are:",repr(myarray_subset.dtype.descr))
        if verbose: print lsafen("(applyto rec): for row %d, I get this as a result:"%j,newrow)
        combined.append(newrow) # add this row to the list
        myarray = myarray[~mask] # mask out everything I have used from the original matrix
        if verbose: print lsafen("(applyto rec): the array is now",repr(myarray))
        j += 1
    #}}}
    combined = concatenate(combined)
    if verbose: print lsafen("(applyto rec): final result",repr(combined),"has length",len(combined))
    return combined
def meanstd_rec(myarray,mylist,verbose = False):
    r'this is something like applyto_rec, except that it applies the mean and creates new rows for the "error," where it puts the standard deviation'
    if type(mylist) is not list and type(mylist) is str:
        mylist = [mylist]
    combined = []
    other_fields = set(mylist)^set(myarray.dtype.names)
    if verbose: print '(meanstd_rec): other fields are',lsafen(other_fields)
    newrow_dtype = [[j,('%s_ERROR'%j[0],)+j[1:]] if j[0] in other_fields else [j] for j in myarray.dtype.descr]
    newrow_dtype = [k for j in newrow_dtype for k in j]
    if verbose: print lsafen('(applyto rec): other fields are:',other_fields)
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
        if verbose: print lsafen('(applyto rec): for row %d, I select these:'%j)
        myarray_subset = myarray[mask]
        if verbose: print lsafen('(applyto rec): ',repr(myarray_subset))
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = mean(myarray_subset[thisfield])
                newrow[thisfield+"_ERROR"] = std(myarray_subset[thisfield])
            except:
                raise CustomError("error in applyto_rec:  You usually get this when one of the fields that you have NOT passed in the second argument is a string.  The fields and types are:",repr(myarray_subset.dtype.descr))
        if verbose: print lsafen("(applyto rec): for row %d, I get this as a result:"%j,newrow)
        combined.append(newrow) # add this row to the list
        myarray = myarray[~mask] # mask out everything I have used from the original matrix
        if verbose: print lsafen("(applyto rec): the array is now",repr(myarray))
        j += 1
    #}}}
    combined = concatenate(combined)
    if verbose: print lsafen("(applyto rec): final result",repr(combined),"has length",len(combined))
    return combined
def make_rec(*args,**kwargs):
    r'input,names or a single argument, which is a dictionary\nstrlen = 100 gives length of the strings (which need to be specified in record arrays)\nyou can also specify (especially useful with the dictionary format) the list order = [str1,str2,...] which orders the output records with the field containing str1 first, then the field containing str2, then any remaining fields'
    strlen = 100
    if 'strlen' in kwargs.keys():
        strlen = kwargs.pop('strlen')
    if 'order' in kwargs.keys():
        order = kwargs.pop('order')
    else:
        order = None
    if len(kwargs)>0:
        raise CustomError("You have kwargs I don't understand!:",kwargs)
    if len(args) == 1:
        names = args[0].keys()
        input = args[0].values()
    elif len(args) == 2:
        input = args[0]
        names = args[1]
    else:
        raise CustomError("I don't understand the arguments you passed to make_rec!!!\nshould be (list of values, list of field names), or a dictionary")
    #{{{ apply the order kwarg
    if order is not None:
        newindeces = []
        for orderitem in order:
            newindeces += [j for j,k in enumerate(names) if (k.find(orderitem)>-1 and j not in newindeces)]
        print 'debug, before adding rest',newindeces
        newindeces += [j for j,k in enumerate(names) if j not in newindeces]
        names = [names[j] for j in newindeces]
        input = [input[j] for j in newindeces]
    #}}}
    if not (type(input) is list and type(names) is list):
        raise CustomError('you must enter a list for both')
    types = map(type,input)
    shapes = map(shape,input)
    for j,k in enumerate(input):
        if type(k) is str:
            types[j] = '|S%d'%strlen
        if type(k) is ndarray:
            types[j] = k.dtype
    try:
        mydtype = dtype(zip(names,types,shapes))
    except:
        raise CustomError('problem trying to make names',names,' types',types,'shapes',shapes)
    try:
        return array([tuple(input)],dtype = mydtype)
    except:
        raise CustomError('problem trying to assign data of type',map(type,input),'\nvalues',input,'\nonto',mydtype,'\ndtype made from tuple:',zip(names,types,shapes))#,'"types" was',types)
#{{{ convert back and forth between lists, etc, and ndarray
def make_ndarray(array_to_conv,name_forprint = 'unknown',verbose = False): 
    if type(array_to_conv) in [int,int32,double,float,complex,complex128,float,bool,bool_]: # if it's a scalar
        pass
    elif type(array_to_conv) is str:
        pass
    elif type(array_to_conv) in [list,ndarray] and len(array_to_conv) > 0:
        array_to_conv = rec.fromarrays([array_to_conv],names = 'LISTELEMENTS') #list(rec.fromarrays([b])['f0']) to convert back
    elif type(array_to_conv) in [list,ndarray] and len(array_to_conv) is 0:
        array_to_conv = None
    elif array_to_conv is  None:
        pass
    else:
        raise CustomError('type of value (',type(array_to_conv),') for attribute name',name_forprint,'passed to make_ndarray is not currently supported')
    return array_to_conv
def unmake_ndarray(array_to_conv,name_forprint = 'unknown',verbose = False): 
    r'Convert this item to an ndarray'
    if (type(array_to_conv) is recarray) or (type(array_to_conv) is ndarray and array_to_conv.dtype.names is not None and len(array_to_conv.dtype.names)>0):
        #{{{ if it's a record/structured array, it should be either a list or dictionary
        if 'LISTELEMENTS' in array_to_conv.dtype.names:
            if array_to_conv.dtype.names == tuple(['LISTELEMENTS']):
                retval = list(array_to_conv['LISTELEMENTS'])
            else:
                raise CustomError('Attribute',name_forprint,'is a recordarray with a LISTELEMENTS field, but it also has other dimensions:',array_to_conv.dtype.names,'not',tuple(['LISTELEMENTS']))
        elif len(array_to_conv)==1:
            thisval = dict(zip(a.dtype.names,a.tolist()[0]))
        else: raise CustomError('You passed a structured array, but it has more than one dimension, which is not yet supported\nLater, this should be supported by returning a dictionary of arrays')
        #}}}
    elif type(array_to_conv) is ndarray and len(array_to_conv)==1:
        #{{{ if it's a length 1 ndarray, then return the element
        retval = array_to_conv.tolist()
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",type(array_to_conv),"is a numpy array of length one"
        #}}}
    elif type(array_to_conv) in [string_,int32,float64,bool_]:
        #{{{ map numpy strings onto normal strings
        retval = array_to_conv.tolist()
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",type(array_to_conv),"is a numpy scalar"
        #}}}
    elif type(array_to_conv) is list:
        #{{{ deal with lists
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"is a list"
        typeofall = map(type,array_to_conv)
        if all(map(lambda x: x is string_,typeofall)):
            if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",typeofall,"are all numpy strings"
            retval = map(str,array_to_conv)
        else:
            if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",typeofall,"are not all numpy string"
            retval = array_to_conv
        #}}}
    else:
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",type(array_to_conv),"is not a numpy string or record array"
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
   if size(x) is 1 and x is None: return True
   if size(x) is 0: return True
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
    if 'spaces' in kwargs.keys():
        spaces = kwargs.pop('spaces')
    if 'wrap' in kwargs.keys():
        wrap = kwargs.pop('wrap')
    else:
        wrap = None
    #}}}
    if type(string) is not str:
        string = repr(string)
    if wrap is True:
        wrap = 60
    if wrap is not None:
        string = '\n'.join(textwrap.wrap(string,wrap))
    string = string.replace('\\','\\textbackslash ')
    if spaces:
        string = string.replace(' ','\\ ')
    string = string.replace('\n\t\t','\n\n\\quad\\quad ')
    string = string.replace('\n\t','\n\n\\quad ')
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
    string = string.replace('|',r'$|$')
    return string
#{{{ errors
class CustomError(Exception):
    def __init__(self, *value, **kwargs):
        if 'figurelist' in kwargs.keys():
            lplotfigures(kwargs.pop('figurelist'),'error_plots.pdf')
        if len(value)>1:
            retval = map(str,value)
        else:
            retval = str(value)
        retval = map(str,value)
        retval = ' '.join(retval)
        retval = '\n'+'\n'.join(textwrap.wrap(retval,90,replace_whitespace = False))
        if traceback.format_exc() != 'None':
            retval += '\n\nOriginal Traceback:\n'+''.join(['V']*40)+'\n'+traceback.format_exc() + '\n' + ''.join(['^']*40) + '\n'
        Exception.__init__(self,retval)
        return
def copy_maybe_none(input):
    if input == None:
        return None
    else:
        if type(input) is list:
            return map(copy,input)
        else:
            return input.copy()
def maprep(*mylist):
    mylist = list(mylist)
    for j in range(0,len(mylist)):
        if type(mylist[j]) is not str:
            mylist[j] = mylist[j].__repr__()
    return ' '.join(mylist)
#}}}
#{{{ HDF5 functions
#{{{ helper function for HDF5 search
def gensearch(labelname,format = '%0.3f',value = None,precision = None):
    'obsolete -- use h5gensearch'
    if value == None:
        raise CustomError('You must pass a value to gensearch')
    if precision == None:
        precision = value*0.01 # the precision is 1% of the value, if we don't give an argument
    searchstring_high = '(%s < %s + (%s))'%tuple([labelname]+[format]*2)
    #print "\n\nDEBUG check format:\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_high = searchstring_high%(value,precision)
    #print "\n\nDEBUG after substitution with",value,precision,":\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_low = '(%s > %s - (%s))'%tuple([labelname]+[format]*2)
    searchstring_low = searchstring_low%(value,precision)
    return searchstring_low + ' & ' + searchstring_high
def h5searchstring(fieldname,value,format = '%g',precision = 0.01):
    'search AROUND a certain value (overcomes some type conversion issues) optional arguments are the format specifier and the fractional precision'
    precision *= value
    searchstring_high = '(%s < %s + (%s))'%tuple([fieldname]+[format]*2)
    #print "\n\nDEBUG check format:\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_high = searchstring_high%(value,precision)
    #print "\n\nDEBUG after substitution with",value,precision,":\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_low = '(%s > %s - (%s))'%tuple([fieldname]+[format]*2)
    searchstring_low = searchstring_low%(value,precision)
    return '(' + searchstring_low + ' & ' + searchstring_high + ')'
#}}}
def h5loaddict(thisnode,verbose = False):
    #{{{ load all attributes of the node
    retval = dict([(x,thisnode._v_attrs.__getattribute__(x))
        for x in thisnode._v_attrs._f_list('user')])
    #}}}
    for k,v in retval.iteritems():#{{{ search for record arrays that represent normal lists
        retval[k]  = unmake_ndarray(v,name_forprint = k,verbose = verbose)
    if type(thisnode) is tables.table.Table:#{{{ load any table data
        if verbose: print "It's a table\n\n"
        if 'data' in retval.keys():
            raise CustomError('There\'s an attribute called data --> this should not happen!')
        retval.update({'data':thisnode.read()})
    elif type(thisnode) is tables.group.Group:
        #{{{ load any sub-nodes as dictionaries
        mychildren = thisnode._v_children
        for thischild in mychildren.keys():
            if thischild in retval.keys():
                raise CustomError('There\'s an attribute called ',thischild,' and also a sub-node called the',thischild,'--> this should not happen!')
            retval.update({thischild:h5loaddict(mychildren[thischild])})
        #}}}
    else:
        raise CustomError("I don't know what to do with this node:",thisnode)
    #}}}
    return retval
def h5child(thisnode,childname,clear = False,create = None,verbose = False):
    r'''grab the child, optionally clearing it and/or (by default) creating it'''
    #{{{ I can't create and clear at the same time
    if create and clear:
        raise CustomError("You can't call clear and create at the same time!\nJust call h5child twice, once with clear, once with create")
    if create is None:
        if clear == True:
            create = False
        else:
            create = True
    #}}}
    h5file = thisnode._v_file
    try:
        childnode = h5file.getNode(thisnode,childname)
        if verbose:
            print lsafe('found',childname)
        if clear:
            childnode._f_remove(recursive = True)
            childnode = None
    except tables.NoSuchNodeError:
        if create is False and not clear:
            raise CustomError('Trying to grab a node that does not exist with create = False')
        elif clear:
            childnode = None
        else:
            childnode = h5file.createGroup(thisnode,childname)
            if verbose:
                print lsafe('created',childname)
    return childnode
def h5remrows(bottomnode,tablename,searchstring):
    try:
        thistable = bottomnode.__getattr__(tablename)
        counter = 0
        try:
            data = thistable.readWhere(searchstring).copy()
        except:
            raise CustomError('Problem trying to remove rows using search string',searchstring,'in',thistable)
        for row in thistable.where(searchstring):
            if len(thistable) == 1:
                thistable.remove()
                counter += 1
            else:
                thistable.removeRows(row.nrow - counter,None) # counter accounts for rows I have already removed.
                counter += 1
        return counter,data
    except tables.NoSuchNodeError:
        return False,None
def h5addrow(bottomnode,tablename,*args,**kwargs):
    'add a row to a table, creating it if necessary, but don\'t add if the data matches the search condition'
    #{{{ process kwargs
    force = False
    if 'force' in kwargs.keys():
        force = kwargs.pop('force')
    match_row = None
    if 'match_row' in kwargs.keys():
        match_row = kwargs.pop('match_row')
    verbose = False
    if 'verbose' in kwargs.keys():
        verbose = kwargs.pop('verbose')
    only_last = True
    if 'only_last' in kwargs.keys():
        only_last = kwargs.pop('only_last')
    if len(kwargs) != 0:
        raise ValueError('kwargs'+repr(kwargs)+'not understood!!')
    #}}}
    try: # see if the table exists
        mytable = h5table(bottomnode,tablename,None)
        #{{{ auto-increment "index"
        newindex = mytable.read()['index'].max() + 1L
        #}}}
        # here is where I would search for the existing data
        if match_row is not None:
            try:
                matches = mytable.readWhere(match_row)
            except NameError:
                raise CustomError('The columns available are',mytable.colnames)
            if len(matches) > 0:
                if only_last:
                    if verbose: print r'\o{',lsafen(len(matches),"rows match your search criterion, returning the last row"),'}'
                    return mytable,matches['index'][-1]
                else:
                    return mytable,matches['index'][:]
            else:
                if add:
                    if verbose: print r'\o{',lsafen("Creating a new row for your data"),'}'
        tableexists = True
    except CustomError: # if table doesn't exist, create it
        newindex = 1L
        tableexists = False
    if len(args) == 1 and (type(args[0]) is dict):
        listofnames,listofdata = map(list,zip(*tuple(args[0].items())))
    elif len(args) == 2 and type(args[0]) is list and type(args[1]) is list:
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
        except ValueError:
            print lsafen("I'm about to flag an error, but it looks like there was an issue appending",myrowdata)
            tabledforerr = mytable.read()
            #raise CustomError('Value of table',tabledforerr.dtype.,'while value of row',myrowdata.dtype)
            raise CustomError('Value of table -- compare names and values table data vs. the row you are trying to add\n','\n'.join(map(repr,zip(mytable.read().dtype.fields.keys(),
                        mytable.read().dtype.fields.values(),
                        myrowdata.dtype.fields.keys(),
                        myrowdata.dtype.fields.values()))))
            #raise CustomError('Value of table',mytable.read().dtype,'while value of row',myrowdata.dtype,'data in row:',myrowdata,
            #        'the one that doesn\'t match is',[mytable.read().dtype.descr[x]
            #            for x in list(myrowdata.dtype)
            #            if mytable.read().dtype.descr[x]==myrowdata.dtype.descr[x]])
        mytable.flush()
    else:
        recorddata = myrowdata
        try:
            mytable = h5table(bottomnode,
                    tablename,
                    recorddata)
        except:
            raise CustomError('Error trying to write record array:',repr(recorddata),'from listofdata',listofdata,'and names',listofnames)
        mytable.flush()
    return mytable,newindex
def h5table(bottomnode,tablename,tabledata):
    'create the table, or if tabledata is None, just check if it exists'
    #{{{ save but don't overwrite the table
    h5file = bottomnode._v_file
    if tablename not in bottomnode._v_children.keys():
        if tabledata is not None:
            datatable = h5file.createTable(bottomnode,tablename,tabledata) # actually write the data to the table
        else:
            raise CustomError('You passed no data, so I can\'t create table',tablename,'but it doesn\'t exist in',bottomnode,'which has children',bottomnode._v_children.keys())
    else:
        if tabledata is not None:
            raise CustomError('You\'re passing data to create the table, but the table already exists!')
        else:
            pass
    return bottomnode._v_children[tablename]
    #}}}
def h5nodebypath(h5path,verbose = False,force = False,only_lowest = False,check_only = False):
    r'''return the node based on an absolute path, including the filename'''
    if verbose: print lsafen("DEBUG: called h5nodebypath on",h5path)
    h5path = h5path.split('/')
    #{{{ open the file / check if it exists
    if verbose: print lsafen('h5path=',h5path)
    try:
        if h5path[0] in listdir('.'):
            if verbose: print 'DEBUG: file exists\n\n'
        else:
            if check_only: raise CustomError("You're checking for a node in a file that does not exist")
            if verbose: print 'DEBUG: file does not exist\n\n'
        mode = 'a'
        #if check_only: mode = 'r'
        h5file = tables.openFile(h5path[0],mode = mode,title = 'test file')
    except IOError:
        raise CustomError('I think the HDF5 file has not been created yet, and there is a bug pytables that makes it freak out, but you can just run again.')
    #}}}
    currentnode = h5file.getNode('/') # open the root node
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
                        verbose = verbose,
                        create = create,
                        clear = clear)
                if verbose: print lsafen("searching for node path: descended to node",currentnode)
            except:
                if verbose: print lsafen("searching for node path: got caught searching for node",h5path[pathlevel])
                h5file.close()
                #print lsafen("DEBUG: Yes, I closed the file")
                raise CustomError('Problem trying to load node ',h5path)
            #}}}
    return h5file,currentnode
def h5attachattributes(node,listofattributes,myvalues):
    #print "DEBUG 5: node passed to h5attachattributes",node
    if node is None:
        raise CustomError('Problem!, node passed to h5attachattributes: ',node,'is None!')
    h5file = node._v_file
    if isinstance(myvalues,nddata):
        attributevalues = map(lambda x: myvalues.__getattribute__(x),listofattributes)
    elif type(myvalues) is list:
        attributevalues = myvalues
    else:
        raise CustomError("I don't understand the type of myvalues, which much be a list or a nddata object, from which the attribute values are retrieved")
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
            print "DEBUG: trying to create dictnode",dictnode
            h5attachattributes(dictnode,
                    thisval.keys(),
                    thisval.values())
            thisval = None
            listout.remove(thisattr)
        else:
            thisval = make_ndarray(thisval,name_forprint = thisattr)
        if thisval is not None:
            node._v_attrs.__setattr__(thisattr,thisval)
            listout.remove(thisattr)
    listofattributes[:] = listout # pointer
def h5inlist(columnname,mylist):
    'returns rows where the column named columnname is in the value of mylist'
    if type(mylist) is not list:
        raise TypeError("the second argument to h5inlist must be a list!!!")
    if any([type(x) in [double,float64] for x in mylist]):
        if all([type(x) in [double,float64,int,int32,int64] for x in mylist]):
            return '('+'|'.join(map(lambda x: "(%s == %g)"%(columnname,x),mylist))+')'
    elif all([type(x) in [int,int32] for x in mylist]):
        return '('+'|'.join(map(lambda x: "(%s == %g)"%(columnname,x),mylist))+')'
    elif all([type(x) is str for x in mylist]):
        return '('+'|'.join(map(lambda x: "(%s == '%s')"%(columnname,x),mylist))+')'
    else:
        raise TypeError("I can't figure out what to do with this list --> I know what to do with a list of numbers or a list of strings, but not this.")
def h5join(firsttuple,secondtuple,
    additional_search = '',
    select_fields = None,
    pop_fields = None,
    verbose = False):
    #{{{ process the first argument as the hdf5 table and indeces, and process the second one as the structured array to join onto
    if not ((type(firsttuple) is tuple) and (type(secondtuple) is tuple)):
        raise ValueError('both the first and second arguments must be tuples!')
    if not ((len(firsttuple) == 2) and (len(secondtuple) == 2)):
        raise ValueError('The length of the first and second arguments must be two!')
    tablenode = firsttuple[0]
    tableindeces = firsttuple[1]
    if verbose: print 'h5join tableindeces looks like this:',tableindeces
    if type(tableindeces) is not list:
        tableindeces = [tableindeces]
    if verbose: print 'h5join tableindeces looks like this:',tableindeces
    mystructarray = secondtuple[0].copy()
    mystructarrayindeces = secondtuple[1]
    if type(mystructarrayindeces) is not list:
        mystructarrayindeces = [mystructarrayindeces]
    #}}}
    #{{{ generate a search string to match potentially more than one key
    search_string = []
    if len(tableindeces) != len(mystructarrayindeces):
        raise ValueError('You must pass either a string or a list for the second element of each tuple!\nIf you pass a list, they must be of the same length, since the field names need to line up!')
    # this can't use h5inlist, because the and needs to be on the inside
    #{{{ this loop creates a list of lists, where the inner lists are a set of conditions that need to be satisfied
    # this is actually not causing  any trouble right now, but needs to be fixed, because of the way that it's doing the type conversion
    for thistableindex,thisstructarrayindex in zip(tableindeces,mystructarrayindeces):
        if thisstructarrayindex not in mystructarray.dtype.names:
            raise ValueError(repr(thisstructarrayindex)+" is not in "+repr(mystructarray.dtype.names))
        if type(mystructarray[thisstructarrayindex][0]) in [str,str_]:
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
    if verbose: print '\n\nh5join generated the search string:',lsafen(search_string)
    retval = tablenode.readWhere(search_string)
    #{{{ then join the data together
    # here I'm debugging the join function, again, and again, and agin
    try:
        retval = decorate_rec((retval,tableindeces),(mystructarray,mystructarrayindeces)) # this must be the problem, since the above looks fine
    except:
        raise CustomError('Some problems trying to decorate the table',retval,'of dtype',retval.dtype,'with the structured array',mystructarray,'of dtype',mystructarray.dtype)
    if pop_fields is not None:
        if select_fields is not None:
            raise ValueError("It doesn't make sense to specify pop_fields and select_fields at the same time!!")
        select_fields = list(set(retval.dtype.names) ^ set(pop_fields))
    if select_fields is not None:
        if verbose: print '\n\nh5join original indeces',lsafen(retval.dtype.names)
        try:
            retval = retval[select_fields]
        except ValueError:
            raise CustomError('One of the fields',select_fields,'is not in',retval.dtype.names)
    #}}}
    return retval
#}}}
#{{{ indeces to slice
#}}}
#{{{ add slashes for dir's
def dirformat(file):
        #{{{ format strings
        if file[-1]!='/':
            file += '/'
        #}}}
        return file
#}}}
#{{{ old grid and tick
def gridandtick(ax,rotation=(0,0),precision=(2,2),labelstring=('',''),gridcolor=r_[0,0,0],formatonly = False,fixed_y_locator = None,logarithmic = False,use_grid = True):
    if not formatonly:
        #{{{x ticks
        # determine the size
        width = abs(diff(ax.get_xlim()))
        if width==0:
            raise CustomError('x axis width is zero')
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
        #{{{ y ticks
        width = abs(diff(ax.get_ylim()))
        if width==0:
            raise CustomError('y axis width is zero')
        widthexp = floor(log(width)/log(10.))-1
        scalefactor = 10**widthexp
        width /= scalefactor
        if fixed_y_locator == None:
            if logarithmic:
                majorLocator = LogLocator(10)
            else:
                majorLocator   = MultipleLocator(5*scalefactor)
        else:
            majorLocator   = MultipleLocator(fixed_y_locator[4::5])
        #majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[1]+'f'+labelstring[1])# labelstring can be used, for instance, for pi
        #ax.yaxis.set_major_formatter(majorFormatter)
        if fixed_y_locator == None:
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
    grid(use_grid,which='major',color=gridcolor,alpha=0.1,linestyle='-')
    grid(use_grid,which='minor',color=gridcolor,alpha=0.05,linestyle='-')
    labels = ax.get_xticklabels()
    setp(labels,rotation=rotation[0],fontsize=10)
    labels = ax.get_yticklabels()
    setp(labels,rotation=rotation[1],fontsize=10)
    fig = gcf()
    fig.autofmt_xdate()
def gridon(gridcolor=r_[0,0,0]):
    grid(True,which='major',color=gridcolor,alpha=0.1,linestyle='-')
    grid(True,which='minor',color=gridcolor,alpha=0.05,linestyle='-')
#}}}
#{{{ a better version?
def othergridandtick(ax,rotation=(0,0),precision=(2,2),labelstring=('',''),gridcolor=r_[0,0,0]):
    #{{{x ticks
    # determine the size
    ax.xaxis.set_major_locator(MaxNLocator(10)) # could use multiplelocator if it keeps try to do multiples of 2
    ax.xaxis.set_minor_locator(MaxNLocator(50))
    #}}}
    #{{{ y ticks
    ax.yaxis.set_major_locator(MaxNLocator(10))
    ax.yaxis.set_minor_locator(MaxNLocator(50))
    #}}}
    grid(True,which='major',color=gridcolor,alpha=0.2,linestyle='-')
    grid(True,which='minor',color=gridcolor,alpha=0.1,linestyle='-')
    labels = ax.get_xticklabels()
    setp(labels,rotation=rotation[0],fontsize=10)
    labels = ax.get_yticklabels()
    setp(labels,rotation=rotation[1],fontsize=10)
#}}}
#{{{ plot wrapper
global OLDplot
OLDplot = plot
global myplotfunc
myplotfunc = OLDplot
def whereblocks(a): # returns contiguous chunks where the condition is true
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
    if 'ax' in kwargs.keys():
        ax_list = [kwargs.pop('ax')]
    else:
        ax_list = [gca()]
    if 'ax2' in kwargs.keys():
        ax_list.append(kwargs.pop('ax2'))
    for ax in ax_list:
        if len(args)==0:
            lg = ax.legend(loc='best')
        elif len(args)==1:
            lg = ax.legend(args[0],'best')
        else:
            lg = ax.legend(args[0],args[1],'best')
        if lg is None:
            print "Warning! you called autolegend, but you don't seem to have anything labeled!!"
        else:
            lg.get_frame().set_alpha(0.45)
    return lg
def autopad_figure(pad = 0.2,centered = False):
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
                   if type(label) is Line2D:
                       pass # just rely on the pad
                       #if any(map(lambda x: x == label.get_transform(),[ax.transData,ax.transAxes,fig.transFigure,None])):
                       #    print 'found it'
                       #else:
                       #    print 'didn not find it'
                       #bbox = label.get_window_extent(fig.canvas).inverse_transformed(ax.transData).inverse_transformed(fig.transFigure)
                   else:
                       try:
                           bbox = label.get_window_extent()
                       except:
                           raise CustomError('type of label = ',type(label))
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
       try:
           if len(spkwargs) > 0:
               if centered and 'left' in spkwargs.keys() and 'right' in spkwargs.keys():
                   big = max(r_[spkwargs['left'],1-spkwargs['right']])
                   spkwargs.update({'left':big,'right':1-big})
               fig.subplots_adjust(**spkwargs) # pad a little
               #print "adjusted to",spkwargs
               fig.canvas.draw()
       except:
           raise CustomError('spwargs = ',spkwargs)
       return False
    fig.canvas.mpl_connect('draw_event', on_draw)
    fig.subplots_adjust(left = 0, right = 1, top = 1, bottom =0)
    fig.canvas.draw()
    #}}}
def expand_x():
    # this is matplotlib code to expand the x axis
    ax = gca()
    xlims = array(ax.get_xlim())
    width = abs(diff(xlims))
    xlims[0] -= width/10
    xlims[1] += width/10
    ax.set_xlim(xlims)
def expand_y():
    # this is matplotlib code to expand the x axis
    ax = gca()
    ylims = array(ax.get_ylim())
    width = abs(diff(ylims))
    ylims[0] -= width/10
    ylims[1] += width/10
    ax.set_ylim(ylims)
def plot_label_points(x,y,labels,**kwargs_passed):
    kwargs = {'alpha':0.5,'color':'g','ha':'left','va':'center','rotation':0}
    kwargs.update(kwargs_passed)
    for j in range(0,len(labels)):
        text(x[j],y[j],labels[j],**kwargs)
def addlabels(labelstring,x,y,labels):
    r'obsolete -- use plot_label_points'
    for j in range(0,len(labels)):
        text(x[j],y[j],labelstring%labels[j],alpha=0.5,color='g',ha='left',va='top',rotation=0)
def plot_color_counter(*args,**kwargs):
    if 'ax' in kwargs.keys():
        ax = kwargs.pop('ax')
    else:
        ax = gca()
    if len(args)>0:
        try:
            ax._get_lines.count = args[0] # set the value of the color counter
        except:
            ax._get_lines.color_cycle = args[0] # set the value of the color counter
    try: # this is different depending on the version of matlablike
        retval = ax._get_lines.count
    except:
        retval = ax._get_lines.color_cycle
    return retval
def contour_plot(xvals,yvals,zvals,color = 'k',alpha = 1.0,npts = 300,**kwargs):
    if 'inline_spacing' in kwargs.keys():
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
    except:
        raise CustomError("Is there something wrong with your levels?:",levels,"min z",zi_min,"max z",zi_max)
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
    if 'force_color' in kwargs.keys() and kwargs['force_color'] == True:
        if hasattr(data,'other_info'):
            if 'plot_color' in data.other_info.keys():
                data.other_info.pop('plot_color')
    plot(data[axis,changemask],color1+symbol,**kwargs)
    if len(kwargs) > 0 and 'label' in kwargs.keys(): kwargs.pop('label') # if I'm doing a legend, I want it on the first
    plot(data[axis,~changemask],color2+symbol,**kwargs)
    return
def nextfigure(figurelist,name):
    'obsolete -- now use class'
    verbose = False # a good way to debug
    if verbose: print lsafe('DEBUG figurelist, called with',name)
    if name in figurelist:
        fig = figure(figurelist.index(name)+1)
        if verbose: print lsafen('in',figurelist,'at figure',figurelist.index(name)+1,'switched figures')
    else:
        fig = figure(len(figurelist)+1)
        fig.add_subplot(111)
        if verbose: print lsafen('added, figure',len(figurelist)+1,'because not in figurelist',figurelist)
        figurelist.append(name)
    return figurelist
def figlistini_old(first_figure):
    verbose = False
    if verbose: print lsafe('DEBUG: initialize figlist')
    if first_figure == None:
        if verbose: print lsafen('empty')
        return []
    else:
        if verbose: print lsafen(first_figure.figurelist)
        return first_figure
class figlist():
    def __init__(self,*arg,**kwargs):
        self.verbose = False
        self.env = ''
        if 'env' in kwargs.keys():
            self.env = kwargs.pop('env')
        if 'verbose' in kwargs.keys():
            self.verbose = kwargs.pop('verbose')
        if self.verbose: print lsafe('DEBUG: initialize figlist')
        if len(arg) == 0:
            if self.verbose: print lsafen('empty')
            self.figurelist = []
        else:
            if self.verbose: print lsafen(arg[0])
            self.figurelist = arg[0]
        if len(kwargs) > 0:
            self.figurelist.append(kwargs)
        return
    def next(self,name,**kwargs):
        if self.verbose: print lsafe('DEBUG figurelist, called with',name)
        if name in self.figurelist:
            fig = figure(self.figurelist.index(name)+1,**kwargs)
            if self.verbose: print lsafen('in',self.figurelist,'at figure',self.figurelist.index(name)+1,'switched figures')
        else:
            fig = figure(len(self.figurelist)+1,**kwargs)
            fig.add_subplot(111)
            if self.verbose: print lsafen('added, figure',len(self.figurelist)+1,'because not in figurelist',self.figurelist)
            self.figurelist.append(name)
        return gca()
    def plot(self,*args,**kwargs):
        plot(*args,**kwargs)#just a placeholder for now, will later keep units + such
    def text(self,mytext):
        self.figurelist.append({'print_string':mytext})
    def show(self,*args):
        if len(args) == 1:
            if (args[0][:-4] == '.pdf') or (args[0][:-4] == '.png') or (args[0][:-4] == '.jpg'):
                print "you passed me a filename, but I'm just burning it"
        show()
def text_on_plot(x,y,thistext,**kwargs):
    ax = gca()
    newkwargs = {'transform':ax.transAxes,'size':'x-large',"horizontalalignment":'center'}
    if 'match_data' in kwargs.keys():
        if kwargs['match_data'].get_plot_color() is not None:
            newkwargs.update({'color':kwargs['match_data'].get_plot_color()})
        else:
            raise CustomError('You passed match_data to text_on_plot, but I can\'t find a color in the object')
    return text(x,y,thistext,**newkwargs)
def plot(*args,**kwargs):
    global myplotfunc
    has_labels = False
    #{{{ deal with axes
    if 'ax' in kwargs.keys():
        ax = kwargs.pop('ax')
    else:
        ax = gca()
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
    elif (len(args)==2) and (type(args[1]) is str):
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
    #{{{ parse nddata
    if isinstance(myy,nddata):
        myy = myy.copy().human_units()
        if myy.get_plot_color() is not None\
            and 'color' not in kwargs.keys():# allow override
            kwargs.update({'color':myy.get_plot_color()})
        #if (len(myy.dimlabels)>1):
        #    myylabel = myy.dimlabels[1]
        #    yunits = myy.get_units(myylabel)
        #    if yunits is not None:
        #        myylabel += ' / $' + yunits + '$'
        if (len(myy.dimlabels)>1):
            yunits = myy.get_units()
            if myy.name() is not None:
                myylabel = myy.name()
            else:
                myylabel = 'data'
            if yunits is not None:
                myylabel = myylabel + ' / ' + yunits
        if (len(myy.dimlabels)>0):
            myxlabel = myy.dimlabels[0]
            xunits = myy.get_units(myxlabel)
            if xunits is not None:
                myxlabel += ' / ' + xunits
        if (myx == None):
            try:
                myx = myy.getaxis(myy.dimlabels[0])
            except:
                myx = r_[0:myy.data.shape[0]]
        if type(myy.data_error) is ndarray and len(myy.data_error)>0: #then this should be an errorbar plot
            def thiserrbarplot(*tebargs,**tebkwargs):
                if type(tebargs[-1]) is str:
                    tebkwargs.update({'fmt':tebargs[-1]})
                    ax.errorbar(*tebargs[:-1],**tebkwargs)
                else:
                    ax.errorbar(*tebargs,**tebkwargs)
            myplotfunc = thiserrbarplot
            #{{{ pop any singleton dims
            myyerror = myy.get_error()
            myyerror = squeeze(myyerror)
            #}}}
            kwargs.update({'yerr':myyerror})
            valueforxerr = myy.get_error(myy.dimlabels[0])
            if valueforxerr != None: # if we have x errorbars too
                #print "DEBUG decided to assign to xerr:",valueforxerr
                kwargs.update({'xerr':valueforxerr})
        #{{{ deal with axis labels along y
        try:
            yaxislabels = myy.getaxis(myy.dimlabels[-1])
            # at this point, if there is no axis label, it will break and go to pass
            if len(yaxislabels) > 0:
                if type(yaxislabels[0]) is string_:
                    has_labels = True
        except:
            pass
        #}}}
        myy = squeeze(myy.data)
    #}}}
    #{{{ semilog where appropriate
    if (myx != None) and (len(myx)>1): # by doing this and making myplotfunc global, we preserve the plot style if we want to tack on one point
        b = diff(log10(myx))
        if (size(b)>3) and all(abs((b-b[0])/b[0])<1e-4) and not ('nosemilog' in kwargs.keys()):
            myplotfunc = ax.semilogx
    if ('nosemilog' in kwargs.keys()):
        #print 'this should pop nosemilog'
        kwargs.pop('nosemilog')
    if 'plottype' in kwargs.keys():
        if kwargs['plottype'] == 'semilogy':
            myplotfunc = ax.semilogy
        elif kwargs['plottype'] == 'semilogx':
            myplotfunc = ax.semilogx
        elif kwargs['plottype'] == 'loglog':
            myplotfunc = ax.loglog
        kwargs.pop('plottype')
    #}}}
    #{{{ take care of manual colors
    if myformat != None:
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
    #{{{ hsv plots when we have multiple lines
    if len(shape(myy))>1 and sum(array(shape(myy))>1):
        #{{{ hsv plots
        hold(True)
        retval = []
        for j in range(0,myy.shape[1]):
            #{{{ this is the way to assign plot arguments
            plotargs = [myx,myy[:,j],myformat]
            while None in plotargs:
                plotargs.remove(None)
            #}}}
            #{{{ here, i update the kwargs to include the specific color for this line
            newkwargs = kwargs.copy() # kwargs is a dict
            newkwargs.update({'color':cm.hsv(double(j)/double(myy.shape[1]))})
            #}}}
            #{{{ here, I update to use the labels
            if has_labels:
                newkwargs.update({'label':yaxislabels[j]})
            #}}}
            myy[isinf(myy)] = NaN # added this to prevent an overflow error
            try:
                retval += [myplotfunc(*tuple(plotargs),**newkwargs)]
                #print "\n\n\\begin{verbatim}DEBUG plot:",plotargs,'\nkwargs:\n',newkwargs,'\\end{verbatim}'
            except: 
                raise CustomError("Error trying to plot using function",myplotfunc,len(plotargs),"arguments",plotargs,"of len",map(len,plotargs),"and",len(newkwargs),"options",newkwargs,"of len",map(len,newkwargs.values()))
        #hold(False)
        #}}}
        #}}}
    else:
        plotargs = [myx,real(myy),myformat]
        while None in plotargs:
            plotargs.remove(None)
        try:
            #print 'DEBUG plotting with args',plotargs,'and kwargs',kwargs,'\n\n'
            retval = myplotfunc(*plotargs,**kwargs)
        except:
            raise CustomError('error trying to plot',myplotfunc,
                    '\nlength of the ndarray arguments:',['shape:'+repr(shape(j)) if type(j) is ndarray else j for j in plotargs],
                    '\nsizes of ndarray kwargs',dict([(j,shape(kwargs[j])) if type(kwargs[j]) is ndarray else (j,kwargs[j]) for j in kwargs.keys()]),
                    '\narguments = ',plotargs,
                    '\nkwargs =',kwargs)
    #{{{ attach labels and such
    if (myxlabel!=None):
        ax.set_xlabel(myxlabel)
    if (myylabel!=None):
        ax.set_ylabel(myylabel)
    try:
        ax.axis('tight')
    except:
        raise CustomError('error trying to set axis tight after plot',myplotfunc,'with arguments',plotargs,'and kwargs',kwargs,'\nsizes of arguments:',[shape(j) for j in plotargs],'\nsizes of ndarray kwargs:',dict([(j,shape(kwargs[j])) for j in kwargs.keys() if type(kwargs[j]) is ndarray]))
    #grid(True)
    #}}}
    return retval
#}}}
#{{{general functions
def box_muller(length):
    r'''algorithm to generate normally distributed noise'''
    s1 = rand(length)
    s2 = rand(length)
    n1 = sqrt(-2*log(s1))*cos(2*pi*s2)
    n2 = sqrt(-2*log(s1))*sin(2*pi*s2)
    return (n1 + 1j * n2)*0.5
#}}}

#{{{nddata
#{{{ shaping and allocating
class ndshape ():
    def __init__(self,*args):
        if len(args) == 2:
            self.shape = list(args[0])
            self.dimlabels = args[1]
        if len(args) == 1: #assum that it's an nddata object
            self.shape = list(args[0].data.shape)
            self.dimlabels = list(args[0].dimlabels)
        return
    def __setitem__(self,reference,setto):
        self.shape[self.dimlabels.index(reference)] = setto
        return
    def copy(self):
        try:
            return deepcopy(self)
        except:
            raise RuntimeError('Some type of error trying to run deepcopy on'+repr(self))
    def matchdims(self,arg):
        r'returns shape with [not in self, len 1] + [overlapping dims between arg + self] + [not in arg] --> this is better accomplished by using sets as I do in the matchdims below'
        for k in set(self.dimlabels) & set(arg.dimlabels):
            a = arg.shape[arg.dimlabels.index(k)]
            b = self.shape[self.dimlabels.index(k)]
            if a != b:
                raise CustomError('the',k,'dimension is not the same for self',self,'and arg',arg)
        if isinstance(arg,nddata):
            arg = ndshape(arg)
        #{{{ add extra 1-len dims
        addeddims = set(self.dimlabels) ^ set(arg.dimlabels) & set(arg.dimlabels)
        self.dimlabels = list(addeddims) + self.dimlabels
        self.shape = [1] * len(addeddims) + list(self.shape)
        #}}}
        return self
    def __add__(self,arg):
        'take list of shape,dimlabels'
        shape = arg[0]
        dimlabels = arg[1]
        self.shape = shape + self.shape
        self.dimlabels = dimlabels + self.dimlabels
        return self
    def __repr__(self): #how it responds to print
        return zip(self.shape,self.dimlabels).__repr__()
    def __getitem__(self,args):
        try:
            mydict = dict(zip(self.dimlabels,self.shape))
        except:
            raise CustomError("either dimlabels=",self.dimlabels,"or shape",self.shape,"not in the correct format")
        try:
            return mydict[args]
        except:
            raise CustomError("one or more of the dimensions named",args,"do not exist in",self.dimlabels)
    def pop(self,label):
        thisindex = self.dimlabels.index(label)
        self.dimlabels.pop(thisindex)
        self.shape.pop(thisindex)
    def alloc(self,dtype='complex128',labels = False,format = 0):
        try:
            if format == 0:
                emptyar = zeros(tuple(self.shape),dtype=dtype)
            elif format == 1:
                emptyar = ones(tuple(self.shape),dtype=dtype)
            elif format is None:
                emptyar = empty(tuple(self.shape),dtype=dtype)
            else:
                emptyar = format*ones(tuple(self.shape),dtype=dtype)
        except TypeError:
            raise CustomError('Wrong type for self.shape',map(type,self.shape))
        retval = nddata(emptyar,self.shape,self.dimlabels)
        if labels:
            retval.labels(self.dimlabels,map(lambda x: zeros(x),self.shape))
        return retval
#}}}
#{{{ format out to a certain decimal place
def dp(number,decimalplaces,scientific=False):
    if scientific:
        tenlog = floor(log(number)/log(10.))
        number /= 10**tenlog
        fstring = '%0.'+'%d'%decimalplaces+r'f\times 10^{%d}'%tenlog
    else:
        fstring = '%0.'+'%d'%decimalplaces+'f'
    return fstring%number
#}}}
#{{{ concatenate datalist along dimname
def concat(datalist,dimname,chop = False,verbose = False):
    #{{{ allocate a new datalist structure  
    newdimsize = 0
    #print 'DEBUG: type(datalist)',type(datalist)
    try:
        shapes = map(ndshape,datalist)
    except:
        if type(datalist) is not list:
            raise CustomError('You didn\'t pass a list, you passed a',type(datalist))
        raise CustomError('Problem with what you passed to concat, list of types,',map(type,datalist))
    other_info_out = datalist[0].other_info
    for j in range(0,len(datalist)):
        #{{{ make list for the shape to check, which contains the dimensions we are NOT concatting along
        if dimname in shapes[j].dimlabels:
            newdimsize += shapes[j][dimname]
            shapetocheck = list(shapes[j].shape)
            shapetocheck.pop(shapes[j].dimlabels.index(dimname))
        else:
            newdimsize += 1
            shapetocheck = list(shapes[j].shape)
        #}}}
        if j is 0:
            shapetocheckagainst = shapetocheck
        else:
            if any(~(array(shapetocheck) == array(shapetocheckagainst))):
                if chop:
                    if verbose:
                        print lsafen(repr(shapetocheck)),lsafen(repr(shapetocheckagainst))
                        raise CustomError('For item ',j,'in concat, ',shapetocheck,'!=',shapetocheckagainst,'where all the shapes of the things you\'re trying to concat are:',shapes)
                else:
                    raise CustomError('For item ',j,'in concat, ',shapetocheck,'!=',shapetocheckagainst,'where all the shapes of the things you\'re trying to concat are:',shapes)
    newdatalist = ndshape(datalist[-1])
    if dimname in newdatalist.dimlabels:
        newdatalist[dimname] = newdimsize
    else:
        newdatalist += ([newdimsize],[dimname])
    #print "DEBUG newdatalist is shaped like",newdatalist
    newdatalist = newdatalist.alloc()
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
        except:
            raise CustomError("trying to attach axes of lengths",map(len,axis_coords),"to",dimlabels)
    #}}}
    newdatalist.other_info = other_info_out
    return newdatalist
#}}}
class nddata (object):
    want_to_prospa_decim_correct = False
    def __init__(self,*args,**kwargs):
        if len(args) > 1:
            if len(args) == 2:
                if len(args[0].shape) == 1 and type(args[1]) is str:
                    self.__my_init__(args[0],[len(args[0])],[args[1]])
                    self.labels(args[1],args[0])
                else:
                    raise ValueError('You can pass two arguments only if you pass a 1d ndarray and a name for the axis') 
            else:
                self.__my_init__(args[0],args[1],args[2],**kwargs)
        else:
            self.__my_init__(args[0],[-1],['value'],**kwargs)
        return
    def __my_init__(self,data,sizes,dimlabels,axis_coords=[],ft_start_time = 0.,data_error = None, axis_coords_error = None,axis_coords_units = None, data_units = None, other_info = {}):
        self.genftpairs = False
        if not (type(data) is ndarray):
            #if (type(data) is float64) or (type(data) is complex128) or (type(data) is list):
            if isscalar(data) or (type(data) is list) or (type(data) is tuple):
                data = array(data)
            else:
                raise CustomError('data is not an array, it\'s',type(data),'!')
        if not (type(dimlabels) is list):
            raise CustomError('labels are not a list')
        try:
            self.data = reshape(data,sizes)
        except:
            raise CustomError("trying to reshape a ",data.shape,"array with list of sizes",sizes)
        self.dimlabels = dimlabels
        self.axis_coords = axis_coords
        #if len(axis_coords) > 0:
        #    testshape = data.shape
        #    if not all([len(axis_coords[j])==testshape[j] if axis_coords[j] is not None else True for j in range(0,len(axis_coords))]):
        #        raise IndexError('The length of your axis labels (axis_coords) (shape %s) and your axis data (shape %s) does not match!!!'%(repr([len(thiscoord) for thiscoord in axis_coords]),repr(data.shape)))
        self.ft_start_time = ft_start_time
        self.data_error = data_error
        self.data_units = data_units
        self.other_info = dict(other_info)
        if axis_coords_error == None:
            self.axis_coords_error = [None]*len(axis_coords)
        else:
            self.axis_coords_error = axis_coords_error
        if axis_coords_units == None:
            self.axis_coords_units = [None]*len(axis_coords)
        else:
            self.axis_coords_units = axis_coords_units 
        return
    def _contains_symbolic(self,string):
        return string[:9] == 'symbolic_' and hasattr(self,string)
    #{{{ for printing
    def __repr__(self):
        return repr(self.data)+'\n\t\t+/-'+repr(self.get_error())+'\n\tdimlabels=['+repr(self.dimlabels)+']\n\taxes='+repr(self.mkd(self.axis_coords))+'\n\t\t+/-'+repr(self.mkd(self.axis_coords_error))+'\n'
    #}}}
    #{{{ for plotting
    def gnuplot_save(self,filename):
        x = self.getaxis(self.dimlabels[0])[:5]
        y = self.getaxis(self.dimlabels[1])[:5]
        z = self.data[:5,:5]
        print "size of x",size(x),"size of y",size(y),"size of z",size(z)
        print "x",x,"y",y,"z",z
        data = empty((z.shape[0]+1,z.shape[1]+1))
        data[1:,1:] = z[:]
        data[0,0] = z.shape[1]
        data[0,1:] = y.flatten()
        data[1:,0] = x.flatten()
        print "data",data
        fp = open('auto_figures/'+filename+'.dat','w')
        fp.write(float32(data).tostring())
        fp.write('\n')
        fp.close()
        return
    #{{{ 3D mesh plot
    def meshplot(self,stride = None,alpha = 0.3,onlycolor = False,light = None,rotation = None,cmap = cm.gray,ax = None,invert = False,**kwargs):
        r'''takes both rotation and light as elevation, azimuth
        only use the light kwarg to generate a black and white shading display'''
        sortedself = self.copy()
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if invert:
            print "trying to invert meshplot"
        if light == True:
            light = [0,0]# I think this is 45 degrees up shining down from the left of the y axis
        if not onlycolor:
            if ax is None: 
                ax = sortedself._init_3d_axis(ax,rotation = rotation)
            else:
                if rotation is not None:
                    raise ValueError("you can only set the rotation once! (you tried"+repr(rotation)+")")
        if len(sortedself.dimlabels) > 2:
            raise CustomError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = sortedself.dimlabels[0]
        y_dim = sortedself.dimlabels[1]
        x_axis = sortedself.retaxis(x_dim).data
        y_axis = sortedself.retaxis(y_dim).data
        #}}}
        rstride = 1
        cstride = 1
        if stride is not None:
            if x_dim in stride.keys():
                rstride = stride[x_dim]
            if y_dim in stride.keys():
                cstride = stride[y_dim]
        X = x_axis*ones(shape(y_axis))
        Y = ones(shape(x_axis))*y_axis
        Z = real(sortedself.data)
        if invert:
            X = X[:,::-1]
            Y = Y[:,::-1]
            Z = Z[:,::-1]
        #mask = isfinite(Z).flatten()
        #X = X.flatten()[mask]
        #Y = Y.flatten()[mask]
        #Z = Z.flatten()[mask]
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
                        shade = False,
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
            ax.set_zlabel(sortedself.name())
        if onlycolor:
            return
        else:
            return ax
    def waterfall(self,alpha = 0.3,ax = None,rotation = None):
        if ax is None: 
            ax = self._init_3d_axis(ax,rotation = rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        if len(self.dimlabels) > 2:
            raise CustomError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        x_axis = self.retaxis(x_dim).data
        y_axis = self.retaxis(y_dim).data
        #}}}
        ax.set_xlabel(x_dim)
        ax.set_ylabel(y_dim)
        ax.set_zlabel(self.name())
        verts = []
        xs = x_axis.flatten()
        xs = r_[xs[0],xs,xs[-1]] # add points for the bottoms of the vertices
        ys = y_axis.flatten()
        for j in range(0,len(ys)):
            zs = self[y_dim,j].data.flatten()
            zs = r_[0,zs,0]
            verts.append(zip(xs,zs)) # one of the faces
        poly = PolyCollection(verts) # the individual facecolors would go here
        poly.set_alpha(0.3)
        fig = gcf()
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
        if ax == None:
            fig = gcf()
            ax = axes3d.Axes3D(fig)
            print "I'm trying to rotate to",rotation
            #ax.view_init(20,-120)
            ax.view_init(elev = 20 + rotation[1],azim = -120 + rotation[0])
        return ax
    def oldtimey(self,alpha = 0.5,ax = None,linewidth = None,sclinewidth = 20.,light = True,rotation = None,invert = False,**kwargs):
        sortedself = self.copy()
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if invert:
            print "trying to invert oldtimey"
        if linewidth == None:
            linewidth = sclinewidth/sortedself.data.shape[1]
            print "setting linewidth to %0.1f"%linewidth
        if ax is None: 
            ax = sortedself._init_3d_axis(ax,rotation = rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        ax = sortedself.meshplot(linewidth = 0,light = light,ax = ax,invert = invert)
        #return
        if len(sortedself.dimlabels) > 2:
            raise CustomError("I don't know how to handle something with more than two dimensions for a surface plot!")
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
        except:
            raise CustomError('Problem getting error because error is',self.get_error())
        return E.toarray()
    #}}}
    #{{{ shortcuts for axes
    def axlen(self,axis):
        return shape(self.data)[self.axn(axis)]
    def axn(self,axis):
        r'return the number for the axis with the name "axis"'
        try:
            return self.dimlabels.index(axis)
        except:
            raise CustomError('there is no axis named',axis,'all axes are named',self.dimlabels)
    #}}}
    #{{{ dictionary functions -- these convert between two formats:
    # dictionary -- stuff labeled according the dimension label.
    # list -- same information, but it's assumed they are listed in the order given by "dimlabels"
    def mkd(self,*arg,**kwargs):
        'make dictionary format'
        #{{{ process kwargs
        give_None = True
        if len(kwargs) > 0:
            if 'give_None' in kwargs.keys():
                give_None = kwargs.pop('give_None')
        if len(kwargs) > 0:
            raise CustomError("you passed mkd kwargs I didn't understand:",kwargs)
        #}}}
        if len(arg) == 1:
            if emptytest(arg[0]):
                return dict(zip(self.dimlabels,
                    [None]*len(self.dimlabels)))
            if len(arg[0]) != len(self.dimlabels):
                print r"{\color{red}WARNING! mkd error (John will fix this later):}"
                print "When making a dictionary with mkd, you must pass a list that has one element for each dimension!  dimlabels is"+repr(self.dimlabels)+"and you passed"+repr(arg)+'\n\n'
                raise ValueError("When making a dictionary with mkd, you must pass a list that has one element for each dimension!  dimlabels is"+repr(self.dimlabels)+"and you passed"+repr(arg))
            for i,v in enumerate(arg[0]):
                if type(v) == ndarray:
                    if v.shape == ():
                        arg[0][i] = None
            if give_None:
                return dict(zip(self.dimlabels,arg[0]))
            else:
                #{{{ don't return values for the things that are None
                mykeys = [self.dimlabels[j] for j in range(0,len(self.dimlabels)) if arg[0][j] is not None]
                myvals = [arg[0][j] for j in range(0,len(self.dimlabels)) if arg[0][j] is not None]
                return dict(zip(mykeys,myvals))
                #}}}
        elif len(arg) == 0:
            if not give_None:
                raise CustomError("You can't tell me not to give none and then not pass me anything!!")
            return dict(zip(self.dimlabels,
                [None]*len(self.dimlabels)))
        else:
            raise CustomError('.mkd() doesn\'t know what to do with %d arguments',len(arg))
    def fld(self,dict_in,noscalar = False):
        'flatten dictionary -- return list'
        return [dict_in[x] for x in self.dimlabels]
    #}}}
    #{{{ set + get the error + units
    #{{{ set units
    def set_units(self,*args):
        if len(args) == 2:
            unitval = args[1] # later, have some type of processing bojive
            if self.axis_coords_units == None or len(self.axis_coords_units) == 0:
                self.axis_coords_units = [None] * len(self.dimlabels)
            self.axis_coords_units[self.axn(args[0])] = unitval
        elif len(args) == 1:
            unitval = args[0] # later, have some type of processing bojive
            self.data_units = unitval
        else:
            raise CustomError(".set_units() takes data units or 'axis' and axis units")
    def human_units(self,verbose = False):
        prev_label = self.get_units()
        oom_names =   ['T' , 'G' , 'M' , 'k' , '' , 'm' , '\\mu' , 'p']
        oom_values = r_[12 , 9   , 6   , 3   , 0  , -3  , -6     , -12]
        if prev_label is not None and len(prev_label)>0:
            average_oom = log10(abs(self.data)).mean()/3. # find the average order of magnitude, rounded down to the nearest power of 3
            average_oom = 3*floor(abs(average_oom))*sign(average_oom)
            oom_index = argmin(abs(average_oom-oom_values))
            self.data /= 10**oom_values[oom_index]
            self.set_units(oom_names[oom_index]+prev_label)
        else:
            if verbose: print 'data does not have a unit label'
        for thisaxis in self.dimlabels:
            prev_label = self.get_units(thisaxis)
            if prev_label is not None and len(prev_label)>0:
                data_to_test = self.getaxis(thisaxis)
                data_to_test = data_to_test[isfinite(data_to_test)]
                #{{{ find the average order of magnitude, rounded down to the nearest power of 3
                average_oom = log10(abs(data_to_test))/3.
                average_oom = average_oom[isfinite(average_oom)].mean()
                #}}}
                if verbose: print "(human units): for axis",thisaxis,"the average oom is",average_oom*3
                average_oom = 3*floor(abs(average_oom))*sign(average_oom)
                if verbose: print "(human units): for axis",thisaxis,"I round this to",average_oom
                oom_index = argmin(abs(average_oom-oom_values))
                if verbose: print "(human units): for axis",thisaxis,"selected an oom value of",oom_values[oom_index]
                x = self.getaxis(thisaxis)
                x[:] /= 10**oom_values[oom_index]
                self.setaxis(thisaxis,x)
                self.set_units(thisaxis,oom_names[oom_index]+prev_label)
            else:
                if verbose: print thisaxis,'does not have a unit label'
        return self
    #}}}
    #{{{ get units
    def get_units(self,*args):
        if len(args) == 1:
            if self.axis_coords_units == None:
                return None
            if len(self.axis_coords_units) == 0:
                return None
            try:
                return self.axis_coords_units[self.axn(args[0])]
            except:
                raise CustomError('problem getting units for',args[0],'dimension',self.dimlabels,self.axis_coords_units)
        elif len(args) == 0:
            return self.data_units
        else:
            raise CustomError(".set_units() takes axis or nothing")
    #}}}
    #{{{ set error
    def set_error(self,*args):
        '''set the errors\neither set_error('axisname',error_for_axis) or set_error(error_for_data)'''
        if (len(args) is 1) and isscalar(args[0]):
            args = (r_[args[0]],)
        if (len(args) is 1) and (type(args[0]) is ndarray):
            self.data_error = reshape(args[0],shape(self.data))
        elif (len(args) is 1) and (type(args[0]) is list):
            self.data_error = reshape(array(args[0]),shape(self.data))
        elif (len(args) is 2) and (type(args[0]) is str) and (type(args[1]) is ndarray):
            self.axis_coords_error[self.axn(args[0])] = args[1]
        elif (len(args) is 1) and args[0] is None:
            self.data_error = None
        else:
            raise CustomError('Not a valid argument to set_error:',map(type,args))
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
        if (len(args) is 0):
            if self.data_error is None:
                return None
            else:
                return real(self.data_error)
        elif (len(args) is 1):
            thearg = args[0]
            if type(thearg) is str_:
                thearg = str(thearg) # like in the other spot, this became necessary with some upgrade, though I'm not sure that I should maybe just change the error functions to treat the numpy string in the same way
            if (type(thearg) is str):
                if len(self.axis_coords_error) == 0: self.axis_coords_error = [None] * len(self.dimlabels) # is we have an empty axis_coords_error, need to fill with None's
                try:
                    errorforthisaxis = self.axis_coords_error[self.axn(thearg)]
                except:
                    raise CustomError('Problem trying to load error',self.axn(thearg),'for axis',thearg,'out of',self.axis_coords_error)
                if errorforthisaxis is None:
                    return None
                else:
                    x = self.axis_coords_error[self.axn(thearg)]
                    if type(x) is ndarray:
                        if x.shape == ():
                            return None
                        else:
                            return real(self.axis_coords_error[self.axn(thearg)])
                    else:
                        return real(self.axis_coords_error[self.axn(thearg)])
        else:
            raise CustomError('Not a valid argument to get_error: *args=',args,'map(type,args)=',map(type,args))
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
    def set_prop(self,propname,val):
        self.other_info.update({propname:val})
        return
    def get_prop(self,propname):
        if propname not in self.other_info.keys():
            return None
        return self.other_info[propname]
    def name(self,*arg):
        r"""args:
           .name(newname) --> Name the object (for storage, etc)
           .name() --> Return the name"""
        if len(arg) == 1:
            self.set_prop('name',arg[0])
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
                raise CustomError('Color',thiscolor,'not in dictionary')
        self.other_info.update({'plot_color':thiscolor})
        return
    def get_plot_color(self):
        if 'plot_color' in self.other_info.keys():
            return self.other_info['plot_color']
        else:
            return None
    #}}}
    #}}}
    #{{{ arithmetic
    def __add__(self,arg):
        a = self.copy()
        if isscalar(arg):
            a.data += arg
            return a
        b = arg.copy()
        if (shape(arg.data)!=shape(a.data)):
            a.matchdims(b)
            b.matchdims(a)
        a.data += b.data
        Aerr = a.get_error()
        Berr = arg.get_error()
        if (Aerr is not None) or (Berr is not None):
            sigma = 0
            if Aerr != None:
                sigma += Aerr**2
            if Berr != None:
                sigma += Berr**2
            a.set_error(sqrt(sigma))
        return a
    def __sub__(self,arg):
        a = self.copy()
        if isscalar(arg):
            a.data -= arg
            return a
        b = arg.copy()
        if (shape(arg.data)!=shape(a.data)):
            a.matchdims(b)
            b.matchdims(a)
        a.data -= b.data
        Aerr = a.get_error()
        Berr = arg.get_error()
        if (Aerr is not None) or (Berr is not None):
            sigma = 0
            if Aerr != None:
                sigma += Aerr**2
            if Berr != None:
                sigma += Berr**2
            a.set_error(sqrt(sigma))
        return a
    def __mul__(self,arg):
        #{{{ do scalar multiplication
        if isscalar(arg):
            A = self.copy()
            A.data *= arg
            if A.get_error() != None:
                error = A.get_error()
                error *= abs(arg)
            return A
        #}}}
        #{{{ shape and multiply
        A,B = self.aligndata(arg)
        retval = A
        retval.data = A.data * B.data
        #}}}
        #{{{ if we have error for both the sets of data, I should propagate that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0 # we can have error on one or both, so we're going to need to add up the variances
        if Aerr != None:
            Rerr += (Aerr * B.data)**2
        if Berr != None:
            Rerr += (Berr * A.data)**2
        Rerr = sqrt(real(Rerr)) # convert back to stdev
        if Aerr == None and Berr == None:
            Rerr = None
        #}}}
        retval.set_error(Rerr)
        return retval
    def __pow__(self,arg):
        if arg == -1:
            x = self.get_error()
            result = self.copy()
            result.data = 1.0/result.data
            if x != None:
                result.set_error(abs(x.copy()/(self.data**2)))
            return result
        else:
            if self.get_error() != None:
                raise CustomError("nothing but -1 supported yet!")
            else:
                result = self.copy()
                result.data = result.data**arg
                return result
    def __div__(self,arg):
        if isscalar(arg):
            A = self.copy()
            A.data /= arg
            if A.get_error() != None:
                error = A.get_error()
                error /= abs(arg)
            return A
        A,B = self.aligndata(arg)
        retval = A
        retval.data = A.data / B.data
        #{{{ if we have error for both the sets of data, I should propagate that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0 # we can have error on one or both, so we're going to need to add up the variances
        if Aerr != None:
            Rerr += (Aerr/B.data)**2
        if Berr != None:
            Rerr += (A.data*Berr/(B.data**2))**2
        try:
            Rerr = sqrt(real(Rerr)) # convert back to stdev
        except:
            raise CustomError("Rerr gave an attribute error when you passed",Rerr)
        #print "Rerr dtype",Rerr.dtype
        if Aerr == None and Berr == None:
            Rerr = None
        #}}}
        retval.set_error(Rerr)
        return retval
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
    def real(self):
        self.data = real(self.data)
        return self
    #}}}
    #{{{ align data
    def aligndata(self,arg):
        r'''This now just returns selfout,argout
        which are aligned to each other, and which contain
        axis labels and axis errors for both'''
        augmentdims = [x for x in arg.dimlabels if x in set(self.dimlabels)^set(arg.dimlabels)] # dims in arg but now self, ordered as they were in arg
        newdims = self.dimlabels + augmentdims # this should return the new dimensions with the order of self preserved, followed by augmentdims
        selfout = self.copy() # copy self
        if len(selfout.data.shape) == 0:
            if len(selfout.dimlabels) == 1:
                selfout.data.reshape(1)
            else:
                raise ValueError("This instance is zero dimensional (It's %s)!!!"%repr(self))
        selfshape = list(selfout.data.shape)+list(ones(len(augmentdims))) # there is no need to transpose self, since its order is preserved
        new_arg_labels = [x for x in newdims if x in arg.dimlabels] # only the labels valid for B, ordered as they are in newdims
        argout = arg.copy()
        if len(argout.data.shape) == 0:
            if len(argout.dimlabels) == 1:
                argout.data = argout.data.reshape(1)
                if argout.get_error() is not None:
                    try:
                        argout.set_error(argout.get_error().reshape(1))
                    except:
                        raise ValueError("error was"+repr(argout.get_error()))
                #print "DEBUG: I reshaped argout, so that it now has shape",argout.data.shape,"dimlabels",argout.dimlabels,"and ndshape",ndshape(argout)
            else:
                raise ValueError("This instance is zero dimensional (It's %s)!!!"%repr(self))
        argshape = list(ones(len(newdims)))
        #{{{ wherever the dimension exists in arg, pull the shape from arg
        for j,k in enumerate(newdims):
            if k in argout.dimlabels:
                try:
                    argshape[j] = argout.data.shape[argout.axn(k)]
                except:
                    raise CustomError("There seems to be a problem because the shape of argout is now len:%d "%len(argout.data.shape),argout.data.shape,"while the dimlabels is len:%d "%len(argout.dimlabels),argout.dimlabels)
        #}}}
        argorder = map(argout.dimlabels.index,new_arg_labels) # for each new dimension, determine the position of the original dimension
        selfout.data = selfout.data.reshape(selfshape) # and reshape to its new shape
        selfout.dimlabels = newdims
        argout.data = argout.data.transpose(argorder).reshape(argshape) # and reshape the data
        argout.dimlabels = newdims
        if selfout.get_error() != None:
            try:
                temp = selfout.get_error().copy().reshape(selfshape)
            except ValueError,Argument:
                raise ValueError("The instance (selfout) has a shape of "+repr(selfout.data.shape)+" but its error has a shape of"+repr(selfout.get_error().shape)+"!!!\n\n(original argument:\n"+repr(Argument)+"\n)")
            selfout.set_error(temp)
        if argout.get_error() != None:
            try:
                temp = argout.get_error().copy().transpose(argorder).reshape(argshape)
            except ValueError,Argument:
                raise ValueError("The argument (argout) has a shape of "+repr(argout.data.shape)+" but its error has a shape of"+repr(argout.get_error().shape)+"(it's "+repr(argout.get_error())+")!!!\n\n(original argument:\n"+repr(Argument)+"\n)")
            argout.set_error(temp)
        if (len(selfout.axis_coords)>0) or (len(argout.axis_coords)>0):
            #{{{ transfer the errors and the axis labels
            #{{{ make dictionaries for both, and update with info from both, giving preference to self
            axesdict = selfout.mkd()
            errordict = selfout.mkd()
            #{{{ add the axes and errors for B
            if type(arg.axis_coords) is list:
                if len(arg.axis_coords) > 0:
                    axesdict.update(arg.mkd(arg.axis_coords))
            if type(arg.axis_coords_error) is list:
                if len(arg.axis_coords_error) > 0 and not all([x is None for x in arg.axis_coords_error]):
                    errordict.update(arg.mkd(arg.axis_coords_error))
            #}}}
            #{{{ add the errors for A
            if type(self.axis_coords) is list:
                if len(self.axis_coords) > 0:
                    axesdict.update(self.mkd(self.axis_coords))
            if type(self.axis_coords_error) is list:
                if len(self.axis_coords_error) > 0 and not all([x is None for x in self.axis_coords_error]):
                    errordict.update(self.mkd(self.axis_coords_error))
            #}}}
            #}}}
            selfout.axis_coords_error = selfout.fld(errordict)
            argout.axis_coords_error = selfout.fld(errordict)
            selfout.axis_coords = selfout.fld(axesdict)
            argout.axis_coords = selfout.fld(axesdict)
            #}}}
        return selfout,argout
    #}}}
    #{{{ integrate, differentiate, and sum
    def integrate(self,thisaxis,backwards = False):
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        self.run_nopop(cumsum,thisaxis)
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        if len(self.axis_coords)>0:
            t = self.getaxis(thisaxis)
            dt = t[1]-t[0]
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
        if (type(axes) is str):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print '|-ERROR FINDING DIMENSION-----'
                print '| dimlabels is: ',self.dimlabels
                print "| doesn't contain: ",axes[j]
                print '|-----------------------------'
                raise
            self.data = sum(self.data,
                    axis=thisindex)
            self._pop_axis_info(thisindex)
        return self
    def sum_nopop(self,axes):
        if (type(axes) is str):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
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
        'return the coefficients and the fit --> later, should probably branch this off as a new type of fit class'
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
        if force_y_intercept != None:
            startingpower = 1
        L =  concatenate([x**j for j in range(startingpower,order+1)],axis=1) # note the totally AWESOME way in which this is done!
        #print 'fitting to matrix',L
        if force_y_intercept != None:
            y -= force_y_intercept
        c = dot(pinv(L),y)
        fity = dot(L,c)
        if force_y_intercept != None:
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
        r'find the max along a particular axis, and get rid of that axis, replacing it with the index number of the max value'
        #{{{ process arguments
        axes = self._possibly_one_axis(*args)
        raw_index = False
        if 'raw_index' in kwargs.keys():
            raw_index = kwargs.pop('raw_index')
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:",repr(kwargs))
        if (type(axes) is str):
            axes = [axes]
        #}}}
        for j in range(0,len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
                raise
            if raw_index:
                self.data = argmax(self.data,
                        axis=thisindex)
            else:
                self.data = self.axis_coords[thisindex][argmax(self.data,
                    axis=thisindex)]
            self._pop_axis_info(thisindex)
        return self
    def argmin(self,axes,raw_index = False):
        r'find the min along a particular axis, and get rid of that axis, replacing it with the index number of the min value'
        if (type(axes) is str):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
                raise
            if raw_index:
                self.data = argmin(self.data,
                        axis=thisindex)
            else:
                self.data = self.axis_coords[thisindex][argmin(self.data,
                    axis=thisindex)]
            self._pop_axis_info(thisindex)
        return self
    def mean_all_but(self,listofdims):
        'take the mean over all dimensions not in the list'
        for dimname in self.dimlabels:
            if not (dimname in listofdims):
                self.mean(dimname)
        return self
    def mean(self,*args,**kwargs):
        r'Take the mean and set the error to the standard deviation'
        #{{{ process arguments
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        axes = self._possibly_one_axis(*args)
        return_error = False
        if return_error in kwargs.keys():
            return_error = kwargs.pop('return_error')
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:",repr(kwargs))
        if (type(axes) is str):
            axes = [axes]
        #}}}
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
                raise
            if self.data_error is not None:
                this_axis_length = self.data.shape[thisindex]
                try:
                    self.data_error = sqrt(sum((self.data*self.data_error)**2,
                            axis=thisindex)/(this_axis_length**2))
                except:
                    raise CustomError('shape of data',shape(self.data),'shape of data error',shape(self.data_error))
            if return_error: # since I think this is causing an error
                thiserror = std(self.data,
                        axis=thisindex)
                if isscalar(thiserror):
                    thiserror = r_[thiserror]
                self.set_error(thiserror) # set the error to the standard deviation
            self.data = mean(self.data,
                    axis=thisindex)
            self._pop_axis_info(thisindex)
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
            try:
                self.axis_coords_error.pop(thisindex)
            except:
                raise CustomError('trying to pop',thisindex,'from',self.axis_coords_error)
        return self
    def popdim(self,dimname):
        thisindex = self.axn(dimname)
        thisshape = list(self.data.shape)
        if thisshape[thisindex]!=1:
            raise CustomError("trying to pop a dim that's not length 1")
        thisshape.pop(thisindex)
        self.data = self.data.reshape(thisshape)
        self._pop_axis_info(thisindex)
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
        func = args[0]
        func = self._wrapaxisfuncs(func)
        if len(args)>1:
            axis = args[1]
            thisindex = self.dimlabels.index(axis)
            self.data = func(self.data,axis=thisindex)
            self._pop_axis_info(thisindex)
        else:
            self.data = func(self.data)
        return self
    def run_nopop(self,func,axis):
        func = self._wrapaxisfuncs(func)
        try:
            thisaxis = self.dimlabels.index(axis)
        except:
            raise CustomError("I couldn't find the dimension",axis,"in the list of axes",self.dimlabels)
        temp = list(self.data.shape)
        temp[thisaxis] = 1
        numnonoptargs = len(getargspec(func)[0])-len(getargspec(func)[3])
        if numnonoptargs == 1:
            self.data = func(self.data,axis=thisaxis)
        elif numnonoptargs == 2:
            self.data = func(self.getaxis(axis),self.data,axis=thisaxis)
        else:
            raise CustomError('you passed a function to run_nopop that doesn\'t have either one or two arguments!')
        #{{{ if the function doesn't rip out the dim, make sure we don't change the dims
        if len(self.data.shape)==len(temp):
            temp[thisaxis] = self.data.shape[thisaxis]
        #}}}
        self.data = self.data.reshape(temp)
        return self
    #}}}
    #{{{ ft-related functions
    def convolve(self,axisname,filterwidth,convfunc = (lambda x,y: exp(-(x**2)/(2.0*(y**2))))):
        r'''perform a normalized convolution'''
        #{{{ make a version of x that is oriented along the correct dimension
        x = self.getaxis(axisname).copy()
        x_centerpoint = (x[-1]+x[0])/2
        x -= x_centerpoint # so that zero is in the center
        x = ifftshift(x) # so that it's convolved about time 0
        thisaxis = self.axn(axisname)
        #}}}
        myfilter = convfunc(x,filterwidth)
        myfilter /= myfilter.sum()
        filtershape = ones_like(self.data.shape)
        filtershape[thisaxis] = len(myfilter)
        myfilter = myfilter.reshape(filtershape)
        self.data = ifft(fft(self.data,axis = thisaxis)*fft(myfilter,axis = thisaxis),axis = thisaxis)
        return self
    def _ft_conj(self,x):
        pairs = [('s','Hz'),('m',r'm^{-1}')]
        a,b = zip(*tuple(pairs))
        if x in a:
            return b[a.index(x)]
        elif x in b:
            return a[b.index(x)]
        else:
            return None
    def ft(self,*args,**kwargs):
        #{{{ process arguments
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        axes = self._possibly_one_axis(*args)
        #kwargs: shiftornot=False,shift=None,pad = False
        shiftornot,shift,pad = self._process_kwargs([('shiftornot',False),('shift',None),('pad',False)],**kwargs)
        if shift != None:
            shiftornot = shift
        if (type(axes) is str):
            axes = [axes]
        if not (type(shiftornot) is list):
            shiftornot = [bool(shiftornot)]*len(axes)
        #}}}
        for j in range(0,len(axes)):
            if self.get_units(axes[j]) is not None:
                self.set_units(axes[j],self._ft_conj(self.get_units(axes[j])))
            try:
                thisaxis = self.dimlabels.index(axes[j])
            except:
                raise CustomError('error, dimlabels is: ',self.dimlabels)
            padded_length = self.data.shape[thisaxis]
            if pad is True:
                padded_length = 2**(ceil(log2(padded_length)))
            elif pad:
                padded_length = pad
            self.data = fft(self.data,n = padded_length,axis=thisaxis)
            if bool(shiftornot[j]):
                self.data = fftshift(self.data,axes=[thisaxis])
            if len(self.axis_coords)>0:
                t = self.axis_coords[thisaxis]
                dt = t[1]-t[0] # the dwell gives the bandwidth, whether or not it has been zero padded
                self.ft_start_time = t[0]
                self.data *= dt
                self.axis_coords[thisaxis] = linspace(0,1./dt,padded_length)
                if bool(shiftornot[j]):
                    mask = self.axis_coords[thisaxis] > 0.5/dt
                    self.axis_coords[thisaxis][mask] -= (1.+1./padded_length)/dt
                    self.axis_coords[thisaxis] = fftshift(self.axis_coords[thisaxis])
        return self
    def ift(self,axes,shiftornot=False,shift=None):
        if shift != None:
            shiftornot = shift
        if (type(axes) is str):
            axes = [axes]
        if not (type(shiftornot) is list):
            shiftornot = [bool(shiftornot)]
        for j in range(0,len(axes)):
            if self.get_units(axes[j]) is not None:
                self.set_units(axes[j],self._ft_conj(self.get_units(axes[j])))
            try:
                thisaxis = self.dimlabels.index(axes[j])
            except:
                raise CustomError('error, dimlabels is: ',self.dimlabels)
            if bool(shiftornot[j]):
                self.data = ifftshift(self.data,axes=[thisaxis])
            self.data = ifft(self.data,axis=thisaxis)
            if len(self.axis_coords)>0:
                t = self.axis_coords[thisaxis]
                dt = t[1]-t[0]
                self.data *= size(t) * dt # here, the algorithm divides by N, so for integration, we need to not do that
                #{{{ shiftornot specifies the shifting of the initial ft, not this result, so we always return a 0->1 time axis
                self.axis_coords[thisaxis] = linspace(0,1./dt,size(t)) + self.ft_start_time # note that I offset by ft_start_time, which I pull from when I ft'd
                #}}}
        return self
    #}}}
    #{{{ interpolation
    def interp(self,axis,axisvalues,past_bounds = None,verbose = False,**kwargs):
        'interpolate data values given axis values'
        oldaxis = self.getaxis(axis)
        if (type(axisvalues) is int) or (type(axisvalues) is int32):
            axisvalues = linspace(oldaxis[0],oldaxis[-1],axisvalues)
        elif (type(axisvalues) not in [ndarray,tuple]):
            raise ValueError("You passed a target axis of type"+repr(type(axisvalues))+"which I don't get")
        if past_bounds == None:
            axisvalues[axisvalues<oldaxis.min()] = oldaxis.min()
            axisvalues[axisvalues>oldaxis.max()] = oldaxis.max()
        elif not (past_bounds == 'fail'):
            if type(past_bounds) is tuple:
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
        if 'kind' in kwargs.keys():
            thiskind = kwargs.pop('kind')
        else:
            thiskind = 'cubic'
            if len(rdata) < 4:
                thiskind = 'quadratic'
                if len(rdata) < 3:
                    thiskind = 'linear'
        thisaxis = self.axn(axis)
        try:
            if verbose: print 'Using %s interpolation'%thiskind
            interpfunc =  interp1d(oldaxis,rdata,kind = thiskind,axis = thisaxis)
            #interpfunc =  UnivariateSpline(oldaxis,rdata,s=1)
            rdata = interpfunc(axisvalues)
            interpfunc =  interp1d(oldaxis,idata,kind = thiskind,axis = thisaxis)
            #interpfunc =  UnivariateSpline(oldaxis,idata,s=1)
            idata = interpfunc(axisvalues)
            self.data = rdata + 1j * idata
        except ValueError:
            raise CustomError('Value error, bounds for data are',self.getaxis(axis)[r_[0,-1]],'bounds for new axis are',axisvalues[r_[0,-1]],
                    'shapes are',
                    shape(axisvalues),shape(oldaxis),shape(self.data),
                    'dtypes are',
                    axisvalues.dtype,self.getaxis(axis).dtype,self.data.dtype)

        except TypeError:
            raise CustomError('Type error, types are',
                    type(axisvalues),type(self.getaxis(axis)),type(self.data),
                    'shapes are',
                    shape(axisvalues),shape(self.getaxis(axis)),shape(self.data),
                    'dtypes are',
                    axisvalues.dtype,self.getaxis(axis).dtype,self.data.dtype)
        self.setaxis(axis,axisvalues)
        return self
    def invinterp(self,axis,values):
        'interpolate axis values given data values'
        origdata = self.data.copy()
        origaxis = self.getaxis(axis).copy()
        self.data = values
        try:
            self.setaxis(axis,interp(values,origdata,origaxis))
        except:
            raise CustomError('problem interpolating for axis',axis,'values',values,'with axis length',shape(self.getaxis(axis)),'and data shape',shape(self.data))
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
    #{{{ functions to manipulate and return the axes
    def reorder(self,*axes):
        if len(axes) == 1:
            axes = axes[0]
        else:
            axes = axes
        try:
            neworder = map(self.dimlabels.index,axes)
        except ValueError:
            raise CustomError('one of',axes,'not in',self.dimlabels)
        self.dimlabels = map(self.dimlabels.__getitem__,neworder)
        if len(self.axis_coords)>0:
            try:
                self.axis_coords = map(self.axis_coords.__getitem__,neworder)
            except:
                raise CustomError('problem mapping',map(len,self.axis_coords),'onto',neworder)
            if len(self.axis_coords_units)>0:
                try:
                    self.axis_coords_units = map(self.axis_coords_units.__getitem__,neworder)
                except:
                    raise CustomError('problem mapping',map(len,self.axis_coords_units),'onto',neworder)
        try:
            self.data = self.data.transpose(neworder)
        except ValueError:
            raise CustomError('you can\'t reorder',self.dimlabels,'as',neworder)
        if self.data_error is not None:
            self.data_error = self.data_error.transpose(neworder)
        return self
    def plot_labels(self,labels,fmt = None):
        r'this only works for one axis now'
        axisname = self.dimlabels[0]
        if format is None:
            plot_label_points(self.getaxis(axisname),self.data,labels)
        else:
            plot_label_points(self.getaxis(axisname),self.data,[fmt%j for j in labels])
        return
    def labels(self,listofstrings,listofaxes):
        r'''label the dimensions, given in listofstrings with the axis labels given in listofaxes -- listofaxes must be a numpy array;
        you can pass either a list or a axis name (string)/axis label (numpy array) pair'''
        if len(self.axis_coords) == 0:
            self.axis_coords = [[]]*len(self.dimlabels)
            self.axis_coords_error = [None]*len(self.dimlabels)
        if type(listofstrings) is str:
            listofstrings = [listofstrings]
            listofaxes = [listofaxes]
        if type(listofstrings) is not list:
            raise TypeError("the arguments passed to the .labels() method must be a list of the axis names followed by the list of the axis arrays")
        elif all(map(( lambda x: type(x) is str_ ),listofstrings)):
            listofstrings = map(str,listofstrings)
        elif not all(map(( lambda x: type(x) in [str,str_] ),listofstrings)):
            raise TypeError("the arguments passed to the .labels() method must be a list of the axis names followed by the list of the axis arrays")
        for j in range(0,len(listofstrings)):
            if listofaxes[j] is None:
                self.axis_coords[self.axn(listofstrings[j])] = None
            else:
                #{{{ test that the axis is the right size
                if isscalar(listofaxes[j]):#interpret as a timestep
                    listofaxes[j] = listofaxes[j] * r_[0:ndshape(self)[listofstrings[j]]]
                if type(listofaxes[j]) not in [ndarray,list]:
                    raise TypeError('You passed an axis label of type '+repr(type(listofaxes[j]))+' for the axis '+listofstrings[j]+' to the labels method, which you can\'t do --> it must be an nddata')
                if (len(listofaxes[j]) != ndshape(self)[listofstrings[j]]) and (len(listofaxes[j])!=0):
                    raise IndexError("You're trying to attach an axis of len %d to the '%s' dimension, which has %d data points"%(len(listofaxes[j]),listofstrings[j],ndshape(self)[listofstrings[j]]))
                #}}}
                try:
                    self.axis_coords[self.dimlabels.index(listofstrings[j])] = listofaxes[j]
                except:
                    try:
                        raise CustomError("Can't assign the coordinates to"+listofstrings[j]+"as"+listofaxes[j].__repr__())
                    except:
                        raise CustomError("listofaxes (",len(listofaxes),") isn't same length as ",listofstrings)
        return self
    def check_axis_coords_errors(self):
        if len(self.axis_coords_error) > len(self.dimlabels):
            raise ValueError("this failed because there are more sets of axis errors than there are axes!\nlen(axis_coords_error) = %s\naxes = %s"%(repr(len(self.axis_coords_error)),repr(self.dimlabels)))
    def sort(self,axisname):
        whichaxis = self.dimlabels.index(axisname)
        order = argsort(self.axis_coords[whichaxis])
        datacopy = self.copy()
        for j in range(0,len(order)): # do it this way, so that it deals with other dimensions correctly
            self.check_axis_coords_errors()
            self[axisname,j] = datacopy[axisname,order[j]]
        self.axis_coords[whichaxis] = self.axis_coords[whichaxis][order]
        return self
    def copyaxes(self,other):
        # in the case that the dimensions match, and we want to copy the labels
        self.axis_coords = other.axis_coords
        return self
    def retaxis(self,axisname):
        newshape = ones(len(self.data.shape),dtype = 'uint')
        thisaxis = self.axn(axisname)
        newshape[thisaxis] = self.data.shape[thisaxis]
        newshape = list(newshape)
        return nddata(self.getaxis(axisname).copy().reshape(newshape),newshape,list(self.dimlabels))
    def getaxis(self,*axisname):
        try:
            return self.axis_coords[self.axn(*axisname)]
        except:
            raise CustomError(axisname,'does not seem to have axis labels')
    def setaxis(self,axis,value):
        if type(value) in [float,int,double,float64]:
           value = linspace(0.,value,self.axlen(axis))
        self.axis_coords[self.axn(axis)] = value
    def getaxisshape(self,axisname):
        thishape = ones(len(self.dimlabels))
        thisaxis = self.dimlabels.index(axisname) 
        thishape[thisaxis] = self.data.shape[thisaxis]
        return thishape
    def circshift(self,axis,amount):
        if amount!=0:
            if abs(amount) > ndshape(self)[axis]:
                CustomError("Trying to circshift by ",amount,"which is bitter than the size of",axis)
            newdata = ndshape(self).alloc(dtype=self.data.dtype)
            newdata[axis,:-amount] = self[axis,amount:]
            newdata[axis,-amount:] = self[axis,:amount]
            self.data = newdata.data
        return self
    #}}}
    #{{{ breaking up and combining axes
    def smash(self,listofhowmany):
        r'''collapse multiple dimensions into one dimension'''
        names = ['error'] * len(listofhowmany)
        newsize = []
        laststop = 0 
        j = 0
        for ndims in listofhowmany:
            thisstop = laststop + ndims
            newsize += [prod(self.data.shape[laststop:thisstop])]
            names[j] = ' x '.join(self.dimlabels[laststop:thisstop])
            laststop = thisstop
            j+=1
        self.data = self.data.reshape(newsize)
        self.dimlabels = names
        return self
    def smashorder(self,listoflists):
        self.reorder(concatenate(listoflists))
        self.smash(map(len,listoflists))
        return self
    def chunk(self,axisin,axesout,shapesout):
        if prod(shapesout) != ndshape(self)[axisin]:
            raise CustomError("The size of the axis (%s) you're trying to split (%s) doesn't match the size of the axes you're trying to split it into (%s = %s)"%(repr(axisin),repr(ndshape(self)[axisin]),repr(axesout),repr(shapesout)))
        thisaxis = self.dimlabels.index(axisin)
        if len(self.axis_coords[thisaxis]) > 0:
            raise CustomError("split not yet supported on axes with labels")
        newaxis_coords = self.axis_coords[0:thisaxis] + [[]]*len(axesout) + self.axis_coords[thisaxis+1:]
        newaxis_coords_error = self.axis_coords_error[0:thisaxis] + [None]*len(axesout) + self.axis_coords_error[thisaxis+1:]
        newshape = list(self.data.shape[0:thisaxis]) + shapesout + list(self.data.shape[thisaxis+1:])
        newnames = list(self.dimlabels[0:thisaxis]) + axesout + list(self.dimlabels[thisaxis+1:])
        self.data = self.data.reshape(newshape)
        self.dimlabels = newnames
        self.axis_coords = newaxis_coords
        self.axis_coords_error = newaxis_coords_error
        return self
    def chunkoff(self,axisin,newaxes,newshapes):
        r'''chunks up axisin, dividing it into newaxes with newshapes on the inside'''
        axesout = [axisin]+newaxes
        shapesout = [ndshape(self)[axisin]/prod(newshapes)]+newshapes
        return self.chunk(axisin,axesout,shapesout)
    #}}}
    #{{{ messing with data -- get, set, and copy
    def __getslice__(self,*args):
        print 'getslice! ',args
    def __setitem__(self,*args):
        righterrors = None
        if isinstance(args[1],nddata):
            #{{{ reorder so the shapes match
            unshared_indeces = list(set(args[1].dimlabels) ^ set(self.dimlabels))
            shared_indeces = list(self.dimlabels)
            for j in unshared_indeces:
                shared_indeces.remove(j)
            if len(args[1].dimlabels) != len(shared_indeces) or (not all([args[1].dimlabels[j] == shared_indeces[j] for j in range(0,len(shared_indeces))])):
                args[1].reorder[shared_indeces]
            #}}}
            rightdata = args[1].data
            righterrors = args[1].get_error()
        else: # assume it's an ndarray
            rightdata = args[1]
            if (type(rightdata) is not ndarray): # in case its a scalar
                rightdata = array([rightdata])
        slicedict,axesdict,errordict = self._parse_slices(args[0]) # pull left index list from parse slices
        leftindex = self.fld(slicedict)
        try:
            self.data[tuple(leftindex)] = rightdata.squeeze() # assign the data
        except:
            raise CustomError('ERROR ASSIGNING NDDATA:','rightdata.shape:',rightdata.shape,'self.data.shape:',self.data.shape,' dimlabels:',self.dimlabels,' leftindex:',tuple(leftindex),'--> shape of left slice: ',self.data[tuple(leftindex)].shape)
            raise
        lefterror = self.get_error()
        if lefterror is not None:
            lefterror[tuple(leftindex)] = righterrors.squeeze()
    def copy(self):
        return deepcopy(self)
    def __getitem__(self,args):
        try:
            slicedict,axesdict,errordict = self._parse_slices(args)
        except:
            raise CustomError('error trying to get slices given by',args)
        if type(args[1]) is list and type(args[0]) is str and len(args) == 2:
            return concat([self[args[0],x] for x in args[1]],args[0])
        indexlist = self.fld(slicedict)
        newlabels = [x for x in self.dimlabels if not isscalar(slicedict[x])] # generate the new list of labels, in order, for all dimensions that are not indexed by a scalar
        #{{{ properly index the data error
        if self.data_error != None:
            try:
                newerror = self.data_error[indexlist]
            except:
                raise ValueError('Problem trying to index data_error'+repr(self.data_error)+' with',repr(indexlist))
        else:
            newerror = None
        #}}}
        if len(self.axis_coords)>0:
            if errordict != None:
                axis_coords_error = [errordict[x] for x in newlabels]
            else:
                axis_coords_error = None
            retval =  nddata(self.data[indexlist],
                    self.data[indexlist].shape,
                    newlabels,
                    axis_coords = [axesdict[x] for x in newlabels],
                    axis_coords_error = axis_coords_error,
                    data_error = newerror,
                    other_info = self.other_info)
            retval.axis_coords_units = self.axis_coords_units
            retval.data_units = self.data_units
            return retval
        else:
            retval = nddata(self.data[indexlist],self.data[indexlist].shape,newlabels,self.other_info)
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
    def _process_kwargs(self,listoftuples,**kwargs):
        kwargnames,kwargdefaultvals = zip(*listoftuples)
        output = []
        for j,val in enumerate(kwargnames):
            output.append(kwargdefaultvals[j])
            if val in kwargs.keys():
                output[-1] = kwargs.pop(val)
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:",repr(kwargs))
        return tuple(output)
    def _parse_slices(self,args):
        #print "DEBUG getitem called with",args
        errordict = None # in case it's not set
        axesdict = None # in case it's not set
        args = list(args)
        for j in range(0,len(args),2):
            if type(args[j]) is str_:
                args[j] = str(args[j]) # on upgrading + using on windows, this became necessary, for some reason I don't understand
        if type(args) in [float,int32,int,double]:
            raise CustomError('You tried to pass just a nddata[',type(args),']')
        if type(args[0]) is str:
            #{{{ create a slicedict and errordict to store the slices
            slicedict = dict(zip(list(self.dimlabels),[slice(None,None,None)]*len(self.dimlabels))) #initialize to all none
            if len(self.axis_coords)>0:
                #print "DEBUG --> trying to make dictionaries from axis coords of len",len(self.axis_coords),"and axis_coords_error of len",len(self.axis_coords_error),"when dimlabels has len",len(self.dimlabels)
                axesdict = self.mkd(self.axis_coords)
                errordict = self.mkd(self.axis_coords_error)
                #print 'DEBUG: made errordict',errordict
            for x,y in zip(args[0::2],args[1::2]):
                slicedict[x] = y
            #}}}
            #{{{ map the slices onto the axis coordinates and errors
            #print "DEBUG slicedict is",slicedict
            testf = lambda x: x+1
            if len(self.axis_coords)>0:
                for x,y in slicedict.iteritems():
                    #print "DEBUG, type of slice",x,"is",type(y)
                    if isscalar(y):
                        axesdict.pop(x) # pop the axes for all scalar dimensions
                    #elif type(y) is slice and type(y.start) in [float,double,float64,float32]:
                    #    mask = logical_and(axesdict[x]>y.start,axesdict[x]<y.stop)
                    #    slicedict[x] = mask
                    #    axesdict[x] = axesdict[x][mask]
                    elif type(y) is type(testf):
                        mask = y(axesdict[x])
                        slicedict[x] = mask
                        axesdict[x] = axesdict[x][mask]
                    else:
                        axesdict[x] = axesdict[x][y]
                if errordict != None and errordict != array(None):
                    for x,y in slicedict.iteritems():
                        if errordict[x] != None:
                            if isscalar(y):
                                errordict.pop(x)
                            elif type(y) is type(emptyfunction):
                                mask = y(axesdict[x])
                                errordict[x] = errordict[x][mask]
                            else:
                                try:
                                    errordict[x] = errordict[x][y] # default
                                except:
                                    raise CustomError('Trying to index',errordict,'-->',x,'=',errordict[x],'with',y,'error started as',self.axis_coords_error)
                    errordict = self.mkd().update(errordict) # make everything none except the ones I've updated
            return slicedict,axesdict,errordict
            #}}}
        else:
            raise CustomError('label your freaking dimensions! (type of args[0] is ',type(args[0]),'and it should be str!)')
    #}}}
    #{{{ hdf5 write
    def hdf5_write(self,h5path,verbose = False):
        #{{{ add the final node based on the name stored in the nddata structure
        if h5path[-1] != '/': h5path += '/' # make sure it ends in a slash first
        try:
            thisname = self.get_prop('name')
        except:
            raise CustomError("You're trying to save an nddata object which does not yet have a name, and you can't do this! Run yourobject.name('setname')")
        if type(thisname) is str:
            h5path += thisname
        else:
            raise CustomError("problem trying to store HDF5 file; you need to set the ``name'' property of the nddata object to a string first!")
        h5file,bottomnode = h5nodebypath(h5path) # open the file and move to the right node
        #print 'DEBUG 1: bottomnode is',bottomnode
        #}}}
        #{{{ print out the attributes of the data
        myattrs = normal_attrs(self)
        #{{{ separate them into data and axes
        mydataattrs = filter((lambda x: x[0:4] == 'data'),myattrs)
        myotherattrs = filter((lambda x: x[0:4] != 'data'),myattrs)
        myaxisattrs = filter((lambda x: x[0:4] == 'axis'),myotherattrs)
        myotherattrs = filter((lambda x: x[0:4] != 'axis'),myotherattrs)
        if verbose: print lsafe('data attributes:',zip(mydataattrs,map(lambda x: type(self.__getattribute__(x)),mydataattrs))),'\n\n'
        if verbose: print lsafe('axis attributes:',zip(myaxisattrs,map(lambda x: type(self.__getattribute__(x)),myaxisattrs))),'\n\n'
        if verbose: print lsafe('other attributes:',zip(myotherattrs,map(lambda x: type(self.__getattribute__(x)),myotherattrs))),'\n\n'
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
            if verbose: print "Writing remaining axis attributes\n\n"
            if len(mydataattrs) > 0:
                h5attachattributes(datatable,mydataattrs,self)
        else:
            raise CustomError("I can't find the data object when trying to save the HDF5 file!!")
        #}}}
        #{{{ write the axes tables
        if 'axis_coords' in myaxisattrs:
            if len(self.axis_coords) > 0:
                #{{{ create an 'axes' node
                axesnode = h5child(bottomnode, # current node
                        'axes', # the child
                        verbose = False,
                        create = True)
                #}}}
                for j,axisname in enumerate(self.dimlabels): # make a table for each different dimension
                    myaxisattrsforthisdim = dict([(x,self.__getattribute__(x)[j])
                        for x in list(myaxisattrs) if len(self.__getattribute__(x)) > 0]) # collect the attributes for this dimension and their values
                    if verbose: print lsafe('for axis',axisname,'myaxisattrsforthisdim=',myaxisattrsforthisdim)
                    if 'axis_coords' in myaxisattrsforthisdim.keys() and myaxisattrsforthisdim['axis_coords'] is not None:
                        if 'axis_coords_error' in myaxisattrsforthisdim.keys() and myaxisattrsforthisdim['axis_coords_error'] is not None and len(myaxisattrsforthisdim['axis_coords_error']) > 0: # this is needed to avoid all errors, though I guess I could use try/except
                            thistable = rec.fromarrays([myaxisattrsforthisdim['axis_coords'],myaxisattrsforthisdim['axis_coords_error']],names='data,error')
                            myaxisattrsforthisdim.pop('axis_coords_error')
                        else:
                            thistable = rec.fromarrays([myaxisattrsforthisdim['axis_coords']],names='data')
                        myaxisattrsforthisdim.pop('axis_coords')
                    datatable = h5table(axesnode,axisname,thistable)
                    #print 'DEBUG 3: axesnode is',axesnode
                    if verbose: print "Writing remaining axis attributes for",axisname,"\n\n"
                    if len(myaxisattrsforthisdim) > 0:
                        h5attachattributes(datatable,myaxisattrsforthisdim.keys(),myaxisattrsforthisdim.values())
        #}}}
        #{{{ Check the remaining attributes.
        if verbose: print lsafe('other attributes:',zip(myotherattrs,map(lambda x: type(self.__getattribute__(x)),myotherattrs))),'\n\n'
        if verbose: print "Writing remaining other attributes\n\n"
        if len(myotherattrs) > 0:
            #print 'DEBUG 4: bottomnode is',bottomnode
            test = repr(bottomnode) # somehow, this prevents it from claiming that the bottomnode is None --> some type of bug?
            try:
                h5attachattributes(bottomnode,
                    [j for j in myotherattrs if not self._contains_symbolic(j)],
                    self)
            except:
                raise CustomError('Problem trying to attach attributes',myotherattrs,'of self to node',bottomnode)
            warnlist = [j for j in myotherattrs if (not self._contains_symbolic(j)) and type(self.__getattribute__(j)) is dict]
            #{{{ to avoid pickling, test that none of the attributes I'm trying to write are dictionaries or lists
            if len(warnlist) > 0:
                print "WARNING!!, attributes",warnlist,"are dictionaries!"
            warnlist = [j for j in myotherattrs if (not self._contains_symbolic(j)) and type(self.__getattribute__(j)) is list]
            if len(warnlist) > 0:
                print "WARNING!!, attributes",warnlist,"are lists!"
            #}}}
            if verbose: print lsafe('other attributes:',zip(myotherattrs,map(lambda x: type(self.__getattribute__(x)),myotherattrs))),'\n\n'
        #}}}
        h5file.close()
    #}}}
class nddata_hdf5 (nddata):
    def __repr__(self):
        if hasattr(self,'_node_children'):
            return repr(self.datanode)
        else:
            return nddata.__repr__(self)
    def __del__(self):
        if hasattr(self,'_node_children'):
            self.h5file.close()
            del self.h5file
            del self.datanode
        return
    def __init__(self,pathstring):
        self.pathstring = pathstring
        try:
            self.h5file,self.datanode = h5nodebypath(pathstring,check_only = True)
        except:
            raise IndexError("I can't find the node"+pathstring)
        self._init_datanode(self.datanode)
    def _init_datanode(self,datanode,verbose = False,**kwargs):
        datadict = h5loaddict(datanode)
        #{{{ load the data, and pop it from datadict
        try:
            datarecordarray = datadict['data']['data'] # the table is called data, and the data of the table is called data
            mydata = datarecordarray['data']
        except:
            raise CustomError("I can't find the nddata.data")
        try:
            kwargs.update({'data_error':datarecordarray['error']})
        except:
            if verbose: print "No error found\n\n"
        datadict.pop('data')
        #}}}
        #{{{ be sure to load the dimlabels
        mydimlabels = datadict['dimlabels']
        if len(mydimlabels) == 1:
            if len(mydimlabels[0]) == 1:
                mydimlabels = list([mydimlabels[0][0]]) # for some reason, think I need to do this for length 1
        #}}}
        #{{{ load the axes and pop them from datadict
        datadict.pop('dimlabels')
        if 'axes' in datadict.keys():
            myaxiscoords = [None]*len(mydimlabels)
            myaxiscoordserror = [None]*len(mydimlabels)
            for axisname in datadict['axes'].keys():
                try:
                    axisnumber = mydimlabels.index(axisname)
                except AttributeError:
                    raise CustomError('mydimlabels is not in the right format!\nit looks like this:\n',mydimlabels,type(mydimlabels))
                recordarrayofaxis = datadict['axes'][axisname]['data']
                myaxiscoords[axisnumber] = recordarrayofaxis['data']
                if 'error' in recordarrayofaxis.dtype.names:
                    myaxiscoordserror[axisnumber] = recordarrayofaxis['error']
                datadict['axes'][axisname].pop('data')
                for k in datadict['axes'][axisname].keys():
                    print lsafen("Warning, attribute",k,"of axis table",axisname,"remains, but the code to load this is not yet supported")
                datadict['axes'].pop(axisname)
            kwargs.update({"axis_coords":myaxiscoords})
            kwargs.update({"axis_coords_error":myaxiscoordserror})
        elif len(mydimlabels)>1:
            raise CustomError("The current version uses the axis labels to figure out the shape of the data\nBecause you stored unlabeled data, I can\'t figure out the shape of the data!!")
            # the reshaping this refers to is done below
        #}}}
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
                if type(temp) is ndarray:
                    temp = len(temp)
                det_shape.append(temp)
            self.data = self.data.reshape(tuple([len(self.getaxis(x)) for x in mydimlabels]))
        #}}}
        for remainingattribute in datadict.keys():
            self.__setattr__(remainingattribute,datadict[remainingattribute])
        self.h5file.close()
        del self.h5file
        del self.datanode
        return
#}}}

#{{{subplot_dim
class subplot_dim():
    def __init__(self,firstdim,seconddim):
        self.num = r_[firstdim,seconddim,0]
    def set(self,args,x='',g=True,y='',t='',a=''):
        if type(args) is int:
            number = args
            ax = subplot(*tuple(self.num+r_[0,0,number]))
            xlabel(x)
            ylabel(y)
            title(t)
            grid(g)
        elif (type(args) is tuple) and (len(args) is 3):
            # the second value passed is 
            whichsmall = args[2]
            break_into = args[1]
            number = args[0]
            mydims = self.num*r_[1,break_into,1]+r_[
                    0,0,break_into*(number-1)+whichsmall]
            try:
                ax = subplot(*tuple(mydims))
            except:
                print 'failed trying subplots: ', mydims
                raise
            xlabel(x)
            ylabel(y)
            title(t)
            grid(g)
        else:
            print "problem, need to pass either 1 or 3 arguments to set"
            print 'type of args: ',type(args)
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
        raise CustomError('pinvr error, U is not finite')
    if any(~isfinite(V)):
        raise CustomError('pinvr error, V is not finite')
    if any(~isfinite(S)):
        raise CustomError('pinvr error, S is not finite')
    S = diag(S / (S**2 + alpha**2))
    if any(~isfinite(S)):
        raise CustomError('pinvr error, problem with S/(S^2+alpha^2) --> set your regularization higher')
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
def image(A,x=[],y=[],**kwargs):
    #{{{ pull out kwargs for imagehsv
    imagehsvkwargs = {}
    for k,v in kwargs.items():
        if k in ['black','logscale']:
            imagehsvkwargs[k] = kwargs.pop(k)
    #}}}
    setlabels = False
    if isinstance(A,nddata):
        setlabels = True
        templabels = list(A.dimlabels)
        x_label = templabels[-1]
        try:
            x = list(A.getaxis(x_label))
        except:
            x = r_[0,ndshape(A)[x_label]]
        templabels.pop(-1)
        y_label = ''
        if len(templabels) == 1:
            try:
                y = list(A.getaxis(templabels[0]))
            except:
                y = r_[0:A.data.shape[A.axn(templabels[0])]]
        while len(templabels)>0:
            y_label += templabels.pop(0)
            if len(templabels)>0:
                y_label += '$\\times$'
        A = A.data
    if type(x) is list:
        x = array(x)
    if type(y) is list:
        y = array(y)
    if len(x)==0:
        x = [1,A.shape[1]]
    else:
        x = x.flatten()
    if len(y)==0:
        y = [1,A.shape[0]]
    else:
        y = y.flatten()
    myext = (x[0],x[-1],y[-1],y[0])
    extralines = 0
    origAndim = A.ndim
    if A.ndim > 2:
        ax = gca()
        setp(ax.get_yticklabels(),visible = False)
        ax.yaxis.set_ticks_position("none")
    while A.ndim>2:# to substitude for imagehsvm, etc., so that we just need a ersion of ft
        # order according to how it's ordered in the memory
        # the innermost two will form the image -- first add a line to the end of the images we're going to join up
        tempsize = array(A.shape) # make a tuple the right shape
        tempsize[-2] = 1 # all dims are the same except the image row, of which there is only one
        A = concatenate((A,nan*zeros(tempsize)),axis=(A.ndim-2)) # concatenate along the rows
        tempsize = r_[A.shape[0:-3],A.shape[-2:]]
        tempsize[-2] *= A.shape[-3]
        A = A.reshape(tempsize) # now join them up
        ++extralines # keep track of the extra lines at the end
    A = A[:A.shape[0]-extralines,:]
    #line_mask = isnan(A)
    #A[line_mask] = A[logical_not(line_mask)].max()
    #A[line_mask] = 0
    if iscomplex(A).any():
        A = imagehsv(A,**imagehsvkwargs)
        imshow(A,extent=myext,**kwargs)
    else:
        imshow(A,extent=myext,**kwargs)
        colorbar()
    if setlabels:
        xlabel(x_label)
        #print y_label
        ylabel(y_label)
    return
def colormap(points,colors,n=256):
    r = interp(linspace(0,1,n),points,colors[:,0].flatten())
    g = interp(linspace(0,1,n),points,colors[:,1].flatten())
    b = interp(linspace(0,1,n),points,colors[:,2].flatten())
    return reshape(r_[r,g,b],(3,n)).T
def imagehsv(A,logscale = False,black = False):
    n = 256
    mask = isnan(A)
    A[mask] = 0
    theta = (n-1.)*mod(angle(A)/pi/2.0,1)# angle in 255*cycles
    hsv = colormap(r_[0.,1./3.,2./3.,1.],double(array([
        [1,0,0],
        [0,1,0],
        [0,0,1],
        [1,0,0]])),n=n)
    hsv_norm = sqrt(sum(hsv*hsv,axis=1))
    hsv_norm = reshape(hsv_norm,(hsv_norm.size,1))
    hsv = hsv/hsv_norm
    colors = hsv[ix_(int32(theta.flatten().round()),[0,1,2])]
    colors = reshape(colors,(A.shape[0],A.shape[1],3))
    mask = mask.reshape(A.shape[0],A.shape[1],1) # reshape the mask into the 3 color shape as well
    mask = tile(mask,(1,3)).reshape(mask.shape[0],mask.shape[1],3) # and copy the mask across all colors
    intensity = abs(A).reshape(A.shape[0],A.shape[1],1)
    intensity /= abs(A).max()
    if logscale:
        intensity = log10(intensity)
    if black:
        colors = intensity*colors
        colors[mask] = 1.0
    else:
        colors = 1.0-intensity*(1.0-colors)
        colors[mask] = 0.0
    return colors
def myfilter(x,center = 250e3,sigma = 100e3):
    x = (x-center)**2
    x /= sigma**2
    return exp(-x)
#}}}

#{{{ fitdata
class fitdata(nddata):
    def __init__(self,*args,**kwargs):
        #{{{ manual kwargs
        fit_axis = None
        if 'fit_axis' in kwargs.keys():
            fit_axis = kwargs.pop('fit_axis')
        #}}}
        if isinstance(args[0],nddata):
            #print "DEBUG trying to transfer",args[0].axis_coords_error
            myattrs = normal_attrs(args[0])
            for j in range(0,len(myattrs)):
                self.__setattr__(myattrs[j],args[0].__getattribute__(myattrs[j]))
            #nddata.__init__(self,
            #        args[0].data,
            #        args[0].data.shape,
            #        args[0].dimlabels,
            #        axis_coords = args[0].axis_coords,
            #        ft_start_time = args[0].ft_start_time,
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
                raise CustomError("I can't figure out the fit axis!")
        self.fit_axis = fit_axis
        #{{{ in the class, only store the forced values and indeces they are set to
        self.set_to = None
        self.set_indeces = None
        self.active_indeces = None
        #}}}
        return
    def parameter_derivatives(self,xvals,set = None,set_to = None,verbose = False):
        r'return a matrix containing derivatives of the parameters, can set dict set, or keys set, vals set_to'
        if verbose: print 'parameter derivatives is called!'
        if iscomplex(self.data.flatten()[0]):
            print lsafen('Warning, taking only real part of fitting data!')
        if type(set) is dict:
            set_to = set.values()
            set = set.keys()
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
            except:
                raise CustomError('error trying to substitute',mydiff_sym[j],'with',solution_list)
            try:
                fprime[j,:] = array([complex(mydiff.subs(x,xvals[k])) for k in range(0,len(xvals))])
            except ValueError:
                raise CustomError('Trying to set index',j, 'shape(fprime)',shape(fprime), 'shape(xvals)',shape(xvals),'the thing I\'m trying to compute looks like this',[mydiff.subs(x,xvals[k]) for k in range(0,len(xvals))])
            except:
                raise CustomError('Trying to set index',j, 'shape(fprime)',shape(fprime), 'shape(xvals)',shape(xvals))
        return fprime
    def parameter_gradient(self,p,x,y,sigma):
        r'this gives the specific format wanted by leastsq'
        # for now, I'm going to assume that it's not using sigma, though this could be wrong
        # and I could need to scale everything by sigma in the same way as errfunc
        return self.parameter_derivatives(x,set = self.symbol_list,set_to = p).T
    def analytical_covariance(self):
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
            except ValueError:
                raise CustomError('shape of fprime_prod',shape(fprime_prod),'shape of inverse',shape(pinv(fprime_prod)),'shape of sigma',shape(sigma))
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
            print 'a'
            minimizer = inv(J.T * Omegainv * J) * J.T * Omegainv
            covarmatrix = minimizer * S * minimizer.T
            #covarmatrix = array(covarmatrix * S * covarmatrix.T)
            #covarmatrix = array((J.T * G.T * G * J)**-1 * J.T * G.T * G * S * G.T * G * J * (J.T * G.T * G * J)**-1)
            #try:
            #    betapremult = (J.T * Omegainv * J)**-1 * J.T * Omegainv
            #except:
            #    print 'from sigma','\n\n',diag(sigma**2),'\n\n','from covarmatrix','\n\n',S,'\n\n'
            #    raise CustomError('problem generating estimator (word?)')
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
        #        #    raise CustomError('Problem multiplying covarmatrix', 'shape(fprime[j,:])',shape(fprime[j,:]), 'shape(fminusE)',shape(fminusE), 'shape(fdprime)',shape(fdprime))
        #        #if j != k:
        #        #    covarmatrix[j,k] *= 2
        return covarmatrix
    def gen_symbolic(self,function_name):
        r'''generates the symbolic representations the function'''
        self.function_name = function_name
        self.symbolic_vars = map(sympy.var,self.symbol_list)
        self.symbolic_x = sympy.var(self.fit_axis)
        #print lsafen('test symbol_list=',self.symbol_list)
        self.symbolic_dict = dict(zip(self.symbol_list,self.symbolic_vars))
        #print lsafen('test symbolic_vars=',self.symbolic_vars)
        #print lsafen('test symbolic_x=',self.symbolic_x)
        if hasattr(self,'fitfunc_raw_symb'):
            self.symbolic_func = self.fitfunc_raw_symb(self.symbolic_vars,self.symbolic_x)
        else:
            self.symbolic_func = self.fitfunc_raw(self.symbolic_vars,self.symbolic_x)
        self.function_string = sympy.latex(self.symbolic_func).replace('$','')
        self.function_string = r'$' + self.function_name + '=' + self.function_string + r'$'
        return self
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
    def gen_indeces(self,set,set_to):
        r'''pass this set and set\_to parameters, and it will return:
        indeces,values,mask
        indeces --> gives the indeces that are forced
        values --> the values they are forced to
        mask --> p[mask] are actually active in the fit'''
        if type(set) is not list:
            set = [set]
        if type(set_to) is not list:
            set_to = [set_to]
        if len(set) != len(set_to):
            raise CustomError('length of set=',set,'and set_to',set_to,'are not the same!')
        set_indeces = map(self.symbol_list.index,set) # calculate indeces once for efficiency
        active_mask = ones(len(self.symbol_list),dtype = bool)
        active_mask[set_indeces] = False # generate the mask of indeces that are actively fit
        return set_indeces,set_to,active_mask
    def remove_inactive_p(self,p):
        return p[self.active_mask]
    def add_inactive_p(self,p):
        if self.set_indeces != None:
            #{{{ uncollapse the function
            temp = p.copy()
            p = zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            #}}}
            p[self.set_indeces] = self.set_to # then just set the forced values to their given values
        return p
    def fitfunc(self,p,x):
        r"this wraps fitfunc_raw (which gives the actual form of the fit function) to take care of forced variables"
        p = self.add_inactive_p(p)
        return self.fitfunc_raw(p,x)
    def errfunc(self,p,x,y,sigma):
        '''just the error function'''
        fit = self.fitfunc(p,x)
        #normalization = sum(1.0/sigma)
        #print 'DEBUG: y=',y,'\nfit=',fit,'\nsigma=',sigma,'\n\n'
        sigma[sigma == 0.0] = 1
        try:
            retval = (y-fit)/sigma #* normalization
            #print 'DEBUG: retval=',retval,'\n\n'
        except ValueError:
            raise CustomError('your error (',shape(sigma),') probably doesn\'t match y (',shape(y),') and fit (',shape(fit),')')
        return retval
    def pinv(self,*args,**kwargs):
        if 'verbose' in kwargs.keys():
            verbose = kwargs.pop('verbose')
        else:
            verbose = False
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
        if verbose:
            print r'\label{fig:pinv_figure_text}y=',y,'yerr=',yerr,'%s='%x_axis,x,'L=',L
            print '\n\n'
            print 'recalc y = ',dot(L,retval)
            print 'recalc E = ',1.0-1.0/dot(L,retval)
            print 'actual E = ',self.data
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
        r'''give the fit value of a particular symbol'''
        if not hasattr(self,'fit_coeff') or self.fit_coeff is None:
            return None
        p = self.fit_coeff.copy()
        if self.set_indeces != None:
            #{{{ uncollapse the function
            temp = p.copy()
            p = zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            #}}}
            p[self.set_indeces] = self.set_to # then just set the forced values to their given values
            #print "DEBUG trying to uncollapse in fitfunc w/ ",self.symbol_list,"; from",temp,"to",p
        # this should also be generic
        if len(name) == 1:
            try:
                return p[self.symbol_list.index(name[0])]
            except:
                raise CustomError("While running output: couldn't find",name,"in",self.symbol_list)
        elif len(name) == 0:
            # return a record array
            return array(tuple(p),{"names":list(self.symbol_list),"formats":['double']*len(p)}).reshape(1)
        else:
            raise CustomError("You can't pass",len(name),"arguments to .output()")
    def _pn(self,name):
        return self.symbol_list.index(name)
    def _active_symbols(self):
        if not hasattr(self,'active_symbols'):
            if self.set_indeces != None:
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
        if (len(names) == 1) and (names[0] is 'recarray'):
            if hasattr(self,'active_mask'):
                active_symbols = [x for x in self.symbol_list if self.active_mask[self._pn(x)]]
            else:
                active_symbols = list(self.symbol_list)
            if len(active_symbols) != self.covariance.shape[0]:
                raise CustomError('length of active symbols',active_symbols,'doesnt match covariance matrix size(',self.covariance.shape[0],')!')
            recnames = ['labels'] + active_symbols
            recdata = []
            for j in range(0,self.covariance.shape[0]): 
                thisdata = [active_symbols[j]] + list(double(self.covariance[j,:].copy())) # the first index is the row
                recdata.append(make_rec(thisdata,recnames))
            return r_[tuple(recdata)]
        if len(names) > 0:
            indeces = map(self._pn_active,names) # slice out only these rows and columns
            return self.covariance[r_[indeces],:][:,r_[indeces]].copy()
        else:
            try:
                return self.covariance.copy()
            except:
                return zeros([len(self.fit_coeff)]*2,dtype = 'double')
    def latex(self,verbose = False):
        r'''show the latex string for the function, with all the symbols substituted by their values'''
        # this should actually be generic to fitdata
        p = self.fit_coeff
        printfstring = self.function_string
        printfargs = []
        allsymb = []
        locations = []
        for j in range(0,len(self.symbol_list)):
            #symbol = self.symbol_list[j]
            symbol = sympy.latex(self.symbolic_vars[j]).replace('$','')
            if verbose: print 'DEBUG: replacing symbol \\verb|',symbol,'|'
            location = printfstring.find(symbol)
            while location != -1:
                if printfstring[location-1] == '-':
                    newstring = printfstring[:location-1]+'+%01.03g'+printfstring[location+len(symbol):] # replace the symbol in the written function with the appropriate number
                    thissign = -1.0
                else:
                    newstring = printfstring[:location]+'%01.03g'+printfstring[location+len(symbol):] # replace the symbol in the written function with the appropriate number
                    thissign = 1.0
                if verbose: print r"\begin{verbatim} trying to replace",printfstring[location:location+len(symbol)],r'\end{verbatim}'
                printfstring = newstring
                printfargs += [thissign*p[j]] # add that number to the printf list
                locations += [location]
                allsymb += [symbol]
                location = printfstring.find(symbol)
        printfargs = [printfargs[x] for x in argsort(locations)]
        if verbose: print r"\begin{verbatim}trying to generate",self.function_string,'\n',printfstring,'\n',[allsymb[x] for x in argsort(locations)],'\n',printfargs,r'\end{verbatim}'
        return printfstring%tuple(printfargs)
    def settoguess(self):
        'a debugging function, to easily plot the initial guess'
        self.fit_coeff = real(self.guess())
        return self
    def _taxis(self,taxis):
        r'You can enter None, to get the fit along the same range as the data, an integer to give the number of points, or a range of data, which will return with 300 points'
        if taxis is None:
            taxis = self.getaxis(self.fit_axis).copy()
        elif type(taxis) is int:
            taxis = linspace(self.getaxis(self.fit_axis).min(),
                    self.getaxis(self.fit_axis).max(),
                    taxis)
        elif not isscalar(taxis) and len(taxis) == 2:
            taxis = linspace(taxis[0],taxis[1],300)
        return taxis
    def eval(self,taxis,set = None,set_to = None):
        r'''after we have fit, evaluate the fit function along the axis taxis
        set and set_to allow you to forcibly set a specific symbol to a specific value --> however, this does not affect the class, but only the return value'''
        if type(set) is dict:
            set_to = set.values()
            set = set.keys()
        taxis = self._taxis(taxis)
        if hasattr(self,'fit_coeff') and self.fit_coeff is not None:
            p = self.fit_coeff.copy()
        else:
            p = array([NaN]*len(self.symbol_list))
        #{{{ LOCALLY apply any forced values
        if set != None:
            if self.set_indeces != None:
                raise CustomError("your'e trying to set indeces in an eval function for a function that was fit constrained; this is not currently supported")
            set_indeces,set_to,active_mask = self.gen_indeces(set,set_to)
            p[set_indeces] = set_to
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
    def fit(self,set = None, set_to = None, force_analytical = False):
        r'''actually run the fit'''
        if type(set) is dict:
            set_to = set.values()
            set = set.keys()
        x = self.getaxis(self.fit_axis)
        if iscomplex(self.data.flatten()[0]):
            print lsafen('Warning, taking only real part of fitting data!')
        y = real(self.data)
        sigma = self.get_error()
        if sigma is None:
            print '{\\bf Warning:} You have no error associated with your plot, and I want to flag this for now\n\n'
            warnings.warn('You have no error associated with your plot, and I want to flag this for now',Warning)
            sigma = ones(shape(y))
        p_ini = real(array(self.guess())) # need the numpy format to allow boolean mask
        if set != None:
            self.set_indeces,self.set_to,self.active_mask = self.gen_indeces(set,set_to)
            p_ini = self.remove_inactive_p(p_ini)
        leastsq_args = (self.errfunc, p_ini)
        leastsq_kwargs = {'args':(x,y,sigma),
                    'full_output':True}# 'maxfev':1000*(len(p_ini)+1)}
        if hasattr(self,'has_grad') and self.has_grad == True:
            leastsq_kwargs.update({'Dfun':self.parameter_gradient})
        try:
            p_out,cov,infodict,mesg,success = leastsq(*leastsq_args,**leastsq_kwargs)
        #{{{ just give various explicit errors
        except TypeError,err:
            if type(x) != ndarray and type(y) != ndarray:
                raise CustomError('leastsq failed because the two arrays aren\'t of the right type','type(x):',type(x),'type(y):',type(y))
            else:
                if any(shape(x) != shape(y)):
                    raise CustomError('leastsq failed because the two arrays do not match in size size','shape(x):',shape(x),'shape(y):',shape(y))
            raise CustomError('leastsq failed because of a type error!','type(x):',showtype(x),'type(y):',showtype(y),'type(sigma)',showtype(sigma),'shape(x):',shape(x),'shape(y):',shape(y),'shape(sigma)',shape(sigma),'p_ini',type(p_ini),p_ini)
        except ValueError,err:
            raise CustomError('leastsq failed with "',err,'", maybe there is something wrong with the input:',self)
        except:
            raise CustomError('leastsq failed; I don\'t know why')
        #}}}
        if success not in [1,2,3,4]:
            #{{{ up maximum number of evals
            if mesg.find('maxfev'):
                leastsq_kwargs.update({ 'maxfev':50000 })
                p_out,cov,infodict,mesg,success = leastsq(*leastsq_args,**leastsq_kwargs)
                if success != 1:
                    if mesg.find('two consecutive iterates'):
                        print r'{\Large\color{red}{\bf Warning data is not fit!!! output shown for debug purposes only!}}','\n\n'
                        print r'{\color{red}{\bf Original message:}',lsafe(mesg),'}','\n\n'
                        infodict_keys = infodict.keys()
                        infodict_vals = infodict.values()
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
                            print r'{\color{red}{\bf %s:}%s}'%(k,v),'\n\n'
                        #self.fit_coeff = None
                        #self.settoguess()
                        #return
                    else:
                        raise CustomError('leastsq finished with an error message:',mesg)
                    #}}}
            else:
                raise CustomError('leastsq finished with an error message:',mesg)
        else:
            print r'{\color{blue}'
            print lsafen("Fit finished successfully with a code of %d and a message ``%s''"%(success,mesg))
            print r'}'
        self.fit_coeff = p_out # note that this is stored in HIDDEN form
        dof = len(x) - len(p_out)
        if hasattr(self,'symbolic_x') and force_analytical:
            self.covariance = self.analytical_covariance()
        else:
            if force_analytical: raise CustomError("I can't take the analytical covariance!  This is problematic.")
            if cov == None:
                #raise CustomError('cov is none! why?!, x=',x,'y=',y,'sigma=',sigma,'p_out=',p_out,'success=',success,'output:',p_out,cov,infodict,mesg,success)
                print r'{\color{red}'+lsafen('cov is none! why?!, x=',x,'y=',y,'sigma=',sigma,'p_out=',p_out,'success=',success,'output:',p_out,cov,infodict,mesg,success),'}\n'
            self.covariance = cov
        if self.covariance is not None:
            try:
                self.covariance *= sum(infodict["fvec"]**2)/dof # scale by chi_v "RMS of residuals"
            except:
                raise CustomError("type(self.covariance)",type(self.covariance),
                        "type(infodict[fvec])",type(infodict["fvec"]),
                        "type(dof)",type(dof))
        #print lsafen("DEBUG: at end of fit covariance is shape",shape(self.covariance),"fit coeff shape",shape(self.fit_coeff))
        return
    def bootstrap(self,points,swap_out = exp(-1.0),seedval = 10347,minbounds = {},maxbounds = {}):
        print r'\begin{verbatim}'
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
                random_indeces = int32((rand(number_to_replace)*(thiscopy.data.size-1.0)).round())
                thiscopy.data = r_[thiscopy.data,thiscopy.data.copy()[random_indeces]]
                thiscopy.labels([thiscopy.fit_axis],[r_[x,x.copy()[random_indeces]]])
                thiscopy.set_error(r_[derr,derr.copy()[random_indeces]])
                #print 'DEBUG: size of data after extension',double(size(thiscopy.data))/origsizecheck
                #}}}
                try:
                    thiscopy.fit()
                    success = True
                    if len(minbounds) > 0:
                        for k,v in minbounds.iteritems():
                            if thiscopy.output(k) < v:
                                success = False
                    if len(maxbounds) > 0:
                        for k,v in maxbounds.iteritems():
                            if thiscopy.output(k) > v:
                                success = False
                except:
                    #print 'WARNING, didn\'t fit'
                    success = False
                # here, use the internal routines, in case there are constraints, etc
                if success is True:
                    for name in thiscopy.symbol_list: # loop over all fit coeff
                        recordlist[runno][name] = thiscopy.output(name)
        print r'\end{verbatim}'
        return recordlist # collect into a single recordlist array
    def guess(self,verbose = False,super_verbose = False):
        r'''provide the guess for our parameters; by default, based on pseudoinverse'''
        self.has_grad = False
        if iscomplex(self.data.flatten()[0]):
            print lsafen('Warning, taking only real part of fitting data!')
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
        guess_dict = dict(zip(self.symbol_list,list(thisguess)))
        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
        f_at_guess = real(self.eval(None,set = guess_dict).data)
        try:
            f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
        except:
            raise CustomError('trying to reshape f_at_ini_guess from',f_at_guess.shape,'to',new_y_shape)
        thisresidual = sqrt((y-f_at_guess)**2).sum()
        #}}}
        lastresidual = thisresidual
        for j in range(0,numguesssteps):
            if super_verbose: print '\n\n(matlablike.guess) '+r'\begin{verbatim} fprime = \n',fprime,'\nf_at_guess\n',f_at_guess,'y=\n',y,'\n',r'\end{verbatim}'
            if super_verbose: print '\n\n(matlablike.guess) shape of parameter derivatives',shape(fprime),'shape of output',shape(y),'\n\n'
            regularization_bad = True
            alpha_max = 100.
            alpha_mult = 2.
            alpha = 0.1 # maybe I can rather estimate this based on the change in the residual, similar to in L-M?
            if verbose: print '\n\n(matlablike.guess) value of residual before regularization %d:'%j,thisresidual
            while regularization_bad:
                newguess = real(array(thisguess) + dot(pinvr(fprime.T,alpha),(y-f_at_guess)).flatten())
                mask = newguess < self.guess_lb
                newguess[mask] = self.guess_lb[mask]
                mask = newguess > self.guess_ub
                newguess[mask] = self.guess_ub[mask]
                if any(isnan(newguess)):
                    if verbose: print '\n\n(matlablike.guess) Regularization blows up $\\rightarrow$ increasing $\\alpha$ to %0.1f\n\n'%alpha
                    alpha *= alpha_mult
                else:
                    #{{{ evaluate f, fprime and residuals
                    guess_dict = dict(zip(self.symbol_list,list(newguess)))
                    # only evaluate fprime once we know this is good, below
                    f_at_guess = real(self.eval(None,set = guess_dict).data)
                    try:
                        f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
                    except:
                        raise CustomError('trying to reshape f_at_ini_guess from',f_at_guess.shape,'to',new_y_shape)
                    thisresidual = sqrt((y-f_at_guess)**2).sum()
                    #}}}
                    if (thisresidual-lastresidual)/lastresidual > 0.10:
                        alpha *= alpha_mult
                        if verbose: print '\n\n(matlablike.guess) Regularized Pinv gave a step uphill $\\rightarrow$ increasing $\\alpha$ to %0.1f\n\n'%alpha
                    else: # accept the step
                        regularization_bad = False
                        thisguess = newguess
                        lastresidual = thisresidual
                        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
                if alpha > alpha_max:
                    print "\n\n(matlablike.guess) I can't find a new guess without increasing the alpha beyond %d\n\n"%alpha_max
                    if which_starting_guess >= len(self.starting_guesses)-1:
                        print "\n\n(matlablike.guess) {\\color{red} Warning!!!} ran out of guesses!!!%d\n\n"%alpha_max
                        return thisguess
                    else:
                        which_starting_guess += 1
                        thisguess = self.starting_guesses[which_starting_guess]
                        print "\n\n(matlablike.guess) try a new starting guess:",lsafen(thisguess)
                        j = 0 # restart the loop
                        #{{{ evaluate f, fprime and residuals for the new starting guess
                        guess_dict = dict(zip(self.symbol_list,list(thisguess)))
                        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
                        f_at_guess = real(self.eval(None,set = guess_dict).data)
                        try:
                            f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
                        except:
                            raise CustomError('trying to reshape f_at_ini_guess from',f_at_guess.shape,'to',new_y_shape)
                        thisresidual = sqrt((y-f_at_guess)**2).sum()
                        #}}}
                        regularization_bad = False # jump out of this loop
            if verbose: print '\n\n(matlablike.guess) new value of guess after regularization:',lsafen(newguess)
            if verbose: print '\n\n(matlablike.guess) value of residual after regularization:',thisresidual
        return thisguess
#}}}
def sqrt(arg):
    if isinstance(arg,nddata):
        return arg**0.5
    else:
        return np_sqrt(arg)
