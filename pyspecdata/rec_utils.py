"utility functions to handle record/structured numpy arrays"

import numpy as np
import logging
import numpy.lib.recfunctions as recf
from .general_functions import (
    strm,
    process_kwargs,
    lsafen,
    autostringconvert,
)


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
    logging.debug(
        strm("(decorate\\_rec):: original list of matching", list_of_matching)
    )
    length_of_matching = np.array([len(j) for j in list_of_matching])
    logging.debug(
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
    logging.debug(
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
    logging.debug(strm("(decorate\\_rec):: new dtypes:", repr(new_dtypes)))
    try:
        retval = newcol_rec(retval, new_dtypes)
    except Exception:
        raise ValueError(
            strm(
                "Problem trying to add new columns with the dtypes", new_dtypes
            )
        )
    # }}}
    logging.debug(strm("(decorate\\_rec):: add data:", repr(add_data)))
    for name in np.dtype(new_dtypes).names:
        logging.debug(
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
        logging.debug(
            strm(
                lsafen(
                    "(applyto np.core.rec): for row %d, I select these:" % j
                )
            )
        )
        myarray_subset = myarray[mask]
        logging.debug(
            strm(lsafen("(applyto np.core.rec): ", repr(myarray_subset)))
        )
        other_fields = set(mylist) ^ set(thisitem.dtype.names)
        logging.debug(
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
        logging.debug(
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
        logging.debug(
            strm(
                lsafen(
                    "(applyto np.core.rec): the np.array is now", repr(myarray)
                )
            )
        )
        j += 1
    # }}}
    combined = np.concatenate(combined)
    logging.debug(
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
    logging.debug(strm("(meanstd_rec): other fields are", lsafen(other_fields)))
    newrow_dtype = [
        [j, ("%s_ERROR" % j[0],) + j[1:]] if j[0] in other_fields else [j]
        for j in myarray.dtype.descr
    ]
    newrow_dtype = [k for j in newrow_dtype for k in j]
    logging.debug(
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
        logging.debug(
            strm(
                lsafen(
                    "(meanstd np.core.rec): for row %d, I select these:" % j
                )
            )
        )
        myarray_subset = myarray[mask]
        logging.debug(
            strm(lsafen("(meanstd np.core.rec): ", repr(myarray_subset)))
        )
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = np.mean(myarray_subset[thisfield])
                if standard_error:
                    newrow[thisfield + "_ERROR"] = np.std(
                        myarray_subset[thisfield]
                    ) / np.sqrt(len(myarray_subset[thisfield]))
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
        logging.debug(
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
        logging.debug(
            strm(
                lsafen(
                    "(meanstd np.core.rec): the np.array is now", repr(myarray)
                )
            )
        )
        j += 1
    # }}}
    combined = np.concatenate(combined)
    logging.debug(
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
