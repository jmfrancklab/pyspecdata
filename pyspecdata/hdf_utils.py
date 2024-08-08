import tables
import numpy as np
from os import listdir
import os
import logging
from .general_functions import strm, lsafen
from .datadir import log_fname, unknown_exp_type_name
from .dict_utils import make_ndarray, unmake_ndarray
logger = logging.getLogger('pyspecdata.hdf_utils')

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
            top_level_node_names = [node._v_name for node in h5file.list_nodes('/')]
            raise RuntimeError('Trying to grab a node that does not exist with create = False.\nHere are the nodes available:'+str(top_level_node_names))
        elif clear:
            childnode = None
        else:
            childnode = h5file.create_group(thisnode,childname)
            logger.debug(strm('created',childname))
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
                searchstring, 'in', thistable))
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
            logger.debug(strm("trying to match row according to",lsafen(match_row)))
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
                list(myrowdata.dtype.fields.values())))))))
        mytable.flush()
    else:
        recorddata = myrowdata
        try:
            mytable = h5table(bottomnode,
                    tablename,
                    recorddata)
        except Exception as e:
            raise RuntimeError(strm('Error trying to write record np.array:',
                repr(recorddata),'from listofdata',listofdata,'and names',listofnames
                ))
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
        log_fname('data_files',h5path[0],directory,unknown_exp_type_name)
    else:
        if check_only:
            errmsg = log_fname('missing_data_files',h5path[0],directory,unknown_exp_type_name)
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
    toplevel_children = ' '.join(currentnode._v_children)
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
                raise IndexError(strm('Problem trying to load node ',h5path,' these are the top level nodes ',toplevel_children))
            #}}}
    return h5file,currentnode
def h5attachattributes(node,listofattributes,myvalues):
    listofattributes = [j for j in listofattributes # need to exclude the properties
            if j not in ['angle','real','imag']]
    if node is None:
        raise IndexError('Problem!, node passed to h5attachattributes: ',node,'is None!')
    h5file = node._v_file
    if hasattr(myvalues,'dimlabels'):# use as a proxy for being nddata
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
                if np.any([isinstance(x,str) for x in thisval]):
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
    if isinstance(mylist, np.ndarray):
        mylist = mylist.tolist()
    if not isinstance(mylist, list):
        raise TypeError("the second argument to h5inlist must be a list!!!")
    if np.any([type(x) in [np.double,np.float64] for x in mylist]):
        if all([type(x) in [np.double,np.float64,int,np.int32,np.int64] for x in mylist]):
            return '('+'|'.join(["(%s == %g)"%(columnname,x) for x in mylist])+')'
    elif all([type(x) in [int,int,np.int32,np.int64] for x in mylist]):
        return '('+'|'.join(["(%s == %g)"%(columnname,x) for x in mylist])+')'
    elif all([isinstance(x, str) for x in mylist]):
        return '('+'|'.join(["(%s == '%s')"%(columnname,x) for x in mylist])+')'
    else:
        raise TypeError("I can't figure out what to do with this list --> I know what to do with a list of numbers or a list of strings, but not a list of type"+repr(list(map(type,mylist))))
def h5join(firsttuple,secondtuple,
    additional_search = '',
    select_fields = None,
    pop_fields = None):
    #{{{ process the first argument as the hdf5 table and indices, and process the second one as the structured np.array to join onto
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
        elif type(mystructarray[thisstructarrayindex][0]) in [int,np.double,float,np.float64,np.float32,np.int32,np.int64]:
            search_string.append(["(%s == %s)"%(thistableindex,str(x)) for x in mystructarray[thisstructarrayindex]])
            #print 'a g mapping for',[x for x in mystructarray[thisstructarrayindex]],'gives',search_string[-1],'\n\n'
        else:
            raise TypeError("I don't know what to do with a structured np.array that has a row of type"+repr(type(mystructarray[thisstructarrayindex][0])))
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
            retval, 'of dtype', retval.dtype, 'with the structured np.array',
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
