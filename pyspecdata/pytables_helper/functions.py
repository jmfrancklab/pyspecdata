from ..general_functions import *
from pylab import * 
logger = logging.getLogger('pyspecdata.load_files.bruker_esr')
import tables
from os import listdir,environ

#{{{ HDF5 functions
#{{{ helper function for HDF5 search
def gensearch(labelname,format = '%0.3f',value = None,precision = None):
    'obsolete -- use h5gensearch'
    if value == None:
        raise ValueError('You must pass a value to gensearch')
    if precision == None:
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
    elif len(args) == 1 and type(args[0]) is dict:
        dict_arg = args[0]
        condlist = []
        for k,v in dict_arg.iteritems():
            condlist.append(h5searchstring(k,v,format = format,precision = precision))
        return ' & '.join(condlist)
    else:
        raise ValueError("pass either field,value pair or a dictionary!")
    if type(value) is str:
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
def h5loaddict(thisnode,verbose = False):
    "Load nddata structure from hdf5 file using pytables and return the result as a dictionary"
    #{{{ load all attributes of the node
    retval = dict([(x,thisnode._v_attrs.__getattribute__(x))
        for x in thisnode._v_attrs._f_list('user')])
    #}}}
    for k,v in retval.iteritems():#{{{ search for record arrays that represent normal lists
        retval[k]  = unmake_ndarray(v,name_forprint = k,verbose = verbose)
    if type(thisnode) is tables.table.Table:#{{{ load any table data
        if verbose: print "It's a table\n\n"
        if 'data' in retval.keys():
            raise AttributeError('There\'s an attribute called data --> this should not happen!')
        retval.update({'data':thisnode.read()})
    elif type(thisnode) is tables.group.Group:
        #{{{ load any sub-nodes as dictionaries
        mychildren = thisnode._v_children
        for thischild in mychildren.keys():
            if thischild in retval.keys():
                raise AttributeError('There\'s an attribute called ',thischild,' and also a sub-node called the',thischild,'--> this should not happen!')
            retval.update({thischild:h5loaddict(mychildren[thischild])})
        #}}}
    else:
        raise AttributeError(strm("I don't know what to do with this node:",thisnode))
    #}}}
    return retval
def h5child(thisnode,childname,clear = False,create = None,verbose = False):
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
        childnode = h5file.getNode(thisnode,childname)
        if verbose:
            print lsafe('found',childname)
        if clear:
            childnode._f_remove(recursive = True)
            childnode = None
    except tables.NoSuchNodeError as e:
        if create is False and not clear:
            raise RuntimeError('Trying to grab a node that does not exist with create = False'+explain_error(e))
        elif clear:
            childnode = None
        else:
            childnode = h5file.createGroup(thisnode,childname)
            if verbose:
                print lsafe('created',childname)
    return childnode
def h5remrows(bottomnode,tablename,searchstring):
    if type(searchstring) is dict:
        searchstring = h5searchstring(searchstring)
    try:
        thistable = bottomnode.__getattr__(tablename)
        counter = 0
        try:
            data = thistable.readWhere(searchstring).copy()
        except Exception as e:
            raise RuntimeError(strm('Problem trying to remove rows using search string',
                searchstring, 'in', thistable, explain_error(e)))
        for row in thistable.where(searchstring):
            if len(thistable) == 1:
                thistable.remove()
                counter += 1
            else:
                try:
                    thistable.removeRows(row.nrow - counter,row.nrow - counter + 1) # counter accounts for rows I have already removed.
                except:
                    print "you passed searchstring",searchstring
                    print "trying to remove row",row
                    print "trying to remove row with number",row.nrow
                    print help(thistable.removeRows)
                    raise RuntimeError("length of thistable is "+repr(len(thistable))+" calling removeRows with "+repr(row.nrow-counter))
                counter += 1
        return counter,data
    except tables.NoSuchNodeError:
        return False,None
def h5addrow(bottomnode,tablename,*args,**kwargs):
    '''add a row to a table, creating it if necessary, but don\'t add if the data matches the search condition indicated by `match_row`
    `match_row` can be either text or a dictionary -- in the latter case it's passed to h5searchstring
    '''
    match_row,verbose,only_last = process_kwargs([('match_row',None),('verbose',False),('only_last',True)],kwargs)
    try: # see if the table exists
        mytable = h5table(bottomnode,tablename,None)
        tableexists = True
    except RuntimeError: # if table doesn't exist, create it
        newindex = 1L
        tableexists = False
    if tableexists:
        #{{{ auto-increment "index"
        newindex = mytable.read()['index'].max() + 1L
        #}}}
        # here is where I would search for the existing data
        if match_row is not None:
            if type(match_row) is dict:
                match_row = h5searchstring(match_row)
            if verbose: obs("trying to match row according to",lsafen(match_row))
            mytable.flush()
            try:
                matches = mytable.readWhere(match_row)
            except NameError as e:
                raise NameError(' '.join(map(str,
                    [e,'\nYou passed',match_row,'\nThe columns available are',mytable.colnames,"condvars are",condvars])))
            except ValueError as e:
                raise NameError(' '.join(map(str,
                    [e,'\nYou passed',match_row,'\nThe columns available are',mytable.colnames])))
            if len(matches) > 0:
                if only_last:
                    if verbose: print r'\o{',lsafen(len(matches),"rows match your search criterion, returning the last row"),'}'
                    return mytable,matches['index'][-1]
                else:
                    return mytable,matches['index'][:]
            else:
                if verbose:
                    obs("I found no matches")
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
        except ValueError as e:
            print lsafen("I'm about to flag an error, but it looks like there was an issue appending",myrowdata)
            tabledforerr = mytable.read()
            raise AttributeError(strm('Compare names and values table data vs. the row you are trying to add\n',
                '\n'.join(map(repr,zip(mytable.read().dtype.fields.keys(),
                mytable.read().dtype.fields.values(),
                myrowdata.dtype.fields.keys(),
                myrowdata.dtype.fields.values())))),explain_error(e))
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
    if tablename not in bottomnode._v_children.keys():
        if tabledata is not None:
            if type(tabledata) is dict:
                tabledata = make_rec(tabledata)
            datatable = h5file.createTable(bottomnode,tablename,tabledata) # actually write the data to the table
        else:
            raise RuntimeError(' '.join(map(str,['You passed no data, so I can\'t create table',tablename,'but it doesn\'t exist in',bottomnode,'which has children',bottomnode._v_children.keys()])))
    else:
        if tabledata is not None:
            raise ValueError('You\'re passing data to create the table, but the table already exists!')
        else:
            pass
    return bottomnode._v_children[tablename]
    #}}}
def h5nodebypath(h5path,verbose = False,force = False,only_lowest = False,check_only = False,directory='.'):
    r'''return the node based on an absolute path, including the filename'''
    logger.debug(strm("DEBUG: called h5nodebypath on",h5path))
    h5path = h5path.split('/')
    #{{{ open the file / check if it exists
    if verbose: print lsafen('h5path=',h5path)
    logger.info(strm('the h5path is',h5path))
    try:
        if h5path[0] in listdir(directory):
            if verbose: print 'DEBUG: file exists\n\n'
        else:
            if check_only: raise AttributeError("You're checking for a node in a file (%s) that does not exist"%h5path[0])
            if verbose: print 'DEBUG: file does not exist\n\n'
        mode = 'a'
        #if check_only: mode = 'r'
        logger.info(strm('so I look for the file',h5path[0],'in directory',directory))
        h5file = tables.openFile(os.path.join(directory,h5path[0]),mode = mode,title = 'test file')
    except IOError as e:
        raise IOError('I think the HDF5 file has not been created yet, and there is a bug pytables that makes it freak out, but you can just run again.'+explain_error(e))
    #}}}
    currentnode = h5file.getNode('/') # open the root node
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
                        verbose = verbose,
                        create = create,
                        clear = clear)
                if verbose: print lsafen("searching for node path: descended to node",currentnode)
                logger.info(strm("searching for node path: descended to node",currentnode))
            except BaseException as e:
                if verbose: print lsafen("searching for node path: got caught searching for node",h5path[pathlevel])
                logger.info(strm("searching for node path: got caught searching for node",h5path[pathlevel]))
                h5file.close()
                #print lsafen("DEBUG: Yes, I closed the file")
                raise IndexError(strm('Problem trying to load node ',h5path,explain_error(e)))
            #}}}
    return h5file,currentnode
def h5attachattributes(node,listofattributes,myvalues,exclusions):
    listofattributes = [j for j in listofattributes # need to exclude the properties
            if j not in exclusions]
    if node is None:
        raise IndexError('Problem!, node passed to h5attachattributes: ',node,'is None!')
    h5file = node._v_file
    if type(myvalues).__name__ == 'nddata':
        attributevalues = map(lambda x: myvalues.__getattribute__(x),listofattributes)
    elif type(myvalues) is list:
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
                    thisval.keys(),
                    thisval.values(),
                    exclusions)
            thisval = None
            listout.remove(thisattr)
        else:
            # {{{ pytables hates <U24 which is created from unicode
            if type(thisval) in [list,tuple]:
                if any(map(lambda x: isinstance(x,unicode),thisval)):
                    logger.info(strm("going to convert",thisval,"to strings"))
                    thisval = [str(x) if isinstance(x,unicode) else x for x in thisval]
                    logger.info(strm("now it looks like this:",thisval))
            thisval = make_ndarray(thisval,name_forprint = thisattr)
            # }}}
        if thisval is not None:
            try:
                node._v_attrs.__setattr__(thisattr,thisval)
            except Exception,e:
                raise RuntimeError("PyTables freaks out when trying to attach attribute "+repr(thisattr)+" with value "+repr(thisval)+"\nOriginal error was:\n"+str(e))
            listout.remove(thisattr)
    listofattributes[:] = listout # pointer
def h5inlist(columnname,mylist):
    'returns rows where the column named columnname is in the value of mylist'
    if type(mylist) is slice:
        if mylist.start is not None and mylist.stop is not None:
            return "(%s >= %g) & (%s < %g)"%(columnname,mylist.start,columnname,mylist.stop)
        elif mylist.stop is not None:
            return "(%s < %g)"%(columnname,mylist.stop)
        elif mylist.start is not None:
            return "(%s > %g)"%(columnname,mylist.start)
        else:
            raise ValueError()
    if type(mylist) is ndarray:
        mylist = mylist.tolist()
    if type(mylist) is not list:
        raise TypeError("the second argument to h5inlist must be a list!!!")
    if any([type(x) in [double,float64] for x in mylist]):
        if all([type(x) in [double,float64,int,int32,int64] for x in mylist]):
            return '('+'|'.join(map(lambda x: "(%s == %g)"%(columnname,x),mylist))+')'
    elif all([type(x) in [int,long,int32,int64] for x in mylist]):
        return '('+'|'.join(map(lambda x: "(%s == %g)"%(columnname,x),mylist))+')'
    elif all([type(x) is str for x in mylist]):
        return '('+'|'.join(map(lambda x: "(%s == '%s')"%(columnname,x),mylist))+')'
    else:
        raise TypeError("I can't figure out what to do with this list --> I know what to do with a list of numbers or a list of strings, but not a list of type"+repr(map(type,mylist)))
def h5join(firsttuple,secondtuple,
    additional_search = '',
    select_fields = None,
    pop_fields = None,
    verbose = False):
    #{{{ process the first argument as the hdf5 table and indices, and process the second one as the structured array to join onto
    if not ((type(firsttuple) is tuple) and (type(secondtuple) is tuple)):
        raise ValueError('both the first and second arguments must be tuples!')
    if not ((len(firsttuple) == 2) and (len(secondtuple) == 2)):
        raise ValueError('The length of the first and second arguments must be two!')
    tablenode = firsttuple[0]
    tableindices = firsttuple[1]
    if verbose: print 'h5join tableindices looks like this:',tableindices
    if type(tableindices) is not list:
        tableindices = [tableindices]
    if verbose: print 'h5join tableindices looks like this:',tableindices
    mystructarray = secondtuple[0].copy()
    mystructarrayindices = secondtuple[1]
    if type(mystructarrayindices) is not list:
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
        if isinstance(mystructarray[thisstructarrayindex][0],basestring):
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
        if verbose: print '\n\nh5join original indices',lsafen(retval.dtype.names)
        try:
            retval = retval[select_fields]
        except ValueError as e:
            raise ValueError(strm('One of the fields', select_fields, 'is not in',
                retval.dtype.names, explain_error(e)))
    #}}}
    return retval
#}}}
