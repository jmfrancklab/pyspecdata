import tables
import numpy as np

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
            raise RuntimeError('Trying to grab a node that does not exist with create = False')
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
#{{{ plot wrapper
global OLDplot
OLDplot = plot
global myplotfunc
myplotfunc = OLDplot
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
    ax,human_units,label_format_string,normalize,noerr,longest_is_x = process_kwargs([('ax',plt.gca()),
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
    if np.isscalar(myx):
        myx = np.array([myx])
    if np.isscalar(myy):
        myy = np.array([myy])
    #}}}
    x_inverted = False
    #{{{ parse nddata
    if isinstance(myy,nddata):
        myy = myy.copy()
        if myy.get_error() is not None:
            logging.debug(strm("shapes at top of function",ndshape(myy), myy.data.shape, myy.data_error.shape))
        # {{{ automatically reduce any singleton dimensions
        if not len(myy.dimlabels) == 1:
            if np.any(np.array(myy.data.shape) == 1):
                for singleton_dim in [lb for j,lb in enumerate(myy.dimlabels) if myy.data.shape[j] == 1]:
                    myy = myy[singleton_dim,0]
        # }}}
        if len(myy.data.shape)>1 and longest_is_x:
            longest_dim = np.argmax(myy.data.shape)
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
        if not noerr and isinstance(myy.data_error, np.ndarray) and len(myy.data_error)>0: #then this should be an errorbar plot
            def thiserrbarplot(*tebargs,**tebkwargs):
                if 'capsize' not in tebkwargs:
                    tebkwargs.update({'capsize':6})
                if isinstance(tebargs[-1], str):
                    tebkwargs.update({'fmt':tebargs[-1]})
                    return ax.errorbar(*tebargs[:-1],**tebkwargs)
                else:
                    return ax.errorbar(*tebargs,**tebkwargs)
            myplotfunc = thiserrbarplot
            logger.debug("this is an errorbar plot")
            #{{{ pop any singleton dims
            myyerror = myy.get_error()
            myyerror = np.squeeze(myyerror)
            #}}}
            kwargs.update({'yerr':None})
            valueforxerr = myy.get_error(myy.dimlabels[longest_dim])
            if valueforxerr is not None: # if we have x errorbars too
                #print "DEBUG decided to assign to xerr:",valueforxerr
                kwargs.update({'xerr':valueforxerr})
            logging.debug(strm("shapes after splitting nddata",myy.data.shape, myyerror.shape))
        #{{{ deal with axis labels along y
        try:
            yaxislabels = myy.getaxis(myy.dimlabels[last_not_longest])
        except:
            yaxislabels = None
        # at this point, if there is no axis label, it will break and go to pass
        if yaxislabels is not None:
            if len(yaxislabels) > 0:
                if isinstance(yaxislabels[0], np.string_):
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
            myy = np.squeeze(myy.data.transpose([longest_dim]+all_but_longest))
            if 'yerr' in kwargs.keys():
                myyerror = np.squeeze(myyerror.transpose([longest_dim]+all_but_longest))
        if 'yerr' in kwargs.keys():
            logging.debug(strm("halfway checkpoint",myy.shape, myyerror.shape))
        #if len(myy.data) == 1 and 'label' not in kwargs.keys() and myy_name is not None:
        #    kwargs.update('label',myy_name)
        # }}}
    #}}}
    # {{{ allow list arguments
    if type(myy) is list:
        myy = np.array(myy)
    if type(myx) is list:
        myx = np.array(myx)
    # }}}
    #{{{ semilog where appropriate
    if (myx is not None) and (len(myx)>1) and all(myx>0.0): # by doing this and making myplotfunc global, we preserve the plot style if we want to tack on one point
        try:
            b = np.diff(np.log10(myx))
        except Exception as e:
            raise Exception(strm('likely a problem with the type of the x label, which is',myx))
        if (np.size(b)>3) and all(abs((b-b[0])/b[0])<1e-4) and not ('nosemilog' in list(kwargs.keys())):
            if 'plottype' not in list(kwargs.keys()):
                myplotfunc = ax.semilogx
    if ('nosemilog' in list(kwargs.keys())):
        #print 'this should pop nosemilog'
        kwargs.pop('nosemilog')
    if 'yerr' in kwargs.keys():
        logging.debug(strm("halfway checkpoint",myy.data.shape, myyerror.shape))
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
    if len(np.shape(myy.squeeze()))>1 and np.sum(np.array(np.shape(myy))>1):
        #{{{ hsv plots
        retval = []
        if 'yerr' in kwargs.keys():
            logger.debug("I see error")
        else:
            logger.debug("I do not see error")
        for j in range(0,myy.shape[1]):
            #{{{ this is the way to assign plot arguments
            plotargs = [k for k in (myx,myy[:,j],myformat) if k is not None]
            #}}}
            if 'yerr' in kwargs.keys():
                kwargs['yerr'] = myyerror[:,j]
            #{{{ here, i update the kwargs to include the specific color for this line
            newkwargs = kwargs.copy() # kwargs is a dict
            newkwargs.update({'color':cm.hsv(np.double(j)/np.double(myy.shape[1]))})
            #}}}
            #{{{ here, I update to use the labels
            if has_labels:
                newkwargs.update({'label':yaxislabels[j]})
            #}}}
            if 'yerr' in newkwargs.keys():
                logging.debug(strm("shapes before plot",plotargs[0].shape,
                    plotargs[1].shape, newkwargs['yerr'].shape))
                #myplotfunc = ax.plot
                #newkwargs.pop('yerr')
            elif len(plotargs)>1 and isinstance(plotargs[1],np.ndarray):
                logging.debug(strm("shapes before plot",plotargs[0].shape,
                    plotargs[1].shape))
            if np.any(np.isinf(myy)):
                myy[np.isinf(myy)] = NaN # added this to prevent an overflow error
            try:
                retval += [myplotfunc(*tuple(plotargs),**newkwargs)]
            except Exception as e:
                raise RuntimeError(strm("Error trying to plot using function",
                    myplotfunc, '\nwith',len(plotargs), "arguments",
                    '\nwhich were\n',plotargs, "\nand had len\n",
                    list(map(len, plotargs)), "and", len(newkwargs),
                    "\noptions", newkwargs, "of len",
                    ', '.join([str(type(j)) + " " + str(j) if np.isscalar(j)
                        else str(len(j)) for j in list(newkwargs.values())])))
            if x_inverted:
                these_xlims = ax.get_xlim()
                ax.set_xlim((max(these_xlims),min(these_xlims)))
        #}}}
        #}}}
    else:
        logger.debug(strm("here are the kwargs",kwargs))
        if 'yerr' in kwargs.keys() and kwargs['yerr'] is None:
            kwargs['yerr'] = myyerror
        plotargs = [j for j in [myx,np.real(myy),myformat] if j is not None]
        try:
            #print 'DEBUG plotting with args',plotargs,'and kwargs',kwargs,'\n\n'
            retval = myplotfunc(*plotargs,**kwargs)
        except Exception as e:
            raise RuntimeError(strm('error trying to plot',type(myplotfunc),'with value',myplotfunc,
                    '\nlength of the np.ndarray arguments:',['shape:'+str(np.shape(j)) if isinstance(j, np.ndarray) else j for j in plotargs],
                    '\nsizes of np.ndarray kwargs',dict([(j,np.shape(kwargs[j])) if isinstance(kwargs[j], np.ndarray) else (j,kwargs[j]) for j in list(kwargs.keys())]),
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
            '\nsizes of arguments:', [np.shape(j) for j in plotargs],
            '\nsizes of np.ndarray kwargs:',
            dict([(j, np.shape(kwargs[j])) for j in
                list(kwargs.keys()) if isinstance(kwargs[j], np.ndarray)])))
    #plt.grid(True)
    #}}}
    return retval
#}}}
