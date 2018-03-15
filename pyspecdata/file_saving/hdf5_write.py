from ..general_functions import *
from pylab import * 

def hdf5_write(self, h5path, directory='.', verbose=False):
    r"""Write the nddata to an HDF5 file.

    `h5path` is the name of the file followed by the node path where
    you want to put it -- it does **not** include the directory where
    the file lives.
    The directory can be passed to the `directory` argument.
    
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
    #{{{ add the final node based on the name stored in the nddata structure
    if h5path[-1] != '/': h5path += '/' # make sure it ends in a slash first
    try:
        thisname = self.get_prop('name')
    except:
        raise ValueError(strm("You're trying to save an nddata object which",
                "does not yet have a name, and you can't do this! Run",
                "yourobject.name('setname')"))
    if type(thisname) is str:
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
        mydataattrs = filter((lambda x: x[0:4] == 'data'),myattrs)
        myotherattrs = filter((lambda x: x[0:4] != 'data'),myattrs)
        myotherattrs = filter(lambda x: x not in self._nosave,myotherattrs)
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
                h5attachattributes(datatable,mydataattrs,self, self._nosave)
        else:
            raise ValueError("I can't find the data object when trying to save the HDF5 file!!")
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
                        h5attachattributes(datatable,myaxisattrsforthisdim.keys(),myaxisattrsforthisdim.values(), self._nosave)
        #}}}
        #{{{ Check the remaining attributes.
        if verbose: print lsafe('other attributes:',zip(myotherattrs,map(lambda x: type(self.__getattribute__(x)),myotherattrs))),'\n\n'
        if verbose: print "Writing remaining other attributes\n\n"
        if len(myotherattrs) > 0:
            #print 'DEBUG 4: bottomnode is',bottomnode
            test = repr(bottomnode) # somehow, this prevents it from claiming that the bottomnode is None --> some type of bug?
            h5attachattributes(bottomnode,
                [j for j in myotherattrs if not self._contains_symbolic(j)],
                self,
                self._nosave)
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
    finally:
        h5file.close()
