"the file-loading subroutines are stored here"
from . import bruker_nmr
from . import prospa
from . import bruker_esr
from . import acert
from ..datadir import getDATADIR
from ..general_functions import process_kwargs
from ..core import *

#{{{ add slashes for dir's
def _dirformat(file):
        #{{{ format strings
        if file[-1] not in ['/',path_sep]:
            file += path_sep
        #}}}
        return file
#}}}
def find_file(searchstring,
            exp_type = None,
            postproc = None,
            print_result = True,
            verbose = False,
            prefilter = None,
            indirect_dimlabels = None,# in case "dimlabels" is not set properly, I can manually pass the value of "indirect_dimlabels"
            **kwargs):
    r'''Find the file  given by the regular expression `searchstring` inside the directory identified by `exp_type`, load the nddata object, and postprocess with the function `postproc`.

    It calls `load_indiv_file`, which finds the specific routine from inside one of the modules (sub-packages) associated with a particular file-type.

    Parameters
    ----------
    searchstring : str
        Most commonly, this is just a fragment of the file name,
        with any literal ``*``, ``.``, or ``?`` characters preceded by
        a backslash.
        More generally, it is a regular expression,
        where ``.*searchstring.*`` matches a filename inside the
        directory appropriate for `exp_type`.
    exp_type : str
        Since the function assumes that you have different types of
        experiments sorted into different directories, this argument
        specifies the type of experiment see :func:`getDATADIR` for
        more info.
    postproc : function, str, or None
        This function is fed the nddata data and the remaining keyword
        arguments (`kwargs`) as arguments.
        It's assumed that each module for each different file type
        provides a dictionary called `postproc_lookup`.
        If `postproc` is a string,
        it looks up the string inside the `postproc_lookup`
        dictionary that's appropriate for the file type.
        If `postproc` is None,
        it checks to see if the any of the loading functions that were
        called set the `postproc_type` property
        -- *i.e.* it checks the value of
        ``data.get_prop('postproc_type')`` --
        if this is set, it uses this as a key
        to pull the corresponding value from `postproc_lookup`.

        For instance, when the acert module loads an ACERT HDF5 file,
        it sets `postproc_type` to the value of
        ``(h5 root).experiment.description['class']``.
        This, in turn, is used to choose the type of post-processing.
        '''
    # {{{ legacy warning
    if 'subdirectory' in kwargs.keys():
        raise ValueError("The `subdirectory` keyword argument is not longer valid -- use `exp_type` instead!")
    # }}}
    #{{{ actually find one file and load it into the h5 object
    directory = getDATADIR(exp_type=exp_type)
    if os.path.isdir(directory):
        files = re.findall('.*' + searchstring + '.*','\n'.join(os.listdir(directory)))
    else:
        raise IOError("I can't find the directory:\n%s\nin order to get a file that matches:\n%s"%(directory,searchstring))
    if len(files) == 0:
        raise IOError("I can't find a file matching the regular expression {:s} in {:s}".format(searchstring,directory))
    else:
        if len(files) > 1:
            warnings.warn('found multiple files:\n'+repr(files)+'\nand opening last')
        elif print_result and verbose:
            obsn("found only one file, and loading it:"+repr(files))
    filename = directory + files[-1]
    #}}}
    # {{{ file loaded here
    loader_module = acert
    data = loader_module.load(filename)
    # }}}
    if hasattr(postproc,'__call__'):
        return postproc(data,**kwargs)
    else:
        if postproc is None:
            postproc_type = data.get_prop('postproc_type')
        if postproc is None:
            return data
        else:
            if postproc_type in loader_module.postproc_lookup.keys():
                data = loader_module.postproc_lookup[desc_class](data)
            else:
                raise ValueError('postprocessing not defined for file with description-->class '+str(data.get_prop('description_class')))
            return data
def format_listofexps(args):
    """This is an auxiliary function that's used to decode the experiment list.

    Parameters
    ----------
    args : list or tuple
        can be in one of two formats
        :``(dirname,[i,j,k,...N])``:  typically used, *e.g.* for
            Bruker NMR experiments.  ``i,j,...N`` are integer numbers
            referring to individual experiments that are stored in
            subdirectories of `dirname` (a string).
        :``([exp_name1,...,exp_nameN])``: just return this list of
            experiments given by the strings
            `exp_name1`...`exp_nameN`.
        :``([exp_name1,...,exp_nameN],[])``: identical to previous
        :``([exp_name1,...,exp_nameN],[])``: identical to previous
        :``(exp_name1,...,exp_nameN)``: identical to previous
        :``(exp_name)`` or ``(exp_name,[])``: works for a single
            experiment
    """
    if type(args[0]) is str: # even if it's just a string, make it a list
        args[0] = [args[0]]
    if len(args) > 1 and ((not isscalar(args[1]))
            and len(args[1]) == 0): args.pop(1)
    if len(args) == 2:
        if isscalar(args[1]):
            if type(args[1]) is str:
                filenames = args
            else:
                args[1] = [args[1]] # if the second argument is a single file number, make it into a list
        if not isscalar(args[1]):
            if len(args[0]) > 1:
                raise ValueError("Invalid arguments to specify a list of experiments (see documentation)")
            filenames = [_dirformat(args[0][0]) + '%d'%x for x in args[1]]
    else:
        filenames = args
    return filenames
def load_file(*args,**kwargs):
    """Load a file or series of files into a single dataset.

    Parameters
    ----------
    args : 
        The files are specified using the format given by
        :func:`format_listofexps`
    """
    args = list(args)
    dimname, calibration, add_sizes, add_dims = process_kwargs(
            [('dimname',''),
                ('calibration',1.0),
                ('add_sizes',[]),
                ('add_dims',[])],kwargs)
    filenames = format_listofexps(args)
    #{{{load all the data into a list
    data = [load_indiv_file(filenames[0],dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    for filename in filenames[1:]:
        data += [load_indiv_file(filename,dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    #}}}
    # for the following, I used to have a condition, but this is incompatible with the pop statement at the end
    newdata = concat(data,dimname) # allocate the size of the indirect array
    newdata_shape = ndshape(newdata)
    # {{{ because the "want_to_prospa_decim_correct" attribute is special --
    #     probably want to get rid of it eventually
    if hasattr(data[0],'want_to_prospa_decim_correct'):
        if all(map(hasattr(x,'want_to_prospa_decim_correct'),data)):
            if data[0].want_to_prospa_decim_correct is True:
                newdata = prospa_decim_correct(newdata)
    # }}}
    if newdata_shape[dimname]==1:
        newdata.popdim(dimname)
    return newdata*calibration
def load_indiv_file(filename, dimname='', return_acq=False,
        add_sizes=[], add_dims=[]):
    """Open the file given by `filename`, use magic (broadly defined)
    to identify the file type, and call the appropriate function to
    open it."""
    #to search for kwargs when separating: \<dimname\>\|\<return_acq\>\|\<add_sizes\>\|\<add_dims\>
    # {{{ first, we search for the file magic to determine the filetype
    file_signatures = {'\x89\x48\x44\x46\x0d\x0a\x1a\x0a':'HDF5'}
    filetype = None
    #{{{ WinEPR
    if os.path.exists(filename+'.spc'):
        return ('winepr',True)
    #}}}
    else:
        filename = dirformat(filename)
        files_in_dir = os.listdir(filename)
        #{{{ Bruker 2D
        if os.path.exists(filename+'ser'):
            return ('bruker',True)
        #}}}
        #{{{ Prospa generic 2D
        elif os.path.exists(filename+'data.2d'):
            return ('prospa',True)
        #}}}
        #{{{ specific Prospa formats
        elif any(map((lambda x:'Delay' in x),files_in_dir)):
            return ('prospa','t1')
        elif os.path.exists(filename+'acqu.par'):
            return ('prospa',False)
        elif os.path.exists(filename+'../acqu.par'):
            return ('prospa','t1_sub')
        #}}}
        #{{{ Bruker 1D
        elif os.path.exists(filename+'acqus'):
            return ('bruker',False)
        #}}}
        else:
            raise CustomError('WARNING! unidentified file type '+filename)
    # }}}
    if filetype == 'winepr':
        data = bruker_esr.winepr(filename, dimname=dimname)
    filename = dirformat(filename)
    if twod and filetype == 'bruker':
        data = bruker_nmr.series(filename, dimname=dimname)
    elif twod and filetype == 'prospa':
        data = prospa.load_2D(filename, dimname=dimname)
    else:
        if filetype == 'bruker':
            #{{{ bruker 1D
            data = bruker_nmr.load_1D(filename, dimname=dimname)
            #}}}
        elif filetype == 'prospa':
            #{{{ prospa 1d
            data = prospa.load_1D(filename)
            #}}}
        else:
            raise CustomError("can't load this file type $\\rightarrow$ \\verb+%s+"%filename)
    #{{{ return, and if necessary, reorganize
    if len(add_sizes)>0:
        data.labels([dimname],[[]]) # remove the axis, so we can reshape
        #print 'DEBUG: data before chunk = ',data
        data.chunkoff(dimname,add_dims,add_sizes)
        #print 'DEBUG: data after chunk = ',data
        data.labels(add_dims,
                [r_[0:x] for x in add_sizes])
    if return_acq:
        raise ValueError('return_acq is deprecated!! All properties are now set directly to the nddata using the set_prop function')
    return data
    #}}}
def load_acqu(filename,whichdim='',return_s = None):
    "This should be deleted -- a placeholder"
    filename = dirformat(filename)
    if det_type(filename)[0] == 'bruker':
        # {{{ this should be OK, since load_acqu from within bruker already calls bruker.load_acqu
        if return_s is not None:
            return bruker.load_acqu(filename,whichdim=whichdim,return_s = return_s)
        else:
            return bruker.load_acqu(filename,whichdim=whichdim)
        # }}}
    elif det_type(filename)[0] == 'prospa':
        # {{{ somehow, I need to deal with the t1_sub within prospa
        if det_type(filename)[1] == 't1_sub':
            filename = dirformat(filename)
            return prospa.load_acqu(filename+'../')
        else:
            return prospa.load_acqu(filename)
        # }}}
    else:
        raise CustomError(det_type(filename),'is not yet supported')
#def load_t1_axis(file):
#    raise RuntimeError("don't use load_t1_axis anymore, the t1 axis should be available as an nddata property called wait_time")
#def bruker_load_t1_axis(file):
#    raise RuntimeError("don't use bruker_load_t1_axis anymore, the t1 axis should be available as an nddata property called wait_time")
#def prospa_t1_info(file):
#    raise RuntimeError("don't use prospa_t1_info anymore, the t1 axis should be available as an nddata property called wait_time")
#def bruker_load_title(file):
#    raise RuntimeError("don't use bruker_load_title -- this should now be loaded as the nddata name")
#def cw(file,**kwargs):
#    raise RuntimeError("don't use the cw method anymore -- just use find_file")

__all__ = ['find_file',
        'load_indiv_file',
        'format_listofexps',
        'bruker_nmr',
        'bruker_esr',
        'acert',
        'prospa',
#        'load_t1_axis',
#        'bruker_load_t1_axis',
#        'prospa_t1_info',
#        'bruker_load_title',
#        'cw',
        ]
