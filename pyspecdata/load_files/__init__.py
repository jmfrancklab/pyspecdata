"the file-loading subroutines are stored here"
from . import bruker_nmr
from . import prospa
from . import bruker_esr
from . import acert
from ..datadir import getDATADIR
from ..general_functions import process_kwargs
from ..core import *
from __builtin__ import any # numpy has an "any" function, which is very annoying
from itertools import tee
import warnings,os,h5py

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
    indirect_dimlabels : 
        passed through to :func:`acert.load_pulse`
    use_sweep : 
        passed through to :func:`acert.load_cw`
        '''
    # {{{ legacy warning
    if 'subdirectory' in kwargs.keys():
        raise ValueError("The `subdirectory` keyword argument is not longer valid -- use `exp_type` instead!")
    # }}}
    # {{{ 
    indirect_dimlabels, use_sweep = process_kwargs([('indirect_dimlabels',None),
        ('use_sweep',False)], kwargs, pass_through=True)
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
                raise ValueError('postprocessing not defined for file with postproc_type'+postproc_type)
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
def _check_extension(filename):
    "Just return the file extension in caps"
    return filename.split('.')[-1].upper()
def _check_signature(filename):
    """Check the filetype by its signature (the leading part of the file).
    If the first several characters are all ASCII, return the string ``TXT``.

    Returns
    -------
    str
        Either
        a short string identifying the filetype
        (currently "HDF5", "DOS Format"  or "TXT")
        OR `None` if the type is unknown.
    """
    file_signatures = {'\x89\x48\x44\x46\x0d\x0a\x1a\x0a':'HDF5',
            'DOS  Format':'DOS Format'}
    max_sig_length = max(map(len,file_signatures.keys()))
    with open(filename,'rb') as fp:
        inistring = fp.read(max_sig_length)
        if any(thiskey in inistring for thiskey in
                file_signatures.keys()):# this will fail without overriding the numpy any( above
            retval = file_signatures[(thiskey for thiskey in
                file_signatures.keys() if thiskey in
                inistring).next()]
            print "Found magic signature, returning",retval
            return retval
        else:
            try:
                inistring.decode('ascii')
                return 'TXT'
            except UnicodeDecodeError:
                # if it failed, it's because the string is not ASCII
                return None
def load_indiv_file(filename, dimname='', return_acq=False,
        add_sizes=[], add_dims=[], use_sweep=None,
        indirect_dimlabels=None):
    """Open the file given by `filename`, use file signature magic and/or
    filename extension(s) to identify the file type, and call the appropriate
    function to open it.
    
    Parameters
    ----------
    dimname : str
        When there is a single indirect dimension composed of several scans,
        call it this.
    return_acq : DEPRECATED
    add_sizes : list
        the sizes associated with the dimensions in add_dims
    add_dims : list
        Can only be used with `dimname`.
        Break the dimension `dimname` into several dimensions,
        with the names given by the list `add_dims` and sizes given by `add_sizes`.
        If the product of the sizes is not the same as the original dimension
        given by `dimname`,
        retain it as the "outermost" (leftmost) dimension. 
        :func:`pyspecdata.core.chunkoff` is used to do this, like so:
        ``data.chunkoff(dimname,add_dims,add_sizes)``
    indirect_dimlabels : str or None
        passed through to `acert.load_pulse` (names an indirect dimension when dimlabels isn't provided)
    """
    #to search for kwargs when separating: \<dimname\>\|\<return_acq\>\|\<add_sizes\>\|\<add_dims\>
    if not os.path.exists(filename):
        if os.path.exists(filename+'.par'):
            # {{{ legacy support for WinEPR without extension
            data = bruker_esr.winepr(filename, dimname=dimname)
            warnings.warn("Don't call load_indiv_file anymore with the"
                    " incomplete filename to try to load the ESR spectrum."
                    " Rather, supply the full name of the .par file.")
            # }}}
        else:
            raise IOError("%s does not exist"%filename)
    # {{{ first, we search for the file magic to determine the filetype
    filetype = None
    if os.path.isdir(filename):# Bruker, prospa, etc, datasets are stored as directories
        files_in_dir = os.listdir(filename)
        if os.path.exists(os.path.join(filename,'ser')):
            #{{{ Bruker 2D
            data = bruker_nmr.series(filename, dimname=dimname)
            #}}}
        elif os.path.exists(os.path.join(filename,'acqus')):
            #{{{ Bruker 1D
            data = bruker_nmr.load_1D(filename, dimname=dimname)
            #}}}
        elif os.path.exists(os.path.join(filename,'data.2d')):
            #{{{ Prospa generic 2D
            data = prospa.load_2D(filename, dimname=dimname)
            #}}}
            #{{{ specific Prospa formats
        elif any(map((lambda x:'Delay' in x),files_in_dir)):
            data = prospa.load_2D(filename, dimname=dimname)
            twod = 't1' # this was set by det_type, not sure what's done with it now
        elif os.path.exists(os.path.join(filename,'acqu.par')):
            data = prospa.load_1D(filename)
        elif os.path.exists(os.path.join(filename,'..','acqu.par')):
            data = prospa.load_2D(filename, dimname=dimname)
            twod = 't1_sub' # this was set by det_type, not sure what's done with it now
            #}}}
        else:
            raise RuntimeError('WARNING! unidentified file type '+filename)
    else:
        type_by_signature = _check_signature(filename)
        type_by_extension = _check_extension(filename)
        if type_by_signature:
            if type_by_signature == 'HDF5':
                # we can have the normal ACERT Pulse experiment or the ACERT CW format
                with h5py.File(filename,'r') as h5:
                    try:
                        description_class = h5['experiment']['description'].attrs['class']
                    except:
                        raise IOError("I am assuming this is an ACERT datafile,"
                                " but can't identify the type, because"
                                " I can't find"
                                " experiment.description['class']")
                if description_class == 'CW':
                    data = acert.load_cw(filename, use_sweep=use_sweep)
                else:
                    data = acert.load_pulse(filename, indirect_dimlabels=indirect_dimlabels)
            elif type_by_signature == 'DOS Format':
                if type_by_extension == 'PAR':
                    # par identifies the old-format WinEPR parameter file, and spc the binary spectrum
                    data = bruker_esr.winepr(filename, dimname=dimname)
                else:
                    raise RuntimeError("I'm not able to figure out what file type %s this is!"%filename)
            elif type_by_signature == 'TXT':
                if type_by_extension == 'DSC':
                    # DSC identifies the new-format XEpr parameter file, and DTA the binary spectrum
                    data = bruker_esr.xepr(filename, dimname=dimname)
                else:
                    raise RuntimeError("I'm not able to figure out what file type %s this is!"%filename)
            else:
                raise RuntimeError("Type %s not yet supported!"%type_by_signature)
        else:
            if type_by_extension == 'SPC':
                pass # ignore SPC, and leave the reading to the PAR file
            elif type_by_extension == 'DTA':
                pass # ignore DTA and load the reading to the DSC file
            else:
                raise RuntimeError("I'm not able to figure out what file type %s this is!"%filename)
    # }}}
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
            return prospa.load_acqu(os.path.join(filename,'..'))
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
#def det_type(file,**kwargs):
#    raise RuntimeError("det_type is deprecated, and should be handled by the file magic inside load_indiv_file.  THE ONE EXCEPTION to this is the fact that det_type would return a second argument that allowed you to classify different types of prospa files.  This is not handled currently")

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
