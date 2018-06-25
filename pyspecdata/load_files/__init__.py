r"""This subpackage holds all the routines for reading raw data in proprietary formats.
It's intended to be accessed entirely through the function :func:`find_file`,
which uses :module:`datadir` to search for the filename, then automatically identifies
the file type and calls the appropriate module to load the data into an nddata.

Currently, Bruker file formats (both ESR and NMR) are supported, as well as
(at least some earlier iteration) of Magritek file formats.

Users/developers are very strongly encouraged to add support for new file types.

.. currentmodule:: pyspecdata.load_files

.. autofunction:: find_file

"""
from . import bruker_nmr
from . import prospa
from . import bruker_esr
from . import acert
from .open_subpath import open_subpath
from ..datadir import getDATADIR
from ..datadir import _my_config
from ..general_functions import process_kwargs,strm
from ..core import *
from __builtin__ import any # numpy has an "any" function, which is very annoying
from itertools import tee
import warnings, os, h5py, re
from zipfile import ZipFile, is_zipfile
logger = logging.getLogger('pyspecdata.load_files')

#{{{ add slashes for dir's
def _dirformat(file):
        #{{{ format strings
        if file[-1] not in ['/',path_sep]:
            file += path_sep
        #}}}
        return file
#}}}
def search_filename(searchstring,exp_type,
        print_result=True):
    r"""Use regular expression `searchstring` to find a file inside the directory indicated by `exp_type`

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
        specifies the type of experiment see :func:`~pyspecdata.datadir.getDATADIR` for
        more info.
    """
    #{{{ actually find the files
    directory = getDATADIR(exp_type=exp_type)
    def look_inside(inp_directory):
        logger.debug(strm("looking inside directory",inp_directory))
        dirlist = os.listdir(inp_directory)
        logger.debug(strm("dirlist inside",inp_directory,"is",dirlist))
        if os.path.isdir(inp_directory):
            files = re.findall('.*' + searchstring + '.*','\n'.join(dirlist))
        else:
            raise IOError("I can't find the directory:\n%s\nin order to get a file that matches:\n%s\nYou might need to change the value associated with this exp_type in %s"%(inp_directory,searchstring,_my_config.config_location))
        logger.debug(strm("after running findall, files is",files))
        if len(files) == 0:
            files = []
            directories_inside = [k for k in dirlist if os.path.isdir(inp_directory+k)]
            logger.debug(strm("I found no matches, but I found the directories",directories_inside))
            for j in directories_inside:
                files += [j+os.path.sep+k for k in look_inside(inp_directory+j)]
            return files
        else:
            return files
    files = look_inside(directory)
    logger.debug(strm("look_inside found the files",files))
    if files is None or len(files) == 0:
        exptype_msg = ""
        if exp_type is None:
            exptype_msg = "\nYou probably need to set exp_type so I know where inside {1:s} to find the file."
        raise IOError(("I can't find a file matching the regular expression {0:s} in {1:s}"+exptype_msg).format(searchstring,directory))
    else:
        if len(files) > 1:
            basenames,exts = map(set,zip(*[j.rsplit('.',1) for j in files if len(j.rsplit('.',1))>1]))
            if len(basenames) == 1 and len(exts) == len(files):
                pass
            else:
                warnings.warn('found multiple files:\n'+'\n\t'.join(files)+'\nand opening last')
        elif print_result:
            logger.info("found only one file, and loading it:"+repr(files))
    #}}}
    return [directory+j for j in files]

def find_file(searchstring,
            exp_type = None,
            postproc = None,
            print_result = True,
            verbose = False,
            prefilter = None,
            expno = None,
            dimname='', return_acq=False,
            add_sizes=[], add_dims=[], use_sweep=None,
            indirect_dimlabels=None,
            **kwargs):
    r'''Find the file  given by the regular expression `searchstring` inside the directory identified by `exp_type`, load the nddata object, and postprocess with the function `postproc`.

    It looks at the top level of the directory first, and if that fails, starts to look recursively.
    Whenever it finds a file in the current directory, it will not return data from files in the directories underneath.

    It calls :func:`~pyspecdata.load_files.load_indiv_file`, which finds the specific routine from inside one of the modules (sub-packages) associated with a particular file-type.

    Parameters
    ----------
    searchstring : str
        Most commonly, this is just a fragment of the file name,
        with any literal ``*``, ``.``, or ``?`` characters preceded by
        a backslash.
        More generally, it is a regular expression,
        where ``.*searchstring.*`` matches a filename inside the
        directory appropriate for `exp_type`.
    expno : int
        For Bruker NMR and Prospa files, where the files are stored in numbered
        subdirectories,
        give the number of the subdirectory that you want.
        Currently, this parameter is needed to load Bruker and Kea files.
        If it finds multiple files that match the regular expression,
        it will try to load this experiment number from all the directories.
    exp_type : str
        Since the function assumes that you have different types of
        experiments sorted into different directories, this argument
        specifies the type of experiment see :func:`~pyspecdata.datadir.getDATADIR` for
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

        :dimname:
            passed to :func:`~pyspecdata.load_files.load_indiv_file`

        :return_acq:
            passed to :func:`~pyspecdata.load_files.load_indiv_file`

        :add_sizes: passed to :func:`~pyspecdata.load_files.load_indiv_file`
        :add_dims: passed to :func:`~pyspecdata.load_files.load_indiv_file`
        :use_sweep: passed to :func:`~pyspecdata.load_files.load_indiv_file`
        :indirect_dimlabels: passed to :func:`~pyspecdata.load_files.load_indiv_file`
        '''
    logger.info(strm("find_file sees indirect_dimlabels",
        indirect_dimlabels))
    # {{{ legacy warning
    if 'subdirectory' in kwargs.keys():
        raise ValueError("The `subdirectory` keyword argument is not longer valid -- use `exp_type` instead!")
    # }}}
    files = search_filename(searchstring, exp_type, print_result=print_result)
    data = None
    while data is None and len(files) > 0:
        filename = files.pop(-1)
        # {{{ file loaded here
        logger.debug(strm("about to call load_indiv_file on",filename))
        data = load_indiv_file(filename,
            dimname=dimname, return_acq=return_acq,
            add_sizes=add_sizes, add_dims=add_dims, use_sweep=use_sweep,
            indirect_dimlabels=indirect_dimlabels,
            expno=expno)
        # }}}
    if data is None:
        raise ValueError(strm(
            "I found no data matching the regexp", searchstring))
    logger.debug("about to look at postproc")
    if hasattr(postproc,'__call__'):
        logger.debug("postproc passed explicitly")
        return postproc(data,**kwargs)
    else:
        if postproc is None:
            postproc_type = data.get_prop('postproc_type')
            logger.debug(strm("found postproc_type",postproc_type))
        else:
            logger.debug("found no postproc_type")
        if postproc_type is None:
            logger.debug("got a postproc_type value of None")
            return data
        else:
            if postproc_type in postproc_lookup.keys():
                data = postproc_lookup[postproc_type](data,**kwargs)
            else:
                raise ValueError('postprocessing not defined for file with postproc_type %s --> it should be defined in the postproc_type dictionary in load_files.__init__.py'+postproc_type)
            return data
def format_listofexps(args):
    """**Phased out**: leaving documentation so we can interpret and update old code

    This is an auxiliary function that's used to decode the experiment list.

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
    raise ValueError("format_listofexps was a legacy function that was used to specify multiple experiments to load -- this should be accomplished manually or (once supported) by passing multiple expno values")
def load_file(*args,**kwargs):
    """**Phased out** -- this was used to concatenate files stored in different experiments

    Parameters
    ----------
    args : 
        The files are specified using the format given by
        :func:`format_listofexps`
    """
    raise ValueError("load_file was a legacy function that was used to concatenate several experiments into a 2D dataset -- this should be accomplished manually or (once supported) by passing multiple expno values")
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
            logger.info(strm("Found magic signature, returning",
                retval))
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
        indirect_dimlabels=None, expno=None):
    """Open the file given by `filename`, use file signature magic and/or
    filename extension(s) to identify the file type, and call the appropriate
    function to open it.
    
    Parameters
    ----------
    dimname : str
        When there is a single indirect dimension composed of several scans,
        call the indirect dimension `dimname`.
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

    Returns
    -------
    nddata or None
        the nddata containing the data,
        or else, `None`, indicating that this is part of a pair of
        files that should be skipped
    """
    logger.debug(strm("load_indiv_file sees indirect_dimlabels",
        indirect_dimlabels))
    #to search for kwargs when separating: \<dimname\>\|\<return_acq\>\|\<add_sizes\>\|\<add_dims\>
    if not os.path.exists(filename):
        if os.path.exists(filename+'.par'):
            # {{{ legacy support for WinEPR without extension
            logger.debug("legacy support for WinEPR without extension")
            data = bruker_esr.winepr(filename, dimname=dimname)
            warnings.warn("Don't call load_indiv_file anymore with the"
                    " incomplete filename to try to load the ESR spectrum."
                    " Rather, supply the full name of the .par file.")
            # }}}
        else:
            raise IOError("%s does not exist"%filename)
    # {{{ first, we search for the file magic to determine the filetype
    filetype = None
    if os.path.isdir(filename) or is_zipfile(filename):# Bruker, prospa, etc, datasets are stored as directories
        if expno is None:
            raise ValueError("currently, you must specify the experiment number when trying to load a directory-style file")
        else:
            expno_as_str = '%d'%expno
        if is_zipfile(filename):
            zf = ZipFile(filename)
            list_of_files = [j.split('/') for j in zf.namelist()]
            basename = os.path.normpath(filename).split(os.path.sep)[-1].split('.')[0]
            assert all([j[0] == basename for j in list_of_files]), strm("I expected that the zip file contains a directory called ",basename,"which contains your NMR data -- this appears not to be the case.  (Note that with future extensions, other formats will be possible.)")
            file_reference = (zf,
                    filename,
                    basename)
        else:
            file_reference = filename
        if open_subpath(file_reference, expno_as_str, 'ser', test_only=True):
            #{{{ Bruker 2D
            logger.debug('Identified a bruker series file')
            data = bruker_nmr.series(file_reference, expno_as_str, 'ser', dimname=dimname)
            #}}}
        elif open_subpath(file_reference, expno_as_str, 'acqus', test_only=True):
            logger.debug('Identified a bruker 1d file')
            #{{{ Bruker 1D
            data = bruker_nmr.load_1D(file_reference, expno_as_str, dimname=dimname)
            #}}}
        else:
            logger.debug('Identified a potential prospa file')
            # the rest are for prospa
            # even though I've made steps towards supporting zipped prospa (by
            # using open_subpath), that's not yet done -- diff against this commit
            # to see what the previous code looked like
            file_reference += '/' + expno_as_str
            files_in_dir = os.listdir(file_reference)
            if open_subpath(file_reference, 'data.2d', test_only=True):
                #{{{ Prospa generic 2D
                data = prospa.load_2D(file_reference, dimname=dimname)
                #}}}
                #{{{ specific Prospa formats
            elif any(map((lambda x:'Delay' in x),files_in_dir)):
                data = prospa.load_2D(file_reference, dimname=dimname)
                twod = 't1' # this was set by det_type, not sure what's done with it now
            elif open_subpath(file_reference, 'acqu.par', test_only=True):
                data = prospa.load_1D(file_reference)
            elif os.path.exists(os.path.join(file_reference,'..','acqu.par')):
                data = prospa.load_2D(file_reference, dimname=dimname)
                twod = 't1_sub' # this was set by det_type, not sure what's done with it now
                #}}}
            else:
                raise RuntimeError('WARNING! unidentified file type '+file_reference)
    else:
        logger.debug(strm("the path",filename,"is a file"))
        type_by_signature = _check_signature(filename)
        type_by_extension = _check_extension(filename)
        logger.debug(strm("signature and extension checks are done"))
        if type_by_signature:
            logger.debug(strm("determining type by signature"))
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
                    logger.debug("old-format WinEPR parameter file")
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
            logger.debug(strm("determining type by extension"))
            if type_by_extension == 'SPC':
                logger.info(strm("skipping SPC file",filename))
                return None # ignore SPC, and leave the reading to the PAR file
            elif type_by_extension == 'DTA':
                logger.info(strm("skipping DTA file",filename))
                return None # ignore DTA and load the reading to the DSC file
            elif type_by_extension == 'YGF':
                logger.info(strm("skipping YGA file",filename))
                return None # ignore YGA and load the reading to the DSC file
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
    logger.debug("done with load_indiv_file")
    if '' in data.dimlabels:
        if ndshape(data)[''] < 2:
            data = data['',0]
        else:
            if 'indirect' in data.dimlabels:
                raise ValueError("Was going to rename unnamed dimension 'indirect', but there's already one called that!")
            else:
                data.rename('','indirect')
    return data
    #}}}
postproc_lookup = {
        'ELDOR':acert.postproc_eldor_old,
        'ELDOR_3D':acert.postproc_eldor_3d,
        'FID':acert.postproc_generic,
        'echo_T2':acert.postproc_echo_T2,
        'B1_se':acert.postproc_B1_se,
        'CW':acert.postproc_cw,
        }
def load_t1_axis(file):
    raise RuntimeError("don't use load_t1_axis anymore, the t1 axis should be available as an nddata property called wait_time")
def bruker_load_t1_axis(file):
    raise RuntimeError("don't use bruker_load_t1_axis anymore, the t1 axis should be available as an nddata property called wait_time")
def prospa_t1_info(file):
    raise RuntimeError("don't use prospa_t1_info anymore, the t1 axis should be available as an nddata property called wait_time")
def bruker_load_title(file):
    raise RuntimeError("don't use bruker_load_title -- this should now be loaded as the nddata name")
def cw(file,**kwargs):
    raise RuntimeError("don't use the cw method anymore -- just use find_file")
def det_type(file,**kwargs):
    raise RuntimeError("det_type is deprecated, and should be handled by the file magic inside load_indiv_file.  THE ONE EXCEPTION to this is the fact that det_type would return a second argument that allowed you to classify different types of prospa files.  This is not handled currently")

__all__ = ['find_file',
        'search_filename',
        'load_indiv_file',
        'format_listofexps',
        'bruker_nmr',
        'bruker_esr',
        'acert',
        'prospa',
        'load_t1_axis',
        'bruker_load_t1_axis',
        'prospa_t1_info',
        'bruker_load_title',
        'cw',
        ]
