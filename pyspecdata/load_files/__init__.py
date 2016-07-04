"the file-loading subroutines are stored here"
from . import bruker_nmr
from . import prospa
from . import bruker_esr
from . import acert
from .datadir import getDATADIR
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
        This function is fed the nddata data and the source text as
        arguments,
        as well as the keyword arguments "kwargs".
        If it is not set, it's chosen automatically.
        For instance, the post-processing for ACERT HDF5 files
        is chosen based on the value of
        ``(h5 root).experiment.description['class']``
        '''
    #{{{ actually find one file and load it into the h5 object
    if 'subdirectory' in kwargs.keys():
        raise ValueError("The `subdirectory` keyword argument is not longer valid -- use `exp_type` instead!")
    filename = search_freed_file(searchstring,exp_type,print_result = print_result,verbose = verbose)
    h5 = h5py.File(filename,'r')
    #}}}
    #{{{ set up the complex number the hard way, for good form
    data = empty_like(h5['experiment']['data.r'],dtype = complex128)
    data_introspect = data.view([('r',double),('i',double)])
    for i_or_r in ['i','r']:
        data_introspect[i_or_r] = h5['experiment']['data.'+i_or_r]
    #}}}
    #{{{ start with the dimensions used in the HDF5 file -- "indirect" is a placeholder (see below)
    if len(data.shape) == 4:
        dimlabels = ['bin','phcyc','indirect','t2']
    elif len(data.shape) == 3:
        dimlabels = ['phcyc','indirect','t2']
    elif len(data.shape) == 2:
        dimlabels = ['phcyc','t2']
    else:
        raise ValueError("I don't know how to interpret data of shape %d"%len(data.shape))
    data = nddata(data,list(data.shape),dimlabels)
    #}}}
    #{{{ put all further information into the nddata in a way such that it can be retrieved with "get_prop"
    data.set_prop(dict([(k,v) for k,v in h5['experiment'].attrs.iteritems() if k[0] != '_' and k[0] != 'dimlabels']))
    data.set_prop(dict([('execution_'+k,v) for k,v in h5['experiment']['execution'].attrs.iteritems() if k[0] != '_']))
    data.set_prop(dict([('description_'+k,v) for k,v in h5['experiment']['description'].attrs.iteritems() if k[0] != '_']))
    #}}}
    #{{{ finish the t2 axis
    t2_steps = data.get_prop('t2_steps')
    t2_inc = data.get_prop('t2_inc')
    if t2_steps is not None:
        if ndshape(data)['t2'] != t2_steps:
            raise ValueError("There is a problem, because the dimensions of your data don't match the value given by t2_steps")
        data.labels('t2',r_[0:t2_steps]*t2_inc)
    else:
        if verbose: print "warning, I couldn't find the t2 steps parameter"
        data.labels('t2',r_[0:ndshape(data)['t2']]*1e-9)
    data.set_units('t2','s')
    if prefilter is not None:
        data.ft('t2',shift = True)
        data = data['t2':prefilter]
    #}}}
    #{{{ now, pull the dimlabels that we want
    if 'dimlabels' in h5['experiment'].attrs.keys():
        dimlabels = h5['experiment'].attrs['dimlabels'].tolist()
    else:
        print "You didn't set dimlabels in your pulse program -- the data is ambiguous and I'm taking my best guess!!!"
    if indirect_dimlabels is not None:# in case it's not stored correctly in the file and I need to hard-code it
        idx = dimlabels.index('indirect')
        dimlabels = dimlabels[:idx] + indirect_dimlabels + dimlabels[idx+1:]
    #}}}
    #{{{ assign all the dimensions
    indirect_names_in_h5 = [k for k,v in h5['experiment'].iteritems() if k not in ['data.i','data.r','bin_switch_times','fields','sweep_currents'] and type(v) == h5py.Dataset]
    if verbose: print "dimlabels are",dimlabels
    expected_dimlabels = set(dimlabels) - {'phcyc','t2','bin'} # we expect to find these axes stored in the HDF5
    if verbose: print "expected dimlabels",expected_dimlabels
    common_dimensions = expected_dimlabels & set(indirect_names_in_h5)
    #{{{ deal with the case where dimlabels and the stored axes don't match up correctly 
    if verbose: print "common dimensions",common_dimensions
    unlabeled_indirect = False
    if len(indirect_names_in_h5) == 0:
        warnings.warn('Warning!! The indirect dimension is not labeled!\nDimensions of the data:%s\nfilename: %s'%(repr(ndshape(data)),filename))
        unlabeled_indirect = True
    if len(common_dimensions) != len(indirect_names_in_h5):
        h5_but_not_dimlabels = common_dimensions ^ set(indirect_names_in_h5)
        dimlabels_but_not_h5 = common_dimensions ^ expected_dimlabels
        if len(h5_but_not_dimlabels) == 1 and len(dimlabels_but_not_h5) == 1:
            real_name = h5_but_not_dimlabels.pop()
            dimlabels_name = dimlabels_but_not_h5.pop()
            warnings.warn("Warning! The dimlabels set in the HDF5 file have a dimension called %s and the arrays given in the HDF5 file have one called %s -- I'm assuming they're the same and using the name from the HDF5 file."%(dimlabels_name,real_name))
            common_dimensions |= {real_name} # so the axis gets labeled below
            dimlabels[dimlabels.index(dimlabels_name)] = real_name
        elif dimlabels_but_not_h5 == {'indirect'}:# I need to split up the indirect dimension
            raise ValueError("You want me to split the indirect dimension, but I don't know how!:\n\tDimensions in the file: %s\n\tDimensions listed in dimlabels %s\n\tDimensions of the data: %s"%(repr(indirect_names_in_h5),repr(dimlabels),repr(ndshape(data))))
        else:
            raise ValueError("This file has dimensions that I can't just automatically sort!:\n\tDimensions in the file: %s\n\tDimensions listed in dimlabels %s\n\tDimensions of the data: %s\n\tDimensions in hdf5 but not dimlabels: %s\n\tDimensions in dimlabels but not HDF5: %s"%(repr(indirect_names_in_h5),repr(dimlabels),repr(ndshape(data)),repr(h5_but_not_dimlabels),repr(dimlabels_but_not_h5)))
    #}}}
    if verbose: print "now dimlabels are",dimlabels,"and indirect_names_in_h5 is",indirect_names_in_h5
    #{{{ chunk the "indirect" dimension up appropriately, and assign the dimensions
    if not unlabeled_indirect:
        if 'indirect' in data.dimlabels:# note that if there was only a single indirect dimension, it's already been used up
            chunk_dict = dict([(j,len(h5['experiment'][j])) for j in dimlabels if j not in ['bin','phcyc','t2']])
            if verbose: print ndshape(data)
            if verbose: print chunk_dict
            # {{{ if an experiment with a single dimension was abandoned early, we know what to do about that
            if (len(chunk_dict) == 1 and chunk_dict.values()[0] >
                    data.data.shape[data.axn('indirect')]):
                chunk_dict[chunk_dict.keys()[0]] = data.data.shape[data.axn('indirect')]
                warnings.warn("Warning: the length of the 'indirect' axis seems to exceed the length of the data!")
            # }}}
            data.chunk('indirect',chunk_dict)
        for this_axis in common_dimensions:
            axis_data = array(h5['experiment'][this_axis])
            temp = ndshape(data)[this_axis]
            if len(axis_data) > temp:
                warnings.warn("Warning: the length of the axis '"+str(this_axis)+"' seems to exceed the length of the data!")
                axis_data = axis_data[0:temp]
            data.labels(this_axis,axis_data)
    #}}}
    #}}}
    #{{{ label the bin dimension, and chunk it up appropriately
    if 'bin' in data.dimlabels:
        #{{{ assign the labels for the bin dimension as a structured array
        forstruct = []
        forstruct_names = []
        for thisaxis in [x for x in 'bin_switch_times','fields','sweep_currents' if x in h5['experiment'].keys() and type(h5['experiment'][x]) is h5py.Dataset]:
            forstruct.append((h5['experiment'][thisaxis])[0])
            forstruct_names.append(thisaxis)
        print "structure is",forstruct,"names",forstruct_names,"zeros like",h5['experiment']['bin_switch_times'].shape
        x = make_rec(forstruct,forstruct_names,zeros_like = h5['experiment']['bin_switch_times'].shape)
        warnings.warn("detected a 'bin' dimension, and associating it with dtype "+repr(x.dtype))
        for thisaxis in forstruct_names:
            x[thisaxis] = h5['experiment'][thisaxis]
        data.labels('bin',x)
        #}}}
        if 'fields' in x.dtype.names:
            fields,indeces = unique(x['fields'],return_index = True)
            if verbose: print "number of fields",len(fields),"number of unique fields",len(x['fields']),'unique fields are',x['fields']
            if len(fields) != len(x['fields']):
                spacing = diff(indeces)
                if all(spacing == spacing[0]):
                    data.chunk('bin',['field','bin'],[spacing[0],len(x)/spacing[0]])
                else:
                    raise ValueError('The fields have some weird/inconsistent spacing')
        if data.getaxis('bin').dtype.names == ('bin_switch_times',):
            data.setaxis('bin',data.getaxis('bin')['bin_switch_times'])
            data.set_units('bin','s')
    #}}}
    #{{{ it always makes sense to put the phase cycle on the outside
    oldorder = list(data.dimlabels)
    data.reorder(['phcyc'])
    #}}}
    h5.close()
    temp = ndshape(data)
    if 'indirect' in temp.dimlabels and temp['indirect'] == 1:
        data = data['indirect',0]
    if postproc is None:
        postproc_dict = {
                'ELDOR':postproc_eldor_old,
                'ELDOR_3D':postproc_eldor_3d,
                'FID':postproc_generic,
                'echo_T2':postproc_echo_T2,
                'B1_se':postproc_B1_se,
                }
        desc_class = data.get_prop('description_class')
        if desc_class in postproc_dict.keys():
            data = postproc_dict[desc_class](data)
        else:
            raise ValueError('postprocessing not defined for file with description-->class '+str(data.get_prop('description_class')))
        return data
    else:
        return postproc(data,**kwargs)
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
    filetype,twod = det_type(filename)
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

__all__ = ['find_file',
        'load_indiv_file',
        'format_listofexps',
        'bruker_nmr',
        'bruker_esr',
        'acert',
        'prospa',
        ]
