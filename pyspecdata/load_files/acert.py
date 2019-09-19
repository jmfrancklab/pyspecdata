'''Open ACERT-format HDF5 files.
Provides post-processing routines for:
'ELDOR'
'ELDOR_3D'
'FID'
'echo_T2'
'B1_se'
(which are typically
experiment names set in
``(h5 root).experiment.description['class']``)
'''        
from ..core import *
from ..general_functions import *
import h5py
logger = logging.getLogger('pyspecdata.load_files.acert')
def load_pulse(filename,
        indirect_dimlabels=None,
        prefilter=None,
        ):
    """Load ACERT pulse data from the 95 GHz.

    Parameters
    ----------
    indirect_dimlabels : str
        In case `dimlabels` is not set properly, I can manually pass the value of `indirect_dimlabels`.
    prefilter : tuple
        If prefilter is set,
        FT the result, and select a specific slice.
        I should think of a more general way of doing this,
        where I pass an ndshape-based slice, instead.
    """
    print("load_pulse sees indirect_dimlabels",indirect_dimlabels)
    with h5py.File(filename,'r') as h5:
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
        data.set_prop(dict([(k,v) for k,v in h5['experiment'].attrs.items() if k[0] != '_' and k[0] != 'dimlabels']))
        data.set_prop(dict([('execution_'+k,v) for k,v in h5['experiment']['execution'].attrs.items() if k[0] != '_']))
        data.set_prop(dict([('description_'+k,v) for k,v in h5['experiment']['description'].attrs.items() if k[0] != '_']))
        #}}}
        #{{{ finish the t2 axis
        t2_steps = data.get_prop('t2_steps')
        t2_inc = data.get_prop('t2_inc')
        if t2_steps is not None:
            if ndshape(data)['t2'] != t2_steps:
                raise ValueError("There is a problem, because the dimensions of your data don't match the value given by t2_steps")
            data.labels('t2',r_[0:t2_steps]*t2_inc)
        else:
            logger.info(strm("warning, I couldn't find the t2 steps parameter"))
            data.labels('t2',r_[0:ndshape(data)['t2']]*1e-9)
        data.set_units('t2','s')
        if prefilter is not None:
            data.ft('t2',shift = True)
            data = data['t2':prefilter]
        #}}}
        #{{{ now, pull the dimlabels that we want
        if 'dimlabels' in list(h5['experiment'].attrs.keys()):
            dimlabels = h5['experiment'].attrs['dimlabels'].tolist()
        else:
            print("You didn't set dimlabels in your pulse program -- the data is ambiguous and I'm taking my best guess!!!")
        if indirect_dimlabels is not None:# in case it's not stored correctly in the file and I need to hard-code it
            idx = dimlabels.index('indirect')
            dimlabels = dimlabels[:idx] + indirect_dimlabels + dimlabels[idx+1:]
        #}}}
        #{{{ assign all the dimensions
        indirect_names_in_h5 = [k for k,v in h5['experiment'].items() if k not in ['data.i','data.r','bin_switch_times','fields','sweep_currents'] and isinstance(v, h5py.Dataset)]
        logger.info(strm("dimlabels are",dimlabels))
        expected_dimlabels = set(dimlabels) - {'phcyc','t2','bin'} # we expect to find these axes stored in the HDF5
        logger.info(strm("expected dimlabels",expected_dimlabels))
        common_dimensions = expected_dimlabels & set(indirect_names_in_h5) # common dimensions are expected, and are also labeled in the HDF5
        #{{{ deal with the case where dimlabels and the stored axes don't match up correctly 
        logger.info(strm("common dimensions",common_dimensions))
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
        logger.info(strm("now dimlabels are",dimlabels,"and indirect_names_in_h5 is",indirect_names_in_h5))
        #{{{ chunk the "indirect" dimension up appropriately, and assign the dimensions
        if not unlabeled_indirect:
            if 'indirect' in data.dimlabels:# note that if there was only a single indirect dimension, it's already been used up
                filtered_dims = [j for j in dimlabels if j not in  ['bin','phcyc','t2']]
                chunk_dict = dict([(j,len(h5['experiment'][j])) if j in list(h5['experiment'].keys())
                    else (j,len(h5['experiment'].attrs[j]))
                    for j in filtered_dims])
                common_dimensions |= set(chunk_dict.keys())
                logger.info(strm("updated common dimensions to include dimensions I pulled from attributes:",common_dimensions))
                logger.info(strm(ndshape(data)))
                logger.info(strm(chunk_dict))
                # {{{ if an experiment with a single dimension was abandoned early, we know what to do about that
                if (len(chunk_dict) == 1 and list(chunk_dict.values())[0] >
                        data.data.shape[data.axn('indirect')]):
                    chunk_dict[list(chunk_dict.keys())[0]] = data.data.shape[data.axn('indirect')]
                    warnings.warn("Warning: the length of the 'indirect' axis seems to exceed the length of the data!")
                # }}}
                data.chunk('indirect',chunk_dict)
            for this_axis in common_dimensions:
                if this_axis in list(h5['experiment'].keys()):
                    axis_data = array(h5['experiment'][this_axis])
                else:
                    axis_data = array(h5['experiment'].attrs[this_axis])
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
            for thisaxis in [x for x in ('bin_switch_times','fields','sweep_currents') if x in list(h5['experiment'].keys()) and isinstance(h5['experiment'][x], h5py.Dataset)]:
                forstruct.append((h5['experiment'][thisaxis])[0])
                forstruct_names.append(thisaxis)
            print("structure is",forstruct,"names",forstruct_names,"zeros like",h5['experiment']['bin_switch_times'].shape)
            x = make_rec(forstruct,forstruct_names,zeros_like = h5['experiment']['bin_switch_times'].shape)
            warnings.warn("detected a 'bin' dimension, and associating it with dtype "+repr(x.dtype))
            for thisaxis in forstruct_names:
                x[thisaxis] = h5['experiment'][thisaxis]
            data.labels('bin',x)
            #}}}
            if 'fields' in x.dtype.names:
                fields,indeces = unique(x['fields'],return_index = True)
                logger.info(strm("number of fields",len(fields),"number of unique fields",len(x['fields']),'unique fields are',x['fields']))
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
    temp = ndshape(data)
    if 'indirect' in temp.dimlabels and temp['indirect'] == 1:
        data = data['indirect',0]
    if data.get_prop('postproc_type') is None:
        data.set_prop('postproc_type',
                data.get_prop('description_class'))
    return data
def load_cw(filename,use_sweep = False):
    """load the cw file given by filename
    
    Parameters
    ----------
    use_sweep : bool
        If true, return the axis labeled by sweep current rather than by field.
    """
    with h5py.File(filename,'r') as h5:
        #{{{ set up the complex number the hard way, for good form
        data = empty_like(h5['experiment']['data.r'],dtype = complex128)
        data_introspect = data.view([('r',double),('i',double)])
        for i_or_r in ['i','r']:
            data_introspect[i_or_r] = h5['experiment']['data.'+i_or_r]
        #}}}
        if len(data.shape) == 1:
            if use_sweep:
                data = nddata(data,data.size,['current']).labels('current',array(h5['experiment']['sweep_currents']))
            else:
                data = nddata(data,data.size,['field']).labels('field',array(h5['experiment']['fields']))
        elif len(data.shape) == 2:
            if use_sweep:
                data = nddata(data,data.shape,['repeats','current']).labels('current',array(h5['experiment']['sweep_currents']))
            else:
                data = nddata(data,data.shape,['repeats','field']).labels('field',array(h5['experiment']['fields']))
        #{{{ put all further information into the nddata in a way such that it can be retrieved with "get_prop"
        data.set_prop(dict([(k,v) for k,v in h5['experiment'].attrs.items() if k[0] != '_' and k[0] != 'dimlabels']))
        data.set_prop(dict([('execution_'+k,v) for k,v in h5['experiment']['execution'].attrs.items() if k[0] != '_']))
        data.set_prop(dict([('description_'+k,v) for k,v in h5['experiment']['description'].attrs.items() if k[0] != '_']))
        #}}}
    if data.get_prop('postproc_type') is None:
        data.set_prop('postproc_type','CW')
        logger.debug(strm("set postproc type for cw to",data.get_prop("postproc_type")))
    return data
def postproc_blank(data):
    return data
def postproc_generic(data):
    data = automagical_phasecycle(data)
    return data
def postproc_B1_se(data):
    data = automagical_phasecycle(data)
    data.rename('indirect','pulse_length')
    data.labels('pulse_length',
            data.get_prop('L_min')+r_[0:data.get_prop('L_step')]*data.get_prop('L_inc'))
    return data
def postproc_eldor_3d(data):
    if data.get_prop('phasecycle') is None:
        warnings.warn("Warning!  There is no phase cycling information -- you should really fix this")
        data.chunk('phcyc',['phcyc1','phcyc2','phcyc3'])
        #data.rename('indirect','t_e')
        data.ift('phcyc3')
        data.ift('phcyc2')
        data.ift('phcyc1')
        data.rename('indirect','t1')
    else:
        data = automagical_phasecycle(data)
    return data
def postproc_echo_T2(data):
    if data.get_prop('phasecycle') is None:
        warnings.warn("Warning!  There is no phase cycling information -- you should really fix this")
        data.chunk('phcyc',['phcyc2','phcyc1']) # the second pulse is on the outside, and the first pulse is on the inside
        data.setaxis('phcyc1',r_[0:1:4j])
        data.setaxis('phcyc2',r_[0:1:4j])
        data.ift('phcyc1') # remember from 2889 that I have to ift
        data.ift('phcyc2')
    else:
        data = automagical_phasecycle(data)
    return data
def postproc_eldor_old(data,**kwargs):
    data = automagical_phasecycle(data)
    data.rename('indirect','t1')
    t1_min = data.get_prop('t1_min')
    t1_inc = data.get_prop('t1_inc')
    t1_steps = data.get_prop('t1_steps')
    data.labels({'t1':r_[0:t1_steps]*t1_inc+t1_min})
    data.set_units('t1','s')
    return data
def postproc_b1_fid_file(data,fl = None,**kwargs):
    if fl is None:
        fl = figlist_var()
    if isinstance(fl, list) and len(fl) == 2:
        figname_append = fl[1]
        fl = fl[0]
        if figname_append[-1] != '_':
            figname_append += '_'
    else:
        figname_append = ''
    data.rename('indirect','plen')
    if 'L_min' in list(kwargs.keys()) and kwargs['L_min'] is not None:
        L_min = kwargs['L_min']
    else:
        L_min = data.get_prop('L_min')
    in_field_sweep = data.get_prop('in_field_sweep_current') 
    obsn("In field sweep current is %0.6f A"%in_field_sweep)
    if 'L_inc' in list(kwargs.keys()) and kwargs['L_inc'] is not None:
        L_inc = kwargs['L_inc']
    else:
        L_inc = data.get_prop('L_inc')
    if 'L_steps' in list(kwargs.keys()) and kwargs['L_steps'] is not None:
        L_steps = kwargs['L_steps']
    else:
        L_steps = data.get_prop('L_steps')
    data.labels({'plen':r_[0:L_steps]*L_inc+L_min,'t2':r_[0:ndshape(data)['t2']]*1e-9,'phcyc':r_[0:4]})
    fl.next(figname_append+'raw_data')
    fl.image(data)
    fl.next(figname_append+'shifted and phase cycled')
    data.ft('t2',shift = True)
    data *= data.fromaxis(['t2','plen'],lambda x,y: exp(2j*pi*x*(y*0.5)))
    data.ift('t2')
    data.ft('phcyc')
    fl.image(data)
    return data['phcyc',1],fl
def postproc_eldor_file(data,fl = None):
    steps = ndshape(data)['indirect']
    if steps != data.get_prop('t1_steps'):
        raise ValueError('I am confused -- the size of your indirect dimension doesn\'t match the t1_steps parameter')
    x = data.get_prop('t1_min')+data.get_prop('t1_inc')*r_[0:steps]
    data.rename('indirect','t_1')
    data.labels('t_1',x)
    data.set_units('t_1','s')
    return data
def postproc_cw(data,phase = True,use_sweep = False):
    r'this opens the cw data, using search_freed_file and open_cw_file, and then autophases it'
    if use_sweep:
        otherdim = 'current'
    else:
        otherdim = 'field'
    if 'repeats' in data.dimlabels:
        logger.info(strm("found",ndshape(data)['repeats'],"averages"))
        data.mean('repeats')
    #{{{ phase the spectrum
    baseline = (data[otherdim,0:2].mean().data + data[otherdim,-3:-1].mean().data)/2.0
    data -= baseline
    if phase:
        testphases = r_[0:2*pi:100j]
        result = zeros_like(testphases)
        for j,ph in enumerate(testphases):
            test = data.copy()
            test *= exp(-1j*ph)
            test_denom = test.runcopy(imag)
            test_denom.mean(otherdim)
            test.run(real)
            test.mean(otherdim)
            result[j] = abs(test.data)**2/abs(test_denom.data)**2
        data *= exp(-1j*testphases[argmin(result)])
        index = argmax(abs(imag(data.data)))
        data *= -1*sign(imag(data.data[index]))
    if use_sweep:
        data.set_units('current','A')
    else:
        data.set_units('field','T')
    #}}}
    return data
def automagical_phasecycle(data,verbose = False):
    "Use the phase cycle list to determine the phase cycles, and then ift them to return coherence skips"
    logger.info(strm("shape of data",ndshape(data)))
    logger.info(strm(data.other_info['phasecycle']))
    phasecyc_origindeces = r_[0:data.other_info['phasecycle'].shape[0]]
    phase_cycles = data.other_info['phasecycle']
    this_dtype = phase_cycles.dtype.descr * phase_cycles.shape[1]
    #{{{ construct the new phase cycle set
    new_fields = ['phcyc%d'%(j+1) for j in range(phase_cycles.shape[1])]
    phase_cycles =  phase_cycles.view(this_dtype).squeeze()
    phase_cycles.dtype.names = new_fields
    #}}}
    prod_of_nuniqueinfield = array([len(unique(phase_cycles[j])) for j in phase_cycles.dtype.names]).prod()
    logger.info(strm('product of unique elements in each field',prod_of_nuniqueinfield))
    logger.info(strm('actual unique tuples',len(unique(phase_cycles))))
    redundancy = prod_of_nuniqueinfield/len(unique(phase_cycles))
    logger.info(strm('redundancy',redundancy))
    final_field = phase_cycles.dtype.names[-1]
    #logger.info(strm("before rotation",data['indirect',0]))
    if redundancy == 2:
        mask = phase_cycles[final_field] == 90
        mask |= phase_cycles[final_field] == 270
        logger.info(strm("these match",phase_cycles[mask]))
        phase_cycles[final_field][mask] -= 90 
        data['phcyc',mask] *= exp(-1j*pi/2) # pulse was rotated 90 forward, which makes signal 90 backwards
        logger.info(strm("after rotation",data['indirect',0]))
    logger.info(strm("argsort:",argsort(phase_cycles,order = phase_cycles.dtype.names)))
    logger.info(strm("argsorted:"))
    for j in phase_cycles[argsort(phase_cycles)]:
        logger.info(strm(j))
    logger.info(strm(ndshape(data)))
    #{{{ actually assign and chunk the phase cycle dimension
    data.setaxis('phcyc',phase_cycles)
    sorted_indeces = argsort(phase_cycles)
    logger.info(strm('sorted indeces are',sorted_indeces))
    data = data['phcyc',sorted_indeces]
    logger.info(strm(data.getaxis('phcyc')))
    #for phcyc_i in [phase_cycles.dtype.names[0]]: # for debugging, just do the first
    for phcyc_i in phase_cycles.dtype.names[:-1]:
        logger.info(strm('chunking out',phcyc_i))
        logger.info(strm('with axis label',data.getaxis('phcyc')))
        data.chunk_auto('phcyc',phcyc_i,dimname = 'phcyc')
        logger.info(strm('just did',phcyc_i,'and got',ndshape(data)))
    data.rename('phcyc',phase_cycles.dtype.names[-1])
    logger.info(strm('finalized and got',ndshape(data)))
    #}}}
    for j in phase_cycles.dtype.names:
        x = data.getaxis(j)
        data.setaxis(j,r_[0:1:1j*len(x)])
        data.ift(j)
    return data
