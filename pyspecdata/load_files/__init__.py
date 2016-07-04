"the file-loading subroutines are stored here"
import .bruker_nmr as bruker_nmr
import .prospa as prospa
import .bruker_esr as bruker_esr
import .acert as acert
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
        data = bruker_nmr.series(filename,dimname=dimname)
    elif twod and filetype == 'prospa':
        data = prospa.load_2D(filename, dimname=dimname)
        #{{{ bruker 1D
    else:
        if filetype == 'bruker':
            v = bruker_load_acqu(filename)
            td2 = int(v['TD'])
            td1 = 1
            td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
            fp = open(filename+'fid','rb')
            data = fp.read()
            data = array(
                    struct.unpack('>%di'%(len(data)/4),data),
                    dtype='complex128')
            data = data[0::2]+1j*data[1::2]
            rg = bruker_det_rg(v['RG'])
            data /= rg
            data = nddata(data,[td1,td2_zf/2],[dimname,'t2'])
            data = data['t2',0:td2/2] # now, chop out their zero filling
            t2axis = 1./v['SW_h']*r_[1:td2/2+1]
            t1axis = r_[1]
            data.labels([dimname,'t2'],[t1axis,t2axis])
            shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the second order phase shift
            #print 'shiftpoints = ',shiftpoints
            data.circshift('t2',shiftpoints)
            # finally, I will probably need to add in the first order phase shift for the decimation --> just translate this
            data.other_info['title'] = bruker_load_title(filename)
        #}}}
        #{{{ prospa 1d
        elif filetype == 'prospa':
            v = prospa_load_acqu(filename)
            data = prospa_load_datafile(filename,dims=1)/v['nrScans']#added this 2/20/13 to allow automatic signal averaging
            data = nddata(data,[v['nrPnts']],['t2'])
            taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
            data.labels(['t2'],[taxis])
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
        return (data,v,v2)
    else:
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
