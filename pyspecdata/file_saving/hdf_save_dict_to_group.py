from ..general_functions import *
from pylab import * 
logger = logging.getLogger('pyspecdata.hdf_save_dict_to_group')

def hdf_save_dict_to_group(group,data):
    '''
    Copied as-is from ACERT hfesr code
    All numpy arrays are datasets.
    '''
    for k,v in data.items():
        if issubclass(type(v),np.ndarray):
            logger.debug('Dataset type %s'%str(v.dtype))
            # Split complex into real and imaginary because Matlab
            # function cannot handle compound data types.
            if np.issubdtype(v.dtype,np.complex):
                logger.debug('Adding %s=%s %s as real and imag datasets' % \
                            (k,str(v.shape),str(v.dtype)))
                group.create_dataset('%s.r'%k,data=v.real)
                group.create_dataset('%s.i'%k,data=v.imag)
            else:
                logger.debug('Adding %s=%s as dataset' % (k,v))
                group.create_dataset(k,data=v,dtype=v.dtype)
        elif issubclass(type(v),dict):
            subgroup=group.create_group(k)
            hdf_save_dict_to_group(subgroup,v)
        else:
            if v is not None and len(v) > 0:
                logger.debug('Adding %s=%s as attribute' % (k,v))
                group.attrs[k]=v
