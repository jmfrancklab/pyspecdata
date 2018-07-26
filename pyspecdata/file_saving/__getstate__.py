from ..general_functions import *
from pylab import * 
from types import MethodType
logger = logging.getLogger('pyspecdata.__getstate__')

def __getstate__(self):
    """This function is called when pickling.  More generally we use it to convert to a dictionary format.
    
    In keeping with the original HDF5 format,
    it returns the following sub-dictionaries:

    :axes:
        A dictionary containing the 
    """
    retval = {}
    all_attributes = dir(self)
    # {{{ process the data and, optionally, the error
    error_processed = False
    if 'data_error' in all_attributes:
        if self.get_error() is not None and len(self.get_error()) > 0:
            if len(self.data.dtype.descr)>1 or len(self.get_error().dtype.descr) > 1:
                raise ValueError("The dtype of your data seems to be"+str(len(self.data.dtype.descr))+", which implies a structured array.\nNot yet supporting data that is a structured array, though this should be fairly easy to implement")
            error_temp = self.get_error()
            # {{{ a simple two-column data structure with data and error
            data_dtype = (
                    [('data',) + self.data.dtype.descr[0][1:]]
                    +[('error',) + self.get_error().dtype.descr[0][1:]])
            # }}}
            retval['data'] = empty(self.data.shape, dtype=data_dtype)
            retval['data']['data'] = self.data
            retval['data']['error'] = self.get_error()
            error_processed = True
        else:
            raise ValueError("found a data_error attribute but it doesn't seem to be sensible")
        all_attributes.remove('data_error')
    elif self.get_error() is not None:
        raise ValueError("get_error is not None, but I didn't find a data_error attribute")
    if not error_processed:
        data_dtype = self.data.dtype.descr
        # just add "data" to the field description, for when the file is saved
        data_dtype = [('data',) + self.data.dtype.descr[0][1:]]
        retval['data'] = empty(self.data.shape, dtype=data_dtype)
        retval['data']['data'] = self.data
        error_processed = True
    all_attributes.remove('data')
    # }}}
    if 'axis_coords' in all_attributes:
        all_attributes.remove('axis_coords')
        retval['axes'] = self.mkd(self.axis_coords)
        logger.debug(strm("the keys of 'axes' at this point are",retval['axes'].keys()))
        if 'axis_coords_error' in all_attributes:
            for this_axis in retval['axes'].keys():
                logger.debug(strm("processing axis",this_axis))
                if self.get_error(this_axis) is not None and len(self.get_error(this_axis)) > 0:
                    if len(retval['axes'][this_axis].dtype.descr)>1 or len(self.get_error(this_axis).dtype.descr) > 1:
                        raise ValueError("The dtype of your data axis for"+this_axis+" seems to be"+str(len(retval['axes'][this_axis].dtype.descr))+", which implies a structured array.\nNot yet supporting data that is a structured array, though this should be fairly easy to implement")
                    error_temp = self.get_error(this_axis)
                    data_temp = retval['axes'][this_axis]
                    # {{{ a simple two-column data structure with data and error
                    data_dtype = (
                            [('data',) + data_temp.dtype.descr[0][1:]]
                            +[('error',) + error_temp.dtype.descr[0][1:]])
                    # }}}
                    retval['axes'][this_axis] = empty(error_temp.shape,
                            dtype=data_dtype)
                    retval['axes'][this_axis]['data'] = data_temp 
                    retval['axes'][this_axis]['error'] = error_temp
            all_attributes.remove('axis_coords_error')
        if 'axis_coords_units' in all_attributes:
            all_attributes.remove('axis_coords_units')
            retval['axes']['units'] = self.mkd(self.axis_coords_units)
    logger.debug(strm("the keys of 'axes' at this point are",retval['axes'].keys()))
    for thisattr in all_attributes:
        if thisattr not in self._nosave and thisattr[0] != '_':
            the_attr = getattr(self,thisattr)
            if type(the_attr) != MethodType:
                logger.debug(strm('converting:',thisattr))
                retval[thisattr] = the_attr
    logger.debug(strm("the keys of 'axes' at this point are",retval['axes'].keys()))
    return retval
