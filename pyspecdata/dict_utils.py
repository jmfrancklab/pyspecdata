import numpy as np
import logging
from numpy.core import rec
from .general_functions import strm
logger = logging.getLogger('pyspecdata.dict_utils')

#{{{ convert back and forth between lists, etc, and np.ndarray
def make_ndarray(array_to_conv,name_forprint = 'unknown'): 
    if type(array_to_conv) in [int,np.int32,np.double,float,complex,np.complex128,float,bool,np.bool_]: # if it's a scalar
        pass
    elif isinstance(array_to_conv, str):
        pass
    elif type(array_to_conv) in [list,np.ndarray] and len(array_to_conv) > 0:
        array_to_conv = rec.fromarrays([array_to_conv],names = 'LISTELEMENTS') #list(rec.fromarrays([b])['f0']) to convert back
    elif type(array_to_conv) in [list,np.ndarray] and len(array_to_conv) == 0:
        array_to_conv = None
    elif array_to_conv is  None:
        pass
    else:
        raise TypeError(strm('type of value (',type(array_to_conv),') for attribute name',
            name_forprint,'passed to make_ndarray is not currently supported'))
    return array_to_conv
def unmake_ndarray(array_to_conv,name_forprint = 'unknown'): 
    r'Convert this item to an np.ndarray'
    if (isinstance(array_to_conv, np.recarray)) or (isinstance(array_to_conv, np.ndarray) and array_to_conv.dtype.names is not None and len(array_to_conv.dtype.names)>0):
        #{{{ if it's a record/structured np.array, it should be either a list or dictionary
        if 'LISTELEMENTS' in array_to_conv.dtype.names:
            if array_to_conv.dtype.names == tuple(['LISTELEMENTS']):
                retval = list(array_to_conv['LISTELEMENTS'])
            else:
                raise ValueError(strm('Attribute',name_forprint,
                    'is a recordarray with a LISTELEMENTS field, but it',
                    'also has other dimensions:',array_to_conv.dtype.names,
                    'not',tuple(['LISTELEMENTS'])))
        elif len(array_to_conv)==1:
            thisval = dict(list(zip(a.dtype.names,a.tolist()[0])))
        else: raise ValueError('You passed a structured np.array, but it has more',
                "than one dimension, which is not yet supported\nLater, this",
                'should be supported by returning a dictionary of arrays')
        #}}}
    elif isinstance(array_to_conv, np.ndarray) and len(array_to_conv)==1:
        #{{{ if it's a length 1 np.ndarray, then return the element
        retval = array_to_conv.tolist()
        logger.debug(strm(name_forprint,"=",type(array_to_conv),"is a numpy np.array of length one"))
        #}}}
    elif type(array_to_conv) in [np.string_,np.int32,np.float64,np.bool_]:
        #{{{ map numpy strings onto normal strings
        retval = array_to_conv.tolist()
        logger.debug(strm("name_forprint","=",type(array_to_conv),"is a numpy scalar"))
        #}}}
    elif isinstance(array_to_conv, list):
        #{{{ deal with lists
        logger.debug(strm(name_forprint,"is a list"))
        typeofall = list(map(type,array_to_conv))
        if all([x is np.string_ for x in typeofall]):
            logger.debug(strm(name_forprint,"=",typeofall,"are all numpy strings"))
            retval = list(map(str,array_to_conv))
        else:
            logger.debug(strm(name_forprint,"=",typeofall,"are not all numpy string"))
            retval = array_to_conv
        #}}}
    else:
        logger.debug(strm(name_forprint,"=",type(array_to_conv),"is not a numpy string or record np.array"))
        retval = array_to_conv
    return retval
#}}}
