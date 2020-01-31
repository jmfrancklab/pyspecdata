from .general_functions import *
from .core import *
from .load_files import *
from .figlist import *
from .nnls import *

#import numpy

# so essentially, __all__ is the namespace that is passed with an import *
#__all__ = ['prop',
#        'nddata',
#        'figlist_var',
#        'plot',
#        'OLDplot',
#        'nddata_hdf5']
#__all__.extend(numpy.__all__)
__all__ = [x for x in dir() if x[0] != '_']
