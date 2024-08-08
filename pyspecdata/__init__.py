from .general_functions import *
from .core import *
from .load_files import *
from .figlist import *
from .nnls import *
from .lmfitdata import lmfitdata
from .DCCT_function import DCCT
from .generate_fake_data import fake_data
from .dict_utils import make_ndarray, unmake_ndarray
from .general_functions import (
    CustomError,
    emptytest,
    balance_clims,
    process_kwargs,
    autostringconvert,
    check_ascending_axis,
    level_str_to_int,
    init_logging,
    inside_sphinx,
    strm,
    reformat_exp,
    complex_str,
    render_matrix,
    redim_F_to_C,
    redim_C_to_F,
    fname_makenice,
    explain_error,
    lsafen,
    lsafe,
    copy_maybe_none,
    whereblocks,
    box_muller,
    dp,
    fa,
    ndgrid,
    pinvr,
    sech,
    myfilter,
)

# import numpy

# so essentially, __all__ is the namespace that is passed with an import *
# __all__ = ['prop',
#        'nddata',
#        'figlist_var',
#        'plot',
#        'OLDplot',
#        'nddata_hdf5']
# __all__.extend(numpy.__all__)
__all__ = [x for x in dir() if x[0] != "_"]
