from .general_functions import (
    CustomError,
    autostringconvert,
    balance_clims,
    box_muller,
    check_ascending_axis,
    complex_str,
    copy_maybe_none,
    det_unit_prefactor,
    dp,
    emptytest,
    fa,
    fname_makenice,
    init_logging,
    inside_sphinx,
    level_str_to_int,
    lsafe,
    lsafen,
    myfilter,
    ndgrid,
    pinvr,
    process_kwargs,
    redim_C_to_F,
    redim_F_to_C,
    reformat_exp,
    render_matrix,
    sech,
    strm,
    whereblocks,
)
from .plot_funcs.DCCT_function import DCCT
from .plot_funcs.image import image
from .core import *
from .load_files import *
from .figlist import *
from .nnls import *
from .lmfitdata import lmfitdata
from .lmfitdataGUI import lmfitdataGUI
from .generate_fake_data import fake_data
from .dict_utils import make_ndarray, unmake_ndarray
from .load_files.zenodo import zenodo_upload, create_deposition
from .datadir import getDATADIR, log_fname, proc_data_target_dir
from .mpl_utils import (
    plot_label_points,
    figlistret,
    figlistini_old,
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
