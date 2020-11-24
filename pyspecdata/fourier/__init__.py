r'''This package contains the Fourier transformation methods belonging to :class:`nddata`'''
from ..general_functions import inside_sphinx
# I name these all with "fourier" since they need to be renamed inside the main class
from ._ft_conj import _ft_conj
from .convolve import convolve
from .ft import ft
from .ft_shift import _ft_shift
from .ftshift import ftshift
from .ft_shift import set_ft_prop
from .ft_shift import get_ft_prop
from .ft_shift import ft_state_to_str
from .ft_shift import ft_clear_startpoints
from .shear import shear
from .shear import extend_for_shear
from .ift import ift
class nddata_ft(object):
    r"""Provides Fourier transformation.

    Deals correctly with units and
    coordinates under zero-filling, choosing to view different aliases
    (including "shifted" vs. "unshifted"), etc.

    All methods here are inherited by (present in) any :class:`nddata` object.

    :meth:`~pyspecdata.nddata.ft` adjusts normalization and units so that the result conforms to

    :math:`\tilde{s}(f)=\int_{x_{min}}^{x_{max}} s(t) e^{-i 2 \pi f t} dt`

    :meth:`~pyspecdata.nddata.ift` adjusts normalization and units so that the result conforms to

    :math:`s(t)=t_{dw} \int_{x_{min}}^{x_{max}} \tilde{s}(f) e^{i 2 \pi f t} df`

    Where :math:`t_{dw}=\frac{1}{\Delta f}`, is the dwell time (with :math:`\Delta f` the spectral width).

    *Why do we do this?* Note that while the analytical integral this corresponds to is normalized, performing
    :meth:`~pyspecdata.nddata.ft` followed by :meth:`~pyspecdata.nddata.ift` on a discrete sequence is NOT completely invertible
    (due to integration of the implied comb function??),
    and would require division by a factor of :math:`\Delta f` (the spectral width) in order
    to retrieve the original function


    """
    _ft_conj = _ft_conj
    ft = ft
    set_ft_prop = set_ft_prop
    get_ft_prop = get_ft_prop
    ft_state_to_str = ft_state_to_str
    ft_clear_startpoints = ft_clear_startpoints
    ift = ift
    _ft_shift = _ft_shift
    ftshift = ftshift
    convolve = convolve
    extend_for_shear = extend_for_shear
    shear = shear
    def copy_ft_props(self, arg, from_axis=None, to_axis=None):
        """Copy the FT properties from one array to another.


        Parameters
        ==========
        arg: nddata
            copy FT properties from this data
        from_axis: str
            copy info from this axis
        from_axis: str
            copy info to this axis
        """
        list_of_ft_props = [j for j in arg.get_props() if j.startswith('FT')]
        if len(list_of_ft_props)==0: raise ValueError("No FT properties in the argument!")
        logger.debug(strm("these are the FT properties",list_of_ft_props))
        assert to_axis in self.dimlabels
        assert from_axis in arg.dimlabels
        for j in list_of_ft_props:
            assert j in ['FT','FT_start_time','FT_start_freq',
                    'FT_time_not_aliased',
                    'FT_freq_not_aliased'],"I don't understand the (ostensibly FT) property "+str(j)
            self.set_prop(j,{to_axis:arg.get_prop(j)[from_axis]})
        return
__all__ = ["_ft_conj",
        "convolve",
        "ft",
        "_ft_shift",
        "ftshift",
        "set_ft_prop",
        "get_ft_prop",
        "ft_state_to_str",
        "ft_clear_startpoints",
        "shear",
        "extend_for_shear",
        "ift",
        ]

