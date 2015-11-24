from ..general_functions import *
from pylab import * 

def skew(self,altered_axis,by_amount,propto_axis,zero_fill = True):
    r'''Use the Fourier shift theorem to skew the data:
    ..math: `f(x,y) \rightarrow f(x+ay,y)`


    Parameters
    ----------
    altered_axis : str
        The coordinate for which data is altered, *i.e.*
        ..math: `x` such that ..math: `f(x+ay,y)`.

    by_amount : double
        The amount of the skew (..math: `a` in the previous)

    propto_axis : str
        The shift along the `altered_axis` dimension is
        proportional to the shift along `propto_axis`.
        The position of data relative to the `propto_axis` is not
        changed.
        Note that by the shift theorem, in the frequency domain,
        an equivalent magnitude, opposite sign, skew is applied
        with the `propto_axis` and `altered_axis` dimensions
        flipped.
    '''
    def calc_double_zero_fill(thisaxis,zero_fill,kwargs):
        if zero_fill:
            n = self.data.shape[self.axn(thisaxis)]
            n = int(2**(ceil(log2(n))))
            n *= 2
            kwargs['pad'] = n
    default_kwargs = {'pad':False}
    kwargs = dict(default_kwargs)
    calc_double_zero_fill(altered_axis,zero_fill,kwargs)
    self.ft(altered_axis,verbose = True,**kwargs)
    #phaseshift = self.fromaxis([altered_axis,propto_axis],
    #        lambda x,y: exp(-2j*pi*by_amount*x*y))
    #self.data *= phaseshift.data
    #kwargs = dict(default_kwargs)
    calc_double_zero_fill(altered_axis,zero_fill,kwargs)
    self.ift(altered_axis,**kwargs)
    return self
