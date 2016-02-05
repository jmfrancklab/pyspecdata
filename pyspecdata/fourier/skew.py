from ..general_functions import *
from pylab import * 

def skew(self,altered_axis,by_amount,propto_axis,zero_fill = True):
    r'''Use the Fourier shift theorem to skew the data:
    ..math: `f(x,y) \rightarrow f(x+ay,y)`

    If `altered_axis` is in the frequency domain, the shift
    will be performed in the time domain, and vice versa.


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
    if self.get_ft_prop(altered_axis) and self.get_ft_prop(propto_axis):
        frequency_domain = True
    elif (self.get_ft_prop(altered_axis) == False) and (self.get_ft_prop(propto_axis) == False):
        frequency_domain = False
    else:
        raise ValueError("In order to skew, both dimensions must be in the same (time vs. frequency) domain")
    if frequency_domain:
        print "entering time domain for",altered_axis
        self.ift(altered_axis)
        self.ift(propto_axis) # before expansion
        self.extend_to_match(altered_axis,propto_axis,by_amount)
        self.ft(propto_axis) # after expansion
        print "applying phase shift"
        phaseshift = self.fromaxis([altered_axis,propto_axis],
                lambda x,y: exp(2j*pi*by_amount*x*y))
        self.data *= phaseshift.data
        print "back to frequency domain"
        self.ft(altered_axis)
    else:
        print "entering time domain for",altered_axis
        self.ft(altered_axis)
        self.ft(propto_axis) # before expansion
        self.extend_to_match(altered_axis,propto_axis,by_amount)
        self.ift(propto_axis) # after expansion
        print "applying phase shift"
        phaseshift = self.fromaxis([altered_axis,propto_axis],
                lambda x,y: exp(-2j*pi*by_amount*x*y))
        self.data *= phaseshift.data
        print "back to frequency domain"
        self.ift(altered_axis)
    return self
