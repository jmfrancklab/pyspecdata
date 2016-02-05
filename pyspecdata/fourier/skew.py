from ..general_functions import *
from pylab import * 

def shear(self,altered_axis,by_amount,propto_axis,zero_fill = False,start_in_conj = False):
    r'''Use the Fourier shift theorem to shear the data:
    ..math: `f(x,y) \rightarrow f(x+ay,y)`

    If `altered_axis` is in the frequency domain, the shift
    will be performed in the time domain, and vice versa.


    Parameters
    ----------

    altered_axis : str

        The coordinate for which data is altered, *i.e.*
        ..math: `x` such that ..math: `f(x+ay,y)`.

    by_amount : double

        The amount of the shear (..math: `a` in the previous)

    propto_axis : str

        The shift along the `altered_axis` dimension is
        proportional to the shift along `propto_axis`.
        The position of data relative to the `propto_axis` is not
        changed.
        Note that by the shift theorem, in the frequency domain,
        an equivalent magnitude, opposite sign, shear is applied
        with the `propto_axis` and `altered_axis` dimensions
        flipped.

    start_in_conj : {False, True}, optional

        Defaults to False

        For efficiency, one can replace a double (I)FT call followed by a
        shear call with a single shear call where `start_in_conj` is set.

        `self` before the call is given in the conjugate domain  (*i.e.*,
        :math:`f` *vs.* :math:`t`) along both dimensions from the one that's
        desired.  This means: (1) `self` after the function call transformed
        into the conjugate domain from that before the call and (2)
        `by_amount`, `altered_axis`, and `propto_axis` all refer to the shear
        in the conjugate domain that the data is in at the end of the
        function call.
    '''
    if zero_fill:
        raise ValueError("zero_fill is put here so that I can avoid aliasing in the domain that I see (by adding another extend), but it's not yet supported")
    #{{{ see if it's in the frequency or time domain
    if self.get_ft_prop(altered_axis) and self.get_ft_prop(propto_axis):
        frequency_domain = True
    elif (self.get_ft_prop(altered_axis) == False) and (self.get_ft_prop(propto_axis) == False):
        frequency_domain = False
    else:
        raise ValueError("In order to shear, both dimensions must be in the same (time vs. frequency) domain")
    if start_in_conj:# then I want to actually specify the shear in the other
        #               domain than the one I start in, and just skip the
        #               beginning
        frequency_domain = not frequency_domain
    #}}}
    if frequency_domain:
        print "entering time domain for",altered_axis
        if not start_in_conj:
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
        if not start_in_conj:
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
