from ..general_functions import *
from pylab import * 

def extend_for_shear(self,b,a,skew_amount,verbose = False):
    "this is a helper function for `.fourier.shear`"
    #{{{ in the time domain, b is the one that's altered (and
    #       needs to be extended), while the shearing is proportional to
    #       -by_amount*a
    if verbose: print "extending to account for the shear along ",b,"by",skew_amount,"which gives lesser and greater expansion amounts of",
    shear_displacement = skew_amount * self.getaxis(a
            )[r_[0,-1]]
    shear_displacement = sort(shear_displacement) # this gives the lesser
    #       and greater, respectively (i.e. le, gt -- not smaller/bigger),
    #       of the two shear displacements, so that I know how to extend.
    #       I need to manually sort, because I don't know if skew_amount is
    #       negative or positive.
    #{{{ actually extend: leave alone if zero.
    if verbose: print " and ".join(map(str,shear_displacement))
    for j in [0,-1]:
        if shear_displacement[j] != 0.:
            if verbose: print ' '.join(map(str,("preparing to extend b",
                    self.getaxis(b)[r_[0,-1]], "to",
                    self.getaxis(b)[j] + shear_displacement[j],
                    "along the", ['lesser','greater'][j],"side of",
                    "b (",self.getaxis(b)[j],")")))
            self.extend(b, self.getaxis(b)[j] +
                    shear_displacement[j])# match greater with greater and
            #       lesser with lesser (doesn't matter to me which side of
            #       a that they came from)
    #}}}
    #}}}
    return self

def shear(self,altered_axis,propto_axis,by_amount,zero_fill = False,start_in_conj = False):
    r'''Use the Fourier shift theorem to shear the data :math:`s`:

    ..math: `s(x',y,z) = s(x+ay,y,z)`

    where :math:`x` is the `altered_axis` and :math:`y` is the
    `propto_axis`.  (Actually typically 2D, but :math:`z` included
    just to illustrate other dimensions that aren't involved)

    This is equivalent to the following in the conjugate domain:

    ..math: `\tilde{s}(f_x,f'_y,z) = \tilde{s}(f_x,f_y-af_x,f_z)`

    Because of this, the algorithm **also** automatically `extend`s the data in `f_y`
    axis.  Equivalently, it increases the resolution (decreases the interval
    between points) in the `propto_axis` dimension.  This prevents aliasing in
    the conjugate domain, which will corrupt the data *w.r.t.* successive
    transformations. *However*, by default (unless you set `zero_fill`), it
    does *not* expand the data in the current domain (*i.e.* frequency *vs.*
    time).  The data in the current domain might alias, but you can see this
    happen.

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
    #{{{ see if it's in the frequency or time domain
    if self.get_ft_prop(altered_axis) and self.get_ft_prop(propto_axis):
        frequency_domain = True
    elif not self.get_ft_prop(altered_axis) and not self.get_ft_prop(propto_axis):
        frequency_domain = False
    else:
        raise ValueError("In order to shear, both dimensions must be in the same (time vs. frequency) domain.  Currently, they are {:s} for {:s} and {:s} for {:s}.".format(
            self.get_ft_prop(altered_axis),altered_axis,
            self.get_ft_prop(propto_axis),propto_axis
            ))
    if start_in_conj:# then I want to actually specify the shear in the other
        #               domain than the one I start in, and just skip the
        #               beginning
        frequency_domain = not frequency_domain
    #}}}
    if frequency_domain:
        print "entering time domain for",altered_axis
        if not start_in_conj:
            if zero_fill:
                self.extend_for_shear(altered_axis,propto_axis,by_amount)
            self.ift(altered_axis)
            self.ift(propto_axis) # before expansion
        else:
            if zero_fill:
                raise ValueError("I can't zero fill  because you chose to start in the conjugate dimension")
        print "conjugate domain extension:"
        self.extend_for_shear(propto_axis,altered_axis,-by_amount) # in
        #       the time domain, propto_axis is the one that's altered
        #       (and needs to be extended), while the shearing is
        #       proportional to -by_amount*altered_axis
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
            if zero_fill:
                self.extend_for_shear(altered_axis,propto_axis,by_amount)
            self.ft(altered_axis)
            self.ft(propto_axis) # before expansion
        else:
            if zero_fill:
                raise ValueError("I can't zero fill  because you chose to start in the conjugate dimension")
        print "conjugate domain extension:"
        self.extend_for_shear(propto_axis,altered_axis,-by_amount) # in
        #       the time domain, propto_axis is the one that's altered
        #       (and needs to be extended), while the shearing is
        #       proportional to -by_amount*altered_axis
        self.ift(propto_axis) # after expansion
        print "applying phase shift"
        phaseshift = self.fromaxis([altered_axis,propto_axis],
                lambda x,y: exp(-2j*pi*by_amount*x*y))
        self.data *= phaseshift.data
        print "back to frequency domain"
        self.ift(altered_axis)
    return self
