from ..general_functions import inside_sphinx
from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like,sqrt,pi
import numpy as np
from pyspecdata import init_logging,strm
import logging

def convolve(self,axisname,filterwidth,convfunc='gaussian',
    enforce_causality=True
    ):
    r'''Perform a convolution.
    
    Parameters
    ==========

    axisname: str
        apply the convolution along `axisname`

    filterwidth: double
        width of the convolution function
        (the units of this value are specified in the
        same domain as that in which the data exists
        when you call this function on said data)

    convfunc: function
        A function that takes two
        arguments -- the first are the axis
        coordinates and the second is
        `filterwidth` (see `filterwidth`).
        Default is a normalized Gaussian of FWHM
        (:math:`\lambda`) `filterwidth`
        For example if you want a complex Lorentzian with `filterwidth` controlled by the rate :math:`R`, 
        *i.e.*
        :math:`\frac{-1}{-i 2 \pi f - R}`
        then ``convfunc = lambda f,R: -1./(-1j*2*pi*f-R)``
    enforce_causality: boolean (default true)
        make sure that the ift of the filter doesn't get aliased to high
        time values.

        "Causal" data here means data derived as the FT of time-domain
        data that starts at time zero -- like an FID -- for which real
        and abs parts are Hermite transform pairs.

        `enforce_causality` should be `True` for frequency-domain data
        whose corresponding time-domain data has a startpoint at or near
        zero, with no negative time values -- like data derived from the
        FT of an IFT.
        In contrast, for example, if you have frequency-domain data that
        is entirely real (like a power spectral density) then you want to
        set enforce_causality to False.

        It is ignored if you call a convolution on time-domain data.
    '''
    if convfunc == 'gaussian':
        def convfunc(x,filterwidth):
            sigma = filterwidth/(2*np.sqrt(2*np.log(2)))
            retval = np.exp(-x**2/2/sigma**2)
            retval /= sigma*np.sqrt(2*pi)
            return retval
    time_domain = True
    x = self.fromaxis(axisname)
    myfilter = convfunc(x,filterwidth)
    if self.get_ft_prop(axisname):
        # detect self in frequency domain
        logging.debug(strm("detect self in frequency domain"))
        self.ift(axisname)
        myfilter.ift(axisname)
        if enforce_causality:
            # this is horribly inefficient, but the most straightforward way I
            # can think of doing this -- would be greatly helped if we could
            # "FT" the "axis object" instead
            orig_start = myfilter.getaxis(axisname)[0]
            tlen = myfilter.getaxis(axisname)[-1]-orig_start
            assert orig_start-2*tlen < 0, "the time origin (%f vs tlen %f) is at too high of a positive value for me to try to enforce causality"%(orig_start,tlen)
            # create a negative part for it to alias into:
            myfilter.extend(axisname,-myfilter.getaxis(axisname)[-1])
            myfilter.ft(axisname)
            # recalculate the filter on the new axis
            x = myfilter.fromaxis(axisname)
            myfilter = convfunc(x,filterwidth)
            myfilter.ift(axisname)
            # throw out the stuff to the left of the original time origin
            myfilter = myfilter[axisname:(orig_start,None)]
        time_domain = False
    elif self.get_ft_prop(axisname,['start','freq']) is None:
        # detect self in time domain, never FT'd
        logging.debug(strm("detect self in time domain, never FT'd"))
        self.ft(axisname, shift=True)
        # here enforcing causality wouldn't have the same meaning, so we don't do it
        myfilter.ft(axisname)
    else:
        # detect self in time domain, already FT'd
        logging.debug(strm("detect self in time domain, already FT'd"))
        self.ft(axisname)
        myfilter.ft(axisname)
    newdata = self*myfilter
    self.data = newdata.data
    if time_domain:
        self.ift(axisname)
    else:
        self.ft(axisname)
    logging.debug(strm("before return",axisname,self.get_ft_prop(axisname)))
    return self
