from ..general_functions import inside_sphinx
from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like
from pyspecdata import init_logging,strm
logger = init_logging("debug")

def convolve(self,axisname,filterwidth,convfunc = (lambda x,y: exp(-(x**2)/(2.0*(y**2))))):
    r'''Perform a convolution on the data in its
    current domain via multiplication in its
    conjugate domain.
    
    Parameters
    ==========

    axisname: str
        apply the convolution along `axisname`

    filterwidth: double
        width of the convolution function
        (this value must be provided in the same
        domain as that in which the data exists
        when you call this function on said data)

    convfunc: function
        A function that takes two
        arguments -- the first are the axis
        coordinates and the second is
        `filterwidth` (see `filterwidth`).
        Default is a normalized Gaussian of width
        (:math:`\sigma`) `filterwidth`
        :math:`\frac{1}{2 \sigma^2}\exp\left( - \frac{x^2}{2 \sigma^2} \right)`
        For example if you want a complex Lorentzian with `filterwidth` controlled by the rate :math:`R`, 
        *i.e.*
        :math:`\frac{-1}{-i 2 \pi f - R}`
        then ``convfunc = lambda f,R: -1./(-1j*2*pi*f-R)``
    '''
    time_domain = True
    if self.get_ft_prop(axisname):
        # detect self in frequency domain
        logger.info(strm("detect self in frequency domain"))
        self.ift(axisname)
        time_domain = False
    elif self.get_ft_prop(axisname,['start','freq']) is None:
        # detect self in time domain, never FT'd
        logger.info(strm("detect self in time domain, never FT'd"))
        self.ft(axisname, shift=True)
    else:
        # detect self in time domain, already FT'd
        logger.info(strm("detect self in time domain, already FT'd"))
        self.ft(axisname)
    x = self.fromaxis(axisname)
    myfilter = convfunc(x,1/filterwidth)
    newdata = self*myfilter
    self.data = newdata.data
    if time_domain:
        self.ift(axisname)
    else:
        self.ft(axisname)
    logger.info(strm("before return",axisname,self.get_ft_prop(axisname)))
    return self
