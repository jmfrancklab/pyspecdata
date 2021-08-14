from ..general_functions import inside_sphinx
from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like,sqrt,pi
from pyspecdata import init_logging,strm
from ..ndshape import ndshape_base as ndshape
logger = init_logging("debug")

def convolve(self,axisname,filterwidth,convfunc = (lambda x,y:
    exp(-(x**2)/(2.0*(y**2)))/(y*sqrt(2*pi))),
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
        Default is a normalized Gaussian of width
        (:math:`\sigma`) `filterwidth`
        :math:`\frac{1}{2 \sigma^2}\exp\left( - \frac{x^2}{2 \sigma^2} \right)`
        For example if you want a complex Lorentzian with `filterwidth` controlled by the rate :math:`R`, 
        *i.e.*
        :math:`\frac{-1}{-i 2 \pi f - R}`
        then ``convfunc = lambda f,R: -1./(-1j*2*pi*f-R)``
    enforce_causality: boolean
        make sure that the ift of the filter doesn't get aliased to high
        time values.
        This only makes sense for frequency-domain data whose
        corresponding time-domain data has a startpoint at or near zero.
        It is ignored if you call a convolution on time-domain data.
    '''
    time_domain = True
    x = self.fromaxis(axisname)
    myfilter = convfunc(x,filterwidth)
    if self.get_ft_prop(axisname):
        # detect self in frequency domain
        logger.debug(strm("detect self in frequency domain"))
        self.ift(axisname)
        myfilter.ift(axisname)
        print("orig ndshape",ndshape(myfilter))
        if enforce_causality:
            # this is horribly inefficient, but the most straightforward way I
            # can think of doing this -- would be greatly helped if we could
            # "FT" the "axis object" instead
            orig_start = myfilter.getaxis(axisname)[0]
            tlen = myfilter.getaxis(axisname)[-1]-orig_start
            assert orig_start-2*tlen < 0, "the time origin (%f vs tlen %f) is at too high of a positive value for me to try to enforce causality"%(orig_start,tlen)
            # create a negative part for it to alias into:
            print("extending to",-myfilter.getaxis(axisname)[-1])
            myfilter.extend(axisname,-myfilter.getaxis(axisname)[-1])
            print("start time",myfilter.get_ft_prop(axisname,
                    ['start','time']))
            myfilter.set_ft_prop(axisname,
                    ['start','time'],
                    -myfilter.getaxis(axisname)[-1])
            print("start time",myfilter.get_ft_prop(axisname,
                    ['start','time']))
            print(myfilter.getaxis(axisname)[r_[0,-1]])
            myfilter.ft(axisname)
            print("start time",myfilter.get_ft_prop(axisname,
                    ['start','time']))
            # recalculate the filter on the new axis
            x = myfilter.fromaxis(axisname)
            print("x start time",x.get_ft_prop(axisname,
                    ['start','time']))
            myfilter = convfunc(x,filterwidth)
            myfilter.ift(axisname)
            print("start time",myfilter.get_ft_prop(axisname,
                    ['start','time']))
            # throw out the stuff to the left of the original time origin
            print("ndshape before slice",ndshape(myfilter),myfilter.getaxis(axisname)[r_[0,-1]])
            myfilter = myfilter[axisname:(orig_start,None)]
            print("ndshape after slice",ndshape(myfilter),myfilter.getaxis(axisname)[r_[0,-1]])
        time_domain = False
    elif self.get_ft_prop(axisname,['start','freq']) is None:
        # detect self in time domain, never FT'd
        logger.debug(strm("detect self in time domain, never FT'd"))
        self.ft(axisname, shift=True)
        # here enforcing causality wouldn't have the same meaning, so we don't do it
        myfilter.ft(axisname)
    else:
        # detect self in time domain, already FT'd
        logger.debug(strm("detect self in time domain, already FT'd"))
        self.ft(axisname)
        myfilter.ft(axisname)
    newdata = self*myfilter
    self.data = newdata.data
    if time_domain:
        self.ift(axisname)
    else:
        self.ft(axisname)
    logger.debug(strm("before return",axisname,self.get_ft_prop(axisname)))
    return self
