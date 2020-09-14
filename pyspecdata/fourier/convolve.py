from ..general_functions import inside_sphinx
if not inside_sphinx():
    from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like

def convolve(self,axisname,filterwidth,convfunc = (lambda x,y: exp(-(x**2)/(2.0*(y**2))))):
    r'''Perform a convolution.
    
    Parameters
    ==========

    axisname: str
        apply the convolution along `axisname`

    filterwidth: double
        width of the convolution function.

    convfunc: function
        A function that takes two arguments -- the first are the axis coordinates and the second is `filterwidth`.
        Default is a normalized Gaussian of width (:math:`\sigma`)
        `filterwidth`
        :math:`\frac{1}{2 \sigma^2}\exp\left( - \frac{x^2}{2 \sigma^2} \right)`
        For example if you want a complex Lorentzian with `filterwidth` controlled by the rate :math:`R`, 
        *i.e.*
        :math:`\frac{-1}{-i 2 \pi f - R}`
        then ``convfunc = lambda f,R: -1./(-1j*2*pi*f-R)``
    '''
    #{{{ make a version of x that is oriented along the correct dimension
    x = self.getaxis(axisname).copy()
    x_centerpoint = (x[-1]+x[0])/2
    x -= x_centerpoint # so that zero is in the center
    x = ifftshift(x) # so that it's convolved about time 0
    thisaxis = self.axn(axisname)
    #}}}
    myfilter = convfunc(x,filterwidth)
    myfilter /= myfilter.sum()
    filtershape = ones_like(self.data.shape)
    filtershape[thisaxis] = len(myfilter)
    myfilter = myfilter.reshape(filtershape)
    #self.data = ifftshift(ifft(fftshift(fft(self.data,axis = thisaxis),axes = thisaxis)*fftshift(fft(myfilter,axis = thisaxis),axes=thisaxis),axis = thisaxis),axes = thisaxis) # for some reason fftconvolve doesn't work!
    self.data = ifft(fft(self.data,axis = thisaxis)*fft(myfilter,axis = thisaxis),axis = thisaxis)
    #self.data = fftconvolve(self.data,myfilter,mode='same') # I need this, so the noise doesn't break up my blocks
    return self
