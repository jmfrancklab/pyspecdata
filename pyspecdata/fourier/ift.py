from ..general_functions import *
from pylab import * 
from ft_conj import _ft_shift

def ift(self,axes,**kwargs):
    ("This performs a fourier transform along the axes identified by the string or list of strings `axes`.\n"
    "   It adjusts normalization and units so that the result conforms to\n"
    "   $$s(t)=\int_{x_min}^{x_max} \tilde{s}(t) e^{i 2 \pi f t} df$$\n"
    "   Note that while the analytical integral this corresponds to is
    "normalized, performing .ft() followed by .ift() on a discrete sequence is "
    "NOT completely invertible (due to integration of the implied comb "
    "function??), and would require division by a factor of $\Delta f$ (the "
    "spectral width) in order to retrieve the original function\n"
    )
    #{{{ process arguments
    if len(args) > 1:
        raise ValueError('you can\'t pass more than one argument!!')
    axes = self._possibly_one_axis(axes)
    if (type(axes) is str):
        axes = [axes]
    #{{{ set the FT property
    x = self.get_prop('FT')
    if x is None:
        x = {}
        self.set_prop('FT',x)
    for j in axes:
        x.update({j:False})
    #}}}
    #kwargs: shiftornot=False,shift=None,pad = False
    if 'shiftornot' in kwargs:
        raise ValueError("shiftornot is obsolete --> use shift instead")
    shift,pad = process_kwargs([
        ('shift',None),
        ('pad',False)],
        kwargs)
    if shift != None:
        shiftornot = shift
    if not (type(shiftornot) is list):
        shiftornot = [bool(shiftornot)]*len(axes)
    #}}}
    for j in range(0,len(axes)):
        if self.get_units(axes[j]) is not None:
            self.set_units(axes[j],self._ft_conj(self.get_units(axes[j])))
        try:
            thisaxis = self.dimlabels.index(axes[j])
        except:
            raise CustomError('error, dimlabels is: ',self.dimlabels)
        padded_length = self.data.shape[thisaxis]
        if pad is True:
            padded_length = int(2**(ceil(log2(padded_length))))
        elif pad:
            padded_length = pad
        #{{{ the pre-IFT shift
        # should become obsolete!!
        if bool(shiftornot[j]):
            if automix:
                raise ValueError("You can't use automix and shift at the same time --> it doesn't make sense")
            n = self.data.shape[thisaxis]
            p2 = n - (n+1) // 2 # this is the size of what starts out as the second half // is floordiv -- copied from scipy -- this whole thing essentially rounds down
            self._ft_shift(p2)
        #}}}
        self.data = ifft(self.data,n = padded_length,axis=thisaxis)
        #{{{ the post-IFT shift
        #}}}
        t = self.getaxis(axes[j])
        if t is not None:
            dt = t[1]-t[0]
            self.data *= padded_length * dt # here, the algorithm divides by padded_length, so for integration, we need to not do that
            #{{{ shiftornot specifies the shifting of the initial ft, not this result, so we always return a 0->1 time axis
            self.axis_coords[thisaxis] = linspace(0,1./dt,padded_length) + self.ft_start_time # note that I offset by ft_start_time, which I pull from when I ft'd
            #}}}
    return self
