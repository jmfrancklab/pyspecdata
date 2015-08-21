from ..general_functions import *
from pylab import * 
from .ft_shift import _find_zero_index,thinkaboutit_message

def ift(self,axes,**kwargs):
    ("This performs a fourier transform along the axes identified by the string or list of strings `axes`.\n"
    "   It adjusts normalization and units so that the result conforms to\n"
    r"   $$s(t)=\int_{x_min}^{x_max} \tilde{s}(t) e^{i 2 \pi f t} df$$"+'\n'
    "   Note that while the analytical integral this corresponds to is "
    "normalized, performing .ft() followed by .ift() on a discrete sequence is "
    "NOT completely invertible (due to integration of the implied comb "
    "function??), and would require division by a factor of $\Delta f$ (the "
    "spectral width) in order to retrieve the original function\n"
    "\tpre-IFT, we use the axis to cyclically permute $f=0$ to the first "
    "index\n"
    "\t post-IFT, we assume that the data has previously been FT'd\n"
    "\t\tIf this is the case, passing `shift`=True will cause an error\n"
    "\t\tIf this is not the case, passing `shift`=True generates a standard ifftshift\n"
    "\t`pad` specifies a zero-filling.  If it's a number, then it gives the"
    " length of the zero-filled dimension.  If it is just `True`, then the size"
    " of the dimension is determined by rounding the dimension size up to the"
    " nearest integral power of 2."
    )
    #{{{ process arguments
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
    if 'shiftornot' in kwargs:
        raise ValueError("shiftornot is obsolete --> use shift instead")
    shift,pad = process_kwargs([
        ('shift',False),
        ('pad',False)],
        kwargs)
    if not (type(shift) is list):
        shift = [bool(shift)]*len(axes)
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
        u = self.getaxis(axes[j]) # here, u is frequency
        #{{{ before anything else, store the start frequency
        startf_dict = self.get_prop("FT_start_freq")
        if startf_dict is None:
            self.set_prop("FT_start_freq",{axes[j]:u[0]})
        else:
            startf_dict.update({axes[j]:u[0]})
        #}}}
        #{{{ need to do the zero-filling manually, so I can properly pre-shift the data
        if not pad is False:
            newdata = list(self.data.shape)
            newdata[thisaxis] = padded_length
            targetslice = [slice(None,None,None)] * len(newdata)
            targetslice[thisaxis] = slice(None,self.data.shape[thisaxis])
            newdata = zeros(newdata,dtype = self.data.dtype)
            newdata[targetslice] = self.data
            self.data = newdata
        #}}}
        #{{{ the pre-IFT shift
        p2 = _find_zero_index(u)
        self._ft_shift(thisaxis,p2)
        #}}}
        self.data = ifft(self.data,
                            n = padded_length,
                            axis=thisaxis)
        if u is not None:
            du = u[1]-u[0] # the dwell gives the bandwidth, whether or not it has been zero padded
            thismsg = "In order to perform FT o IFT, the axis must be equally spaced and ascending"
            assert allclose(diff(u),du,atol = 0), thismsg# absolute tolerance can be large relative to a du of ns
            assert du > 0, thismsg
            self.data *= padded_length * du # here, the algorithm divides by padded_length, so for integration, we need to not do that
            self.axis_coords[thisaxis] = linspace(0,1./du,padded_length)
            u = self.axis_coords[thisaxis]
        #{{{ the post-IFT shift
        startt_dict = self.get_prop("FT_start_time")
        if startt_dict is not None and axes[j] in startt_dict.keys():
            if shift[j]:
                raise ValueError("you are not allowed to shift an array for which the index for $t=0$ has already been determined!")
            #{{{ the starting time is <0 and aliased over, and I want to shift it to 0
            assert startt_dict[axes[j]] <= 0 , ("Trying to reset to a time value greater than"
                        " zero ("+repr(startt_dict[axes[j]])+") which is not"
                        " supported.  "+thinkaboutit_message)
            p2 = argmin(abs(u-(
                        1/du + startt_dict[axes[j]])))
            self._ft_shift(thisaxis,p2,shift_axis = True)
            #}}}
        elif shift[j]:
            n = self.data.shape[thisaxis]
            p2 = n - (n+1) // 2 # this is the size of what starts out as the second half // is floordiv -- copied from scipy -- this whole thing essentially rounds down
            self._ft_shift(thisaxis,p2,shift_axis = True)
        #}}}
    return self
