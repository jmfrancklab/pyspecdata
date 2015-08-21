from ..general_functions import *
from pylab import * 
from .ft_shift import _find_zero_index,thinkaboutit_message

def ft(self,axes,**kwargs):
    ("This performs a fourier transform along the axes identified by the string or list of strings `axes`.\n"
    "   It adjusts normalization and units so that the result conforms to\n"
    r"   $$\tilde{s}(f)=\int_{x_min}^{x_max} s(t) e^{-i 2 \pi f t} dt$$"+'\n'
    "   Note that while the analytical integral this corresponds to is "
    "normalized, performing .ft() followed by .ift() on a discrete sequence is "
    "NOT completely invertible (due to integration of the implied comb "
    "function??), and would require division by a factor of $\Delta f$ (the "
    "spectral width) in order to retrieve the original function\n"
    "\tpre-FT, we use the axis to cyclically permute $t=0$ to the first "
    "index\n"
    "\t post-FT, we assume that the data has previously been IFT'd\n"
    "\t\tIf this is the case, passing `shift`=True will cause an error\n"
    "\t\tIf this is not the case, passing `shift`=True generates a standard fftshift\n"
    "\t`pad` specifies a zero-filling.  If it's a number, then it gives the"
    " length of the zero-filled dimension.  If it is just `True`, then the size"
    " of the dimension is determined by rounding the dimension size up to the"
    " nearest integral power of 2.\n"
    "\t`automix` can be set to the approximate frequency value.  This is useful"
    " for the specific case where the data has been captured on a sampling scope,"
    " and it's severely aliased over.\n"
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
        x.update({j:True})
    #}}}
    if 'shiftornot' in kwargs:
        raise ValueError("shiftornot is obsolete --> use shift instead")
    shift,pad,automix = process_kwargs([
        ('shift',False),
        ('pad',False),
        ('automix',False)],
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
        u = self.getaxis(axes[j]) # here, u is time
        #{{{ before anything else, store the start time
        startt_dict = self.get_prop("FT_start_time")
        if startt_dict is None:
            self.set_prop("FT_start_time",{axes[j]:u[0]})
        else:
            startt_dict.update({axes[j]:u[0]})
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
        #{{{ the pre-FT shift
        p2 = _find_zero_index(u)
        self._ft_shift(thisaxis,p2)
        #}}}
        self.data = fft(self.data,
                            n = padded_length,
                            axis=thisaxis)
        if u is not None:
            du = u[1]-u[0] # the dwell gives the bandwidth, whether or not it has been zero padded
            thismsg = "In order to perform FT o IFT, the axis must be equally spaced and ascending"
            assert allclose(diff(u),du,atol = 0), thismsg# absolute tolerance can be large relative to a du of ns
            assert du > 0, thismsg
            self.data *= du # this gives the units in the integral noted in the docstring
            self.axis_coords[thisaxis] = linspace(0,1./du,padded_length)
            u = self.axis_coords[thisaxis]
        #{{{ the post-FT shift
        startf_dict = self.get_prop("FT_start_freq")
        if startf_dict is not None and axes[j] in startf_dict.keys():
            if shift[j]:
                raise ValueError("you are not allowed to shift an array for which the index for $f=0$ has already been determined!")
            #{{{ the starting frequency is <0 and aliased over, and I want to shift it to 0
            assert startf_dict[axes[j]] <= 0 , ("Trying to reset to a frequency value greater than"
                        " zero ("+repr(startf_dict[axes[j]])+") which is not"
                        " supported.  "+thinkaboutit_message)
            p2 = argmin(abs(u-(
                        1/du + startf_dict[axes[j]])))
            self._ft_shift(thisaxis,p2,shift_axis = True)
            #}}}
        elif shift[j]:
            if automix:
                raise ValueError("You can't use automix and shift at the same time --> it doesn't make sense")
            n = self.data.shape[thisaxis]
            p2 = (n+1) // 2 # this is the starting index of what starts out as the second half (// is floordiv) -- copied from scipy -- this essentially rounds up (by default assigning more negative frequencies than positive ones)
            self._ft_shift(thisaxis,p2,shift_axis = True)
        #}}}
        if automix:
            sw = 1.0/du
            carrier = abs(self).mean_all_but(axes[j]).argmax(axes[j]).data
            print "I find carrier at",carrier
            add_to_axis = (automix - carrier) / sw
            print "which is",add_to_axis,"times the sw of",sw,"off from the automix value of",automix
            x = self.getaxis(axes[j])
            x += round(add_to_axis)*sw
    return self
