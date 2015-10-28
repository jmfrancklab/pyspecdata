from ..general_functions import *
from pylab import * 
from .ft_shift import _find_index,thinkaboutit_message

def ft(self,axes,tolerance = 1e-5,verbose = False,**kwargs):
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
        if verbose: print "check the startf_dict",self.get_prop('FT_start_freq')
        do_post_shift = False
        p2_post_discrepancy = None
        p2_pre_discrepancy = None
        #{{{ grab the axes, set the units, and determine padded_length
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
        #}}}
        #{{{ before anything else, store the start time
        startt_dict = self.get_prop("FT_start_time")
        if startt_dict is None:
            self.set_prop("FT_start_time",{axes[j]:u[0]})
        else:
            startt_dict.update({axes[j]:u[0]})
        #}}}
        #{{{ calculate the post-FT shift -- calculate it first, in case it's non-integral
        if u is None:
            raise ValueError("seems to be no axis for"+repr(axes[j])+"set an axis before you try to FT")
        else:
            #{{{ need to calculate du and all checks here so I can calculate new u
            du = u[1]-u[0] # the dwell gives the bandwidth, whether or not it has been zero padded
            thismsg = "In order to perform FT or IFT, the axis must be equally spaced and ascending"
            assert all(abs(diff(u) - du) < tolerance), thismsg# absolute
            #   tolerance can be large relative to a du of ns -- don't use
            #   allclose/isclose, since they are more recent numpy additions
            assert du > 0, thismsg
            #}}}
            dv = double(1) / du / double(padded_length) # so padded length gives the SW
            v = r_[0:padded_length] * dv # v is the name of the *new* axis.  Note
            #   that we stop one index before the SW, which is what we want
        startf_dict = self.get_prop("FT_start_freq")
        if verbose: print "pulled up the startf_dict",startf_dict,repr(startf_dict)
        if startf_dict is not None and axes[j] in startf_dict.keys():
            if shift[j]:
                raise ValueError("you are not allowed to shift an array for"
                        " which the index for $f=0$ has already been"
                        " determined!")
            if verbose: print "check for p2_post_discrepancy"
            desired_startpoint = startt_dict[axes[j]]
            if verbose: print "desired startpoint",desired_startpoint
            if desired_startpoint < 0: # if it's negative, just grab the replicate starting 1 SW away
                desired_startpoint += 1/du # + SW
            if verbose: print "aliased startpoint",desired_startpoint
            p2_post,p2_post_discrepancy = _find_index(v,origin = desired_startpoint)
            if verbose: print "p2_post,p2_post_discrepancy,v at p2_post, and v at p2_post-1:", p2_post, p2_post_discrepancy, v[p2_post], v[p2_post - 1]
            do_post_shift = True
        elif shift[j]: # a default fftshift
            if automix:
                raise ValueError("You can't use automix and shift at the same time --> it doesn't make sense")
            n = self.data.shape[thisaxis]
            p2_post = (n+1) // 2 # this is the starting index of what starts out as the second half (// is floordiv) -- copied from scipy -- this essentially rounds up (by default assigning more negative frequencies than positive ones)
            do_post_shift = True
        #}}}
        #{{{ I might need to perform a phase-shift now...
        #          in order to adjust for a final u-axis that doesn't pass
        #          exactly through zero
        if p2_post_discrepancy is not None:
            phaseshift =  self.fromaxis(axes[j],
                    lambda q: exp(1j*2*pi*q*p2_post_discrepancy))
            self.data *= phaseshift.data
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
        p2_pre,p2_pre_discrepancy = _find_index(u)
        self._ft_shift(thisaxis,p2_pre)
        #}}}
        self.data = fft(self.data,
                            axis=thisaxis)
        self.axis_coords[thisaxis] = v
        #{{{ actually run the post-FT shift
        if do_post_shift:
            self._ft_shift(thisaxis,p2_post,shift_axis = True)
        #}}}
        #{{{ finally, I must allow for the possibility that "p2_post" in the pre-shift was not
        #    actually at zero, but at some other value, and I must apply
        #    a phase shift to reflect the fact that I need to add back
        #    that time
        if p2_post_discrepancy is not None:
            if verbose: print "adjusting axis by",p2_post_discrepancy,"where du is",u[1]-u[0]
            self.axis_coords[thisaxis][:] += p2_post_discrepancy # reflect the
            #   p2_post_discrepancy that we have already incorporated via a
            #   phase-shift above
        #}}}
        #{{{ adjust the normalization appropriately
        if u is not None:
            self.data *= du # this gives the units in the integral noted in the docstring
        #}}}
        #{{{ finally, if "p2_pre" for the pre-shift didn't correspond exactly to
        #       zero, then the pre-ift data was shifted, and I must reflect
        #       that by performing a post-ift phase shift
        if p2_pre_discrepancy is not None:
            result = self * self.fromaxis(axes[j],
                    lambda f: exp(1j*2*pi*f*p2_pre_discrepancy))
            self.data = result.data
        #}}}
        if automix:
            sw = 1.0/du
            carrier = abs(self).mean_all_but(axes[j]).argmax(axes[j]).data
            if verbose: print "I find carrier at",carrier
            add_to_axis = (automix - carrier) / sw
            if verbose: print "which is",add_to_axis,"times the sw of",sw,"off from the automix value of",automix
            x = self.getaxis(axes[j])
            x += round(add_to_axis)*sw
    return self
