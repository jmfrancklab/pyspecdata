from ..general_functions import *
from pylab import * 
from .ft_shift import _find_index,thinkaboutit_message

def ift(self,axes,tolerance = 1e-5,verbose = False,**kwargs):
    ("This performs a fourier transform along the axes identified by the string"
    " or list of strings `axes`.\n"
    "It adjusts normalization and units so that the result conforms to\n\t"
    r"$$s(t)=\int_{x_min}^{x_max} \tilde{s}(t) e^{i 2 \pi f t} df$$"+'\n'
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
    if verbose: print "check the startt_dict",self.get_prop('FT_start_time')
    if verbose: print "check 1",self.data.dtype
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
        if verbose: print "check the startt_dict",self.get_prop('FT_start_time')
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
        u = self.getaxis(axes[j]) # here, u is frequency
        #}}}
        #{{{ before anything else, store the start frequency
        startf_dict = self.get_prop("FT_start_freq")
        if startf_dict is None:
            self.set_prop("FT_start_freq",{axes[j]:u[0]})
        else:
            startf_dict.update({axes[j]:u[0]})
        #}}}
        #{{{ calculate new axis and post-IFT shift..
        #       Calculate it first, in case it's non-integral.  Also note that
        #       I calculate the post-discrepancy here
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
            v = linspace(0,1./du,padded_length) # v is the name of the *new*
            #   axis, which is only assigned right before the shift, below
        startt_dict = self.get_prop("FT_start_time")
        if verbose: print "pulled up the startt_dict",startt_dict,repr(startt_dict)
        if startt_dict is not None and axes[j] in startt_dict.keys():# FT_start_time is set
            if shift[j]:
                raise ValueError("you are not allowed to shift an array for"
                        " which the index for $t=0$ has already been"
                        " determined!")
            if verbose: print "check for p2_post_discrepancy"
            desired_startpoint = startt_dict[axes[j]]
            if verbose: print "desired startpoint",desired_startpoint
            if desired_startpoint < 0:
                desired_startpoint += 1/du
            if verbose: print "aliased startpoint",desired_startpoint
            p2,p2_post_discrepancy = _find_index(v,origin = desired_startpoint)
            if verbose: print p2,p2_post_discrepancy
            do_post_shift = True
        elif shift[j]: # a default fftshift
            n = self.data.shape[thisaxis]
            p2 = n - (n+1) // 2 # this is the size of what starts out as the second half // is floordiv -- copied from scipy -- this whole thing essentially rounds down
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
        #{{{ do zero-filling manually and first, so I can properly pre-shift the data
        if not pad is False:
            newdata = list(self.data.shape)
            newdata[thisaxis] = padded_length
            targetslice = [slice(None,None,None)] * len(newdata)
            targetslice[thisaxis] = slice(None,self.data.shape[thisaxis])
            newdata = zeros(newdata,dtype = self.data.dtype)
            newdata[targetslice] = self.data
            self.data = newdata
        #}}}
        #{{{ pre-IFT shift so that we start at u=0
        if verbose: print "check for p2_pre_discrepancy"
        p2,p2_pre_discrepancy = _find_index(u)
        self._ft_shift(thisaxis,p2)
        #}}}
        self.data = ifft(self.data,
                            n = padded_length,
                            axis=thisaxis)
        self.axis_coords[thisaxis] = v
        #{{{ actually run the post-IFT shift
        if do_post_shift:
            self._ft_shift(thisaxis,p2,shift_axis = True)
        if p2_post_discrepancy is not None:
            if verbose: print "adjusting axis by",p2_post_discrepancy,"where du is",u[1]-u[0]
            self.axis_coords[thisaxis][:] += p2_post_discrepancy # reflect the
            #   p2_post_discrepancy that we have already incorporated via a
            #   phase-shift above
        #}}}
        #{{{ adjust the normalization appropriately
        self.data *= padded_length * du # here, the algorithm divides by
        #       padded_length, so for integration, we need to not do that
        #}}}
        #{{{ finally, if "p2" for the pre-shift didn't correspond exactly to
        #       zero, then the pre-ift data was shifted, and I must reflect
        #       that by performing a post-ift phase shift
        if p2_pre_discrepancy is not None:
            result = self * self.fromaxis(axes[j],
                    lambda f: exp(1j*2*pi*f*p2_pre_discrepancy))
            self.data = result.data
        #}}}
    return self
