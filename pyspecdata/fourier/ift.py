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
    "\t\t`shift`=None will choose True, if and only if this is not the case\n"
    "\t`pad` specifies a zero-filling.  If it's a number, then it gives the"
    " length of the zero-filled dimension.  If it is just `True`, then the size"
    " of the dimension is determined by rounding the dimension size up to the"
    " nearest integral power of 2."
    "\tIt uses the `start_time` ft property to determine the start of the axis.  To"
    " do this, it assumes that it is a stationary signal (convolved with infinite comb"
    " function).  The value of `start_time` can differ from by a non-integral multiple"
    " of $\Delta t$, though the routine will check whether or not it is safe to do"
    " this."
    )
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
        shift = [shift]*len(axes)
    #}}}
    for j in range(0,len(axes)):
        do_post_shift = False
        p2_post_discrepancy = None
        p2_pre_discrepancy = None
        #{{{ if this is NOT the source data, I need to mark it as not alias-safe!
        if self.get_ft_prop(axes[j],['start','time']) is None:#  this is the same
            #           as saying that I have NOT run .ft() on this data yet,
            #           meaning that I must have started with frequency as the
            #           source data, and am now constructing an artificial and
            #           possibly aliased time-domain
            self.set_ft_prop(axes[j],['time','not','aliased'],False)
            #{{{ but on the other hand, I am taking the frequency as the
            #    not-aliased "source", so as long as I haven't explicitly set it as
            #    unsafe, declare it "safe"
            if self.get_ft_prop(axes[j],['freq','not','aliased']) is not False:
                self.set_ft_prop(axes[j],['freq','not','aliased'])
            #}}}
        #}}}
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
        self.set_ft_prop(axes[j],['start','freq'],u[0]) # before anything else, store the start frequency
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
            dv = double(1) / du / double(padded_length) # so padded length gives the SW
            v = r_[0:padded_length] * dv # v is the name of the *new* axis.  Note
            #   that we stop one index before the SW, which is what we want
        desired_startpoint = self.get_ft_prop(axes[j],['start','time'])
        if desired_startpoint is not None:# FT_start_time is set
            if shift[j]:
                raise ValueError("you are not allowed to shift an array for"
                        " which the index for $t=0$ has already been"
                        " determined!")
            if verbose: print "check for p2_post_discrepancy"
            if verbose: print "desired startpoint",desired_startpoint
            p2_post,p2_post_discrepancy,alias_shift_post = _find_index(v,origin = desired_startpoint)
            if verbose: print "p2_post,p2_post_discrepancy,v at p2_post, and v at p2_post-1:", p2_post, p2_post_discrepancy, v[p2_post], v[p2_post - 1]
            if p2_post != 0 or p2_post_discrepancy is not None:
                do_post_shift = True
            else:
                do_post_shift = False
        elif shift[j] or shift[j] is None: # a default fftshift
            n = self.data.shape[thisaxis]
            p2_post = n - (n+1) // 2 # this is the size of what starts out as the second half // is floordiv -- copied from scipy -- this whole thing essentially rounds down
            alias_shift_post = 0
            do_post_shift = True
            #{{{ if I start with an alias-safe axis, and perform a
            #    traditional shift, I get an alias-safe axis
            if self.get_ft_prop(axes[j],'freq_not_aliased'):
                self.set_ft_prop(axes[j],'time_not_aliased')
            #}}}
        #}}}
        #{{{ I might need to perform a phase-shift now...
        #          in order to adjust for a final u-axis that doesn't pass
        #          exactly through zero
        if p2_post_discrepancy is not None:
            assert self.get_ft_prop(axes[j],['freq','not','aliased']),("in order to"
                " shift by a time that is not integral w.r.t. the dwell time, you need"
                " to be sure that the frequency-domain spectrum is not aliased.  This"
                " is typically achieved by starting from a time domain spectrum and"
                " generating the frequency domain by an FT.  If you **know** by other"
                " means that the frequency-domain spectrum is not aliased, you can also"
                " set the `freq_not_aliased` FT property to `True`")
            assert abs(p2_post_discrepancy)<abs(dv),("I expect the discrepancy to be"
                    " smaller than dv ({:0.2f}), but it's {:0.2f} -- what's going"
                    " on??").format(dv,p2_post_discrepancy)
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
            u = r_[0:padded_length] * du + u[0]
        #}}}
        #{{{ pre-IFT shift so that we start at u=0
        p2_pre,p2_pre_discrepancy,alias_shift_pre = _find_index(u)
        self._ft_shift(thisaxis,p2_pre)
        #}}}
        #{{{ the actual (I)FFT portion of the routine
        self.data = ifft(self.data,
                            axis=thisaxis)
        self.axis_coords[thisaxis] = v
        #}}}
        #{{{ actually run the post-IFT shift
        if do_post_shift:
            self._ft_shift(thisaxis,p2_post,shift_axis = True)
            if alias_shift_post != 0:
                self.axis_coords[thisaxis] += alias_shift_post
        #}}}
        #{{{ finally, I must allow for the possibility that "p2_post" in the
        #    pre-shift was not actually at zero, but at some other value, and I
        #    must apply a phase shift to reflect the fact that I need to add
        #    back that frequency
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
        #{{{ finally, if "p2_pre" for the pre-shift didn't correspond exactly to
        #       zero, then the pre-ift data was shifted, and I must reflect
        #       that by performing a post-ift phase shift
        if p2_pre_discrepancy is not None:
            assert abs(p2_pre_discrepancy)<abs(du),("I expect the discrepancy to be"
                    " smaller than du ({:0.2f}), but it's {:0.2f} -- what's going"
                    " on??").format(du,p2_pre_discrepancy)
            result = self * self.fromaxis(axes[j],
                    lambda f: exp(-1j*2*pi*f*p2_pre_discrepancy))
            self.data = result.data
        #}}}
    return self
