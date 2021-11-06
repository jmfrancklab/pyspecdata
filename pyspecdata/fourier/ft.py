from ..general_functions import *
from numpy import r_
import numpy as np
from .ft_shift import _find_index,thinkaboutit_message

def ft(self,axes,tolerance = 1e-5,cosine=False,verbose = False,unitary=None,**kwargs):
    r"""This performs a Fourier transform along the axes identified by the string or list of strings `axes`.

    It adjusts normalization and units so that the result conforms to
            :math:`\tilde{s}(f)=\int_{x_{min}}^{x_{max}} s(t) e^{-i 2 \pi f t} dt`

    **pre-FT**, we use the axis to cyclically permute :math:`t=0` to the first index

    **post-FT**, we assume that the data has previously been IFT'd
    If this is the case, passing ``shift=True`` will cause an error
    If this is not the case, passing ``shift=True`` generates a standard fftshift
    ``shift=None`` will choose True, if and only if this is not the case

    Parameters
    ----------
    pad : int or boolean
        `pad` specifies a zero-filling.  If it's a number, then it gives
        the length of the zero-filled dimension.  If it is just `True`,
        then the size of the dimension is determined by rounding the
        dimension size up to the nearest integral power of 2.
    automix : double
        `automix` can be set to the approximate frequency value.  This is
        useful for the specific case where the data has been captured on a
        sampling scope, and it's severely aliased over.
    cosine : boolean
        yields a sum of the fft and ifft, for a cosine transform
    unitary : boolean (None)
        return a result that is vector-unitary
    """
    if self.data.dtype == np.float64:
        self.data = np.complex128(self.data) # everything is done assuming complex data
    #{{{ process arguments
    axes = self._possibly_one_axis(axes)
    if (isinstance(axes, str)):
        axes = [axes]
    #{{{ check and set the FT property
    for j in axes:
        if j not in self.dimlabels: raise ValueError("the axis "+j+" doesn't exist")
        if self.get_ft_prop(j):
            errmsg = "This data has been FT'd along "+str(j)
            raise ValueError(errmsg + "-- you can't FT"
                    " again unless you explicitly"
                    " .set_prop('FT',None), which is"
                    " probably not what you want to do")
        self.set_ft_prop(j) # sets the "FT" property to "true"
    #}}}
    if 'shiftornot' in kwargs:
        raise ValueError("shiftornot is obsolete --> use shift instead")
    shift,pad,automix = process_kwargs([
        ('shift',False),
        ('pad',False),
        ('automix',False),
        ],
        kwargs)
    if not (isinstance(shift, list)):
        shift = [shift]*len(axes)
    if not (isinstance(unitary, list)):
        unitary = [unitary]*len(axes)
    for j in range(0,len(axes)):
        #print("called FT on",axes[j],", unitary",unitary[j],"and property",
        #        self.get_ft_prop(axes[j],"unitary"))
        if self.get_ft_prop(axes[j],"unitary") is None: # has not been called
            if unitary[j] is None:
                unitary[j]=False
            self.set_ft_prop(axes[j],"unitary",unitary[j])
        else:
            if unitary[j] is None:
                unitary[j] = self.get_ft_prop(axes[j],"unitary")
            else:
                raise ValueError("Call ft or ift with unitary only the first time, and it will be set thereafter.\nOR if you really want to override mid-way use self.set_ft_prop(axisname,\"unitary\",True/False) before calling ft or ift")
        #print("for",axes[j],"set to",unitary[j])
    #}}}
    for j in range(0,len(axes)):
        do_post_shift = False
        p2_post_discrepancy = None
        p2_pre_discrepancy = None
        #{{{ if this is NOT the source data, I need to mark it as not alias-safe!
        if self.get_ft_prop(axes[j],['start','freq']) is None:#  this is the same
            #           as saying that I have NOT run .ift() on this data yet,
            #           meaning that I must have started with time as the
            #           source data, and am now constructing an artificial and
            #           possibly aliased frequency-domain
            if self.get_ft_prop(axes[j],['freq','not','aliased']) is not True:#
                #                              has been manually set/overridden
                self.set_ft_prop(axes[j],['freq','not','aliased'],False)
            #{{{ but on the other hand, I am taking the time as the
            #    not-aliased "source", so as long as I haven't explicitly set it as
            #    unsafe, declare it "safe"
            if self.get_ft_prop(axes[j],['time','not','aliased']) is not False:
                self.set_ft_prop(axes[j],['time','not','aliased'])
            #}}}
        #}}}
        #{{{ grab the axes, set the units, and determine padded_length
        if self.get_units(axes[j]) is not None:
            self.set_units(axes[j],self._ft_conj(self.get_units(axes[j])))
        try:
            thisaxis = self.dimlabels.index(axes[j])
        except:
            raise RuntimeError(strm("I can't find",axes[j],
                "dimlabels is: ",self.dimlabels))
        padded_length = self.data.shape[thisaxis]
        if pad is True:
            padded_length = int(2**(np.ceil(np.log2(padded_length))))
        elif pad:
            padded_length = pad
        u = self.getaxis(axes[j]) # here, u is time
        if u is None:
            raise ValueError("seems to be no axis for"+repr(axes[j])+"set an axis before you try to FT")
        #}}}
        self.set_ft_prop(axes[j],['start','time'],u[0]) # before anything else, store the start time
        #{{{ calculate new axis and post-IFT shift..
        #       Calculate it first, in case it's non-integral.  Also note that
        #       I calculate the post-discrepancy here
        #{{{ need to calculate du and all checks here so I can calculate new u
        du = check_ascending_axis(u,tolerance,"In order to perform FT or IFT")
        self.set_ft_prop(axes[j],['dt'],du)
        #}}}
        dv = np.double(1) / du / np.double(padded_length) # so padded length gives the SW
        self.set_ft_prop(axes[j],['df'],dv)
        v = r_[0:padded_length] * dv # v is the name of the *new* axis.  Note
        #   that we stop one index before the SW, which is what we want
        desired_startpoint = self.get_ft_prop(axes[j],['start','freq'])
        if desired_startpoint is not None:# FT_start_freq is set
            if shift[j]:
                raise ValueError("you are not allowed to shift an array for"
                        " which the index for $f=0$ has already been"
                        " determined!")
            if verbose: print("check for p2_post_discrepancy")
            if verbose: print("desired startpoint",desired_startpoint)
            p2_post,p2_post_discrepancy,alias_shift_post = _find_index(v,origin = desired_startpoint,verbose = verbose)
            if verbose: print("p2_post,p2_post_discrepancy,alias_shift_post,v at p2_post, and v at p2_post-1:", p2_post, p2_post_discrepancy, alias_shift_post, v[p2_post], v[p2_post - 1])
            if p2_post != 0 or p2_post_discrepancy is not None:
                do_post_shift = True
            else:
                do_post_shift = False
        elif shift[j] or shift[j] is None: # a default fftshift
            if automix:
                raise ValueError("You can't use automix and shift at the same time --> it doesn't make sense")
            n = padded_length
            p2_post = (n+1) // 2 # this is the starting index of what starts out as the second half (// is floordiv) -- copied from scipy -- this essentially rounds up (by default assigning more negative frequencies than positive ones)
            alias_shift_post = 0
            do_post_shift = True
            #{{{ if I start with an alias-safe axis, and perform a
            #    traditional shift, I get an alias-safe axis
            if self.get_ft_prop(axes[j],'time_not_aliased'):
                self.set_ft_prop(axes[j],'freq_not_aliased')
            #}}}
        #}}}
        #{{{ I might need to perform a phase-shift now...
        #          in order to adjust for a final u-axis that doesn't pass
        #          exactly through zero
        if p2_post_discrepancy is not None:
            asrt_msg = r"""You are trying to shift the frequency axis by (%d+%g) du (%g).

            In order to shift by a frequency that is not
            integral w.r.t. the frequency resolution step, you need to be sure
            that the time-domain spectrum is not aliased.
            This is typically achieved by starting from a time domain spectrum and
            generating the frequency domain by an FT.
            If you **know** by other means that the time-domain spectrum
            is not aliased, you can also set the `time_not_aliased` FT property
            to `True`."""%(p2_post, p2_post_discrepancy, du)
            assert self.get_ft_prop(axes[j],['time','not','aliased']),(asrt_msg)
            assert abs(p2_post_discrepancy)<abs(dv),("I expect the discrepancy to be"
                    " smaller than dv ({:0.2f}), but it's {:0.2f} -- what's going"
                    " on??").format(dv,p2_post_discrepancy)
            phaseshift =  self.fromaxis(axes[j],
                    lambda q: np.exp(-1j*2*pi*q*p2_post_discrepancy))
            try:
                self.data *= phaseshift.data
            except TypeError as e:
                if self.data.dtype != 'complex128':
                    raise TypeError("You tried to ft nddata that was of type "+str(self.data.dtype))
                else:
                    raise TypeError(e)
        #}}}
        #{{{ do zero-filling manually and first, so I can properly pre-shift the data
        if not pad is False:
            newdata = list(self.data.shape)
            newdata[thisaxis] = padded_length
            targetslice = [slice(None,None,None)] * len(newdata)
            targetslice[thisaxis] = slice(None,self.data.shape[thisaxis])
            targetslice = tuple(targetslice)
            newdata = np.zeros(newdata,dtype = self.data.dtype)
            newdata[targetslice] = self.data
            self.data = newdata
            u = r_[0:padded_length] * du + u[0]
        #}}}
        #{{{ pre-FT shift so that we start at u=0
        p2_pre,p2_pre_discrepancy,alias_shift_pre = _find_index(u,verbose = verbose)
        self._ft_shift(thisaxis,p2_pre)
        #}}}
        #{{{ the actual (I)FFT portion of the routine
        if cosine:
            self.data = np.fft.fft(self.data,
                    axis=thisaxis) + np.ifft(self.data,
                            axis=thisaxis)
            self.data *= 0.5
        else:
            self.data = np.fft.fft(self.data,
                                axis=thisaxis)
        self.axis_coords[thisaxis] = v
        #}}}
        #{{{ actually run the post-FT shift
        if do_post_shift:
            self._ft_shift(thisaxis,p2_post,shift_axis = True)
            if alias_shift_post != 0:
                self.axis_coords[thisaxis] += alias_shift_post
        #}}}
        #{{{ finally, I must allow for the possibility that "p2_post" in the
        #    pre-shift was not actually at zero, but at some other value, and I
        #    must apply a phase shift to reflect the fact that I need to add
        #    back that time
        if p2_post_discrepancy is not None:
            if verbose: print("adjusting axis by",p2_post_discrepancy,"where du is",u[1]-u[0])
            self.axis_coords[thisaxis][:] += p2_post_discrepancy # reflect the
            #   p2_post_discrepancy that we have already incorporated via a
            #   phase-shift above
        #}}}
        #{{{ adjust the normalization appropriately
        if unitary[j]:
            self.data /= np.sqrt(padded_length)
        else:
            self.data *= du # this gives the units in the integral noted in the docstring
        #}}}
        #{{{ finally, if "p2_pre" for the pre-shift didn't correspond exactly to
        #       zero, then the pre-ft data was shifted, and I must reflect
        #       that by performing a post-ft phase shift
        if p2_pre_discrepancy is not None:
            assert abs(p2_pre_discrepancy)<abs(du) or np.isclose(
                    abs(p2_pre_discrepancy),abs(du)),("I expect the discrepancy to be"
                    " smaller than du ({:0.5g}), but it's {:0.5g} -- what's going"
                    " on??").format(du,p2_pre_discrepancy)
            result = self * self.fromaxis(axes[j],
                    lambda f: np.exp(1j*2*pi*f*p2_pre_discrepancy))
            self.data = result.data
        #}}}
        if automix:
            sw = 1.0/du
            carrier = abs(self).mean_all_but(axes[j]).argmax(axes[j]).data
            if verbose: print("I find carrier at",carrier)
            add_to_axis = (automix - carrier) / sw
            if verbose: print("which is",add_to_axis,"times the sw of",sw,"off from the automix value of",automix)
            x = self.getaxis(axes[j])
            x += np.round(add_to_axis)*sw
    return self
