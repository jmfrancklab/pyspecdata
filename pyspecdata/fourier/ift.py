from ..general_functions import *
from numpy import r_
import numpy as np
from .ft_shift import _find_index,thinkaboutit_message

def ift(self,axes,n=False,tolerance = 1e-5,verbose = False,unitary=None,**kwargs):
    r"""This performs an inverse Fourier transform along the axes identified by the string or list of strings `axes`.

    It adjusts normalization and units so that the result conforms to
            :math:`s(t)=\int_{x_{min}}^{x_{max}} \tilde{s}(f) e^{i 2 \pi f t} df`

    **pre-IFT**, we use the axis to cyclically permute :math:`f=0` to the first index

    **post-IFT**, we assume that the data has previously been FT'd
    If this is the case, passing ``shift=True`` will cause an error
    If this is not the case, passing ``shift=True`` generates a standard ifftshift
    ``shift=None`` will choose True, if and only if this is not the case

    Parameters
    ----------
    pad : int or boolean
        `pad` specifies a zero-filling.  If it's a number, then it gives
        the length of the zero-filled dimension.  If it is just `True`,
        then the size of the dimension is determined by rounding the
        dimension size up to the nearest integral power of 2.   It uses the
        `start_time` ft property to determine the start of the axis.  To
        do this, it assumes that it is a stationary signal
        (convolved with infinite comb function).
        The value of `start_time` can differ from by a non-integral multiple of
        :math:`\Delta t`, though the routine will check whether or not it is safe to
        do this.

        ..note ::
            In the code, this is controlled by `p2_post` (the integral
            :math:`\Delta t` and `p2_post_discrepancy` -- the non-integral.
    unitary : boolean (None)
        return a result that is vector-unitary
    """
    if verbose: print("check 1",self.data.dtype)
    if self.data.dtype == np.float64:
        self.data = np.complex128(self.data) # everything is done assuming complex data
    #{{{ process arguments
    axes = self._possibly_one_axis(axes)
    if (isinstance(axes, str)):
        axes = [axes]
    #{{{ check and set the FT property
    for j in axes:
        if self.get_ft_prop(j) == False:
            errmsg = "This data has been IFT'd along "+str(j)
            raise ValueError(errmsg + "-- you can't IFT"
                    " again unless you explicitly"
                    " .set_ft_prop('"+str(j)+"',None), which is"
                    " probably not what you want to do")
        self.set_ft_prop(j,False) # sets the "FT" property to "false"
    #}}}
    if 'shiftornot' in kwargs:
        raise ValueError("shiftornot is obsolete --> use shift instead")
    shift,pad = process_kwargs([
        ('shift',False),
        ('pad',False),
        ],
        kwargs)
    if not (isinstance(shift, list)):
        shift = [shift]*len(axes)
    if not (isinstance(unitary, list)):
        unitary = [unitary]*len(axes)
    for j in range(0,len(axes)):
        #print("called IFT on",axes[j],", unitary",unitary[j],"and property",
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
        if self.get_ft_prop(axes[j],['start','time']) is None:#  this is the same
            #           as saying that I have NOT run .ft() on this data yet,
            #           meaning that I must have started with frequency as the
            #           source data, and am now constructing an artificial and
            #           possibly aliased time-domain
            if self.get_ft_prop(axes[j],['time','not','aliased']) is not True:#
                #                              has been manually set/overridden
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
            raise RuntimeError(strm("I can't find",axes[j],
                "dimlabels is: ",self.dimlabels))
        padded_length = self.data.shape[thisaxis]
        if pad is True:
            padded_length = int(2**(np.ceil(np.log2(padded_length))))
        elif pad:
            padded_length = pad
        u = self.getaxis(axes[j]) # here, u is frequency
        if u is None:
            raise ValueError("seems to be no axis for"+repr(axes[j])+"set an axis before you try to FT")
        #}}}
        self.set_ft_prop(axes[j],['start','freq'],u[0]) # before anything else, store the start frequency
        #{{{ calculate new axis and post-IFT shift..
        #       Calculate it first, in case it's non-integral.  Also note that
        #       I calculate the post-discrepancy here
        #{{{ need to calculate du and all checks here so I can calculate new u
        du = check_ascending_axis(u,tolerance,"In order to perform FT or IFT")
        self.set_ft_prop(axes[j],['df'],du)
        #}}}
        dv = np.double(1) / du / np.double(padded_length) # so padded length gives the SW
        self.set_ft_prop(axes[j],['dt'],dv)
        v = r_[0:padded_length] * dv # v is the name of the *new* axis.  Note
        #   that we stop one index before the SW, which is what we want
        desired_startpoint = self.get_ft_prop(axes[j],['start','time'])
        if desired_startpoint is not None:# FT_start_time is set
            if shift[j]:
                raise ValueError("you are not allowed to shift an array for"
                        " which the index for $t=0$ has already been"
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
            n = padded_length
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
            asrt_msg = r"""You are trying to shift the time axis by (%d+%g) du (%g).

            In order to shift by a time that is not
            integral w.r.t. the dwell time, you need to be sure
            that the frequency-domain spectrum is not aliased.
            If you **know** that the frequency-domain spectrum
            is not aliased, you can also set the `freq_not_aliased` FT property
            to `True`."""%(p2_post, p2_post_discrepancy, du)
            assert self.get_ft_prop(axes[j],['freq','not','aliased']),(asrt_msg)
            assert abs(p2_post_discrepancy)<abs(dv),("I expect the discrepancy to be"
                    " smaller than dv ({:0.2f}), but it's {:0.2f} -- what's going"
                    " on??").format(dv,p2_post_discrepancy)
            phaseshift =  self.fromaxis(axes[j],
                    lambda q: np.exp(1j*2*pi*q*p2_post_discrepancy))
            try:
                self.data *= phaseshift.data
            except TypeError as e:
                if self.data.dtype != 'complex128':
                    raise TypeError("You tried to ift nddata that was of type "+str(self.data.dtype))
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
        #{{{ pre-IFT shift so that we start at u=0
        p2_pre,p2_pre_discrepancy,alias_shift_pre = _find_index(u,verbose = verbose)
        self._ft_shift(thisaxis,p2_pre)
        #}}}
        #{{{ the actual (I)FFT portion of the routine
        self.data = np.fft.ifft(self.data,
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
            if verbose: print("adjusting axis by",p2_post_discrepancy,"where du is",u[1]-u[0])
            self.axis_coords[thisaxis][:] += p2_post_discrepancy # reflect the
            #   p2_post_discrepancy that we have already incorporated via a
            #   phase-shift above
        #}}}
        #{{{ adjust the normalization appropriately
        if unitary[j]:
            self.data *= np.sqrt(padded_length)
        else:
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
                    lambda f: np.exp(-1j*2*pi*f*p2_pre_discrepancy))
            self.data = result.data
        #}}}
    return self
