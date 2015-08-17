from ..general_functions import *
from pylab import * 

def ift(self,*args,**kwargs):
    #{{{ process arguments
    if len(args) > 1:
        raise ValueError('you can\'t pass more than one argument!!')
    axes = self._possibly_one_axis(*args)
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
    shiftornot,shift,pad = process_kwargs([
        ('shiftornot',False),
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
        if bool(shiftornot[j]):
            newdata = list(shape(self.data))
            newdata[thisaxis] = padded_length
            newdata = zeros(tuple(newdata),dtype = self.data.dtype)
            n = self.data.shape[thisaxis]
            p2 = n - (n+1) // 2 # floordiv -- copied from scipy -- this essentially rounds up
            sourceslice = [slice(None,None,None)] * len(self.data.shape)
            targetslice = [slice(None,None,None)] * len(self.data.shape)
            # move second half first -- the following are analogous to the numpy function, but uses slices instead
            sourceslice[thisaxis] = slice(p2,n)
            targetslice[thisaxis] = slice(None,n-p2)
            newdata[targetslice]  = self.data[sourceslice]
            # move first half second (the negative frequencies)
            sourceslice[thisaxis] = slice(None,p2)
            targetslice[thisaxis] = slice(-p2,None)
            newdata[targetslice]  = self.data[sourceslice]
            self.data = newdata
            #self.data = ifftshift(self.data,axes=[thisaxis])
        self.data = ifft(self.data,n = padded_length,axis=thisaxis)
        t = self.getaxis(axes[j])
        if t is not None:
            dt = t[1]-t[0]
            self.data *= size(t) * dt # here, the algorithm divides by N, so for integration, we need to not do that
            #{{{ shiftornot specifies the shifting of the initial ft, not this result, so we always return a 0->1 time axis
            self.axis_coords[thisaxis] = linspace(0,1./dt,padded_length) + self.ft_start_time # note that I offset by ft_start_time, which I pull from when I ft'd
            #}}}
    return self
