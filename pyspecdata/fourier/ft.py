from ..general_functions import *
from pylab import * 

def ft(self,*args,**kwargs):
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
        x.update({j:True})
    #}}}
    #kwargs: shiftornot=False,shift=None,pad = False
    shiftornot,shift,pad,automix = process_kwargs([
        ('shiftornot',False),
        ('shift',None),
        ('pad',False),
        ('automix',False)],
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
            padded_length = 2**(ceil(log2(padded_length)))
        elif pad:
            padded_length = pad
        self.data = fft(self.data,n = padded_length,axis=thisaxis)
        if bool(shiftornot[j]):
            if automix:
                raise ValueError("You can't use automix and shift at the same time --> it doesn't make sense")
            self.data = fftshift(self.data,axes=[thisaxis])
        t = self.getaxis(axes[j])
        if t is not None:
            dt = t[1]-t[0] # the dwell gives the bandwidth, whether or not it has been zero padded
            self.ft_start_time = t[0]
            self.data *= dt
            self.axis_coords[thisaxis] = linspace(0,1./dt,padded_length)
            if bool(shiftornot[j]):
                mask = self.axis_coords[thisaxis] > 0.5/dt
                #{{{ just the axis part of ftshift
                x = self.axis_coords[thisaxis]
                x[:] = fftshift(x)
                j = len(x)/2 # given the floor, this works out to be the central index
                x_subset = x[:j]
                x_subset -= x_subset[-1] + x[j+1] # zero and set to this
                #}}}
        if automix:
            sw = 1.0/dt
            carrier = abs(self).mean_all_but(axes[j]).argmax(axes[j]).data
            print "I find carrier at",carrier
            add_to_axis = (automix - carrier) / sw
            print "which is",add_to_axis,"times the sw of",sw,"off from the automix value of",automix
            x = self.getaxis(axes[j])
            x += round(add_to_axis)*sw
    return self
