from ..general_functions import *
from pylab import * 

def register_axis(self,arg):
    r'''Interpolate the data so that the given axes are in register with a set of specified values. Does not change the spacing of the axis labels.
    
    It finds the axis label position that is closest to the values given in `arg`, then interpolates (Fourier/sinc method) the data onto a new, slightly shifted, axis that passes exactly through the value given.
    To do this, it uses 
    :func:`.ft_clear_startpoints() <pyspecdata.core.nddata.ft_clear_startpoints>`
    and uses
    :func:`.set_ft_prop() <pyspecdata.core.nddata.set_ft_prop>`
    to override the "not aliased" flag.

    Parameters
    ----------
    arg : dict (key,value = str,double)
        A list of the dimensions that you want to place in register, and the values you want them registered to.
    '''
    for k,v in arg.items():
        x = self.getaxis(k)
        idx = argmin(abs(x - v))
        offset = (v # where I want to be
                - x[idx]) # where I actually am
        offset += x[0] # since the following takes the startpoint
        if self.get_ft_prop(k):
            self.ift(k).ft_clear_startpoints(k,t='current',f=offset)
            self.set_ft_prop(k,'freq_not_aliased').ft(k)
        elif self.get_ft_prop(k) is False:
            self.ft(k).ft_clear_startpoints(k,t=offset,f='current')
            self.set_ft_prop(k,'time_not_aliased').ift(k)
        else: raise ValueError("???")
        return self
