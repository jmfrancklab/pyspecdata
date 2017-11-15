def ftshift(self,axis,value):
    """FT-based shift.  Currently only works in time domain.
    
    This was previously made obsolete, but is now a demo of how to use the ft properties.
    It is not the most efficient way to do this.
    """
    self.ft(axis, shift=True)
    t_start = self.get_ft_prop(axis,['start_time'])
    #print "original start time",t_start
    t_start -= value
    #print "new start time",t_start
    self.ft_clear_startpoints(axis,
            t=t_start,
            f='current')
    #self.set_ft_prop(axis,['freq_not_aliased'])
    self.ift(axis)
    self.setaxis('t',lambda x: x+value)
    return self
