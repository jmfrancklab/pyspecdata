def register_axis(self,arg):
    for k,v in arg.iteritems():
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
