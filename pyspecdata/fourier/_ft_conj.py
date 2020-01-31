def _ft_conj(self,x):
    pairs = [('s','Hz'),('m',r'm^{-1}')]
    a,b = list(zip(*tuple(pairs)))
    if x in a:
        return b[a.index(x)]
    elif x in b:
        return a[b.index(x)]
    else:
        return None
