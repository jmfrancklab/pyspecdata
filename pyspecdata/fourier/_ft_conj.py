import re
def _ft_conj(self,x):
    pairs = [('s','Hz')]
    a,b = list(zip(*tuple(pairs)))
    if x in a:
        return b[a.index(x)]
    elif x in b:
        return a[b.index(x)]
    else:
        if 'cyc' in x:
            isfourier = re.compile('cyc · \((.*)\)\$\^\{-1\}\$')
            m = isfourier.match(x)
            if m:
                retval, = m.groups()
                return retval
        else:
            return f'cyc · ({x})$^{{-1}}$'
