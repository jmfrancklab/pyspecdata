"These are general functions that need to be accessible to everything inside pyspecdata.core.  I can't just put these inside pyspecdata.core, because that would lead to cyclic imports, and e.g. submodules of pyspecdata can't find them."

from pylab import *

def process_kwargs(listoftuples,kwargs,pass_through = False):
    '''In order, return the value of keyword arguments `kwargs` named with key, value pairs in listoftuples
    Note that having kwargs as an explicit argument avoids errors where the user forgets to pass the **kwargs.'''
    kwargnames,kwargdefaultvals = zip(*listoftuples)
    output = []
    for j,val in enumerate(kwargnames):
        output.append(kwargdefaultvals[j])
        if val in kwargs.keys():
            output[-1] = kwargs.pop(val)
    if not pass_through and len(kwargs) > 0:
        raise ValueError("I didn't understand the kwargs:",repr(kwargs))
    return tuple(output)
def autostringconvert(arg):
    if type(arg) in [unicode,str_]:
        return str(arg)
    else:
        return arg
