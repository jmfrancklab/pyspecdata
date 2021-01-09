from pyspecdata import *
init_logging('debug')
# a demonstration of chunk and smoosh
# to do this, build some fake CPMG-like data
# (series of echoes one after the over, which are decaying overall)
t = nddata(r_[0:1:1024j],'t')
T2 = 0.5
T2star = 0.02
sig = exp(-t/T2) # create a t2 decay
sig.add_noise(0.1)
# in the following, we specify that there are 16 echoes on the "outer" dimension
# and let it figure out the number of points inside
sig.chunk('t',['echo','t'],[16,-1])
# change the time axis within each echo so that it's centered about t=0
t_max = sig.getaxis('t')[-1]
sig.setaxis('t', lambda x: x-t_max/2)
sig *= exp(-abs(sig.fromaxis('t'))/T2star)
with figlist_var() as fl:
    fl.next('echo decay')
    fl.image(sig)
    sig.smoosh(['echo','t'],'t')
    # the following shows a structured
    # array -- in this case (where
    # dimensions are nested times) not
    # sure if this is the behavior we
    # want or not, but at least it works
    # generally
    print(repr(sig.getaxis('t')))
    fl.next('smooshed time domain')
    fl.plot(sig.C.setaxis('t','#').set_units('t','datapoint'))
    fl.next('reconstitute after smoosh')
    sig.chunk_auto('t')
    print('t axis',sig.getaxis('t'))
    print('echo',sig.getaxis('echo'))
    fl.image(sig)
