from pylab import *
from pyspecdata import *
from numpy.random import normal, seed
from matplotlib.ticker import NullLocator, MultipleLocator, FormatStrFormatter
seed(7919)
d = nddata(normal(size=10000)+1j*normal(size=10000), [100,100], ['y','x']).labels(
        {'x':r_[0:1:100j],
        'y':r_[0:0.1:100j]})
with figlist_var() as fl:
    fl.next('random data')
    fl.image(d)
    d.extend('x',extent=10)
    fl.next('extend along $x$')
    fl.image(d)
    fl.show();quit()
    d.extend('y',-0.05,fill_with=1)
    fl.next('extend along $y$')
    fl.image(d)
