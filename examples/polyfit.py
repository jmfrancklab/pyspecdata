from pylab import *
from pyspecdata import *
# {{{ generate fake data
x = nddata(r_[0:10:7j], 'x')
y = (x-2.0)**2
y.add_noise(5)
# }}}
plot(y,'o')
c = y.polyfit('x', order=2)
assert len(c)==3
# generate a polynomial that's more finely spaced
x = nddata(r_[0:10:100j], 'x')
plot(x.apply_poly(c,'x'))
show()
