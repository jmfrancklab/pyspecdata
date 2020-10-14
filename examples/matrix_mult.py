# -*- coding: utf-8 -*-
from pyspecdata import *
from numpy.random import random
import time
init_logging('debug')

# Note that in all these examples, the pyspecdata version *appears* more
# complicated.
# But, that's because these are toy examples, where we have no need for the
# dimension names or axes.
# Nonetheless, we wanted to give the simplest working example possible.

# First, we demonstrate matrix multiplication

a_nd = nddata(random(10*2048),[10,2048],['x','y'])
a = a_nd.data
# note how only the dimension that goes away is named the same!
#
# if think about the matrix a transforming from one vector space (labeled y) to
# another (labeled x) this makes sense
a2_nd = nddata(random(10*2048),[2048,10],['y','z'])
a2 = a2_nd.data

# multiply two different matrices

time1 = time.time()
b = a @ a2
time2 = time.time()
# note that here, I have to rename the column space
b_nd = a_nd.along('y') @ a2_nd
time3 = time.time()
assert b_nd.dimlabels == ['x','z'], b_nd.dimlabels
assert all(isclose(b,b_nd.data))
print("total time",(time2-time1),"vs raw",((time3-time2)/(time2-time1)))
assert ((time3-time2)/(time2-time1))<1

# calculate a projection matrix

time1 = time.time()
b = a @ a.T
time2 = time.time()
# note that here, I have to rename the column space
b_nd = a_nd.along('y',('x','x_new')) @ a_nd
time3 = time.time()
assert b_nd.dimlabels == ['x_new','x'], b_nd.dimlabels
assert all(isclose(b,b_nd.data))
if time2-time1>0:
    print("total time",(time2-time1),"vs raw",((time3-time2)/(time2-time1)))
    assert ((time3-time2)/(time2-time1))<1.1

# now, a standard dot product

a = random(10)
b = random(10)

a_nd = nddata(a,[10],['myaxis'])
b_nd = nddata(b,[10],['myaxis'])
assert all(isclose(a.dot(b),(a_nd @ b_nd).data))

# finally, let's show what happens when we *don't* rename x and multiply the
# same matrix
# 
# This will take the dot product of our 10 2048-long vectors, and present them
# 10-long array

a_nd = nddata(random(10*2048),[10,2048],['x','y'])
a = a_nd.data
b_nd = a_nd.along('y') @ a_nd
b = matmul(a_nd.data.reshape(10,1,2048),
        a_nd.data.reshape(10,2048,1))
assert all(isclose(b,b_nd.data))
assert len(b.data) == 10
