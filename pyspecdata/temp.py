# TODO: Delete before merging.
from pyspecdata import *
import numpy as np

a = a = ndshape([10], ["t"])
print(a)
print("Allocated", a.alloc())
print(
    "Adding labels",
    a.alloc(dtype=np.dtype([("field", "f8"), ("power", "f8")])),
)
a = a.alloc(dtype=np.dtype([("field", "f8"), ("power", "f8")]))
a.setaxis("t", 0.1 * np.r_[0:10])
print("Axis is set", a)
a.data["power"] = np.r_[11:21]
a.data["field"] = np.r_[1:11]
print(a)
a = a.mean("t")
print(a)
