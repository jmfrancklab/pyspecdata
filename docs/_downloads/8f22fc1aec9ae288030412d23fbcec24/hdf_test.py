"""Save data to HDF5
=================

Save nddata to standard HDF5 format.
"""
from pyspecdata import nddata, nddata_hdf5
from numpy import r_
import os

a = nddata(r_[0:9], [3, 3], ["a", "b"]).labels(
    dict(
        a=r_[0:3],
        b=r_[1:4],
    )
)
a.name("test_data")
a.set_units("a", "s").set_units("W")
a.hdf5_write("test.h5")
a_reload = nddata_hdf5("test.h5/test_data")
print("Units set for the a axis are:", a_reload.get_units("a"))
assert (a.data == a_reload.data).all()
assert (a.getaxis("a") == a_reload.getaxis("a")).all()
assert (a.getaxis("b") == a_reload.getaxis("b")).all()
assert a.get_units("a") == "s"
assert a.get_units() == "W"
print(a_reload)
os.remove("test.h5")
