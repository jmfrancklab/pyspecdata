from pyspecdata import nddata, nddata_hdf5
from np import r_
import os

a = nddata(r_[0:9], [3, 3], ["a", "b"]).labels(
    dict(
        a=r_[0:3],
        b=r_[1:4],
    )
)
a.name("test_data")
a.set_units("a", "s")
a.hdf5_write("test.h5")
a_reload = nddata_hdf5("test.h5/test_data")
print("Units set for the a axis are:", a_reload.get_units("a"))
assert all(a.data == a_reload.data)
assert all(a.getaxis("a") == a_reload.getaxis("a"))
assert all(a.getaxis("b") == a_reload.getaxis("b"))
print(a_reload)
os.remove("test.h5")
