from pyspecdata import nddata, r_

a = nddata(r_[0:1e-6:500j], "t").set_units("t", "s")
a.human_units()
print("this should be 1e-9", a.div_units("t", "s"))
try:
    print(a.div_units("t", "T"))
except Exception:
    print("yay! I failed when I should have")
