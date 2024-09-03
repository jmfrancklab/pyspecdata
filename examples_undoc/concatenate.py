import pyspecdata as psd
from numpy import r_

# {{{ When the dimname exists
coarse = psd.nddata(r_[0:4:4j], ["t"])
coarse += 0.1  # so data and axis coords are not the same
coarse.setaxis("t", r_[0:4:4j])
fine = psd.nddata(r_[2:3:5j], ["t"])
fine += 0.1
fine.setaxis("t", r_[2:3:5j])
combined = psd.concat([coarse, fine], dimname="t")  # not sure if it does --
#                                                     thiskwargs should match the
#                                                     existing kwarg! (dimname, dimlabel,
#                                                     etc -- may also not be a kwarg, but
#                                                     positional)
print("concatenate without sort:")
print(combined)
combined.sort(
    "t"
)  # sort is already applied in function but included here for clarity
print("concatenate after sort:")
print(combined)
# }}}
# {{{ when we make a new dimension - demos the existing functionality -- see above about kwarg
data1 = psd.nddata(r_[0:4:4j], "t")
data2 = psd.nddata(r_[0:4:4j], "t") * 1.3
combined = psd.concat([data1, data2], dimname="indirect")
print("concatenate along new dimension:")
print(combined)  # should be 2D with an indirect dimension
# }}}
