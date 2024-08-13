""" Test units are transferred from HDF5 file
=============================================
Test that the units saved in the HDF5 file are transferred upon opening
"""

import pyspecdata as psd

filename = "240805_amp1_nutation"
nodename = "nutation_4"
data = psd.find_file(
    filename, exp_type="ODNP_NMR_comp/nutation", expno=nodename
)
print("units of direct axis saved in HDF5 file are", data.get_units("t2"))
