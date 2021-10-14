""" HDF5 error message
==================

Trying to load an HDF5 file without an expno generates an error with a list of possible node names
"""
from pyspecdata import *
try:
    d = find_file("210603_EtOH_cap_probe_COSY", exp_type="ODNP_NMR_comp/COSY")
except ValueError as e:
    # NOTE: we use a try-except here so that the example doesn't show up as
    # "broken" -- typically, you would just use the command inside the "try"
    # command as-is
    print("this will print the error message:\n\n")
    print(e)
