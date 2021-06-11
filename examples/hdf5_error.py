""" HDF5 error message
==================

Trying to load an HDF5 file without an expno generates an error with a list of possible node names
"""
from pyspecdata import *
d = find_file("210603_EtOH_cap_probe_COSY", exp_type="ODNP_NMR_comp/COSY")
