"""
Saving an HDF5 file
===================
Example of how to save 
data to an HDF5 file.
"""
from pyspecdata import *
init_logging('debug')
a = nddata(r_[0:5:10j], 'x')
a.name('test_data')
try:
    a.hdf5_write('example.h5',getDATADIR(exp_type='francklab_esr/Sam'))
except:
    print("file already exists, not creating again -- delete the file or node if wanted")
# read the file by the "raw method"
b = nddata_hdf5('example.h5/test_data',
        getDATADIR(exp_type='francklab_esr/Sam'))
print("found data:",b)
# or use the find file method
c = find_file('example.h5', exp_type='francklab_esr/Sam',
        expno='test_data')
print("found data:",c)
