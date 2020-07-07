""" written based on igor code
(under GPL license)
taken from `here <https://github.com/Yohko/importtool/blob/master/Igor%20Procedures/import%20tool/import_VARIAN_CARY_BIN.ipf>`_.

"""
from numpy import *
import os

def load_data(fp,header):
    x_mode_rep = "nm;Å;cm-1;°".split(';')
    y_mode_rep = "Abs;%T;Absorptivity;%R;4?;5?;Log(Abs);Absolute %R;F(R);Log(F(R));Log(1/R)".split(';')
    factoreV_list = [1239.8424,1239.8424,12398.424,8065.54429,0]
    x_unit = x_mode_rep[header['x_mode']]
    factoreV = factoreV_list[header['x_mode']]
    y_unit = y_mode_rep[header['y_mode']]
    tmps = str(header['spec_counter'])+"_"+header['spectrum_name']+"_"+x_unit
def load_header(fp):
    retval = fromfile(fp, dtype=[
        ("softversion","<u4"),
        ("tmp1","<u4"),
        ("end_x","<i4"),
        ("start_x","<i4"),
        ("V_min","<i4"),
        ("V_max","<i4"),
        ("points","<u4"),
        ("tmp2","<i4"),
        ("tmp3","<u2"),
        ("tmp4","<u2"),
        ("tmp5","<i4"),
        ("spectrum_number","<i4"),
        ("tmp6","<f4"),
        ("tmp7","<i4"),
        ("tmp8","<f4"),
        ("tmp9","<i4"),
        ("x_mode","<f4"),
        ("y_mode","<f4")], count=1)
    return dict(zip(retval.dtype.fields,retval.item()))

def load_cary(filename):
    filesize = os.stat(filename).st_size
    print("filesize",filesize)
    with open(filename,'rb') as fp:
        # line 92
        thislen = fromfile(fp,dtype='<u1', count=1).item()
        print("thislen is",thislen)
        magic = fromfile(fp,dtype=f'a{thislen}', count=1).item()
        print("magic is",magic)
        _ = fromfile(fp,dtype=f'a{61-thislen}', count=1)
        print("position",fp.tell())
        # line 100
        thislen = fromfile(fp,dtype='<u4', count=1).item()
        Tstore_type = fromfile(fp,dtype=f'a{thislen}', count=1).item().decode('ascii')
        print("Tstore_type",Tstore_type)
        # line 102
        blocklen = fromfile(fp,dtype='<u4', count=1).item()
        print("blocklen",blocklen)
        if Tstore_type == "TContinuumStore":
            # line 106
            header = load_header(fp)
            print("header:",header)
            data = load_bindata(fp)
            footer = load_footer(fp)
