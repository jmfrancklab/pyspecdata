""" written based on igor code
(under GPL license)
taken from `here <https://github.com/Yohko/importtool/blob/master/Igor%20Procedures/import%20tool/import_VARIAN_CARY_BIN.ipf>`_.

"""
from numpy import *
import os

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
            data = load_bindata(fp)
            footer = load_footer(fp)
