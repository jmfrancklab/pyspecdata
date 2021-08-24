"""
Simple Cary UV-Vis loading
==========================

A simple demo of loading Cary UV-Vis data.
This example just loads a file and
plots all the spectra in a file, without embellishment

Here we have a file called Pure_T177R1a_pR_210615.BSW on our computer.
There are three requirements for where this file must be stored:

-   It **must** be stored in a folder called "proteorhodopsin" that's itself
    inside a folder called "UV_Vis" (as indicated by the ``exp_type`` argument).  Typically, this will be achieved by just
    cloning/syncing the entire "UV_Vis" directory of data shared by your lab on
    google drive, etc, etc.
-   Our pyspecdata config file (``~/.pyspecdata`` on Linux/Mac or ``~/_pyspecdata``
    on Windows) must know about this "UV_Vis" directory.
    If not, you can use the ``pyspecdata_register_dir`` command on the command line
    (see :func:`~pyspecdata.datadir.register_directory`).
-   The name of the file itself must contain the string "T177R1a_pR_210615" →
    note that you don't need to specify the whole file name, just enough for it
    to be unique.

"""
from pylab import *
from pyspecdata import *
data = find_file('T177R1a_pR_210615',
        exp_type='UV_Vis/proteorhodopsin')
print("the experiments present in this file are:",data.keys())
with figlist_var() as fl:
    fl.next("UV data")
    for j in data.keys():
        fl.plot(data[j], label=j, alpha=0.5)
