from pyspecdata import *
from pyspecdata.load_files.load_cary import load_cary
data = load_cary(search_filename('200703_Ellman_before_SL.DSW',
    exp_type='UV_Vis/Ellmans_Assay',
    unique=True))
with figlist_var() as fl:
    fl.next("UV data")
    for thisspectrum in data:
        fl.plot(thisspectrum,
                alpha=0.5,
                label=thisspectrum.name())
    ylabel(thisspectrum.get_units())
    ylim((-0.05,1))
