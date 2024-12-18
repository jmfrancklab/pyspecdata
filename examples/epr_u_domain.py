"""
EPR u-domain
============
Show the :math:`u`-domain in EPR.

"""
import pyspecdata as psd
import re

Bname = "$B_0$"
d = psd.find_file(re.escape("220307_S175_KCl.DSC"), exp_type="francklab_esr/Farhana")
d.chunk_auto("harmonic")
d = d["harmonic", 0]["phase", 0]
d[Bname] *= 1e-4
d.set_units(Bname, "T")
d.set_ft_initial(Bname, "f")
with psd.figlist_var() as fl:
    fl.next("initial spectrum", figsize=(3 * 1.618, 3))
    fl.plot(d)
    d.ift(Bname)
    fl.next("u-domain", figsize=(3 * 1.618, 3))
    fl.plot(d)
