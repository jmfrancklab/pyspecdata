"""Qualitative Shape-Based FT Tests
================================
Here, we perform various FTs whose results can be
predicted from quatlitative rules.
We show that they still behave as expected when we
do mean things like frequency filtering them, etc."""

import numpy as np
from numpy import r_, pi
import pyspecdata as psd

with psd.figlist_var(file_name="ft_demo_150820.pdf") as fl:
    t = r_[-3:10:1024j]
    t -= t[np.argmin(abs(t))]  # to be sure that an index exactly equals zero
    data = psd.nddata(np.empty_like(t), [-1], ["t"]).setaxis("t", t)
    data.set_units(
        "t", "s"
    )  # set the units to s, which are automatically converted to Hz upon FT
    sinc_width = 0.1
    data = data.fromaxis(
        "t",
        lambda x: np.where(
            x == 0, 1, np.sin(x / sinc_width) / (x / sinc_width)
        ),
    )  # `np.where` catches nan exception
    # {{{ add a simple Lorentzian peak
    R = 0.1
    f_center = 10.0
    A = 0.25
    data += A * data.fromaxis(
        "t",
        lambda x: np.where(
            x >= 0, np.exp(-2 * pi * x * R + 1j * 2 * pi * x * 10.0), 0
        ),
    )  # `np.where` acts as heaviside
    # }}}

    default_plot_kwargs = dict(alpha=0.5, linewidth=2)

    fl.next("time domain")
    fl.plot(data, **default_plot_kwargs)
    fl.plot(data.imag, **default_plot_kwargs)

    fl.next("A) frequency domain")
    data.ft("t", shift=True)
    fl.plot(data, **default_plot_kwargs)
    fl.plot(data.runcopy(np.imag), **default_plot_kwargs)

    fl.next("B) ift (A) without any funny business")
    data_save = data.copy()
    data.ift("t")
    fl.plot(data, **default_plot_kwargs)
    fl.plot(data.runcopy(np.imag), **default_plot_kwargs)

    fl.next("C) frequency filter (A)\nand generate (downsampled) ift")
    data = data_save["t":(-10, 25)]
    data_save = data.copy()
    data.ift("t")
    fl.plot(data, **default_plot_kwargs)
    fl.plot(data.runcopy(np.imag), **default_plot_kwargs)

    fl.next("D) zero-fill (1024) and ft (C)")
    data.ft("t", pad=1024)
    data_zf_ft = data.copy()
    fl.plot(data, **default_plot_kwargs)
    fl.plot(data.runcopy(np.imag), **default_plot_kwargs)

    fl.next("E) pad and ift (D), then ft back")
    data = data_zf_ft
    data.ift("t", pad=2048)
    data.ft("t")
    fl.plot(data, **default_plot_kwargs)
    fl.plot(data.runcopy(np.imag), **default_plot_kwargs)

    fl.next("F) ift (C before ift),\nzero-filling up to 1024")
    data = data_save
    data.ift("t", pad=1024)
    fl.plot(data, **default_plot_kwargs)
    fl.plot(data.imag, **default_plot_kwargs)
