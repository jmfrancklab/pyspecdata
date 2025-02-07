"""Using Different FT "Viewports" with "Start Point"
=================================================
Here, I prepare a Gaussian in the frequency domain.
In one case, the Gaussian is not aliased, but
centered at zero.
In one case, the Gaussian is aliased, so that half
is at negative frequencies, and the other half at
very high positive frequencies.

I know this should give a Gaussian in the time
domain,
but I see what happens when I start from different
combinations of:

    -   (non-)/aliased signal
    -   shifting my startpoint by (non-)/integer multiples of dt
"""
import numpy as np
from numpy import r_, pi
import pyspecdata as psd
from collections import ChainMap

# https://matplotlib.org/3.2.1/api/_as_gen/matplotlib.pyplot.plot.html
default_plot_kwargs = dict(alpha=0.5, lw=2, mew=2, ms=4, marker="o", ls="none")

t = psd.nddata(
    r_[-256:257] * 10 / 256,  # to be sure that an index exactly equals zero
    "t",
).set_units(
    "t", "s"
)  # set the units to s, which are automatically converted to Hz upon FT
sigma = 1.0
with psd.figlist_var(file_name="interpolation_test_150824.pdf") as fl:
    for test_non_aliased in [True, False]:
        fl.basename = "not aliased f-domain" if test_non_aliased else "aliased f-domain"
        data = np.exp(-(t**2) / 2.0 / sigma**2)
        data.ft(
            "t", shift=test_non_aliased
        )  # this is required for the non-integral shift!
        print("is it safe?", data.get_ft_prop("t", ["freq", "not", "aliased"]))
        fl.next("ft")
        f_start = data["t"][0]
        df = data.get_ft_prop("t", "df")
        fl.plot(
            data,
            alpha=0.5,
            label=(
                f"$f_{{start}}={f_start:#0.4g}$, $f_{{start}}/\\Delta"
                f" f={f_start/df:#0.6g}$"
            ),
        )
        fl.plot(data.imag, alpha=0.5)
        psd.expand_x()
        psd.expand_y()
        print(
            "what is the initial desired startpoint?",
            data.get_prop("FT_start_time"),
        )
        print("-----------------------")
        print("starting standard")
        forplot = data.C  # keep and re-use the gaussian
        forplot.set_plot_color_next()
        print(
            "what is the initial desired startpoint?",
            forplot.get_prop("FT_start_time"),
        )
        fl.next("ift")
        forplot.ift("t")
        t_start = forplot["t"][0]
        dt = forplot.get_ft_prop("t", "dt")
        fl.plot(
            forplot,
            label=(
                f"$t_{{start}}={t_start:#0.4g}$, $t_{{start}}/\\Delta"
                f" t={t_start/dt:#0.6g}$"
            ),
            **default_plot_kwargs,
        )
        fl.plot(
            forplot.imag,
            **ChainMap({'alpha':0.03},default_plot_kwargs),
        )
        symbols = iter(["d", "x", "s", "o"])
        for this_integer in [2, -250, 1000]:
            print("-----------------------")
            print("starting integral shift for", this_integer)
            forplot = data.C  # keep and re-use the gaussian
            forplot.set_plot_color_next()
            print(
                "what is the initial desired startpoint?",
                forplot.get_ft_prop("t", "start_time"),
            )
            new_startpoint = t_start + this_integer * dt
            print("now, I try to reset the startpoint to", new_startpoint)
            print("my dt", dt, data.get_ft_prop("t", "dt"))
            forplot.ft_new_startpoint("t", "time", new_startpoint)
            print("is it safe?", data.get_ft_prop("t", ["freq", "not", "aliased"]))
            fl.next("ift")
            forplot.ift("t")
            print("And the actual t startpoint after ift? ", forplot.getaxis("t")[0])
            print(
                "the difference between the two?",
                forplot.getaxis("t")[0] - forplot.get_ft_prop("t", "start_time"),
            )
            default_plot_kwargs["marker"] = next(symbols)
            fl.plot(
                forplot,
                label=(
                    f"$t_{{start}}={t_start:#0.4g}$, $t_{{start}}/\\Delta"
                    f" t={t_start/dt:#0.6g}$"
                ),
                **default_plot_kwargs,
            )
            fl.plot(
                forplot.imag,
                **ChainMap({'alpha':0.03},default_plot_kwargs),
            )
        psd.expand_x()
        psd.expand_y()
        #
        # if test_non_aliased:
        #    symbols = iter(["d", "x", "s", "o"])
        #    for this_float in [0.5, 0.25, 10.75]:
        #        print("-----------------------")
        #        print("starting non-integral shift for", this_float)
        #        forplot = data.copy()  # keep and re-use the gaussian
        #        print(
        #            "what is the initial desired startpoint?",
        #            forplot.get_ft_prop("t", "start_time"),
        #        )
        #        print("is it safe?", data.get_ft_prop("t", ["freq", "not", "aliased"]))
        #        new_startpoint = t_start + this_float * dt
        #        print("now, I try to reset the startpoint to", new_startpoint)
        #        forplot.ft_clear_startpoints("t", t=new_startpoint, f="current")
        #        fl.next("ift -- non-integral")
        #        print("is it safe?", data.get_ft_prop("t", ["freq", "not", "aliased"]))
        #        forplot.ift("t")
        #        print(
        #            "And the actual t startpoint after ift? ", forplot.getaxis("t")[0]
        #        )
        #        print(
        #            "the difference between the two?",
        #            forplot.getaxis("t")[0] - forplot.get_ft_prop("t", "start_time"),
        #        )
        #        default_plot_kwargs["marker"] = next(symbols)
        #        default_plot_kwargs["markersize"] = 10.0
        #        fl.plot(
        #            forplot,
        #            label="$t_{start}$: shifted by %0.0fpts $\\rightarrow$ %0.2fs"
        #            % (this_float, new_startpoint),
        #            **default_plot_kwargs,
        #        )
        #        # fl.plot(forplot.runcopy(np.imag),label = 'I: integral shifted',**default_plot_kwargs)
        #    # {{{ these are manually set for a nice view of the peak of the gaussian
        #    xlim(-1, 1)
        #    ylim(0.9, 1.04)
        #    # }}}
        #
