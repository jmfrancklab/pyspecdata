from .fl_dummy import fl_dummy
from sympy import Symbol as sympy_symbol
from sympy import lambdify
from pyspecdata import nddata, ndshape
import numpy as np
from numpy.random import normal

def fake_data(
    expression,
    axis_coords,
    signal_pathway,
    direct="t2",
    SD_sigma=[0.05, 0.003],
    SD_amp=[1, 5],
    scale=100,
    fake_data_noise_std=1.0,
    fl=fl_dummy
):
    """Generate fake data subject to noise and frequency variation.

    This includes a variation of the resonance frequency.  The user can adjust the scale and the timescale of the frequency variation, which is modeled by way of spectral density function that describes the random fluctuation of the resonance frequency.  (To avoid confusion, note that this spectral density function does NOT control the noise voltage, which is given by a standard normal distribution of constant variation.)

    Parameters
    ==========
    expression: sympy expression
        Gives the functional form of the data.
    axis_coords: OrderedDict
        Gives nddata objects providing all the axis coordinates.
        **Very importantly**, these must be listed in the loop nesting order
        (outside in) in which they occur in the pulse program,
        or the frequency drift will not be modeled correctly.

        To enable *simulating echo-like data*, you can specify a direct axis
        that starts at a negative number.
        If you do this, the beginning of the axis will be re-set to 0 before returning.
    SD_sigma: list of floats
        Gives the Gaussian σ for the spectral density of the
        time-dependent variation of the resonance frequency.
        Typically, there are more than one gaussian terms used to
        create the spectral density.

        A small σ gives rise to slow fluctuations while a large σ
        gives rise to fast fluctuations.
        The units of σ are in cycles per scan.
    SD_amp: list of floats
        Amplitudes associated with SD_sigma terms -- must be the
        same length.
    signal_pathway: dict
        Gives the signal pathway, with keys being phase cycling dimensions, and
        values being the corresponding Δp.
    scale: float (default 100)
        amplitude of frequency variation
    """
    mysymbols = list(expression.atoms(sympy_symbol))
    missing = set(str(j) for j in mysymbols) - set(axis_coords.keys())
    assert len(missing) == 0, (
        "all non-phase cycling symbols in your expression must have matching axis coordinates in the dictionary -- you are missing %s!"
        % str(missing)
    )
    thefunction = lambdify(mysymbols, expression, "numpy")
    clean_data = thefunction(*tuple(axis_coords[str(j)] for j in mysymbols))
    for j in signal_pathway.keys():
        clean_data *= np.exp(signal_pathway[j] * 1j * 2 * np.pi * axis_coords[j])
    ## {{{ model frequency drift
    indirect_size = np.prod(
        [ndshape(clean_data)[j] for j in axis_coords.keys() if j != direct]
    )
    frq_noise = normal(scale=scale, size=indirect_size)
    frq_noise = frq_noise + 1j * normal(scale=scale, size=frq_noise.size)
    # {{{ control the spectral density of the shifts to be gaussian
    frq_noise = nddata(frq_noise, [-1], ["temp"])
    N = ndshape(frq_noise)["temp"]
    frq_noise.setaxis("temp", -0.5 + np.r_[0:N] / N).set_units("temp", "cycperscan")
    SD_gen = zip(SD_sigma, SD_amp)
    sigma, A = next(SD_gen)
    frq_noise_dens = A * np.exp(-frq_noise.fromaxis("temp") ** 2 / 2 / sigma ** 2)
    for sigma, A in SD_gen:
        frq_noise_dens += A * np.exp(-frq_noise.fromaxis("temp") ** 2 / 2 / sigma ** 2)
    frq_noise *= frq_noise_dens
    fl.push_marker()
    fl.next("frq-noise density")
    fl.plot(frq_noise)
    frq_noise.ift("temp")
    frq_noise /= np.sqrt(ndshape(frq_noise)["temp"]) * frq_noise.get_ft_prop(
        "temp", "df"
    )  # normalization
    fl.next("frq-noise time domain")
    fl.plot(frq_noise)
    # }}}
    frq_noise = nddata(
        frq_noise.data.real,
        [ndshape(clean_data)[j] for j in axis_coords.keys() if j != direct],
        [j for j in axis_coords.keys() if j != direct],
    )
    ## }}}
    data = clean_data.C
    data.add_noise(fake_data_noise_std)
    data *= np.exp(
        1j * 2 * np.pi * frq_noise * (data.fromaxis(direct))
    )  # the frequency shift
    # at this point, the fake data has been generated
    for j in signal_pathway.keys():
        data.ft(j, unitary=True)
    fl.pop_marker()
    data.setaxis(direct, lambda x: x - data.getaxis(direct)[0])
    data.ft(direct, shift=True)
    data.ift(direct)
    data.register_axis({direct: 0})
    data.set_units(direct, "s")
    return data
