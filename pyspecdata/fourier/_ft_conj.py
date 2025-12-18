from numpy import pi, log10, round

from ..general_functions import Q_


def _ft_conj(self, x):
    """Return the Fourier-conjugate unit label for *x*.

    The pyspecdata ``Q_`` helper is used to interpret the incoming unit label,
    convert it to base units, and then form the conjugate by dividing the
    ``cyc`` unit by that base quantity.  If the base-unit form introduces
    radians, the quantity is scaled by the 2π ratio between radians and cycles
    before compacting so the final result reports in Hz (or its reciprocal)
    with unit magnitude.
    """

    if not isinstance(x, str) and not hasattr(x, "to_base_units"):
        return x

    quantity = x if hasattr(x, "to_base_units") else Q_(1, x)
    base_quantity = quantity.to_base_units()
    conjugate = (Q_(1, "cyc") / base_quantity).to_base_units()

    rad_power = conjugate.units._units.get("radian", 0)
    if rad_power:
        rad_to_hz = Q_(2 * pi, "Hz") / Q_(1, "rad/s")
        if rad_power > 0:
            conjugate *= rad_to_hz**rad_power
        else:
            conjugate /= rad_to_hz ** (-rad_power)

    conjugate = conjugate.to_nice()
    to_round = (
        10 ** (round(log10(conjugate.magnitude) * 1e5) / 1e5)
        / conjugate.magnitude
    )
    conjugate *= to_round
    conjugate = conjugate.to_nice()
    assert conjugate.magnitude == 1, (
        "conjugate magnitude of"
        f" {conjugate.magnitude} {log10(conjugate.magnitude)} is not 1 when"
        f" trying to find conjugate units for {base_quantity}"
    )
    return f"{conjugate.units:~P}"
