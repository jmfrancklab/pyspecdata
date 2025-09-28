from ..general_functions import Q_


def _ft_conj(self, x):
    """Return the Fourier-conjugate unit label for *x*.

    The pyspecdata ``Q_`` helper is used to interpret the incoming unit label,
    convert it to base units, and then form the conjugate by dividing the
    ``cyc`` unit by that base quantity.  The resulting quantity is compacted and
    re-normalized to ensure the magnitude remains one before formatting the
    LaTeX-ready unit string.
    """

    if not isinstance(x, str) and not hasattr(x, "to_base_units"):
        return x

    quantity = x if hasattr(x, "to_base_units") else Q_(1, x)
    base_quantity = quantity.to_base_units()
    conjugate = (Q_(1, "cyc") / base_quantity).to_base_units().to_compact()
    assert conjugate.magnitude == 1
    return f"{conjugate.units:~L}"
