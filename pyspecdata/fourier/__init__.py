# I name these all with "fourier" since they need to be renamed inside the main class
__all__ = [ "_ft_conj",
            "convolve",
            "ft",
            "ftshift",
            "ift",
            ]

from . import * # needed so that pyspecdata.fourier contains all the previous names
