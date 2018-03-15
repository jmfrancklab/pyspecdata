r'''This package contains the functions involved in saving files'''

# import all the functions listed in "functions"
from .functions import gensearch, h5searchstring, h5loaddict, h5child, h5remrows, h5addrow, h5table, h5nodebypath, h5attachattributes, h5inlist

# and make them all accessible
__all__ = [ 
        "gensearch",
        "h5searchstring",
        "h5loaddict",
        "h5child",
        "h5remrows",
        "h5addrow",
        "h5table",
        "h5nodebypath",
        "h5attachattributes",
        "h5inlist",
        ]
