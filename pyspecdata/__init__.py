r'''
The :class:`figlist` is the base class for "Figure lists."
Figure lists allows you to organize plots and text and to refer to plots
by name, rather than number.
They are designed so that same code can be used seamlessly from within
ipython, jupyter, a python script, or a python environment within latex
(JMF can also distribute latex code for this -- nice python based
installer is planned).
The user does not initialize the figlist class directly,
but rather initializes ``figlist_var``.
At the end of this file,
there is a snippet of code that sets
``figlist_var`` to choice that's appropriate for the working environment
(*i.e.*, python, latex environment, *etc.)
'''

from .general_functions import *
from .core import *
from .load_files import *
from .figlist import *
from .nnls import *
from .lmfitdata import lmfitdata
from .DCCT_function import DCCT
from .generate_fake_data import fake_data
from .plot_funcs.image import image
#import numpy

# {{{ determine the figure style, and load the appropriate modules
_figure_mode_setting = pyspec_config.get_setting('figures', section='mode', environ='pyspecdata_figures')
if _figure_mode_setting is None:
    print("Warning!  Figure mode is not set, so I'm going to set it to standard by default!!!")
    _figure_mode_setting = 'standard'
    pyspec_config.set_setting('mode','figures','standard')
if _figure_mode_setting == 'latex':
    environ['ETS_TOOLKIT'] = 'qt4'
    import matplotlib; matplotlib.use('Agg')
# }}} -- continued below
# {{{ determine the figure style, and load the appropriate modules
if _figure_mode_setting == 'latex':
    from .fornotebook import *
    figlist_var = figlistl
elif _figure_mode_setting == 'standard':
    def obsn(*x): #because this is used in fornotebook, and I want it defined
        print(''.join(x),'\n')
    def obs(*x): #because this is used in fornotebook, and I want it defined
        print(''.join(map(repr,x)))
    def lrecordarray(*x,**kwargs):
        return repr(x) # if I'm not using tex, it's easier to not use the formatting
    def lsafe(*string,**kwargs):
        "replacement for normal lsafe -- no escaping"
        if len(string) > 1:
            lsafewkargs = lambda x: lsafe(x,**kwargs)
            return ' '.join(list(map(lsafewkargs,string)))
        else:
            string = string[0]
        #{{{ kwargs
        spaces = False
        if 'spaces' in list(kwargs.keys()):
            spaces = kwargs.pop('spaces')
        if 'wrap' in list(kwargs.keys()):
            wrap = kwargs.pop('wrap')
        else:
            wrap = None
        #}}}
        if not isinstance(string, str):
            string = repr(string)
        if wrap is True:
            wrap = 60
        if wrap is not None:
            string = '\n'.join(textwrap.wrap(string,wrap))
        return string
    figlist_var = figlist
else:
    raise ValueError("I don't understand the figures mode "+_figure_mode_setting)
# }}}

# so essentially, __all__ is the namespace that is passed with an import *
#__all__ = ['prop',
#        'nddata',
#        'figlist_var',
#        'plot',
#        'OLDplot',
#        'nddata_hdf5']
#__all__.extend(numpy.__all__)
__all__ = [x for x in dir() if x[0] != '_']
