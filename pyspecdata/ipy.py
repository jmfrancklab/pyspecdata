r'''Provides the jupyter extension:

    %load_ext pyspecdata.ipy

That allows for fancy representation nddata instances -- *i.e.* you can type the name of
an instance and hit shift-Enter, and a plot will appear rather than some text
representation.

Also overrides plain text representation
of numpy arrays with latex representation that we build ourselves
or pull from sympy.

Also known as "generalized jupyter awesomeness" in only ~150 lines of code!

See [O'Reilly Book](https://www.safaribooksonline.com/blog/2014/02/11/altering-display-existing-classes-ipython/)
for minimal guidance if you're interested.'''

# I should implement this as an extension module
# https://mindtrove.info/4-ways-to-extend-jupyter-notebook/
# {{{ this allows me to skip if I called from sphinx
from .general_functions import inside_sphinx
if inside_sphinx():
    pass # called from sphinx
else:
    get_ipython().magic('pylab inline')
# }}}
import numpy
import sympy
from sympy.interactive import printing
printing.init_printing(use_latex=True,
                       wrap_line=False,
                       num_columns=10000)
if not inside_sphinx():
    from .core import image as pyspec_image
    from .core import plot as pyspec_plot
    from .core import nddata as pyspec_nddata
    from .core import gca

import re
from IPython.display import Math
def load_ipython_extension(ip):
    list(ip.display_formatter.formatters.keys())

    tex_formatters = ip.display_formatter.formatters['text/latex']
    plain_formatters = ip.display_formatter.formatters['text/plain']
    sympy_formatter = tex_formatters.for_type(sympy.Matrix)
    exp_re = re.compile(r'(.*)e([+\-])0*([0-9]+)')
    def reformat_exp(arg):
        m = exp_re.match(arg)
        if 'i' not in arg and float(arg) == 0:
            return ''
        if m:
            retstr,pm,fin_numb = m.groups()
            retstr += r'\mathrm{\scriptsize E}\!'
            retstr += pm
            retstr += r'\!\!'
            #retstr += pm if pm == '-' else ''
            retstr += fin_numb
            return retstr
        else:
            return arg
    def complex_str(arg, fancy_format=False,format_code = '%.4g'):
        "render a complex string -- leaving out imaginary if it's real"
        retval = [format_code%arg.real]
        if arg.imag != 0.0:
            retval.append((format_code+"i")%arg.imag)
        retval = [reformat_exp(j) for j in retval]
        if len(retval)>1 and retval[1][0] not in '+-':
            retval[1] = '+'+retval[1]
        return ''.join(retval)
    def render_matrix(arg):
        "return latex string representing 2D matrix"
        import IPython.display as d
        try:
            math_str = r'\begin{bmatrix}'
            math_str += '\n'
            if hasattr(arg.dtype,'fields') and arg.dtype.fields is not None:
                math_str += '\\\\\n'.join([' & '.join([', '.join([r'\text{'+f[0]+r'}\!=\!\text{"'+elem[f[0]]+'"}'
                                                                  if isinstance(elem[f[0]],str)
                                                                  else r'\text{%s}\!=\!%g'%(f[0],elem[f[0]])
                                                                  for f in arg.dtype.descr])# f[0] is the name (vs. size)
                                                       for elem in arg[k,:]]) for k in range(arg.shape[0])])
            else:
                math_str += '\\\\\n'.join([' & '.join([complex_str(j) for j in arg[k,:]]) for k in range(arg.shape[0])])
            math_str += '\n'
            math_str += r'\end{bmatrix}'
        except:
            math_str = "\\text{(unavailable pretty print)}"
            print(repr(arg))
        return math_str
    def _print_plain_override_for_ndarray(arg,p,cycle):
        """caller for pretty, for use in IPython 0.11"""
        import IPython.display as d
        if isinstance(arg, numpy.ndarray):
            if hasattr(arg.dtype,'fields') and arg.dtype.fields is not None:
                mkd = 'structured array with fields: '+', '.join([j[0] for j in arg.dtype.descr])
            else:
                mkd = "array"
            if any(numpy.array(arg.shape) > 20):
                print("Matrix is too large ("+'x'.join([str(j) for j in arg.shape])+") -- using text numpy print")
                print(arg)
            elif len(arg.shape) == 2:
                math_str = render_matrix(arg)
                d.display(d.Markdown("**numpy 2D "+mkd+"** represented as a matrix:"))
                d.display(d.Math(math_str))
            elif len(arg.shape) == 1:
                d.display(d.Markdown("**numpy 1D "+mkd+"** represented as a row vector:"))
                d.display(d.Math(render_matrix(arg.reshape(1,-1))))
            elif len(arg.shape) == 3:
                d.display(d.Markdown("***numpy 3D "+mkd+",*** represented as a series of matrices:"))
                math_str = r'\begin{bmatrix}'+'\n'
                for j in range(arg.shape[0]):
                    math_str += '\\text{Select slice with outermost dimension set to %d}'%j
                    math_str += r'\\' + '\n'
                    math_str += render_matrix(arg[j,:,:])
                    math_str += r'\\' + '\n'
                math_str += r'\end{bmatrix}'+'\n'
                d.display(d.Math(math_str))
            elif len(arg.shape) > 3:
                d.display(d.Markdown("***numpy ND array*** $N>3$ dimensions ($"+r'\times'.join(map(str,arg.shape))+"$), so I'm just giving the text representation:"))
                d.display(str(arg))
    from IPython.display import display
    def _print_plain_override_for_nddata(arg,p,cycle):
        """caller for pretty, for use in IPython 0.11"""
        import IPython.display as d
        arg_copy = arg.copy()
        arg_copy.human_units()
        if len(arg_copy.dimlabels) == 1:
            if arg_copy.data.dtype == numpy.complex128:
                pyspec_plot(arg_copy.real,'g',alpha=0.5)
                pyspec_plot(arg_copy.imag,'y',alpha=0.5)
            else:
                pyspec_plot(arg_copy)
        elif len(arg_copy.dimlabels) == 2 and any([j == 1 for j in arg_copy.data.shape]):
            arg_copy.reorder(arg_copy.dimlabels[numpy.array(arg_copy.data.shape).argmin()],
                    first=False)
            if arg_copy.data.dtype == numpy.complex128:
                pyspec_plot(arg_copy.real,'g',alpha=0.5)
                pyspec_plot(arg_copy.imag,'y',alpha=0.5)
            else:
                pyspec_plot(arg_copy)
        elif len(arg_copy.dimlabels) == 2 and any([j < 10 for j in arg_copy.data.shape]):
            arg_copy.reorder(arg_copy.dimlabels[numpy.array(arg_copy.data.shape).argmin()],
                    first=False)
            pyspec_plot(arg_copy)
        else:
            pyspec_image(arg_copy)
            if arg_copy.name() is not None:
                gca().set_title(arg_copy.name())
    plain_formatters.for_type(numpy.ndarray,_print_plain_override_for_ndarray)
    plain_formatters.for_type(pyspec_nddata,_print_plain_override_for_nddata)
    ip.ex("fancy_legend = lambda: legend(**dict(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.))")
    ip.ex("from pyspecdata import *")
def unload_ipython_extension(ip):
    print("I will not not go gentle into that good night!!!")
