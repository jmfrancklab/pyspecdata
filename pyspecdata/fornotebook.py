r'''This provides :class:`figlistl`, the Latex figure list.
Any other functions here are helper functions for the class.
:class:`figlist` is generally **not chosen manually**,
but ``figlist_var`` will be assigned to :class:`figlistl` when
python code is embedded in a python environment inside latex.
'''
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from pylab import *
import warnings
# sympy doesn't like to be imported from fornotebook as part of a *
warnings.filterwarnings("ignore")
from .core import *
from .figlist import figlist
from .general_functions import fname_makenice
from scipy.io import savemat,loadmat
from os.path import exists as path_exists
from os import name as os_name
from scipy.optimize import leastsq
from scipy.signal import fftconvolve
from datetime import datetime
from time import mktime
from PIL import Image
import numpy as np
import re
import sys

golden_ratio = (1.0 + np.sqrt(5))/2.0

def dprint(*stuff):
    print('\n\nDEBUG',' '.join(map(repr(stuff))),'\n\n')
def thisjobname():
    if path_exists('pythonjobname.txt'):
        with open('pythonjobname.txt') as fp:
            retval = fp.read().rstrip()
    else:
        retval = 'JobnameUndefined'
    return retval
def show_matrix(B,plot_name,first_figure = None):
    fl = figlistini(first_figure)
    #nextfigure(fl,{'width':1}) # only do this once I am able to pop the command later
    nextfigure(fl,'matrix'+plot_name)
    B[~isfinite(B)] = max(B.flatten())
    image(B)
    colorbar()
    #nextfigure(fl,'matrix'+plot_name+'_finite')
    #image(~isfinite(B))
    return figlistret(first_figure,fl)
class figlistl (figlist):
    def __init__(self,*args,**kwargs):
        super(figlistl,self).__init__(*args,**kwargs)
        self.black = False
        self._print_at_end = False
        return
    def par_break(self):
        self.text("\\par")
        return
    def show(self,string,line_spacing=True,**kwargs):
        '''latexify the series of figures, where "string" gives the base file name
        
        Parameters
        ----------
        line_spacing : bool
            if false, suppress empty lines between output
        '''
        self.basename = None # must be turned off, so it can cycle through lists, etc, on its own
        #{{{ process kwargs
        verbose = False
        if 'verbose' in list(kwargs.keys()):
            verbose = kwargs.pop('verbose')
        #}}}
        if line_spacing: print('\n\n')
        self.show_prep()
        if not len(kwargs):
            kwargs = {}
            mlab = False
            if hasattr(self,'mlab') and self.mlab:
                kwargs.update({'mlab':self.mlab})
                mlab = True
        if hasattr(self,'lplot_kwargs'):
            kwargs.update(self.lplot_kwargs)
        for figname in self.figurelist:
            if verbose: print("showing figure"+lsafen(figname))
            if isinstance(figname, dict):
                kwargs.update(figname)
                if 'print_string' in kwargs:
                    if line_spacing: sys.stdout.buffer.write(b'\n\n')
                    outstr = kwargs.pop('print_string')
                    # {{{ to "print" with correct encoding
                    sys.stdout.buffer.write(outstr.encode('utf-8'))
                    # }}}
                    if line_spacing: sys.stdout.buffer.write(b'\n\n')
                    sys.stdout.buffer.flush()
            else:
                j = self.get_fig_number(figname)
                if mlab:
                    #self.mlab.options.offscreen = True
                    thefigure = self.figdict[figname] # an object
                    fig = self.mlab.figure(thefigure,bgcolor = (1,0,0),fgcolor = (1,1,1)) # set background to red (no antialiasing) to be set to transparent later
                    fig.scene.render_window.aa_frames = 0
                    fig.scene.anti_aliasing_frames = 0
                    #fig.scene.off_screen_rendering = True
                else:
                    figure(j)
                sep = ''
                if len(string)>0:
                    sep = '_'
                if self.env == 'test':
                    print('not running lplot')
                else:
                    if mlab:
                        lplot(figname.replace('.','_')+sep+string,fig = thefigure,**kwargs)
                    else:
                        if figname in list(self.twinx_list.keys()):
                            kwargs.update(autopad = False)# because the autopad is done manually, since it freaks out with twinx
                        lplot(figname.replace('.','_')+sep+string,**kwargs)
        while len(self.figurelist) > 0:
            self.figurelist.pop(-1)
        return
    def obs(self,*args):
        self.text(obs_repr(*args))
        return
    def obsn(self,*args):
        self.text(obs_repr(*args)+'\n\n')
        return
    def __enter__(self):
        if not hasattr(self,'file_name'):
            raise ValueError("You can't use a with statement if you don't pass the initialization a file name")
        return self
figlist = figlistl
def lplotfigures(figurelist,string,**kwargs):
    'obsolete, use the class!'
    print('\n\n')
    if not len(kwargs):
        kwargs = {}
    for j,figname in enumerate(figurelist):
        if isinstance(figname, dict):
            kwargs.update(figname)
            if 'print_string' in kwargs:
                print('\n\n'+kwargs.pop('print_string')+'\n\n')
        else:
            figure(j+1)
            try:
                lplot(figname+string,**kwargs)
            except:
                raise CustomError('figurelist = ',figurelist,'figname = ',figname,'while string = ',string)
    while len(figurelist) > 0:
        figurelist.pop(-1)
    return
def figlisterr(figurelist,*args,**kwargs):
    if 'basename' in list(kwargs.keys()):
        basename = kwargs['basename']
    else:
        basename = thisjobname()
    print(lsafen("Trying to plot the figurelist",figurelist))
    lplotfigures(figurelist,basename+'errplot.pdf')
    return args
def see_if_math(recnames):
    'return latex formatted strings --> ultimately, have this actually check to see if the individual strings are math or not, rather than just assuming that they are.  Do this simply by seeing whether or not the thing starts with a full word (more than two characters) or not.'
    def ismath(myin):
        out = myin.split('_')
        myin = out
        out = []
        for j in range(0,len(myin)):
            out.extend(myin[j].split('('))
        myin = out
        out = []
        for j in range(0,len(myin)):
            out.extend(myin[j].split(')'))
        return not any([len(x)>1 and not (x[0] == '{' or x[0] == '\\') for x in out])
    return [r'\ensuremath{'+v.replace('\\\\','\\')+'}' if ismath(v) else v.replace('_',' ') for v in recnames]
def lrecordarray_broken(recordlist,rows=30,numwide= 5):
    width = 1.0/numwide
    retval = ''
    for j in range(0,len(recordlist),rows):
        retval += r'\begin{minipage}[t]{%0.3f\linewidth}'%width
        retval += lrecordarray(recordlist[j:j+rows],resizebox = True)+'\n'
        retval += r'\end{minipage}'
    if j+rows<len(recordlist):
        retval += r'\begin{minipage}[t]{%0.3\linewidth}'%width
        retval += lrecordarray(recordlist[j+rows:],resizebox = True)+'\n'
        retval += r'\end{minipage}'
    return retval
def lrecordarray(recordlist,columnformat = True,smoosh = True,multi = True,resizebox = False,showpipe = True,return_only=False,format = '%0.3f',std_sf = 2,scientific_notation = True):
    r'''generate latex representation of a structured array
    if set to True, resizebox will automatically scale down the table so it fits on the page (but it will also scale up a small table to fit the width of the page)
    resizebox can also be a fractional number, so that it is resized to a fraction of the page'''
    error_names = [j for j in recordlist.dtype.names if j[-6:] == '_ERROR']
    data_names = [j for j in recordlist.dtype.names if j[-6:] != '_ERROR']
    data_has_error_names = [j for j in data_names if j+'_ERROR' in error_names]
    retval = ''
    if len(data_has_error_names)<len(error_names):
        retval += r'{\color{red}Warning! unassigned error names!'+'\n'
        retval += lsafe(set(error_names)-set([x+'_ERROR' for x in data_has_error_names]))+'\n'
        retval += '}\n\n'
    recordlist_errors = recordlist[error_names]
    recordlist = recordlist[data_names]
    final_retval = []
    if len(recordlist) == 0:
        retval += r'{\color{red}This array is empty!!}'+'\n'
        return retval
    # previously documented
    def this_format_function(x,error_val = None,scientific_notation = scientific_notation):
        if isinstance(x, np.str_):
            if error_val is None:
                return str(x)
            else:
                return str(x)+'\\pm '+str(error_val)
        elif std_sf is not None:
            #just return two significant figures
            number_sf = std_sf
            if error_val is None or error_val == 0.0:
                if x == 0.0:
                    return '0.0'
                powerof = np.floor(np.log10(abs(x)))
                if scientific_notation and abs(powerof)>3:
                    if abs(powerof)>3:
                        return ('%%0.%df' % int(number_sf))%float(x/(10**powerof))+r'\magn{%d}'%int(powerof)
                else:
                    highest_significant = np.floor(np.log10(x))
                    lowest_significant = highest_significant-(number_sf-1) # the -1 is for number_sf significant figures, not just one
                    x /= 10**lowest_significant
                    x = round(x)
                    x *= 10**lowest_significant
                    if lowest_significant<0:
                        thisformat = '%%0.%df'%(-1*lowest_significant)
                    else:
                        thisformat = '%0.0f'
                    return (thisformat)%(x)
            else:
                highest_significant = np.floor(np.log10(error_val))
                lowest_significant = highest_significant-(number_sf-1) # the -1 is for number_sf significant figures, not just one
                error_val /= 10**lowest_significant
                x /= 10**lowest_significant
                error_val = round(error_val)
                x = round(x)
                error_val *= 10**lowest_significant
                x *= 10**lowest_significant
                if lowest_significant<0:
                    thisformat = '%%0.%df'%int(-1*lowest_significant)
                else:
                    thisformat = '%0.0f'
                return (thisformat+'\\pm '+thisformat)%(x,error_val)
        else:
            if error_val is None:
                return (format)%(x)
            else:
                return (format+'\\pm '+format)%(x,error_val)
    def maybe_matrix(x,error = None):
        if error is not None:
            if isinstance(x, np.ndarray):
                if x.shape != error.shape:
                    raise ValueError("shape of the error does not match the shape of the value/matrix")
        if isinstance(x, np.ndarray):
            if len(x.shape) == 2:
                if error is None:
                    #retval = r'\begin{tabular}{'+'c'*x.shape[1]+'}'
                    retval = r'$\displaystyle\begin{pmatrix}'
                    matrixlist = []
                    for j in range(0,x.shape[0]):
                        matrixlist.append(' & '.join(map(this_format_function,list(x[j,:]))))
                    retval += (r'\\' +'\n').join(matrixlist)
                    retval += r'\end{pmatrix}$'
                else:
                    # values first
                    retval = r'$\displaystyle\begin{pmatrix}'
                    matrixlist = []
                    for j in range(0,x.shape[0]):
                        matrixlist.append(' & '.join(map(this_format_function,list(x[j,:]),list(error[j,:]))))
                    retval += (r'\\' +'\n').join(matrixlist)
                    retval += r'\end{pmatrix}$'
                return retval
            else:
                return "an array:"+lsafe(str(x))
        else:
            if isinstance(x,str):
                return lsafe(this_format_function(x,error))
            else:
                return '$'+this_format_function(x,error)+'$'
    if resizebox not in [True,False]:
        final_retval.append(r'\resizebox{%0.3f\linewidth}{!}{'%resizebox)
    elif resizebox:
        final_retval.append(r'\resizebox{\linewidth}{!}{')
    reclen = len(recordlist)
    recnames = recordlist.dtype.names
    alltext = ''
    if columnformat:
        colheadings = [r'c']*len(recnames)
        recordstrings = see_if_math(recnames)
        datastrings = []
        clinerow = [True]*len(recnames)
        clinelist = []
        def gencline(inputar):
            ind = 0
            listoftuples = []
            while ind < len(inputar):
                if inputar[ind]:
                    startval = ind
                    while ind < len(inputar) and inputar[ind]:
                        ind += 1
                    endval = ind - 1
                    listoftuples.append((startval,endval))
                else:
                    ind += 1
            final_retval = ''
            for thistuple in listoftuples:
                if thistuple[0] == thistuple[1]:
                    final_retval = final_retval + r'\cline{%d-%d}'%tuple([thistuple[0]+1]*2)
                elif thistuple[0] == 0 and thistuple[1] == len(recnames)-1:
                    final_retval = final_retval + r'\hline'
                else:
                    final_retval = final_retval + r'\cline{%d-%d}'%tuple(map((lambda x: x+1),thistuple))
            return final_retval
        #{{{ first, store the strings for the data
        for j in range(0,len(recordlist)):
            datastrings.append([maybe_matrix(recordlist[v][j],error = recordlist_errors[v+'_ERROR'][j])
                if v in data_has_error_names
                else maybe_matrix(recordlist[v][j])
                for v in recnames])
            clinelist.append(list(clinerow))
        #}}}
        if multi:
            #{{{ now, do a multirow replacement
            for j in range(0,len(datastrings[0])): # loop over the columns
                k = 0
                while k<len(datastrings): # loop over rows (outer)
                    # see if the next row is the same
                    numrows = 1
                    while ((numrows+k < len(datastrings))
                            and
                            (datastrings[k+numrows][j] == datastrings[k][j])):
                        datastrings[k+numrows][j] = '' # clear the data of the later row
                        clinelist[k+numrows-1][j] = False # remove this cline
                        numrows += 1
                    #clinelist[k+numrows-1][j] = True # add the cline for the last one
                    # if there are matching rows, replace the first with a multirow statement
                    if numrows > 1:
                        datastrings[k][j] = r'\multirow{%d}{*}{%s}'%(numrows,datastrings[k][j])
                    k += numrows # jump ahead past the end of the processed rows
            #}}}
        for j in range(0,len(recordlist)):
            if smoosh:
                for k,l in enumerate(map(len,datastrings)):
                    if l>40:
                        colheadings[k] = r'p{%0.2f\textwidth}'%(0.8/len(recnames))
            alltext+=' & '.join(datastrings[j])+r'\\ '+gencline(clinelist[j])+'\n'
        alltext += r'\end{tabular}'
        if showpipe:
            final_retval.append(r'\begin{tabular}{'+('|'.join(colheadings))+'}')
        else:
            final_retval.append(r'\begin{tabular}{'+(''.join(colheadings))+'}')
        final_retval.append((' & '.join(recordstrings))+r'\\ \hline\hline'+'\n')
        final_retval.append(alltext)
    else:
        final_retval.append(r'\begin{tabular}{'+(''.join(['r']+['l']*reclen))+'}')
        for v in recnames:
            final_retval.append(r'\ensuremath{'+v+'}$=$ &'+(' & '.join(map(lambda x: lsafe(str(x))+list(recordlist[v]))))+r'\\')
        final_retval.append(r'\end{tabular}')
    if resizebox:
        final_retval.append(r'}')
    if return_only:
        return '\n'.join(final_retval)
    else:
        retval += '\n'.join(final_retval)
        return retval
def lplot(fname, width=0.33, figure=False, dpi=72, grid=False,
        alsosave=None, gensvg=False, print_string=None, centered=False,
        equal_aspect=False, autopad=True, bytextwidth=None,
        showbox=True, boundaries=True, genpng=False, mlab=False,
        fig=None, verbose=False):
    '''
    used with python.sty instead of savefig
    
    by default, if width is less than 1, 
    it's interpreted as bytextwidth = True (i.e. width given as a fraction of the linewidth)
    if it's greater than, the width is interpreted in inches.
    '''
    if width <= 1.0:
        bytextwidth = True
    else:
        bytextwidth = False
    if print_string is not None:
        print(print_string)
    fname = r'auto_figures/'+fname_makenice(fname)
    if alsosave != None:
        alsosave = r'auto_figures/'+alsosave
    if gensvg == True:
        alsosave = fname.replace('.pdf','.svg')
    if genpng == True:
        alsosave = fname.replace('.pdf','.png')
    if grid:
        gridandtick(plt.gca())
    if fig is None:
        fig = plt.gcf()
    if mlab:
        temp = fig.scene.anti_aliasing_frames
        fig.scene.disable_render = True
        fig.scene.anti_aliasing_frames = 0 # mayavi antialiasing is terrible, so just acquire at a high dpi setting
        mlab.savefig(fname, magnification=int(dpi/72+0.5),
                bbox_inches='tight')
        #{{{ convert to transparent, find red pixels, and change them to transparent
        data = imread(fname)
        data_shape = list(data.shape)
        data_shape[-1] = 4
        new_data = ones(data_shape,dtype = data.dtype)
        new_data[:,:,:3] = data
        data = new_data
        new_data = data.view(dtype = [('',data.dtype)]*4)
        mask = new_data == array([(1,0,0,1)],dtype = new_data.dtype)
        new_data[mask] = array([(1,1,1,0)],dtype = new_data.dtype)
        #data = imresize(data,data.shape[:2]/4,interp = 'bilinear') # this takes huge amounts of memory
        imsave(fname,data)
        if os_name == 'posix':# seems like windows can't handle the resize
            imshape = data.shape[0:2]
            resize_factor = 4
            data = uint8((data*255).round())
            img = Image.fromarray(data)
            if verbose: print("old size",img.size)
            newimgsize = tuple(uint((np.r_[img.size[0],img.size[1]]/resize_factor).round()))
            if verbose: print("target new size",newimgsize)
            img = img.resize(newimgsize,Image.ANTIALIAS)
            if verbose: print("new size, after resize",img.size)
            del data
            img.save(fname)
            del img
        #}}}
        fig.scene.disable_render = False
        fig.scene.anti_aliasing_frames = temp
    else:
        ax = plt.gca()
        if equal_aspect:
            ax.set_aspect('equal')
        fig.autofmt_xdate()
        if autopad:
            autopad_figure(centered = centered,figname = fname)
    # replaced outer_legend with appropriate modification to the "legend" option of figlist.show_prep(), same with legend option
        if not boundaries:
            ax = plt.gca()
            for j in list(ax.spines.keys()):
                ax.spines[j].set_visible(False)
            setp(ax.get_xticklabels(),visible = False)
            setp(ax.get_yticklabels(),visible = False)
            setp(ax.get_xticklines(),visible = False)
            setp(ax.get_yticklines(),visible = False)
            this_xlabel = ax.get_xlabel()
            if len(this_xlabel) > 0:
                ax.set_xlabel(this_xlabel + r" $\rightarrow$")
            this_ylabel = ax.get_ylabel()
            if len(this_ylabel) > 0:
                ax.set_ylabel(this_ylabel + r" $\rightarrow$")
        try:
            plt.savefig(fname, dpi=dpi,
                    facecolor=(1,1,1,0),
                    bbox_inches='tight')
        except ValueError as exc_string:
            if exc_string.find('finite numbers') > -1:
                raise ValueError("It gives this error because you're trying to do a bar graph with zero width")
            else:
                raise ValueError(exc_string)
        except IOError as e:
            raise IOError("This is giving an IOError -- check that you haven't maybe changed directories" + os.getcwd().replace('\\','/'))
        if alsosave != None:
            plt.savefig(alsosave,
                    dpi=dpi,
                    facecolor=(1,1,1,0))
    if figure:
        print(r"""
        \begin{figure}[h]
        \end{figure}
        """)
    if bytextwidth:
        if showbox:
            mpwidth = r'%0.2f\linewidth'%width
            figwidth = r'\linewidth'
        else:
            figwidth = r'%0.2f\linewidth'%width
    else:
        mpwidth = r'%0.2fin'%width
        figwidth = mpwidth
    def bufwr(inpstr, end='\n'):
        endstr = end.encode('utf-8')
        sys.stdout.buffer.write(inpstr.encode('utf-8')+endstr)
    if showbox:
        bufwr(r'''\mbox{\begin{minipage}{%s}'''%mpwidth)
        if alsosave != None:
            bufwr(r'also saved \fn{%s}'%alsosave+'\n\n')
    bufwr(r'\includegraphics[width=%s]{%s}'%(figwidth,fname.replace(r'auto_figures',r'\autofiguredir')))
    if showbox:
        bufwr('\n\n'+r'\hrulefill'+'\n\n')
        bufwr(r'''{\color{red}{\tiny %s}:}\begin{tiny}\fn{%s}\end{tiny}'''%('file:',fname))
        bufwr('\n\n'+r'\hrulefill'+'\n\n')
        bufwr(r'''\end{minipage} }''', end=' ')
    plt.clf()
    sys.stdout.buffer.flush()
    return
def ordplot(x,y,labels,formatstring):
        order = argsort(x)
        plot(x[order],y[order],'o-')
        newlabels=[]
        for j in range(0,len(order)): 
            newlabels += [labels[order[j]]]
        addlabels(formatstring,x[order],y[order],newlabels)
def obs_repr(*arg):
    thisstr = ''
    thisstr += r'\o{'
    num = 0
    for j in arg:
        if type(j) in [np.double,float,np.float32]:
            if np.log10(j) > 3 or np.log10(j) < -3:# greater than a thousand or less than a thousandth
                power = np.floor(np.log10(j))
                j /= 10**power
                thisstr += r'$%0.3f\times 10^{%d}$'%(j,power)
            else:
                if np.log10(j)<-1:
                    temp = '%%0.%df'%(np.floor(-np.log10(j))+3)
                    thisstr += temp%j
                else:
                    thisstr += r'%0.3f'%j
        else:
            thisstr += j
        if num < len(arg)-1:
            thisstr += " "
            num += 1
    thisstr += r'}'
    return thisstr
def obs(*arg):
    print(obs_repr(*arg))
    return
obsndef = True
def obsn(*arg):
    obs(*arg)
    print('\n\n')
    return
def txt_to_dict(file='data.txt'):
    fp = open(file,'r')
    retval = {}
    for j in fp.readlines():
        j = j.replace('\\n','\n')
        k,v = tuple(j.split('::'))
        retval.update({k:eval(v)})
    fp.close()
    return retval
def dict_to_txt(mydict,file='data.txt'):
    set_printoptions(precision = 16)
    fp = open(file,'w')
    for k,v in mydict.items():
        fp.write('%s::%s\n'%(k,repr(v).replace('\n','\\n')))
    fp.close()
def save_data(inputdata={},file = 'data.txt'):
    # concatenate data to data already in data.mat file
    if path_exists(file):
        #data = loadmat(mat_file,struct_as_record=True)
        data = txt_to_dict(file = file)
    else:
        data = {}
    if not(inputdata=={}):
        try:
            data.update(inputdata)
        except:
            raise CustomError('trying to update',data,'with',inputdata)
        try:
            dict_to_txt(data,file = file)
            #print "DEBUG, saved",data,'to',mat_file
        except:
            raise CustomError('trying to write',data,'to',file)
    return data
def save_local(inputdata={}):
    data = save_data(inputdata,file = 'local.txt')
    return data
def clear_local(inputdata=[]):
    obs(r'{\it clear local}'+'\n\n')
    if path_exists('local.txt'):
        os.unlink('local.txt')
def save_variable(variable,content,disp=True,file = 'data.txt'):
    if path_exists(file):
        #data = loadmat('data.mat',struct_as_record=True)
        data = txt_to_dict(file = file)
    else:
        data = {}
    data.update({variable:content})
    dict_to_txt(data,file = file)
    if disp:
        obs(variable.replace('_',r'\_'),'=',dp(content,5))
    return data
#{{{ specific convenience functions
def calcdistance(freq,optdistance):
    print(r'$%0.1f\;GHz$ @ '%(freq/1e9),dp((optdistance+119e-3)*1e3,1),r'$mm$')
def calcfield(elfreq,elratio = 28.113,nmrelratio = 1.5167):
    obs('$%0.6f'%elfreq,r'\;\text{GHz}$ ${\tiny ',dp(elratio,5),r'\;\text{GHz/T}=}',dp(1e4/elratio,5),r'\;\text{G/GHz}$, $',dp(nmrelratio,6),r'\permil\;\Rightarrow$ $',dp(elfreq/elratio*1e4,3),r'\;\text{G}$, $',dp(elfreq*nmrelratio,5),r'\;\text{MHz}$')
def qcalc(freq1,freq2):
    centerfreq = (freq2+freq1)/2
    q = centerfreq/abs(freq2-freq1)
    obs('Q=',dp(q,2))
def ernstangle(pulsetime=None,pulsetip=None,Tr=None,T1=None):
    ninetytime = pulsetime/pulsetip*pi/2.
    if not(Tr==None):
        alpha = arccos(exp(-Tr/T1))
        obs(r'optimal $T_r=%0.2f\;s$ with $\alpha=%0.2f^o$'%(Tr,alpha*180/pi))
    obs(r'$\tau_{90}=%0.2f\;\mu s$ $\tau_\alpha=%0.2f\;\mu s$'%(ninetytime,ninetytime*alpha/pi*2.))
def calcfielddata(freq,substance,spec=''):
    data = save_data()
    calcfield(freq,elratio=data[substance+'_elratio'+spec],nmrelratio=data[substance+'_nmrelratio'])
    obs('(for radical',substance,')')
    save_data({'current_frequency':freq})
    save_data({'current_ppt':data[substance+'_nmrelratio']})
def cpmgseries(filename,plotlabel,tau=None,alpha=None,alphaselect=None):
    data = load_file.prospa.load_datafile(filename,dims=2)
    plot(data)
    lplot(plotlabel+'.pdf')
    #data = load_2d(filename)
    #data.ft('t2')
    #image(data)
    #lplot(plotlabel+'ft.pdf',grid=False)
    data = process_cpmg(filename)
    if (tau!=None):
        coeff,fit,rms = regularize1d(data.data,data.getaxis('echo'),tau,alpha)
        plot(alpha,rms)
        if (alphaselect!=None):
            coeff,fit,rms = regularize1d(data.data,data.getaxis('echo'),tau,[alphaselect])
            plot([alphaselect],rms,'rx')
        axis('tight')
        lplot(plotlabel+'_reg.pdf')
        print('\n\n')
    plot(data)
    if (alphaselect!=None):
        plot(data.getaxis('echo'),fit.flatten(),'r',alpha=0.5,linewidth=2)
    lplot(plotlabel+'_filt.pdf')
    if (alphaselect!=None):
        plot(tau,coeff.flatten(),'r')
    lplot(plotlabel+'_coeff.pdf')
def cpmgs(exp,number,tau=None,alpha=None,alphaselect=None,first=False):
    #{{{ carry over stored data
    if first:
        clear_local()
    if (tau!=None):
        local_data = save_local({'tau':tau})
    if (alpha!=None):
        local_data = save_local({'alpha':alpha})
    if (alphaselect!=None):
        local_data = save_local({'alphaselect':alphaselect})
    local_data = save_local()
    tau = local_data['tau']
    alpha = local_data['alpha']
    alphaselect = local_data['alphaselect']
    #}}}
    filename = getDATADIR()+'franck_hanlabMagritek/'+exp+'/%d/'%number
    print(r'\fn{%s %d}'%(exp,number)+'\n\n')
    cpmgseries(filename,exp+thisjobname(),tau,alpha,alphaselect)
#{{{ esr_saturation
def esr_saturation(file,powerseries,smoothing=0.2,threshold=0.8,figname = None,hn23adjustment = 1.,show_avg = False):
    if figname == None:
        figname = thisjobname()
    print('\n\n\\fn{',file,'}\n\n')
    data = load_indiv_file(file,dimname='power')
    #plot(data,'.-')
    x = data.getaxis('$B_0$').flatten()
    k = exp(-(x-x.mean())**2/2./smoothing**2)
    nslices = ndshape(data)['power']
    allpeaks_top = []
    allpeaks_bottom = []
    allpeaks_top_x = []
    allpeaks_bottom_x = []
    imageformat = False
    for j in range(0,nslices):
        #{{{ make a smoothed version of this slice
        thisslice = data['power',j].data.flatten()
        #curvature = diff(fftconvolve(thisslice,k,mode='same'),n=2)
        smoothed = fftconvolve(thisslice,k,mode='same') # I need this, so the noise doesn't break up my blocks
        #}}}
        #{{{ make lists to put the peaks in
        top_peak_x = []
        bottom_peak_x = []
        top_peak = []
        bottom_peak = []
        #}}}
        minind = 1
        #{{{ find the peaks for the segments above threshold
        peakmask = whereblocks(smoothed>smoothed.max()*threshold)
        indeces = np.r_[0:len(thisslice)]
        for peakset in peakmask: # peakset gives the indeces for a given slice
            if len(peakset)>minind:
                peak_ind = peakset[argmax(thisslice[peakset])]
                top_peak_x += [x[peak_ind]]
                top_peak += [thisslice[peak_ind]]
        #}}}
        #{{{ find the peaks for the segments below lower threshold
        peakmask = whereblocks(smoothed<smoothed.min()*threshold)
        for peakset in peakmask:
            if len(peakset)>minind:
                peak_ind = peakset[argmin(thisslice[peakset])]
                bottom_peak_x += [x[peak_ind]]
                bottom_peak += [thisslice[peak_ind]]
        #}}}
        if (not imageformat):
            plot(x,thisslice,color=cm.hsv(np.double(j)/np.double(nslices)),alpha=0.5)
            plot(bottom_peak_x,bottom_peak,'o',color=cm.hsv(np.double(j)/np.double(nslices)),alpha=0.5)
            plot(top_peak_x,top_peak,'o',color=cm.hsv(np.double(j)/np.double(nslices)),alpha=0.5)
        allpeaks_top += [top_peak]
        allpeaks_top_x += [top_peak_x]
        allpeaks_bottom += [bottom_peak]
        allpeaks_bottom_x += [bottom_peak_x]
    num_peaks = len(allpeaks_top_x[0])
    try:
        allpeaks_top_x = nddata(allpeaks_top_x,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    except:
        print(r'\begin{verbatim} If you have an error here, probably change smoothing (%0.2f) or threshold (%0.2f)\end{verbatim}'%(smoothing,threshold),'\n\n')
        plt.clf()
        for j in range(0,nslices):
            thisslice = data['power',j].data
            #curvature = diff(fftconvolve(thisslice,k,mode='same'),n=2)
            smoothed = fftconvolve(thisslice,k,mode='same') # I need this, so the noise doesn't break up my blocks
            plot(x,smoothed,alpha=0.1)
            peakmask = whereblocks(smoothed>smoothed.max()*threshold)
            for peakset in peakmask:
                plot(x[peakset],smoothed[peakset])
        lplot('error_plot'+figname+'.png',width=6)
        print(r'lengths: ',list(map(len,allpeaks_top_x)),'')
        return
    try:
        allpeaks_bottom_x = nddata(allpeaks_bottom_x,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    except:
        print(r'\begin{verbatim} If you have an error here, probably change smoothing (%0.2f) or threshold (%0.2f)\end{verbatim}'%(smoothing,threshold),'\n\n')
        plt.clf()
        for j in range(0,nslices):
            thisslice = data['power',j].data
            #curvature = diff(fftconvolve(thisslice,k,mode='same'),n=2)
            smoothed = fftconvolve(thisslice,k,mode='same') # I need this, so the noise doesn't break up my blocks
            plot(x,smoothed,alpha=0.1)
            peakmask = whereblocks(smoothed<smoothed.min()*threshold)
            for peakset in peakmask:
                plot(x[peakset],smoothed[peakset])
        lplot('error_plot'+figname+'.png',width=6)
        print(r'\begin{verbatim}lengths: ',list(map(len,allpeaks_top_x)),r'\end{verbatim}')
        return
    allpeaks_top = nddata(allpeaks_top,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    allpeaks_bottom = nddata(allpeaks_bottom,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    if imageformat:
        image(data.data,x=x,y=np.r_[0:len(powerseries)])
        plot(np.r_[0:len(powerseries)],allpeaks_top_x.data)
        #plot(np.r_[0:shape(data.data)[1]],allpeaks_bottom_x.data)
        lplot('esr_dataset'+figname+'.png',width=6,grid=False)
    else:
        lplot('esr_dataset'+figname+'.png',width=6)
    print('\n\n')
    #{{{ peak to peak
    peaktopeak = allpeaks_bottom_x - allpeaks_top_x
    peaktopeak_squared = peaktopeak.copy()
    peaktopeak.labels(['power'],[np.sqrt(powerseries)])
    peaktopeak.rename('power','$B_1$ / arb')
    plot(peaktopeak,'.-',nosemilog=True)
    plt.ylabel(r'$\Delta B_{pp}$')
    lplot('esr_dataset'+figname+'_pp.pdf')
    #{{{ linearity test
    peaktopeak_squared.data = peaktopeak_squared.data**2
    peaktopeak_squared.labels(['power'],[powerseries])
    peaktopeak_squared.rename('power',r'$p$ / $mW$')
    plot(peaktopeak_squared,'.-',nosemilog=True)
    plt.ylabel(r'$\Delta B_{pp}^2\propto s^{-1}$')
    lplot('esr_dataset'+figname+'_pp2.pdf')
    #}}}
    #}}}
    print('\n\n')
    #{{{ height
    height = (allpeaks_top - allpeaks_bottom)/2
    height_n23 = height.copy()
    height.labels(['power'],[np.sqrt(powerseries)])
    height.rename('power','$B_1$ / arb')
    plot(height,'.-',nosemilog=True)
    plt.ylabel(r"$y'_m$")
    lplot('esr_dataset'+figname+'_height.pdf')
    #{{{linearity test
    b1 = ndshape(height_n23)
    b1['peak'] = 1
    b1 = b1.alloc()
    #b1['power',:] = powerseries.copy().reshape(-1,1)
    b1.data = np.sqrt(powerseries).copy().reshape(b1.data.shape)
    height_n23 = height_n23/b1
    height_n23.data = height_n23.data**(-2./3.)
    height_n23.labels(['power'],[powerseries])
    newpname = r'$p$ / $mW$'
    height_n23.rename('power',newpname)
    height_n23_avg = height_n23.copy()
    height_n23_avg.mean('peak')
    if show_avg == True:
        plot(height_n23_avg/hn23adjustment,'.',nosemilog=True)
    else:
        plot(height_n23/hn23adjustment,'.',nosemilog=True)
    maxp = lambda x: x == x.max()
    plot(np.r_[0.0,height_n23_avg[newpname,maxp].getaxis(newpname)],np.r_[1.0,height_n23_avg[newpname,maxp].data[0]/hn23adjustment],'k-',linewidth=2,alpha=0.3,nosemilog=True)
    plt.ylabel(r"$\propto\frac{y'_m}{B_1}^{-2/3}\propto \frac{1}{1-s_{ESR}}$")
    lplot('esr_dataset'+figname+'_hn23.pdf',width = 3.5)
    #}}}
    #}}}
#}}}
#{{{ dnp
def standard_noise_comparison(name,path = 'franck_cnsi/nmr/', data_subdir = 'reference_data',expnos = [3]):
    print('\n\n')
    # noise tests
    close(1)
    figure(1,figsize=(16,8))
    v = save_data();our_calibration = np.double(v['our_calibration']);cnsi_calibration = np.double(v['cnsi_calibration'])
    calibration = cnsi_calibration*np.sqrt(50.0/10.0)*np.sqrt(50.0/40.0)
    path_list = []
    explabel = []
    noiseexpno = []
    signalexpno = []
    plotlabel = name+'_noise'
    #
    path_list += [getDATADIR()+'%s/nmr/popem_4mM_5p_pct_110610/'%data_subdir]
    explabel += ['control without shield']
    noiseexpno += [3] # 3 is the noise scan 2 is the reference
    path_list += [getDATADIR()+'%s/nmr/noisetest100916/'%data_subdir] + [getDATADIR()+path+name+'/']*len(expnos)
    explabel += ['']+[r'$\mathbf{this experiment}$']*len(expnos)
    noiseexpno += [2]+expnos # 3 is the noise scan 2 is the reference
    #
    mask_start = -1e6
    mask_stop = 1e6
    ind = 0
    smoothing = 5e3
    for j in range(0,1): # for multiple plots $\Rightarrow$ add in j index below if this is what i want
       figure(1)
       ind += 1
       legendstr = []
       linelist = []
       plt.subplot(121) # so that legend will fit
       for k in range(0,len(noiseexpno)):
          retval = plot_noise(path_list[k],noiseexpno[k],calibration,mask_start,mask_stop,smoothing = smoothing, both = False,retplot = True)
          linelist += retval[0]
          legendstr.append('\n'.join(textwrap.wrap(explabel[k]+':'+retval[1][0],50))+'\n')
       plt.ylabel(r'$\Omega$')
       titlestr = 'Noise scans (smoothed %0.2f $kHz$) for CNSI spectrometer\n'%(smoothing/1e3)
       plt.title(titlestr+r'$n V$ RG/ disk units = %0.3f, mask (%0.3f,%0.3f)'%(calibration*1e9,mask_start,mask_stop))
       ax = plt.gca()
       ylims = list(ax.get_ylim())
       #gridandtick(plt.gca(),formatonly = True)
       gridandtick(plt.gca(),logarithmic = True)
       plt.subplot(122)
       grid(False)
       lg = autolegend(linelist,legendstr)
       ax = plt.gca()
       ax.get_xaxis().set_visible(False)
       ax.get_yaxis().set_visible(False)
       for x in list(ax.spines.values()):
           x.set_visible(False)
       lplot('noise'+plotlabel+'_%d.pdf'%ind,grid=False,width=5,gensvg=True)
       print('\n\n')
       figure(2)
       legendstr = []
       for k in range(0,len(signalexpno)):
          data = load_file(dirformat(path_list[k])+'%d'%noiseexpno[k],calibration=calibration)
          data.ft('t2',shift = True)
          x = data.getaxis('t2')
          data['t2',abs(x)>1e3] = 0
          data.ift('t2',shift = True)
          plot(abs(data['t2',0:300])*1e9)
          xlabel('signal / $nV$')
          legendstr += [explabel[k]]
       if len(signalexpno)>0:
           autolegend(legendstr)
           lplot('signal'+plotlabel+'_%d.pdf'%ind,grid=False)
       if (ind % 2) ==  0:
          print('\n\n')
#}}}
