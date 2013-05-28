import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matlablike import *
from string import rstrip
from scipy.io import savemat,loadmat
from os.path import exists as path_exists
from scipy.optimize import leastsq
from scipy.signal import fftconvolve
from tables import openFile
from datetime import datetime
from time import mktime
from datadir import grab_data_directory
from nmr import *
import re

golden_ratio = (1.0 + sqrt(5))/2.0

def dprint(*stuff):
    print '\n\nDEBUG',' '.join(map(repr(stuff))),'\n\n'
def thisjobname():
    if path_exists('pythonjobname.txt'):
        fp = open('pythonjobname.txt')
        retval = rstrip(fp.read())
        fp.close()
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
    def show(self,string,**kwargs):
        'latexify the series of figures, where "string" gives the base file name'
        print '\n\n'
        if not len(kwargs):
            kwargs = {}
        for j,figname in enumerate(self.figurelist):
            if type(figname) is dict:
                kwargs.update(figname)
                if 'print_string' in kwargs:
                    print '\n\n'
                    print kwargs.pop('print_string')
                    print '\n\n'
            else:
                figure(j+1)
                sep = ''
                if len(string)>0:
                    sep = '_'
                if self.env == 'test':
                    print 'not running lplot'
                else:
                    try:
                        lplot(figname+sep+string,**kwargs)
                    except:
                        raise CustomError('self.figurelist = ',self.figurelist,'figname = ',figname,'while string = ',string)
        while len(self.figurelist) > 0:
            self.figurelist.pop(-1)
        return
    def obs(self,*args):
        self.text(obs_repr(*args))
        return
    def obsn(self,*args):
        self.text(obs_repr(*args)+'\n\n')
        return
def lplotfigures(figurelist,string,**kwargs):
    'obsolete, use the class!'
    print '\n\n'
    if not len(kwargs):
        kwargs = {}
    for j,figname in enumerate(figurelist):
        if type(figname) is dict:
            kwargs.update(figname)
            if 'print_string' in kwargs:
                print '\n\n'+kwargs.pop('print_string')+'\n\n'
        else:
            figure(j+1)
            try:
                lplot(figname+string,**kwargs)
            except:
                raise CustomError('figurelist = ',figurelist,'figname = ',figname,'while string = ',string)
    while len(figurelist) > 0:
        figurelist.pop(-1)
    return
def figlistini(first_figure):
    r"""processes a figure list argument:
    typically, you want to have a figure_list keyword argument for every function, which is by default set to None, then call this on the argument -- it always returns a figure list, creating a new one if required
    similarly, somewhere I have another guy that processes the output, so that if it's set to None, it will by default dump and show the figure list,
    and not return a figure list in the output"""
    verbose = False
    if verbose: print lsafe('DEBUG: initialize figlist')
    if first_figure == None:
        if verbose: print lsafen('empty')
        return figlistl() # note that this is changed, and it used to be an empty list
    else:
        if verbose: print lsafen(first_figure.figurelist)
        return first_figure
def figlisterr(figurelist,*args,**kwargs):
    if 'basename' in kwargs.keys():
        basename = kwargs['basename']
    else:
        basename = thisjobname()
    print lsafen("Trying to plot the figurelist",figurelist)
    lplotfigures(figurelist,basename+'errplot.pdf')
    return args
def figlistret(first_figure,figurelist,*args,**kwargs):
    if 'basename' in kwargs.keys():
        basename = kwargs['basename']
    else:
        basename = thisjobname()
    if first_figure is None:
        lplotfigures(figurelist,basename+'.pdf')
        return args
    else:
        args += (figurelist,)
        if len(args) == 1:
            return args[0]
        else:
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
        return not any(map(lambda x: len(x)>1 and not (x[0] == '{' or x[0] == '\\'),out))
    return [r'\ensuremath{'+v.replace('\\\\','\\')+'}' if ismath(v) else v.replace('_',' ') for v in recnames]
def lrecordarray_broken(recordlist,rows=30,numwide= 5):
    width = 1.0/numwide
    for j in range(0,len(recordlist),rows):
        print r'\begin{minipage}[t]{%0.3f\linewidth}'%width,
        lrecordarray(recordlist[j:j+rows],resizebox = True)
        print r'\end{minipage}',
    if j+rows<len(recordlist):
        print r'\begin{minipage}[t]{%0.3\linewidth}'%width,
        lrecordarray(recordlist[j+rows:],resizebox = True)
        print r'\end{minipage}',
def lrecordarray(recordlist,columnformat = True,smoosh = True,multi = True,resizebox = False,showpipe = True,return_only=False,format = '%0.3f',std_sf = None):
    r'''generate latex representation of a structured array
    if set to True, resizebox will automatically scale down the table so it fits on the page (but it will also scale up a small table to fit the width of the page)
    resizebox can also be a fractional number, so that it is resized to a fraction of the page'''
    error_names = [j for j in recordlist.dtype.names if j[-6:] == '_ERROR']
    data_names = [j for j in recordlist.dtype.names if j[-6:] != '_ERROR']
    data_has_error_names = [j for j in data_names if j+'_ERROR' in error_names]
    if len(data_has_error_names)<len(error_names):
        print r'{\color{red}Warning! unassigned error names!',lsafe(set(error_names)-set(map(lambda x: x+'_ERROR',data_has_error_names))),'}\n\n'
    recordlist_errors = recordlist[error_names]
    recordlist = recordlist[data_names]
    final_retval = []
    if len(recordlist) == 0:
        print r'{\color{red}This array is empty!!}'
        return
    # previously documented
    def this_format_function(x,error_val = None):
        if type(x) is str_:
            if error_val is None:
                return str(x)
            else:
                return str(x)+'\\pm '+str(error_val)
        elif std_sf is not None:
            #just return two significant figures
            number_sf = std_sf
            if error_val is None or error_val == 0.0:
                highest_significant = floor(log10(x))
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
                highest_significant = floor(log10(error_val))
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
            if type(x) is ndarray:
                if x.shape != error.shape:
                    raise ValueError("shape of the error does not match the shape of the value/matrix")
        if type(x) is ndarray:
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
            if type(x) in [str_,str]:
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
        print '\n'.join(final_retval)
def lplot(fname,width=0.33,figure=False,dpi=72,grid=False,alsosave=None,gensvg = False,print_string = None,centered = False,legend = False,equal_aspect = False,autopad = True,bytextwidth = None,showbox = True,outer_legend = False,boundaries = True):
    '''
    used with python.sty instead of savefig
    
    by default, if width is less than 1, 
    it's interpreted as bytextwidth = True (i.e. width given as a fraction of the linewidth)
    if it's greater than, the width is interpreted in inches.
    '''
    if width < 1.0:
        bytextwidth = True
    else:
        bytextwidth = False
    if print_string is not None:
        print print_string
    fname = r'auto_figures/'+fname
    if alsosave != None:
        alsosave = r'auto_figures/'+alsosave
    if gensvg == True:
        alsosave = fname.replace('.pdf','.svg')
    if grid:
        gridandtick(gca())
    fig = gcf()
    ax = gca()
    if equal_aspect:
        ax.set_aspect('equal')
    fig.autofmt_xdate()
    if autopad: autopad_figure(centered = centered)
    if legend:
        autolegend()
    if outer_legend:
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if not boundaries:
        ax = gca()
        for j in ax.spines.keys():
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
        savefig(fname,dpi=dpi)
    except ValueError,exc_string:
        if exc_string.find('finite numbers') > -1:
            raise ValueError("It gives this error because you're trying to do a bar graph with zero width")
        else:
            raise ValueError(exc_string)
    if alsosave != None:
        savefig(alsosave,dpi=dpi)
    if figure:
        print r"""
        \begin{figure}[h]
        \end{figure}
        """
    if bytextwidth:
        if showbox:
            mpwidth = r'%0.2f\linewidth'%width
            figwidth = r'\linewidth'
        else:
            figwidth = r'%0.2f\linewidth'%width
    else:
        mpwidth = r'%0.2fin'%width
        figwidth = mpwidth
    if showbox:
        print r'''\mbox{\begin{minipage}{%s}'''%mpwidth
        if alsosave != None:
            print r'also saved \fn{%s}'%alsosave+'\n\n'
    print r'\includegraphics[width=%s]{%s}'%(figwidth,fname.replace(r'auto_figures',r'\autofiguredir'))
    if showbox:
        print '\n\n'+r'\hrulefill'+'\n\n'
        print r'''{\color{red}{\tiny %s}:}\begin{tiny}\fn{%s}\end{tiny}'''%('file:',fname)
        print '\n\n'+r'\hrulefill'+'\n\n'
        print r'''\end{minipage}
}''',
    clf()
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
        if type(j) in [double,float,float32]:
            if log10(j) > 3 or log10(j) < -3:# greater than a thousand or less than a thousandth
                power = floor(log10(j))
                j /= 10**power
                thisstr += r'$%0.3f\times 10^{%d}$'%(j,power)
            else:
                if log10(j)<-1:
                    temp = '%%0.%df'%(floor(-log10(j))+3)
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
    print obs_repr(*arg)
    return
def obsn(*arg):
    obs(*arg)
    print '\n\n'
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
    for k,v in mydict.iteritems():
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
    print r'$%0.1f\;GHz$ @ '%(freq/1e9),dp((optdistance+119e-3)*1e3,1),r'$mm$'
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
    data = prospa_load_datafile(filename,dims=2)
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
        print '\n\n'
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
    print r'\fn{%s %d}'%(exp,number)+'\n\n'
    cpmgseries(filename,exp+thisjobname(),tau,alpha,alphaselect)
#{{{ esr_saturation
def esr_saturation(file,powerseries,smoothing=0.2,threshold=0.8,figname = None,hn23adjustment = 1.,show_avg = False):
    if figname == None:
        figname = thisjobname()
    print '\n\n\\fn{',file,'}\n\n'
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
        indeces = r_[0:len(thisslice)]
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
            plot(x,thisslice,color=cm.hsv(double(j)/double(nslices)),alpha=0.5)
            plot(bottom_peak_x,bottom_peak,'o',color=cm.hsv(double(j)/double(nslices)),alpha=0.5)
            plot(top_peak_x,top_peak,'o',color=cm.hsv(double(j)/double(nslices)),alpha=0.5)
        allpeaks_top += [top_peak]
        allpeaks_top_x += [top_peak_x]
        allpeaks_bottom += [bottom_peak]
        allpeaks_bottom_x += [bottom_peak_x]
    num_peaks = len(allpeaks_top_x[0])
    try:
        allpeaks_top_x = nddata(allpeaks_top_x,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    except:
        print r'\begin{verbatim} If you have an error here, probably change smoothing (%0.2f) or threshold (%0.2f)\end{verbatim}'%(smoothing,threshold),'\n\n'
        clf()
        for j in range(0,nslices):
            thisslice = data['power',j].data
            #curvature = diff(fftconvolve(thisslice,k,mode='same'),n=2)
            smoothed = fftconvolve(thisslice,k,mode='same') # I need this, so the noise doesn't break up my blocks
            plot(x,smoothed,alpha=0.1)
            peakmask = whereblocks(smoothed>smoothed.max()*threshold)
            for peakset in peakmask:
                plot(x[peakset],smoothed[peakset])
        lplot('error_plot'+figname+'.png',width=6)
        print r'lengths: ',map(len,allpeaks_top_x),''
        return
    try:
        allpeaks_bottom_x = nddata(allpeaks_bottom_x,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    except:
        print r'\begin{verbatim} If you have an error here, probably change smoothing (%0.2f) or threshold (%0.2f)\end{verbatim}'%(smoothing,threshold),'\n\n'
        clf()
        for j in range(0,nslices):
            thisslice = data['power',j].data
            #curvature = diff(fftconvolve(thisslice,k,mode='same'),n=2)
            smoothed = fftconvolve(thisslice,k,mode='same') # I need this, so the noise doesn't break up my blocks
            plot(x,smoothed,alpha=0.1)
            peakmask = whereblocks(smoothed<smoothed.min()*threshold)
            for peakset in peakmask:
                plot(x[peakset],smoothed[peakset])
        lplot('error_plot'+figname+'.png',width=6)
        print r'\begin{verbatim}lengths: ',map(len,allpeaks_top_x),'\end{verbatim}'
        return
    allpeaks_top = nddata(allpeaks_top,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    allpeaks_bottom = nddata(allpeaks_bottom,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    if imageformat:
        image(data.data,x=x,y=r_[0:len(powerseries)])
        plot(r_[0:len(powerseries)],allpeaks_top_x.data)
        #plot(r_[0:shape(data.data)[1]],allpeaks_bottom_x.data)
        lplot('esr_dataset'+figname+'.png',width=6,grid=False)
    else:
        lplot('esr_dataset'+figname+'.png',width=6)
    print '\n\n'
    #{{{ peak to peak
    peaktopeak = allpeaks_bottom_x - allpeaks_top_x
    peaktopeak_squared = peaktopeak.copy()
    peaktopeak.labels(['power'],[sqrt(powerseries)])
    peaktopeak.rename('power','$B_1$ / arb')
    plot(peaktopeak,'.-',nosemilog=True)
    ylabel(r'$\Delta B_{pp}$')
    lplot('esr_dataset'+figname+'_pp.pdf')
    #{{{ linearity test
    peaktopeak_squared.data = peaktopeak_squared.data**2
    peaktopeak_squared.labels(['power'],[powerseries])
    peaktopeak_squared.rename('power',r'$p$ / $mW$')
    plot(peaktopeak_squared,'.-',nosemilog=True)
    ylabel(r'$\Delta B_{pp}^2\propto s^{-1}$')
    lplot('esr_dataset'+figname+'_pp2.pdf')
    #}}}
    #}}}
    print '\n\n'
    #{{{ height
    height = (allpeaks_top - allpeaks_bottom)/2
    height_n23 = height.copy()
    height.labels(['power'],[sqrt(powerseries)])
    height.rename('power','$B_1$ / arb')
    plot(height,'.-',nosemilog=True)
    ylabel(r"$y'_m$")
    lplot('esr_dataset'+figname+'_height.pdf')
    #{{{linearity test
    b1 = ndshape(height_n23)
    b1['peak'] = 1
    b1 = b1.alloc()
    #b1['power',:] = powerseries.copy().reshape(-1,1)
    b1.data = sqrt(powerseries).copy().reshape(b1.data.shape)
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
    plot(r_[0.0,height_n23_avg[newpname,maxp].getaxis(newpname)],r_[1.0,height_n23_avg[newpname,maxp].data[0]/hn23adjustment],'k-',linewidth=2,alpha=0.3,nosemilog=True)
    ylabel(r"$\propto\frac{y'_m}{B_1}^{-2/3}\propto \frac{1}{1-s_{ESR}}$")
    lplot('esr_dataset'+figname+'_hn23.pdf',width = 3.5)
    #}}}
    #}}}
#}}}
#{{{ dnp
def standard_noise_comparison(name,path = 'franck_cnsi/nmr/', data_subdir = 'reference_data',expnos = [3]):
    print '\n\n'
    # noise tests
    close(1)
    figure(1,figsize=(16,8))
    v = save_data();our_calibration = double(v['our_calibration']);cnsi_calibration = double(v['cnsi_calibration'])
    calibration = cnsi_calibration*sqrt(50.0/10.0)*sqrt(50.0/40.0)
    path_list = []
    explabel = []
    noiseexpno = []
    signalexpno = []
    plotlabel = name+'_noise'
    #
    path_list += [DATADIR+'%s/nmr/popem_4mM_5p_pct_110610/'%data_subdir]
    explabel += ['control without shield']
    noiseexpno += [3] # 3 is the noise scan 2 is the reference
    path_list += [DATADIR+'%s/nmr/noisetest100916/'%data_subdir] + [DATADIR+path+name+'/']*len(expnos)
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
       subplot(121) # so that legend will fit
       for k in range(0,len(noiseexpno)):
          retval = plot_noise(path_list[k],noiseexpno[k],calibration,mask_start,mask_stop,smoothing = smoothing, both = False,retplot = True)
          linelist += retval[0]
          legendstr.append('\n'.join(textwrap.wrap(explabel[k]+':'+retval[1][0],50))+'\n')
       ylabel(r'$\Omega$')
       titlestr = 'Noise scans (smoothed %0.2f $kHz$) for CNSI spectrometer\n'%(smoothing/1e3)
       title(titlestr+r'$n V$ RG/ disk units = %0.3f, mask (%0.3f,%0.3f)'%(calibration*1e9,mask_start,mask_stop))
       ax = gca()
       ylims = list(ax.get_ylim())
       #gridandtick(gca(),formatonly = True)
       gridandtick(gca(),logarithmic = True)
       subplot(122)
       grid(False)
       lg = autolegend(linelist,legendstr)
       ax = gca()
       ax.get_xaxis().set_visible(False)
       ax.get_yaxis().set_visible(False)
       map((lambda x: x.set_visible(False)),ax.spines.values())
       lplot('noise'+plotlabel+'_%d.pdf'%ind,grid=False,width=5,gensvg=True)
       print '\n\n'
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
          print '\n\n'
#}}}
