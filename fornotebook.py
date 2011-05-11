import matplotlib; matplotlib.use('Agg')
from matlablike import *
from string import rstrip
from scipy.io import savemat,loadmat
from os.path import exists as path_exists
from scipy.optimize import leastsq
from scipy.signal import fftconvolve
from tables import openFile
from datetime import datetime
from time import mktime
from nmr import *
import re

golden_ratio = (1.0 + sqrt(5))/2.0

def dprint(*stuff):
    print '\n\nDEBUG',' '.join(map(repr(stuff))),'\n\n'
def thisjobname():
    fp = open('pythonjobname.txt')
    retval = rstrip(fp.read())
    fp.close()
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
    def show(self,string,numwide = 2,**kwargs):
        'latexify the series of figures, where "string" gives the base file name'
        print '\n\n'
        if not len(kwargs):
            kwargs = {}
        for j,figname in enumerate(self.figurelist):
            if type(figname) is dict:
                kwargs.update(figname)
                if 'print_string' in kwargs:
                    print kwargs.pop('print_string')
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
                if (j+1) % numwide == 0:
                    print '\n\n'
        while len(self.figurelist) > 0:
            self.figurelist.pop(-1)
        return
def lplotfigures(figurelist,string,numwide = 2,**kwargs):
    'obsolete, use the class!'
    print '\n\n'
    if not len(kwargs):
        kwargs = {}
    for j,figname in enumerate(figurelist):
        if type(figname) is dict:
            kwargs.update(figname)
            if 'print_string' in kwargs:
                print kwargs.pop('print_string')
        else:
            figure(j+1)
            try:
                lplot(figname+string,**kwargs)
            except:
                raise CustomError('figurelist = ',figurelist,'figname = ',figname,'while string = ',string)
            if (j+1) % numwide == 0:
                print '\n\n'
    while len(figurelist) > 0:
        figurelist.pop(-1)
    return
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
def lrecordarray(recordlist,columnformat = True):
    reclen = len(recordlist)
    recnames = recordlist.dtype.names
    if columnformat:
        print r'\begin{tabular}{',''.join(['c']*len(recnames)),'}'
        recordstrings = [r'\ensuremath{'+v.replace('\\\\','\\')+'}' for v in recnames]
        print ' & '.join(recordstrings),r'\\'# \hline'
        for j in range(0,len(recordlist)):
            datastrings = [lsafe(str(recordlist[v][j])) for v in recnames]
            print ' & '.join(datastrings),r'\\'# \hline'
        print r'\end{tabular}'
    else:
        print r'\begin{tabular}{',''.join(['r']+['l']*reclen),'}'
        for v in recnames:
            print r'\ensuremath{',v,'}$=$ &',' & '.join(map(lambda x: lsafe(str(x)),list(recordlist[v]))),r'\\'
        print r'\end{tabular}'
def lplot(fname,width=3,figure=False,dpi=300,grid=False,alsosave=None,gensvg = False,print_string = None,centered = False,legend = False,equal_aspect = False,autopad = True):
    '''
    used with python.sty instead of savefig
    '''
    if print_string is not None:
        print print_string
    fname = 'auto_figures/'+fname
    if alsosave != None:
        alsosave = 'auto_figures/'+alsosave
    if gensvg == True:
        alsosave = fname.replace('.pdf','.svg')
    if grid:
        gridandtick(gca())
    fig = gcf()
    ax = gca()
    if equal_aspect:
        ax.set_aspect('equal')
    fig.autofmt_xdate()
    if legend:
        autolegend()
    if autopad: autopad_figure(centered = centered)
    savefig(fname,dpi=dpi)
    if alsosave != None:
        savefig(alsosave,dpi=dpi)
    if figure:
        print r"""
        \begin{figure}[h]
        \end{figure}
        """
    print r'''\begin{minipage}{%0.2fin}
    \fn{%s: %s}'''%(width,thisjobname(),fname)
    if alsosave != None:
        print r'also saved \fn{%s}'%alsosave
    print r'''

    \includegraphics[width=%0.2fin]{%s}
    \end{minipage}'''%(width,fname)
    clf()
    return
def addlabels(labelstring,x,y,labels):
    for j in range(0,len(labels)):
        text(x[j],y[j],labelstring%labels[j],alpha=0.5,color='g',ha='left',va='top',rotation=0)
def ordplot(x,y,labels,formatstring):
        order = argsort(x)
        plot(x[order],y[order],'o-')
        newlabels=[]
        for j in range(0,len(order)): 
            newlabels += [labels[order[j]]]
        addlabels(formatstring,x[order],y[order],newlabels)
def obs(*arg):
    print r'\o{'
    for j in arg:
        print j,
    print r'}'
def save_data(inputdata={},mat_file = 'data.mat'):
    # concatenate data to data already in data.mat file
    if path_exists(mat_file):
        #data = loadmat(mat_file,struct_as_record=True)
        data = loadmat(mat_file)
    else:
        data = {}
    if not(inputdata=={}):
        try:
            data.update(inputdata)
        except:
            raise CustomError('trying to update',data,'with',inputdata)
        try:
            savemat(mat_file,data)
            #print "DEBUG, saved",data,'to',mat_file
        except:
            raise CustomError('trying to write',data,'to',mat_file)
    return data
def save_local(inputdata={},verb = True):
    data = save_data(inputdata,mat_file = 'local.mat')
    if verb:
        for k,v in inputdata.items():
            obs(k.replace('_',r'\_'),r'$\Rightarrow$',v.replace('_',r'\_'))
    return data
def clear_local(inputdata=[]):
    obs(r'{\it clear local}'+'\n\n')
    if path_exists('local.mat'):
        os.unlink('local.mat')
def save_variable(variable,content,disp=True):
    if path_exists('data.mat'):
        #data = loadmat('data.mat',struct_as_record=True)
        data = loadmat('data.mat')
    else:
        data = {}
    data.update({variable:content})
    savemat('data.mat',data)
    if disp:
        obs(variable.replace('_',r'\_'),'=',dp(content,5))
    return data
#{{{ specific convenience functions
def calcdistance(freq,optdistance):
    print r'$%0.1f\;GHz$ @ '%(freq/1e9),dp((optdistance+119e-3)*1e3,1),r'$mm$'
def calcfield(elfreq,elratio = 28.113,nmrelratio = 1.5167):
    obs(elfreq,'$GHz$,',dp(elratio,5),'$GHz/T\;',dp(nmrelratio,6),'\;ppt\;\Rightarrow$',dp(elfreq/elratio*1e4,3),'$G$',dp(elfreq*nmrelratio,5),'$MHz$')
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
    filename = '/mnt/esr/Magritek/john/'+exp+'/%d/'%number
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
#}}}
