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
                    print '\n\n'+kwargs.pop('print_string')+'\n\n'
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
                print '\n\n'+kwargs.pop('print_string')+'\n\n'
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
def lrecordarray(recordlist,columnformat = True,smoosh = True,multi = True):
    r'generate latex representation of a structured array'
    # previously documented
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
            retval = ''
            for thistuple in listoftuples:
                if thistuple[0] == thistuple[1]:
                    retval = retval + r'\cline{%d-%d}'%tuple([thistuple[0]+1]*2)
                elif thistuple[0] == 0 and thistuple[1] == len(recnames)-1:
                    retval = retval + r'\hline'
                else:
                    retval = retval + r'\cline{%d-%d}'%tuple(map((lambda x: x+1),thistuple))
            return retval
        #{{{ first, store the strings for the data
        for j in range(0,len(recordlist)):
            datastrings.append([lsafe(str(recordlist[v][j])) for v in recnames])
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
        print r'\begin{tabular}{','|'.join(colheadings),'}'
        print ' & '.join(recordstrings),r'\\ \hline\hline'+'\n'
        print alltext
    else:
        print r'\begin{tabular}{',''.join(['r']+['l']*reclen),'}'
        for v in recnames:
            print r'\ensuremath{',v,'}$=$ &',' & '.join(map(lambda x: lsafe(str(x)),list(recordlist[v]))),r'\\'
        print r'\end{tabular}'
def lplot(fname,width=3,figure=False,dpi=300,grid=False,alsosave=None,gensvg = False,print_string = None,centered = False,legend = False,equal_aspect = False,autopad = True,bytextwidth = False):
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
    if bytextwidth:
        mpwidth = r'%0.2f\textwidth'%width
        figwidth = r'\textwidth'
    else:
        mpwidth = r'%0.2fin'%width
        figwidth = mpwidth
    print r'''\fbox{
    \begin{minipage}{%s}
    {\color{red}{\tiny %s}:}\begin{tiny}\fn{%s}'''%(mpwidth,thisjobname(),fname)
    if alsosave != None:
        print r'also saved \fn{%s}'%alsosave
    print '\n\n'+r'\hrulefill'+'\n\n'
    print r'''\includegraphics[width=%s]{%s}
    \end{tiny}\end{minipage}
    }'''%(figwidth,fname)
    clf()
    return
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
def standard_noise_comparison(name):
    print '\n\n'
    # noise tests
    close(1)
    figure(1,figsize=(16,8))
    v = save_data();our_calibration = double(v['our_calibration']);cnsi_calibration = double(v['cnsi_calibration'])
    calibration = cnsi_calibration*sqrt(50.0/10.0)*sqrt(50.0/40.0)
    path = []
    explabel = []
    noiseexpno = []
    signalexpno = []
    plotlabel = name+'_noise'
    #
    path += [DATADIR+'cnsi_data/popem_4mM_5p_pct_110610/']
    explabel += ['control without shield']
    noiseexpno += [3] # 3 is the noise scan 2 is the reference
    path += [DATADIR+'cnsi_data/noisetest100916/'] + [DATADIR+'cnsi_data/'+name+'/']
    explabel += ['',r'$\mathbf{this experiment}$']
    noiseexpno += [2,3] # 3 is the noise scan 2 is the reference
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
          retval = plot_noise(path[k],noiseexpno[k],calibration,mask_start,mask_stop,smoothing = smoothing, both = False,retplot = True)
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
          data = load_file(dirformat(path[k])+'%d'%noiseexpno[k],calibration=calibration)
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
