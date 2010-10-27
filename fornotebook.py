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

def thisjobname():
    fp = open('pythonjobname.txt')
    retval = rstrip(fp.read())
    fp.close()
    return retval
def lplotfigures(figurelist,string,numwide = 2,**kwargs):
    print '\n\n'
    for j,figname in enumerate(figurelist):
        figure(j+1)
        lplot(figname+string,grid = False,**kwargs)
        if (j+1) % numwide == 0:
            print '\n\n'
    return []
def lplot(fname,width=3,figure=False,dpi=100,grid=False,alsosave=None,gensvg = False):
    '''
    used with python.sty instead of savefig
    '''
    fname = 'auto_figures/'+fname
    if alsosave != None:
        alsosave = 'auto_figures/'+alsosave
    if gensvg == True:
        alsosave = fname.replace('.pdf','.svg')
    if grid:
        gridandtick(gca())
    fig = gcf()
    fig.autofmt_xdate()
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
#{{{ functions to pull info out of file names
def pull_chemical_from_filename(name):
    a = re.compile('(.*?)_?([0-9]+)([mu]M)')
    b = re.compile('(.*?)(_|[0-9])')
    m = a.search(name)
    if m:
        g = m.groups()
        chemical_name = g[0]
        print "\n\nfound chemical name \\verb|",g[0],'|\n\n'
    elif b.search(name):
        chemical_name = b.search(name).groups()[0]
        print "\n\nfound chemical name \\verb|",chemical_name,'|\n\n'
    else:
        chemical_name = 'not_found'
        print "\n\nfound no chemical name \n\n"
    return chemical_name
def pull_concentration_from_filename(name):
    a = re.compile('([0-9]+)_?([mu]M)')
    m = a.search(name)
    if m:
        g = m.groups()
        if g[1]=='uM':
            concentration = 1e-6
        elif g[1]=='mM':
            concentration = 1e-3
        concentration *= double(int(g[0]))
    else:
        concentration = 0
    return concentration
def pull_date_from_filename(name):
    a = re.compile('([0-9]{2})([0-9]{2})([0-9]{2})')
    m = a.search(name)
    g = m.groups()
    date_of_file = mktime(datetime(int(g[0])+2000,int(g[1]),int(g[2])).timetuple())
    return date_of_file
#}}}
#{{{ helper function for HDF5 search
def gensearch(labelname,format,value,precision):
    searchstring_high = '(%s < %s + (%s))'%tuple([labelname]+[format]*2)
    #print "\n\nDEBUG check format:\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_high = searchstring_high%(value,precision)
    #print "\n\nDEBUG after substitution with",value,precision,":\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_low = '(%s > %s - (%s))'%tuple([labelname]+[format]*2)
    searchstring_low = searchstring_low%(value,precision)
    return searchstring_low + ' & ' + searchstring_high
#}}}
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
def esr_saturation(file,powerseries,smoothing=0.2,threshold=0.8):
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
        lplot('error_plot'+thisjobname()+'.png',width=6)
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
        lplot('error_plot'+thisjobname()+'.png',width=6)
        print r'\begin{verbatim}lengths: ',map(len,allpeaks_top_x),'\end{verbatim}'
        return
    allpeaks_top = nddata(allpeaks_top,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    allpeaks_bottom = nddata(allpeaks_bottom,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
    if imageformat:
        image(data.data,x=x,y=r_[0:len(powerseries)])
        plot(r_[0:len(powerseries)],allpeaks_top_x.data)
        #plot(r_[0:shape(data.data)[1]],allpeaks_bottom_x.data)
        lplot('esr_dataset'+thisjobname()+'.png',width=6,grid=False)
    else:
        lplot('esr_dataset'+thisjobname()+'.png',width=6)
    print '\n\n'
    #{{{ peak to peak
    peaktopeak = allpeaks_bottom_x - allpeaks_top_x
    peaktopeak_squared = peaktopeak.copy()
    peaktopeak.labels(['power'],[sqrt(powerseries)])
    peaktopeak.rename('power','$B_1$ / arb')
    plot(peaktopeak,'.-',nosemilog=True)
    ylabel(r'$\Delta B_{pp}$')
    lplot('esr_dataset'+thisjobname()+'_pp.pdf')
    #{{{ linearity test
    peaktopeak_squared.data = peaktopeak_squared.data**2
    peaktopeak_squared.labels(['power'],[powerseries])
    peaktopeak_squared.rename('power',r'$p$ / $mW$')
    plot(peaktopeak_squared,'.-',nosemilog=True)
    ylabel(r'$\Delta B_{pp}^2\propto s^{-1}$')
    lplot('esr_dataset'+thisjobname()+'_pp2.pdf')
    #}}}
    #}}}
    print '\n\n'
    #{{{ height
    height = allpeaks_top# - allpeaks_bottom
    height_n23 = height.copy()
    height.labels(['power'],[sqrt(powerseries)])
    height.rename('power','$B_1$ / arb')
    plot(height,'.-',nosemilog=True)
    ylabel(r"$y'_m$")
    lplot('esr_dataset'+thisjobname()+'_height.pdf')
    #{{{linearity test
    b1 = ndshape(height_n23)
    b1['peak'] = 1
    b1 = b1.alloc()
    #b1['power',:] = powerseries.copy().reshape(-1,1)
    b1.data = powerseries.copy().reshape(b1.data.shape)
    height_n23 = height_n23/b1
    height_n23.data = height_n23.data**(-2./3.)
    height_n23.labels(['power'],[powerseries])
    height_n23.rename('power',r'$p$ / $mW$')
    plot(height_n23,'.-',nosemilog=True)
    height_n23_avg = height_n23.copy()
    height_n23_avg.mean('peak')
    plot(powerseries[r_[0,-1]],height_n23_avg.data[r_[0,-1]],'k-',linewidth=2,alpha=0.3,nosemilog=True)
    ylabel(r"$\propto(y'_m/B_1)^{-2/3}\propto 1/s$")
    lplot('esr_dataset'+thisjobname()+'_hn23.pdf')
    #}}}
    #}}}
#}}}
#{{{ dnp
#{{{ base function for dealing with dnp data
def dnp_for_rho(path,name,powerseries,
        concentration = None,
        chemical_name = None,
        expno = r_[5:32],
        t1expnos = [4,36,37],
        t1powers = None,
        first_figure = None,
        show_t1_raw = False,
        show_plots = True,
        verbose = False,
        integration_width = 150,
        basename = '',
        h5file = 'temperature_paper.h5',
        gensvg = False,
        coarse_power = False,
        t1_offset_corr = 0,
        t1_phnum = [],
        t1_phchannel = [],
        **kwargs):
    r'This is essentially a function that matches the jf_dnp au program\nit will process the resulting data\n(using a standard run for that program as a default) and pull out all the info needed to calculate rho,\nplacing them in a HDF5 file'
    if first_figure == None:
        figurelist = []
    else:
        figurelist = first_figure
    #{{{ pull all the relevant info out of the file name
    if h5file != None:
        if chemical_name == None:
            chemical_name = pull_chemical_from_filename(name)
        if concentration == None:
            concentration = pull_concentration_from_filename(name)
        date_of_file = pull_date_from_filename(name)
    #}}}
    #{{{ initialize HDF5 info
    if h5file != None:
        h5file = openFile(h5file,'r+')
        if coarse_power:
            Emax_table = h5file.root.concentration_series.Emax_coarse
        else:
            Emax_table = h5file.root.concentration_series.Emax
        T1_table = h5file.root.concentration_series.T1
        search_string = '(chemical_name == \'%s\') &'%chemical_name + gensearch('concentration','%0.5g',concentration,1e-6)+' & '+gensearch('date','%0.3f',date_of_file,0.1)
        search_string_lowt1 = search_string + " & (power < -99)"
        search_string_powert1 = search_string + " & (power > -99)"
    #}}}
    #{{{ process the Emax data
    if len(expno)>0:
        #{{{ if desired, mask out some subset of the powers
        if coarse_power:
            desired_powers = r_[-999,powerseries.max()-r_[0:10]] # pull the zero power and 9 integral attenuation settings down from the max power
            #{{{ a BEAUTIFUL way to pull the closest value to what I want -- just shuffle out of old lists and into new lists
            powerlist = list(powerseries) # convert to a list so we can pop
            explist = list(expno)
            nearest = lambda value: abs(array(powerlist)-value).argmin() # return the index nearest to this value
            powerseries = []
            expno = []
            for j in desired_powers:
                expno += [explist.pop(nearest(j))]
                powerseries += [powerlist.pop(nearest(j))]
            powerseries = array(powerseries)
            expno = array(expno)
            print '\n\ngenerated new coarse power list:',powerseries
            expno = expno[argsort(powerseries)]
            powerseries = sort(powerseries)
            print ' $\Rightarrow$ then I sorted it','\n\n'
            #}}}
        #}}}
        #{{{ calculate the powers from the powers
        powers = dbm_to_power(powerseries)
        max_power = max(powers)
        #}}}
        #{{{ find the normalized integrals
        # this should be pretty portable -- maybe stick in fornotebook.py
        integral,figurelist = integrate(path+name+'/',expno,integration_width=integration_width,return_noise=True,first_figure = figurelist,**kwargs)
        #print '\n\nDEBUG test figurelist',figurelist,'\n\n'
        figurelist = lplotfigures(figurelist,name+basename+'.pdf')
        #print '\n\nDEBUG 1 integral data ',zip(powerseries,integral.data),'\n\n'
        normalizer = integral['power',0:1].copy()
        integral /= normalizer
        #print '\n\nDEBUG 2 integral data ',zip(powerseries,integral.data),'\n\n'
        #}}}
        if (h5file != None):
            #{{{ find the power scans
            lowt1 = T1_table.readWhere(search_string_lowt1)
            hight1 = T1_table.readWhere(search_string_powert1)
            if (size(hight1)>0):
                hight1 = hight1[argmax(abs(hight1['power']))]
                testf = open('test_output.txt','a')
                T1fit_table = h5file.root.concentration_series.T1fit
                fit_data = T1fit_table.readWhere("(chemical_name == 'water')")
                if fit_data['c_len']>0:
                    fit_coeff = fit_data['c'][0:fit_data['c_len']].flatten()
                    testf.write("found low power %g and high power %g for %g and T10 at max %g\n"%(lowt1['power'],hight1['power'],concentration,fit_coeff[0] + max_power * fit_coeff[1]))
            #testf.write("found low power %g and high power %g for %g and T10 fit_coeff %s\n"%(lowt1['power'],hight1['power'],concentration,repr(fit_coeff).replace('\t',' ')))
                testf.close()
            #}}}
        #{{{ convert integral to Emax fit class, and plot the fit
        figurelist = nextfigure(figurelist,'emax')
        integral = emax(integral,fit_axis = 'power')
        integral.labels(['power'],[powers])
        integral.makereal()
        change = r_[1,diff(integral.getaxis('power'))]
        changemask = change > 0
        plot(integral['power',changemask],'ko')
        plot(integral['power',~changemask],'ro')
        integral.fit()
        power_axis = integral.getaxis('power').copy()
        power_axis = sort(power_axis)
        power_axis_forplot = linspace(power_axis[0],power_axis[-1],100)
        print '\n\n$c_0$ is',integral.output('c_0'),'\n\n'
        plot(integral.eval(power_axis_forplot))
        gridandtick(gca())
        Emax = integral.output('c_0') - integral.output('A')/integral.output('B')
        ax = gca()
        text(0.7,0.5,r'$E_{max}=$ %0.02f'%Emax,transform = ax.transAxes,size = 'x-large', horizontalalignment = 'center',color = 'b')
        #lplot('Evp'+name+basename+'.pdf',gensvg = gensvg)
        #}}}
        #{{{ show the linear plot
        figurelist = nextfigure(figurelist,'linearemax')
        #{{{ find the plot limits
        lims = integral.linear().data.copy()
        lims = lims[isfinite(lims)]
        lims = r_[lims.min(),lims.max()]
        #}}}
        ax = gca()
        plot(integral.linear()['1 / power',changemask],'k.')
        plot(integral.linear()['1 / power',~changemask],'r.')
        power_axis_forplot = linspace(power_axis[1],power_axis[-1],100)
        Evipplot = integral.linear(power_axis_forplot,set='c_0',set_to=1.0)
        plot(Evipplot,'b-')
        ax.set_ylim(lims)
        #lplot('Evip'+name+basename+'.pdf',grid=False,gensvg = gensvg)
        #}}}
        #{{{ do the pinv fit
        print '\n\n'
        coeff = integral.pinv(verbose = verbose)
        print '\n\n'
        if verbose:
            print 'coeff = ',coeff
        print 'coeff[1]',coeff[1]
        print r'(from linear) \Emax = ',1.0-1.0/coeff[1]
        #}}}
        if h5file != None:
            #{{{ save the Emax data
            foundrow = False
            for row in Emax_table.where(search_string):
                foundrow = True
                row['chemical_name'] = chemical_name
                row['concentration'] = concentration
                row['value'] = Emax
                row['date'] = date_of_file
                row.update()
                print "found row with ",concentration,"M ",chemical_name.replace('_',r'\_'),"date",datetime.fromtimestamp(date_of_file).strftime('%m/%d/%y'),"and overwrote with Emax=",Emax
            if not foundrow:
                row = Emax_table.row
                row['chemical_name'] = chemical_name
                row['concentration'] = concentration
                row['value'] = Emax
                row['date'] = date_of_file
                row.append()
                print "found no row with ",concentration,"M",chemical_name.replace('_',r'\_'),"date",datetime.fromtimestamp(date_of_file).strftime('%m/%d/%y'),"so appended Emax=",Emax
            #}}}
        print '\n\ntest figurelist',figurelist,'\n\n'
        figurelist = lplotfigures(figurelist,name+basename+'.pdf')
    #}}}
    if not coarse_power: # since this would just waste space, since I'm going to do it both ways anyways
        #{{{ now, go ahead and process the T_1 data
        t1expnos = map(str,t1expnos)
        if len(t1_phchannel)>0:
            kwargs['phchannel'] = t1_phchannel
            kwargs['phnum'] = t1_phnum
        for j,t1expno in enumerate(t1expnos):
            t1file = [path+name+'/'+t1expno]
            fit,figurelist = process_t1(t1file,[],
               show_image = True,
               integration_width = integration_width,
               usebaseline = False,
               center_peak = True,
               return_noise = True,
               show_integral = True,
               plotcheckbaseline = False,
               abs_image = False,
               first_figure = figurelist,
               offset_corr = t1_offset_corr,
               **kwargs)
            print '\n\n'
            #{{{ show the raw data
            #}}}
            #{{{ show the T1 plots and integral peaks
            if show_plots:
                t1lplotkwargs = {'numwide':4,'width':2}
                figurelist = lplotfigures(figurelist,name+basename+'_'+t1expno+'.pdf',**t1lplotkwargs)
            #}}}
            if t1powers == None:
                #{{{ search for the power in the title string
                titlestr = load_title(t1file[0])
                a = re.compile('(?:^| |\$)([0-9\.]+).*dB')
                m = a.search(titlestr)
                if m:
                    g = m.groups()
                    power = 10.0-double(g[0])
                else:
                    print "\n\nfound no power in the title titlestr"
                    power = -999.
                #}}}
            else:
                power = t1powers[j]
            if h5file != None:
                #{{{ add the data into the HDF5 file
                search_string_w_power = search_string + ' & '+gensearch('power',r'%0.2f',power,0.01)# because we have multiple T1's per power
                foundrow = False
                for row in T1_table.where(search_string_w_power):
                    foundrow = True
                    row['chemical_name'] = chemical_name
                    row['concentration'] = concentration
                    row['value'] = fit.output(r'T_1')
                    row[r'Minf'] = fit.output(r'M(\infty)')
                    row['power'] = power
                    row['date'] = date_of_file # also same
                    row.update()
                    print "\n\nfound row with ",concentration,'M',chemical_name.replace('_',r'\_'),"date",datetime.fromtimestamp(date_of_file).strftime('%m/%d/%y'),'power',power,"and overwrote with T1=",fit.output(r'T_1'),"Minf=",fit.output(r'M(\infty)')
                if not foundrow:
                    row = T1_table.row
                    row['chemical_name'] = chemical_name
                    row['concentration'] = concentration
                    row['value'] = fit.output(r'T_1')
                    row[r'Minf'] = fit.output(r'M(\infty)')
                    row['power'] = power
                    row['date'] = date_of_file # also same
                    row.append()
                    print "\n\nfound no row with ",concentration,'M',chemical_name.replace('_',r'\_'),"date",datetime.fromtimestamp(date_of_file).strftime('%m/%d/%y'),'power',power,"so appended T1=",fit.output(r'T_1'),"Minf=",fit.output(r'M(\infty)')
                    #print r'search string:\begin{verbatim}',search_string_w_power,'\end{verbatim}'
                print '\n\n'
                #}}}
        #}}}
    if h5file != None:
        h5file.close()
#}}}
#{{{ retrieve the t1 power series results from an HDF file
def plot_t1series(path,name,t1expnos,dbm,h5file = 'temperature_paper.h5',gensvg = False):
    dnp_for_rho(path,name,[],expno = [],t1expnos = t1expnos,t1powers = dbm,show_plots = False)
    h5file = openFile(h5file,'r+')
    T1_table = h5file.root.concentration_series.T1
    print r'\begin{verbatim}'
    search_string = '(chemical_name == "%s")'%pull_chemical_from_filename(name)
    search_string += ' & '+gensearch('concentration','%0.8g',pull_concentration_from_filename(name),1e-6)
    search_string += ' & '+gensearch('date','%0.3f',pull_date_from_filename(name),0.1)
    results = T1_table.readWhere(search_string)
    print results
    print r'\end{verbatim}'
    figure(100)
    clf()
    t1data = nddata(results['value'],[len(results),1],['power',r'$T_1$'])
    #print 'powers are ',results['power']
    t1data.labels(['power'],[dbm_to_power(results['power'])])
    t1data.sort('power')
    c,straightline = t1data.polyfit('power')
    print r'found coeff \begin{verbatim}',c,r'\end{verbatim}'
    plot(straightline,'-')
    plot(t1data,'o')
    lplot('T1vp'+name+'_%d_%d.pdf'%(t1expnos[0],t1expnos[-1]),gensvg = gensvg)
    figure(101)
    clf()
    m0data = nddata(results['Minf'],[len(results),1],['power',r'$M(\infty)$'])
    m0data.labels(['power'],[dbm_to_power(results['power'])])
    m0data.sort('power')
    m0data /= m0data['power',0]
    plot(m0data,'o')
    lplot('T1vpMinf'+name+'_%d_%d.pdf'%(t1expnos[0],t1expnos[-1]))
    c_forstore = zeros(10)
    c_forstore[0:size(c)] = c.flatten()
    foundrow = False
    T1fit_table = h5file.root.concentration_series.T1fit
    for row in T1fit_table.where(search_string):
        foundrow = True
        row['chemical_name'] = pull_chemical_from_filename(name)
        row['concentration'] = pull_concentration_from_filename(name)
        row['date'] = pull_date_from_filename(name)
        row['c_len'] = len(c)
        row['c'] = c_forstore
        row.update()
    if not foundrow:
        row = T1fit_table.row
        row['chemical_name'] = pull_chemical_from_filename(name)
        row['concentration'] = pull_concentration_from_filename(name)
        row['date'] = pull_date_from_filename(name)
        row['c_len'] = len(c)
        row['c'] = c_forstore
        row.append()
    h5file.close()
    return m0data,c,straightline
#}}}
#{{{ standard power parsing
def parse_powers(path,name,power_file,expno = None,
        scanlength_estimation_error = 0.2,
        t_minlength = 10.0, # in seconds
        t_start = 4.9, # in minutes
        t_stop_extra = 0.0, # in minutes
        extra_time = 7.86, # amount extra per scan
        tolerance = 2.0):
    fullpath = dirformat(dirformat(path)+name)+'%d/'%1
    fileinfo = load_acqu(fullpath)
    if expno == None:
        expno = r_[5:5+3+1+int(fileinfo['CNST'][6])] # get the number of experiments stored by the au program
        obs(r'automatically determined that you\'re using expno %d $\Rightarrow$ %d'%(expno[0],expno[-1]))
    try:
        fullpath = dirformat(dirformat(path)+name)+'%d/'%expno[0]
    except:
        raise CustomError('type of expno is',expno.dtype)
    fileinfo = load_acqu(fullpath)
    if det_type(fullpath)[1] == True:
        td1 = load_acqu(fullpath,whichdim='2')['TD']
    else:
        td1 = 1
    scan_length = td1*fileinfo['NS']*(fileinfo['D'][1]+fileinfo['TD']/fileinfo['SW_h']/2.0+fileinfo['D'][11]) + 11.0 + extra_time
    t_stop = t_start*60.0+t_stop_extra*60.0+(len(expno)-1)*scan_length
    obs(r'scan\_length is',scan_length,'and there are %d scans for t$_{stop}$=%0.3f\n\n'%(len(expno)-1,t_stop/60.0))
    scanlength_estimation_error *= scan_length # convert from percentage to time
    templen = scan_length-scanlength_estimation_error
    if templen > t_minlength:
       t_minlength = templen
    t_maxlen = scan_length+scanlength_estimation_error
    dbm = auto_steps(path+power_file,t_start = t_start*60,t_stop = t_stop,t_minlength = t_minlength,t_maxlen = t_maxlen,tolerance = tolerance, threshold = -36)
    figure(1)
    gridandtick(gca())
    lplot('powerlog'+name+'.pdf')
    figure(2)
    gridandtick(gca())
    lplot('powerlog_raw'+name+'.pdf')
    return expno,dbm
#}}}
#}}}
