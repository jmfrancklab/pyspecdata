from pyspecdata import *
from string import rstrip
from scipy.io import savemat,loadmat
from os.path import exists as path_exists
from scipy.optimize import leastsq
from scipy.signal import fftconvolve
from pyspecdata.nmr import *
from pyspecdata.load_files import *
# this file should be like fornotebook, except that thisjobname and lplot are stripped of the appropriate lines, so that rather than dumping a file and latex code, it just plots stuff

def thisjobname():
	return ''
def lplot(fname,width=3,dpi=200,grid=True,alsosave=None):
	'''
	used with python.sty instead of savefig
	'''
	if grid:
		gridandtick(gca())
	figure()
	if alsosave != None:
		savefig(alsosave,dpi=dpi)
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
	print(r'\o{')
	for j in arg:
		print(j, end=' ')
	print(r'}')
def save_data(inputdata=[]):
	# concatenate data to data already in data.mat file
	if path_exists('data.mat'):
		data = loadmat('data.mat',struct_as_record=True)
	else:
		data = {}
	if not(data==[]):
		data.update(inputdata)
		savemat('data.mat',data)
	return data
def save_local(inputdata=[]):
	# concatenate data to data already in data.mat file
	if path_exists('local.mat'):
		data = loadmat('local.mat',struct_as_record=True)
	else:
		data = {}
	if not(data==[]):
		data.update(inputdata)
		savemat('local.mat',data)
	return data
def clear_local(inputdata=[]):
	obs(r'{\it clear local}'+'\n\n')
	if path_exists('local.mat'):
		os.unlink('local.mat')
def save_variable(variable,content,disp=True):
	if path_exists('data.mat'):
		data = loadmat('data.mat',struct_as_record=True)
	else:
		data = {}
	data.update({variable:content})
	savemat('data.mat',data)
	if disp:
		obs(variable.replace('_',r'\_'),'=',dp(content,5))
	return data
#{{{ Emax fitting
def eplots(filename,setting,dbm,expno,integral,invalid=[],plot_slope=False,normalize = False,order=True,forceone=False,inverse=True):
	print('\n\n')
	if plot_slope:
		Emax,E,relative_slope = emax_plot(setting,dbm,expno,integral,invalid=invalid,return_slope=True,normalize=normalize,forceone=forceone)
	else:
		Emax = emax_plot(setting,dbm,expno,integral,invalid=invalid,normalize=normalize,order=order,forceone=forceone)
	lplot(filename+r'.pdf')
	emax_T(setting,dbm,expno,integral,invalid=invalid,order=order,showemax=Emax,inverse=inverse)
	title(filename)
	lplot(filename+r'_t.pdf')
	if plot_slope:
		print('\n\n')
		order = argsort(E)[::-1]
		plot(E[order],relative_slope[order])
		lplot(filename+r'_slope.pdf')
		xlabel('E')
		ylabel('$\frac{dE}{dp}/\left.\frac{dE}{dp}\right|_{p=0}$')
def emax_fit(power,integral,fit_axis_size = 100,emax=[],return_slope=False,normalize=False,order=False,forceone=False):
	zeropower = 1e-3*10.0**((-99+40.0)/10.0) #20 db for each atten
	third_integral_with_power= where(power>zeropower)[0][2]
	initial_slope = (integral[third_integral_with_power]-integral[0])/(power[third_integral_with_power]-power[0])
	if emax==[]:
		approx_emax = integral[-1]
	else:
		approx_emax = emax
	p0_ini = [1,-initial_slope,-initial_slope/(1-approx_emax)]
	#p0_ini = [1,1,1]
	if forceone:
		p_out,success = leastsq(emax_forceone_errfunc, p0_ini[1:], args = (
			power # the x axis
			,
			integral # the data we're fitting to
			),
			ftol = 1e-4, xtol = 1e-6)
		p_out = r_[1,p_out]
	else:
		p_out,success = leastsq(emax_errfunc, p0_ini, args = (
			power # the x axis
			,
			integral # the data we're fitting to
			),
			ftol = 1e-4, xtol = 1e-6)
	Emax = p_out[0]-p_out[1]/p_out[2]
	#print 'p\_out is ',p_out
	x = linspace(power.min(),power.max(),100)
	if normalize:
		normalization = p_out[0] # normalize against the "1" of the fit
	else:
		normalization = 1.0 # don't normalize
	if return_slope:
		# just for clarity, to match to the notebook
		A = p_out[0]
		B = p_out[1]
		C = p_out[2]
		relative_slope = (B*(1.+C*x)-C*B*x)/((1.+C*x)**2)/B
		print('A is ',A,'\n\n')
		return (x,emax_fit_fitfunc(p_out,x),Emax,relative_slope,normalization)
	else:
		return (x,emax_fit_fitfunc(p_out,x),Emax,normalization)
def emax_fit_fitfunc(p,x):
	# initial slope is p[1]
	# emax is p[0]-p[1]/p[2]
	return (p[0]-(p[1]*x/(1.+p[2]*x)))
def emax_forceone_errfunc(p,x,y):
	fit = emax_fit_fitfunc(r_[1.,p],x)
	return fit-y
def emax_errfunc(p,x,y):
	fit = emax_fit_fitfunc(p,x)
	return fit-y
def emax_plot(setting,dbm,expno,integral,invalid=[],emax=[],return_slope=False,normalize=False,order=True,forceone=False):
	setting = setting[0:len(expno)]
	power = 1e-3*10.0**((dbm+40.0)/10.0) #20 db for each atten
	labelstring = []
	for j in range(0,len(expno)):
		labelstring += [r'#%d,set %d %0.2f$dBm$'%(expno[j],setting[j],dbm[j])]
	inipower = power.copy()
	iniintegral = integral.copy()
	thisorder = argsort(power)
	power = power[thisorder]
	setting = setting[thisorder]
	integral = integral[thisorder]
	expno = expno[thisorder]
	dbm = dbm[thisorder]
	#{{{ switch the mask to one for all "invalid" entries	
	invalid_mask = zeros(len(expno),dtype='bool')
	for j in range(0,len(invalid)):
		invalid_mask |= (expno == invalid[j])
	#}}}
	if return_slope:
		fitpower,fitintegral,Emax,relative_slope,normalization = emax_fit(
				power[logical_not(invalid_mask)],
				integral[logical_not(invalid_mask)],emax=emax,return_slope=True,normalize=normalize,forceone=forceone)
	else:
		fitpower,fitintegral,Emax,normalization = emax_fit(
				power[logical_not(invalid_mask)],
				integral[logical_not(invalid_mask)],emax=emax,normalize=normalize,forceone=forceone)
	if order:
		ordplot(inipower,iniintegral/normalization,labelstring,'%s')
	else:
		plot(inipower,iniintegral/normalization)
		addlabels('%s',inipower,iniintegral/normalization,labelstring)
	if normalize:
		title(r'$E_{max} = %0.2f$ (normalized by %0.2f)'%(Emax/normalization,normalization))
	else:
		if forceone:
			title(r'$E_{max} = %0.2f$ (one forced)'%Emax)
		else:
			title(r'$E_{max} = %0.2f$'%Emax)
	plot(fitpower,fitintegral/normalization,'r',alpha=0.5,linewidth=2)
	plot(power[invalid_mask],integral[invalid_mask]/normalization,'rx') # put an x over all the "invalid" entries we've excluded from the fit
	xlabel('power / $W$')
	ylabel('integral')
	if return_slope:
		return Emax,fitintegral,relative_slope
	else:
		return Emax
def emax_T(setting,dbm,expno,integral,invalid=[],order=False,showemax = None,inverse=True):
	setting = setting[0:len(expno)]
	zeropower = 1e-3*10.0**((-99+40.0)/10.0) #20 db for each atten
	nonzeropower = dbm>-99
	power = 1e-3*10.0**((dbm+40.0)/10.0) # 20 db for each atten
	labelstring = []
	for j in range(0,len(expno)):
		labelstring += [r'#%d,set %d %0.2f$dBm$'%(expno[j],setting[j],dbm[j])]
	#wa = where(integral!=1.)[0]
	wa = where(power>zeropower)[0]
	if inverse:
		x = 1./power[wa]
		y = 1./(1.-integral[wa])
	else:
		x = power[wa]
		y = 1.-integral[wa]
	if order:
		#ordplot(1./power[wa],1./(1.-integral[wa]),map((lambda x: labelstring[x]),wa),'%s')
		ordplot(x,y,list(map((lambda x: labelstring[x]),wa)),'%s')
	else:
		plot(x,y,'x-')
		addlabels('%s',x,y,list(map((lambda z: labelstring[z]),wa)))
	invalid_mask = zeros(len(expno),dtype='bool')
	for j in range(0,len(invalid)):
		invalid_mask |= (expno == invalid[j])
	if inverse:
		x = 1./power[invalid_mask]
		y = 1./(1.-integral[invalid_mask])
	else:
		x = power[invalid_mask]
		y = 1.-integral[invalid_mask]
	plot(x,y,'rx') # put an x over all the "invalid" entries we've excluded from the fit
	##{{{ add extra lines to check for temperature effects
	if (showemax!=None) and inverse:
		minpower_idx = argmax(1./power*nonzeropower) # zero if false
		maxpower_idx = argmin(1./power*(nonzeropower*max(power)))
		minpower_used_idx = argmax(1./power*logical_and(nonzeropower,logical_not(invalid_mask))) # zero if either is false
		maxpower_used_idx = argmin(1./power*(max(power)*logical_and(nonzeropower,logical_not(invalid_mask))))
		#{{{ make line to test for fit to Emax
		plot(r_[1./power[minpower_used_idx],0],
				r_[1./(1.-integral[minpower_used_idx]),1./(1.-showemax)],
				'r',linewidth=2,alpha=0.5)
		#}}}
		#{{{ make line to show direction of 1
		#}}}
	#}}}
	#plot((1./power[wa]).flatten(),(1./(1.-integral[wa])).flatten())
	if inverse:
		xlabel('1/ power / $W^{-1}$')
		ylabel('1/1-integral')
	else:
		xlabel('power / $W$')
		ylabel('1-integral')
	#redo
#}}}
#{{{ specific convenience functions
def calcdistance(freq,optdistance):
	print(r'$%0.1f\;GHz$ @ '%(freq/1e9),dp((optdistance+119e-3)*1e3,1),r'$mm$')
def calcfield(elfreq,elratio = 28.113,nmrelratio = 1.5167):
	obs(elfreq,'$GHz$,',dp(elratio,5),'$GHz/T\;',dp(nmrelratio,4),'\;ppt\;\Rightarrow$',dp(elfreq/elratio*1e4,2),'$G$',dp(elfreq*nmrelratio,5),'$MHz$')
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
def calcfielddata(freq,substance):
	data = save_data()
	calcfield(freq,elratio=data[substance+'_elratio'],nmrelratio=data[substance+'_nmrelratio'])
	save_data({'current_frequency':freq})
	save_data({'current_ppt':data[substance+'_nmrelratio']})
def cpmgseries(filename,plotlabel,tau=None,alpha=None,alphaselect=None):
	data = prospa.load_datafile(filename,dims=2)
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
	filename = '/mnt/esr/Magritek/john/'+exp+'/%d/'%number
	print(r'\fn{%s %d}'%(exp,number)+'\n\n')
	cpmgseries(filename,exp+thisjobname(),tau,alpha,alphaselect)
#{{{ esr_saturation
def esr_saturation(file,powerseries,smoothing=0.2,threshold=0.8):
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
		print(r'\begin{verbatim} If you have an error here, probably change smoothing (%0.2f) or threshold (%0.2f)\end{verbatim}'%(smoothing,threshold),'\n\n')
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
		print(r'\begin{verbatim}lengths: ',list(map(len,allpeaks_top_x)),'\end{verbatim}')
		return
	try:
		allpeaks_bottom_x = nddata(allpeaks_bottom_x,[nslices,num_peaks],['power','peak']).reorder(['power','peak'])
	except:
		print(r'\begin{verbatim} If you have an error here, probably change smoothing (%0.2f) or threshold (%0.2f)\end{verbatim}'%(smoothing,threshold),'\n\n')
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
		print(r'\begin{verbatim}lengths: ',list(map(len,allpeaks_top_x)),'\end{verbatim}')
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
	print('\n\n')
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
	print('\n\n')
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
#}}}
