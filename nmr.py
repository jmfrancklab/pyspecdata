# just make this a library of all NMR reading software
from matlablike import *
import re
import string
import struct
import os
import fornotebook
from scipy.optimize import leastsq
from scipy.io import loadmat

def OUTPUT_notebook():
	return True

def auto_steps(filename,threshold = -35, t_minlength = 0.5*60,minstdev = 0.1,showplots = True, showdebug = False,t_start=0):
	v = loadmat(filename)
	p_ini = v['powerlist']
	t_ini = v['timelist']
	a = whereblocks(logical_and(t_ini>t_start,logical_or(p_ini>threshold,~(isfinite(p_ini)))))
	a = a[argmax(map(len,a))]
	minsteps = sum(t_ini<t_minlength)
	t = t_ini[a]
	p = p_ini[a]
	flattened = NaN*zeros(len(p))
	plotdev = zeros(len(p))
	plotstd = zeros(len(p))
	powerlist = []
	nextpos = 0
	while nextpos < len(t)-1:
		# grab the stdev and average for the minimum number of steps
		blockstart = nextpos
		subset= p[nextpos:minsteps+nextpos]
		curmean = mean(subset[isfinite(subset)])
		stdev = std(subset[isfinite(subset)])
		if stdev < minstdev:
			stdev = minstdev
		nextpos += minsteps
		# now that we have a decent stdev, just test to see that every step is less than the stdev
		while (nextpos < len(t)-1) and (abs(p[nextpos]-curmean)<2*stdev or ~isfinite(p[nextpos])):
			subset= p[blockstart:nextpos]
			curmean = mean(subset[isfinite(subset)])
			stdev = std(subset[isfinite(subset)])
			if stdev < minstdev:
				stdev = minstdev
			if showdebug:
				plotdev[nextpos] = p[nextpos]-curmean
				plotstd[nextpos] = 2*stdev 
			nextpos += 1
		#uptohere = flattened[0:nextpos]
		#uptohere[uptohere==0] = cursum/curnum
		#flattened[nextpos-1] = cursum/curnum
		subset = flattened[0:nextpos]
		subset[isnan(subset)] = curmean
		powerlist += [curmean]
		#flattened[nextpos-1] = curmean
	if showplots:
		plot(t_ini/60,p_ini.flatten(),'b')
		plot(t/60,flattened,'g',linewidth=5,alpha=0.5)
		title(filename)
	if showdebug:
		if showplots:
			twinx()
		plot(t/60,plotdev,'k')
		plot(t/60,plotstd,'r')
		title('Power meter log')
	return array(powerlist)
def error_plot(*arg):
	width=6
	dpi=200
	fname = 'error_plot.png'
	fname = 'auto_figures/'+fname
	if grid:
		gridandtick(gca())
	savefig(fname,dpi=dpi)
	if figure:
		print r"""
		\begin{figure}[h]
		\end{figure}
		"""
	print r'\includegraphics[width=%0.2fin]{%s}'%(width,fname)
	clf()
	raise CustomError(*arg)
#{{{ lower level functions
#{{{ wrappers/generic functions to load acq and data files
b0 = r'$B_0$'
def show_acqu(vars):
	print '\\begin{verbatim}',vars.__repr__().replace(',','\n'),'\\end{verbatim}\n\n'
#{{{ is this (bruker or prospa, 1d or nd)?
def det_type(filename):
	filetype = None
	if os.path.exists(filename+'.spc'):
		return ('winepr',True)
	else:
		filename = dirformat(filename)
		files_in_dir = os.listdir(filename)
		if os.path.exists(filename+'ser'):
			return ('bruker',True)
		elif os.path.exists(filename+'data.2d'):
			return ('prospa',True)
		elif any(map((lambda x:'Delay' in x),files_in_dir)):
			return ('prospa','t1')
		elif os.path.exists(filename+'acqu.par'):
			return ('prospa',False)
		elif os.path.exists(filename+'acqus'):
			return ('bruker',False)
		elif os.path.exists(filename+'../acqu.par'):
			return ('prospa','t1_sub')
		else:
			raise CustomError('WARNING! unidentified file type '+filename)
#}}}
#{{{ load an nddata structure for a 2d set -- give the data needed to load
def load_file(filenames,dimname='',calibration=1.0):
	if type(filenames) is not list:
		filenames = [filenames]
	filenames = map(dirformat,filenames)
	#{{{load all the data into a list
	data = [load_indiv_file(filenames[0],dimname=dimname)]
	for filename in filenames[1:]:
		data += [load_indiv_file(filename,dimname=dimname)]
	#}}}
	newdata = concat(data,dimname) # allocate the size of the indirect array
	newdata_shape = ndshape(newdata)
	if all(map((lambda x:det_type(x)[0]=='prospa'),filenames)):
		newdata = prospa_decim_correct(newdata,dimname)
	if newdata_shape[dimname]==1:
		newdata.popdim(dimname)
	return newdata*calibration
def bruker_det_rg(a):
	'''determine the actual voltage correction from the value of rg for a bruker NMR file'''
	return a
def load_indiv_file(filename,dimname='',return_acq=False):
	filetype,twod = det_type(filename)
	if filetype == 'winepr':
			fp = open(filename+'.spc','rb')
			data = fp.read()
			data = array(
					struct.unpack('<%df'%(len(data)/4),data),
					dtype='double')
			v = winepr_load_acqu(filename)
			xpoints = v['RES']
			rg = v['RRG']
			data /= rg
			ypoints = len(data)/xpoints
			if ypoints>1:
				if ypoints != v['REY']:
					raise CustomError('I thought REY was the indirect dim, guess not')
				if dimname=='':
					dimname = v['JEY']
				data = nddata(data,[ypoints,xpoints],[dimname,b0])
			else:
				data = nddata(data,[xpoints],[b0])
			xlabels = linspace(v['HCF']-v['HSW']/2.,v['HCF']+v['HSW']/2.,xpoints)
			if len(data.dimlabels)>1:
				data.labels([dimname,b0],[linspace(0,1,ypoints),xlabels])
				data.reorder([b0,dimname])
			else:
				data.labels([b0],[xlabels])
			return data
	filename = dirformat(filename)
	if twod:
		if filetype == 'bruker':
			v = bruker_load_acqu(filename)
			v2 = bruker_load_acqu(filename,whichdim='2')
			td2 = int(v['TD'])
			rg = bruker_det_rg(float(v['RG']))
			td1 = int(v2['TD'])
			td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
			#print 'loading bruker:',td2,'x',td1
			#print 'file should be:',td2_zf,'x',td1
			fp = open(filename+'ser','rb')
			data = fp.read()
			#print 'data size',len(data),'should be',4*td2_zf*td1
			data = array(
					struct.unpack('>%di'%(len(data)/4),data),
					dtype='complex128')
			data = data[0::2]+1j*data[1::2]
			data /= rg
			data = nddata(data,[td1,td2_zf/2],[dimname,'t2'])
			data = data['t2',0:td2/2] # now, chop out their zero filling
			t2axis = 1./v['SW_h']*r_[1:td2/2+1]
			t1axis = r_[0:td1]
			data.labels([dimname,'t2'],[t1axis,t2axis])
			shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the second order phase shift
			#print 'shiftpoints = ',shiftpoints
			data.circshift('t2',shiftpoints)
			# finally, I will probably need to add in the first order phase shift for the decimation --> just translate this
			if return_acq:
				return (data,v,v2)
			else:
				return data
		if (filetype == 'prospa') and (twod == 't1_sub'):
			v = prospa_load_acqu(filename+'../')
			data = prospa_load_datafile(filename,dims=1)
			data = nddata(data,[1,v['nrPnts']],[dimname,'t2'])
			taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
			data.labels([dimname,'t2'],[r_[1],taxis])
			return data
		elif filetype == 'prospa':
			v = prospa_load_acqu(filename)
			if v['experiment'][0:4]=='cpmg':
				data = prospa_load_datafile(filename,dims=2)
				data = nddata(data,[v['nrEchoes'],v['nrPnts']],['echo','t2'])
				taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
				echotime = (r_[0:v['nrEchoes']]+0.5)*v['echoTime']/1e6
				data.labels(['echo','t2'],[echotime,taxis])
				return data
			else:
				print "ERROR, load_file, doesn't understand the type of pulse sequenced used, so I can't properly sort the 2d data"
	else:
		if filetype == 'bruker':
			v = bruker_load_acqu(filename)
			td2 = int(v['TD'])
			td1 = 1
			td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
			fp = open(filename+'fid','rb')
			data = fp.read()
			data = array(
					struct.unpack('>%di'%(len(data)/4),data),
					dtype='complex128')
			data = data[0::2]+1j*data[1::2]
			rg = bruker_det_rg(v['RG'])
			data /= rg
			data = nddata(data,[td1,td2_zf/2],[dimname,'t2'])
			data = data['t2',0:td2/2] # now, chop out their zero filling
			t2axis = 1./v['SW_h']*r_[1:td2/2+1]
			t1axis = r_[1]
			data.labels([dimname,'t2'],[t1axis,t2axis])
			shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the second order phase shift
			#print 'shiftpoints = ',shiftpoints
			data.circshift('t2',shiftpoints)
			# finally, I will probably need to add in the first order phase shift for the decimation --> just translate this
			if return_acq:
				return (data,v,v2)
			else:
				return data
		elif filetype == 'prospa':
			v = prospa_load_acqu(filename)
			data = prospa_load_datafile(filename,dims=1)
			data = nddata(data,[v['nrPnts']],['t2'])
			taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
			data.labels(['t2'],[taxis])
			return data
		else:
			raise CustomError("can't load this file type")
#}}}
#{{{ t1 axis
def bruker_load_t1_axis(files):
	wait_time = []
	if type(files) is str:
		files = [files]
	#print 'DEBUG: trying to load files: ',files
	for thisfile in files:
		thisfile = dirformat(thisfile)
		thisfiletype = det_type(thisfile)
		if thisfiletype[0] == 'prospa':
			print 'need to copy over code from prospa'
		if thisfiletype[0] == 'bruker':
			wait_time += [bruker_load_vdlist(thisfile)]
		else:
			print 'couldn\'t determine thisfile type'
	return array(wait_time).flatten()
#}}}
#}}}
#{{{ routines specific to prospa
#{{{ load acqu
def prospa_decim_correct(data,indirect_dim_name):
	#{{{ get rid of the finite rise time	
	data_abs = abs(data)
	data_abs.mean(indirect_dim_name)
	data_abs = data_abs.run(argmax,'t2')
	top = int(data_abs.data)
	data.circshift('t2',top)
	#}}}
	print 'Applied prospa decimation correction'
	return data
def prospa_load_acqu(file):
	file = dirformat(file)
	fp = open(file+'acqu.par')
	lines = fp.readlines()
	line_re = re.compile(r'([^ \t]+) *= *(.+)')
	vars = {}
	for j in range(0,len(lines)):
		lines[j] = string.rstrip(lines[j])
		m = line_re.match(lines[j])
		if m:
			exec 'temp = %s'%m.groups()[1]
			vars.update({m.groups()[0]:temp})
		else:
			print "error, acqu.par line not parsed: ",lines[j]
	fp.close()
	return vars
#}}}
#{{{ load_datafile
def prospa_load_datafile(file,dims=1):
	file = dirformat(file)
	if dims == 1:
		fp = open(file+'data.1d','rb')
	elif dims == 2:
		fp = open(file+'data.2d','rb')
	else:
		print 'ERROR: wrong number of dims'
	data = fp.read()
	data = array(struct.unpack('%df'%(len(data)/4),data))
	data = data[7:]
	# the following is junk!!!
	#elif precision=='b':
	#	data = array(struct.unpack('%db'%(len(data)/1),data))
	#	data = data[7*4:]
	#else:
	#	print 'error, precision wrong'
	data = reshape(data,(-1,2))
	data = data[:,0]+1j*data[:,1]
	fp.close()
	return data
#}}}
#{{{ OBSOLETE -- load_prospa1d -- just replace with a higher-level function
def load_prospa1d(file):
	file = dirformat(file)
	data = prospa_load_datafile(file)
	vars = prospa_load_acqu(file)
	return (vars,data)
#}}}
#{{{ load the data from a t1 file based on file names
def prospa_t1_info(file):
	file = dirformat(file)
	files = os.listdir(file)
	file_re = re.compile(r'([0-9]+)msDelay$')
	datafiles = []
	wait_time = []
	for j in range(0,len(files)):
		m = file_re.match(files[j])
		if m:
			datafiles += [file+files[j]]
			wait_time += [int(m.groups()[0])]
	return datafiles,array(wait_time)*1e-3
#}}}
#}}}
#{{{ routines specific to Bruker
#{{{ calculate decimation correction
def bruker_det_phcorr(v):
	if v['DIGMOD']==1:
		gdparray=array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[179,201,533,709,1097,1449,2225,2929,4481,5889,8993,11809,18017,23649,36065,47329,72161,94689,144353,189409,288737],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[184,219,384,602,852,1668,2292,3368,4616,6768,9264,13568,18560,27392,36992,55040,73856,110336,147584,220928,295040]])
		decimarray=array([2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024])
		dspfvs = v['DSPFVS']
		decim = v['DECIM']
		return gdparray[dspfvs,where(decimarray==decim)[0]]/2/decim;
	else:
		return array([0])
#}}}
#{{{ load acqu
def bruker_match_line(line,number_re,string_re,array_re):
	m = number_re.match(line)
	if m:
		retval = (0,m.groups()[0],double(m.groups()[1]))
	else:
		m = string_re.match(line)
		if m:
			retstring = m.groups()[1]
			if retstring[-1]=='>':
				retstring = retstring[:-1]
			retval = (1,m.groups()[0],0,retstring)
		else:
			m = array_re.match(line)
			if m:
				retval = (2,m.groups()[0],(double(m.groups()[1]),double(m.groups()[2])),m.groups()[3])
			else:
				retval = (3,line)
	return retval
#{{{ bruker_load_vdlist
def bruker_load_vdlist(file):
	fp = open(file+'vdlist')
	lines = fp.readlines()
	lines = map(string.rstrip,lines)
	lines = map((lambda x: x.replace('m','e-3')),lines)
	lines = map((lambda x: x.replace('s','')),lines)
	lines = map((lambda x: x.replace('u','e-6')),lines)
	lines = map(double,lines)
	return array(lines)
#}}}
#{{{ bruker_load_acqu
def load_title(file):
	return bruker_load_title(file)
def bruker_load_title(file):
	file = dirformat(file)
	fp = open(file+'pdata/1/title')
	lines = fp.readlines()
	emptystring = '\r\n'
	while emptystring in lines:
		lines.pop(lines.index(emptystring))
	emptystring = '\n'
	while emptystring in lines:
		lines.pop(lines.index(emptystring))
	return ''.join(lines)
def bruker_load_acqu(file,whichdim=''):
	fp = open(file+'acqu'+whichdim+'s')
	lines = fp.readlines()
	vars = {}
	number_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *([0-9\-\.]+)')
	string_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *<(.*)')
	array_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *\(([0-9]+)\.\.([0-9]+)\)(.*)')
	lines = map(string.rstrip,lines)
	j=0
	retval =  bruker_match_line(lines[j],number_re,string_re,array_re)
	j = j+1
	retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re) #always grab the second line
	while j < len(lines):
		isdata = False
		if retval[0]==1 or retval[0]==2:
			name = retval[1]
			thislen = retval[2]
			data = retval[3]
			while (retval2[0] == 3) and (j<len(lines)): # eat up the following lines
				data += ' '+retval2[1]
				j = j+1
				retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re)
			isdata = True
		elif retval[0]==0:
			name = retval[1]
			data = retval[2]
			isdata = True
		#else:
		#	print 'not a data line:',retval[1]
		if(isdata):
			if retval[0]==2: #if it's an array
				data = data.split(' ')
				if len(data)>0:
					while '' in data:
						data.remove('')
					data = map(double,data)
					if len(data)-1!= thislen[1]:
						print 'error:',len(data)-1,'!=',thislen[1]
			vars.update({name:data})
		# at this point, the string or array data is loaded into data and we have something in retval2 which is definitely a new line
		retval = retval2
		j = j+1
		if j<len(lines):
			retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re)
	fp.close()
	return vars
#}}}
#{{{ winepr_load_acqu
def winepr_load_acqu(file):
	fp = open(file+'.par','rU') # the U automatically converts dos format
	lines = fp.readlines()
	vars = {}
	line_re = re.compile(r'([_A-Za-z0-9]+) +(.*)')
	lines = map(string.rstrip,lines)
	v = {'DRS':4096,'RES':1024,'HSW':50}
	for line in lines:
		m = line_re.match(line)
		name = m.groups()[0]
		value = m.groups()[1]
		try:
			value = int(value)
		except:
			try:
				value = double(value)
			except:
				pass
		v[name]=value
	jss = long(v['JSS'])
	parameters = [ 'DUAL', '2D', 'FT', 'MAN0', 'MAN1', 'PROT', 'VEPR', 'POW', 'ABS', 'FTX', 'FTY', 'POW2', 'ABS2']
	parameters = map((lambda x: 's_'+x),parameters)
	masks = [ 0x00000001L, 0x00000002L, 0x00000004L, 0x00000008L, 0x00000010L, 0x00000020L, 0x00000040L, 0x00000080L, 0x00000100L, 0x00000200L, 0x00000400L, 0x00000800L, 0x00001000L]
	values = map((lambda x: x&jss),masks)
	values = map(bool,values)
	values = map(bool,values)
	v.update(dict(zip(parameters,values)))
	return v
#}}}
#}}}
#}}}
#{{{ higher level functions
#{{{ optimize the first-order phase the right way
def phaseopt(curve):
	curve = curve.copy()
	#{{{ find bestindex once
	phases = linspace(-pi/2,pi/2,100).reshape(1,-1) # this should work w/out altering the sign
	rotated_data = (curve.reshape(-1,1))*exp(-1j*phases)
	success = (real(rotated_data)**2).sum(axis=0)/((imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
	bestindex = where(success==max(success))[0][0]
	#}}}
	#{{{ find the bestindex based on that 
	if (bestindex>phases.shape[1]-2):
		bestindex = phases.shape[1]-2
	phases = linspace(
			phases[0,bestindex-1],
			phases[0,bestindex+1],
			100).reshape(1,-1)
	rotated_data = (curve.reshape(-1,1))*exp(-1j*phases)
	success = (real(rotated_data)**2).sum(axis=0)/((imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
	bestindex = where(success==max(success))[0][0]
	#}}}
	return exp(-1j*phases[0,bestindex])
#}}}
#{{{ fitting functions
#{{{ filter fitting 
def exp_fit(x,y):
	'''fit an exponential function, disregarding the first point'''
	testpoint = int(len(x)/10)
	startpoint = y[1:6].mean() #because of the filter, sometimes, the very first point can start at 0
	endpoint = y[-testpoint-6:-testpoint].mean()
	initial_slope = (y[testpoint]-startpoint)/(x[testpoint]-x[3])
	# initial slope = -startpoint/t1
	p_ini = [startpoint,-startpoint/initial_slope,endpoint]
	p_out,success = leastsq(exp_errfunc, p_ini,
			args = (x[1:],y[1:]),maxfev=3000,ftol = 1e-4, xtol = 1e-6)
	if (success<0) or (success>4):
		p_out = r_[1.,len(x)*4.,0]
		clf()
		plot(x,exp_fitfunc(p_ini,x))
		plot(x,exp_fitfunc(p_out,x))
		plot(x,y)
		legend(['initial guess','after fit','data'])
		error_plot(success,'There was an error attempting to fit')
	else:
		return p_out
def exp_fitfunc(p,x):
	return p[0]*exp(-x/p[1])+p[2]
def exp_errfunc(p,x,y):
	fit = exp_fitfunc(p,x)
	return fit-y
#}}}
#{{{ T1 fitting
# M(t) = M0+(M_ini-M0)exp(-t/T1)
# M'(t) = -(M_ini-M0)exp(-t/T1)/T1 = -(M(t)-M0)/T1
def t1_fit(x,y):
	testpoint = len(y)/3
	initial_slope = (y[testpoint]-y[0])/(x[testpoint]-x[0])
	A = y[-1]
	B = y[testpoint]-x[testpoint]*initial_slope
	C = (A-B)/initial_slope
	if (C < 0):
		raise CustomError(maprep('Negative T1!!! A-B=',A-B,'initial_slope=',initial_slope,x,y))
	#else:
	#	print 'C is ',C
	p_ini = r_[A,B,C]
	#p_ini = [y[-1],y[0],1/initial_slope*(y[-1]-y[0])]
	try:
		p_out,success = leastsq(t1_errfunc, p_ini,
				args = (x,y),ftol = 1e-4, xtol = 1e-6)
	except:
		print 'type(x):', type(x), 'type(y):', type(y)
		raise CustomError('leastsq failed')
	#print 'success:',success
	#print p_ini
	return p_out
def t1_fitfunc(p,x):
	return p[0]+(p[1]-p[0])*exp(-x/p[2])
def t1_errfunc(p,x,y):
	fit = t1_fitfunc(p,x)
	return fit-y
def t1_sat_fitfunc(p,x):
	return p[0]-p[0]*exp(-x/p[1])
def t1_sat_errfunc(p,x,y):
	fit = t1_sat_errfunc(p,x)
	return fit-y
#}}}
#}}}
#{{{ process_emax, which just takes a list of expno and returns the list of integrals
def rg_check(file,expno):
	#print 'diagnose: before emax'
	data = load_emax(file,expno,printinfo = False) # load the data
	data.reorder(['t2','power'])
	plot(data['t2',0:100],'k')
	plot(data['t2',0:100]*(-1j),'b',alpha=0.3)
	title('receiver gain check')
	axis('tight')
	return
def integrate_emax(file,expno,integration_width=1e3,intpoints=None,showimage=False,usephase=True,filteredge = 10,center_peak=False,usebaseline=False,plotcheckbaseline=False,filter_direct = False,return_noise=False,show_integral = False,indiv_phase = False,abs_image = False):
	#print 'diagnose: before emax'
	data = load_emax(file,expno) # load the data
	see_fid = False
	if see_fid:
		#{{{ for debug
		data.reorder(['t2','power'])
		plot(abs(data['t2',0:500]))
		return r_[0,1,2],r_[0,1,2]
		#}}}
	debug_ini_data_size = ndshape(data)
	data.ft('t2',shiftornot=True) # ft along t2
	#{{{baseline correction
	if usebaseline:
		for j in range(0,ndshape(data)['power']):
			if plotcheckbaseline:
				hold(True)
				baseline_data = baseline_spectrum(data['power',j].copy(),showplots = True)
			else:
				baseline_data = baseline_spectrum(data['power',j].copy())
				data['power',j].data[:] -= baseline_data.data.flatten()
				if any(isnan(data['power',j].data[:])):
					print 'isnan!!'
				if any(isinf(data['power',j].data[:])):
					print 'isinf!!'
		if plotcheckbaseline:
			error_plot('Wanted to check baseline')
	#}}}
	debug_data_afterbaseline = data.data.copy()
	if filter_direct:
		data.ift('t2',shiftornot=True) # ft along t2
		filter = matched_filter(data,'t2')
		data *= filter
		data.ft('t2',shiftornot=True) # ft along t2
	#print 't2:',data.getaxis('t2')
	#print 'f2:',data.getaxis('t2')
	#{{{ find the top of the peak	
	data_shape = ndshape(data)
	#{{{ get data_abs, which is for the peak-picking from matched-filtered data
	data_abs = data.copy()
	data_abs.ift('t2',shiftornot=True) # ft along t2
	filter = matched_filter(data_abs,'t2')
	data_abs *= filter
	data_abs.ft('t2',shiftornot=True) # ft along t2
	data_abs = abs(data_abs)
	#}}}
	if not(center_peak):
		data_abs.mean('power')
	data_topvals = data_abs.copy()
	data_topvals.run(amax,'t2')
	data_abs.run(argmax,'t2')
	if not(center_peak):
		data_abs.data = ones(data_shape['power'])*data_abs.data
	top = int32(data_abs.data)
	topvals = data_topvals.data
	#}}}
	#{{{ if integration points unspec, pull out the integration width
	if intpoints==None:
		df = data.getaxis('t2').copy()
		df = df[1]-df[0]
		intpoints = floor(integration_width/(df))
	#}}}
	#{{{ pull out the data + and - intpoints from the peak	
	data_shape['t2'] = intpoints*2+1
	#newdata = data_shape.alloc()
	newdata = []
	newnoise = []
	top[top<intpoints] = intpoints # this is non-ideal but prevents a bug where it gives a bad range
	for j in range(0,data_shape['power']):
		#newdata['power',j,'t2',:] = data['power',j,'t2',top[j]-intpoints:top[j]+intpoints+1]
		newdata += [data['power',j,'t2',top[j]-intpoints:top[j]+intpoints+1]]
		#print 'trying slice',top[j]-intpoints,'to',top[j]+intpoints+1
		if return_noise:
			newnoise += [data['power',j,'t2',10:10+intpoints]]
	debug_top = top
	debug_slice = [top[j]-intpoints,top[j]+intpoints+1]
	debug_newdata_len = len(newdata)
	newdata = concat(newdata,'power')
	if return_noise:
		newnoise = concat(newnoise,'power')
	debug_newdata_shape = ndshape(newdata)
	debug_newdata_data_shape = newdata.data.shape
	#}}}
	#{{{ if we're using absolute value, calculate the noiselevel
	if not usephase:
		noiselevel = data['t2',-(2*intpoints+1)-filteredge:-filteredge]
		noiselevel = abs(noiselevel).mean('t2').mean_nopop('power')
		newdata = abs(newdata)
	#}}}
	#{{{ show what we're integrating
	if showimage:
		clf()
		figure(1)
		try:
			#image(newdata)
			center_for_image = int32(sum(top*topvals)/sum(topvals))
			if abs_image:
				image(abs(data['t2',center_for_image-3*intpoints:center_for_image+3*intpoints+1]))
				title('2D plot of spectra')
			else:
				image(data['t2',center_for_image-3*intpoints:center_for_image+3*intpoints+1])
				title('2D plot of spectra')
		except:
			print 'top = ',top
			print 'topvals = ',topvals
			if(any(isnan(newdata.data))):
				raise CustomError('image isnan')
			elif(any(isinf(newdata.data))):
				raise CustomError('image isinf')
			else:
				clf()
				plot(debug_data_afterbaseline.T)
				error_plot()
				raise CustomError('cant make image from type',type(newdata.data),newdata.data.dtype,'data_shape:',data_shape,'debug_newdata_len:',debug_newdata_len,'debug_newdata_shape:',debug_newdata_shape,'debug_newdata_data_shape:',debug_newdata_data_shape,'debug_top:',debug_top,'debug_slice:',debug_slice,'debug_ini_data_size:',debug_ini_data_size,'debug_data_afterbaseline(',shape(debug_data_afterbaseline),'):',debug_data_afterbaseline)
	newdatacopy = newdata.copy()
	#}}}
	newdata.mean('t2') # integrate
	#{{{ either autophase, or subtract out the abs baseline	
	if usephase:
		if not indiv_phase:
			#newdata.data *= exp(-1j*angle(newdata.data.sum()))
			phaseoptval = phaseopt(newdata.data)
			newdata.data *= phaseoptval
			newdatacopy.data *= phaseoptval
			if return_noise:
				newnoise.data *= phaseoptval
		else:
			for j in range(0,len(newdata.data)):
				phcorr =  newdata['power',j]
				phcorr /= abs(phcorr)
				try:
					newdatacopy['power',j] *= phcorr
					newdata['power',j] *= phcorr
				except:
					print 'shape of newdatacopy',ndshape(newdatacopy)
					print 'shape of newdata',ndshape(newdata)
	else:
		newdata -= noiselevel
		newdatacopy.data -= noiselevel
	if showimage:
		figure(2)
		clf()
		plot((newdatacopy).reorder(['t2','power']))
		title('Peaks, zoomed in to integration region')
		if show_integral:
			#{{{this does work to plot the integral
			newdatacopy.integrate('t2') #newdatacopy.run_nopop(cumsum,'t2')
			gridandtick(gca())
			twinx()
			plot(newdatacopy,alpha=0.5)
			#}}}
		figure(1)
	#}}}
	if return_noise:
		newnoise.data = abs(newnoise.data)
		newnoise.data **= 2
		newnoise.mean('t2') # which should give us the sigma --> note that if there is no baseline correction, this will not be a sigma
		return real(newdata.data.flatten()),sqrt(newnoise.data.flatten()*sqrt(2.*intpoints+1.)) # we have to account for the signal averaging from the integral
	else:
		return real(newdata.data.flatten())
#}}}
#{{{ load the data from a emax series based on input array
def print_info(filename):
	filetype,twod = det_type(filename)
	if filetype == 'bruker':
		v = bruker_load_acqu(dirformat(filename))
		if OUTPUT_notebook():
			print r'\fn{%s}: sfo1:%0.5f aq:%0.3f swh:%0.3f ns:%d ds: %d rg:%0.1f d1:%0.1f p1:%0.2f pl1:%0.1f'%(filename,v['SFO1'],v['TD']/v['SW_h']/2.0,v['SW_h'],v['NS'],v['DS'],v['RG'],v['D'][1],v['P'][1],v['PL'][1])
			data = fornotebook.save_data()
			pptvalue = v['SFO1']/data['current_frequency']
			if abs(pptvalue-data['current_ppt'])>1e-4:
				print ('WARNING! ppt value is off!! (',pptvalue,' ppt)')
			else:
				print '%0.4f $ppt$'%pptvalue
			print '\n\n'
def load_emax(file,datafiles,printinfo=True):
	datafiles = list(datafiles)
	if len(datafiles)==0:
		data = load_file(file,dimname='power') # got rid of a [0] here --> don't know why it was there --> load_file loads lists of files
		if printinfo:
			print_info(file[0])
	else:
		file = dirformat(file)
		datafiles = map(str,datafiles)
		spectype,dim = det_type(file+datafiles[0])
		if spectype == 'bruker':
			#{{{ concatenate all the files into one	
			data = []
			for j in range(0,len(datafiles)):
				data += [load_file(file+datafiles[j],dimname='power')]
			data = concat(data,'power')
			#}}}
			if printinfo:
				print_info(file+datafiles[0])
				print_info(file+datafiles[-1])
		elif spectype == 'prospa':
			#{{{ kea version
			datafiles = map((lambda x: file+str(x)), datafiles)
			vars = prospa_load_acqu(datafiles[0])
			for j in range(0,len(datafiles)):
				vars_check = prospa_load_acqu(datafiles[j])
				if (vars_check['acqTime']!=vars['acqTime']) or (vars_check['nrPnts']!=vars['nrPnts']):
					raise CustomError('Problem! doesn\'t match first file:\n')
					show_acqu(vars)
					show_acqu(vars_check)
					return
			npoints = vars['nrPnts']
			taxis = linspace(0,1,npoints)*vars['acqTime']/1e3
			#npoints = load_datafile(datafiles[0]).shape[0] #hack
			data_shape = ndshape([npoints,len(datafiles)],['t2','power'])
			data = data_shape.alloc()
			data.labels(['power','t2'],[range(0,len(datafiles)),taxis])
			for j in range(0,data_shape['power']):
					data['power',j] = prospa_load_datafile(datafiles[j])
			data = prospa_decim_correct(data,'power')
			#}}}
	return data
#}}}
#{{{ process_t1
def process_t1(file,expno,integration_width = 0.5e3,usebaseline = False,filter_direct = False,return_fit=False,center_peak = True,showimage = True,show_integral = True,plotcheckbaseline = False,abs_image = False,saturation = False):
	if type(file) is str:
		file = [file]
	titlestr = load_title(file[0])
	if det_type(file[0]) == ('prospa','t1'):
		print 'Found a prospa format T1 file'
		file,wait_time = prospa_t1_info(file[0])
	elif det_type(file[0])[0] == 'bruker':
		print 'Found a bruker format T1 file'
		wait_time = bruker_load_t1_axis(file[0])
	else:
		CustomError("You are trying to process a T1 format file, but I'm not finding the right types of files!")
	order = argsort(wait_time)
	integral,noiselevel = integrate_emax(file,expno,
	   showimage = showimage,
	   integration_width = integration_width,
	   usebaseline = usebaseline,
	   filter_direct = filter_direct,
	   center_peak = center_peak,
	   return_noise = True,
	   show_integral = show_integral, plotcheckbaseline = plotcheckbaseline,abs_image = abs_image)
	#return # to DEBUG
	wait_time = wait_time[order]
	integral = integral[order]
	noiselevel = noiselevel[order]
	#{{{ finally, show the fit	
	figure(3)
	taxis = wait_time
	imagpart = imag(integral*phaseopt(integral))
	integral = real(integral*phaseopt(integral))
	plot(taxis,integral,'o')
	plot(taxis,imagpart,'o')
	p = t1_fit(taxis,integral)
	if saturation:
		plot(taxis,t1_sat_fitfunc(p,taxis))
		plot(taxis,t1_sat_fitfunc(r_[p[0:2],p[2]*1.2],taxis),'y')
		plot(taxis,t1_sat_fitfunc(r_[p[0:2],p[2]*0.8],taxis),'y')
		p = r_[p[0],0,p[2]]
		xlabel(r'$M(t)=%0.3g+(%0.3g-%0.3g)e^{-t/%0.4g}$'%(p[0],p[1],p[0],p[2]))
	else:
		plot(taxis,t1_fitfunc(p,taxis))
		plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*1.2],taxis),'y')
		plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*0.8],taxis),'y')
		xlabel(r'$M(t)=%0.3g+(%0.3g-%0.3g)e^{-t/%0.4g}$'%(p[0],p[1],p[0],p[2]))
	title(titlestr)
	if return_fit:
		return p
	#}}}
#}}}
#{{{ process cpmg 
def process_cpmg(file):
	data = load_file(file)
	data.ft('t2')
	findmax = abs(data)
	findmax.sum('echo')
	findmax = findmax.run(argmax,'t2').data
	data = data['t2',findmax]
	data.data *= phaseopt(data.data) # I added this in, not sure why it was gone!
	return data
#}}}
#{{{regularization
def regularize1d(b,t,tau,alpha):
	# for Lx=b
	if size(b) != size(t):
		print "ERROR, size of b doesn't match size of t"
	tau = tau.reshape(1,-1)
	t = t.reshape(-1,1)
	L = exp(-t/tau)
	U,s,V = svd(L,full_matrices=0)
	rms = zeros(len(alpha),dtype='double')
	coeff = zeros((size(tau),size(alpha)),dtype='double')
	fit = zeros((size(t),size(alpha)),dtype='double')
	for j in range(0,len(alpha)):
		S = diag(s / (s**2 + alpha[j]**2))
		x = dot(
				dot(
					conj(transpose(V)),
					dot(S,conj(transpose(U))))
				,b)# was b
		fit[:,j] = dot(L,x)
		try:
			coeff[:,j] = x.flatten()
		except:
			print 'shape(coeff)',shape(coeff),'shape(x)',shape(x)
			print 'first',shape(coeff[:,j]),'second',shape(x.reshape(-1,1))
			raise
		rms[j] = linalg.norm(fit[:,j]-b)
	return (coeff,fit,rms)
#}}}
#{{{ matched filter
def matched_filter(data,along_dim,decay_rate = 1,return_fit=False):
	#{{{ actually find the filter	
	data_abs = abs(data)
	timeaxis = data_abs.getaxis('t2')
	labels = list(data.dimlabels)
	labels.pop(labels.index(along_dim))
	for thisdim in labels:
		data_abs.mean(thisdim)
	p = exp_fit(timeaxis,data_abs.data)
	#}}}
	#{{{ actually apply the filter
	filter = ndshape(data)
	for thisdim in labels:
		filter.pop(thisdim)
	filter = filter.alloc()
	if (not return_fit): # don't mess with it if we want to check the fit
		p[1] /= decay_rate
	filter.data = exp_fitfunc(p,data.getaxis(along_dim).copy())
	#print 'for matched filter, the x axis is ',data.getaxis(along_dim).copy()
	if not return_fit:
		filter.data[:] -= filter.data.flatten()[-1] # drop so that end is at zero (since we have a noise baseline)
		filter.data[:] /= filter.data.flatten()[0] # normalize
	filter.labels(['t2'],[data_abs.getaxis('t2')])
	return filter
	#}}}
#}}}
def __baseline_gen_L(data_filt1):
	x = data_filt1.getaxis('t2').copy()
	x_norm = x[-1].copy()
	x = x/x_norm # normalize, otherwise we get ridiculously large higher order terms
	#L = array([ones(shape(x)),x,x**2/2,x**3/6,x**4/24,x**5/120]).T
	L = array([ones(shape(x)),x,x**2/2]).T
	return x,x_norm,L
def baseline_spectrum(data,showplots=False,threshold=10,check_filter=False,set_errorlevel=False):
	#print 'diagnose: start baseline'
	data = data.copy()
	prelim_offset = mean(r_[data.data[0],data.data[1],data.data[-2],data.data[-1]])
	data.data[:] -= prelim_offset
	data.ift('t2',shiftornot=True)
	if check_filter:
		clf()
		#plot(abs(data))
		plot(matched_filter(data,'t2',return_fit=True))
		legend(['fits to']) #legend(['data','fits to'])
		twinx()
		plot(matched_filter(data,'t2'))
		legend(['filter with'])
		return
	#{{{ make the abs of a broadened spectrum	
	data_widefilt = data.copy() * matched_filter(data,'t2',decay_rate=10)
	data_widefilt.ft('t2',shiftornot=True)
	#data_widefilt.data *= phaseopt(data_widefilt.data)
	#data_widefilt.data = real(data_widefilt.data)
	#data_widefilt.data /= sign(data_widefilt.data[argmax(abs(data_widefilt.data))])
	data_widefilt.data -= (data_widefilt.data[0]+data_widefilt.data[-1])/2
	data_widefilt = abs(data_widefilt)
	mask = (data_widefilt.data<(data_widefilt.data.max()/threshold)) # generate the mask according to the threshold
	#}}}
	#{{{ ft the data
	data.ft('t2',shiftornot=True)
	#}}}
	if sum(mask)==0:
		erroronnopeak = False
		if erroronnopeak:
			legendstring = []
			plot(abs(data))
			legendstring += ['mf data']
			plot(data_widefilt)
			legendstring += ['wide filter']
			legend(legendstring)
			error_plot("Error -- not able to find any non-baseline data")
		else:
			print "Warning, fit entire spectrum to baseline.\n\n"
		mask = (data_widefilt.data<(data_widefilt.data.max()/1.5))
		mask[:] = True
	data_baseline = data.copy()
	data_baseline.data = data_baseline.data[mask]
	data_baseline.axis_coords[data_baseline.dimlabels.index('t2')] = data_baseline.getaxis('t2')[mask]
	legendstring = []
	x,x_norm,L = __baseline_gen_L(data_baseline)
	#print 'diagnose: ready to invert'
	try:
		fit_coeff = dot(pinv(L,rcond=1e-5),data_baseline.data) # L * fit_coeff = data
	except:
		raise CustomError(maprep('problem inverting:',shape(L),shape(data_baseline.data)))
	#print 'diagnose: inverted'
	if set_errorlevel:
		if any(abs(dot(L,fit_coeff))>set_errorlevel):
			showplots = True
	if showplots:
		#plot(abs(data))
		#legendstring += ['data']
		plot(data_widefilt)
		legendstring += ['wide filter']
		plot(abs(data_baseline))
		legendstring += ['baseline portion']
		show_L = False
		if show_L:
			plot(x*x_norm,L)
			legendstring += map(
					(lambda x:'L'+str(x)),
					range(1,1+L.shape[1])
					)
		plot(x*x_norm,abs(dot(L,fit_coeff)))
		legendstring += ['fit to baseline']
	x,x_norm,L = __baseline_gen_L(data)
	if showplots:
		plot(x*x_norm,dot(L,fit_coeff))
		legendstring += ['entire fit']
		legend(legendstring,'best')
	baseline_data = nddata(prelim_offset+dot(L,fit_coeff),[size(x),1],['t2','power'])
	baseline_data.labels(['t2'],[x])
	#print 'diagnose: shape of baseline ',ndshape(baseline_data)
	return baseline_data

#{{{ plot_noise
def plot_noise(path,j,calibration,mask_start,mask_stop,rgmin=0,k_B = None,smoothing = False, both = False, T = 293.0,plottype = 'semilogy'):
	'''plot noise scan as resistance'''
	data = load_file(r'%s%d'%(path,j),calibration=calibration)
	k_B = 1.3806504e-23
	data.ft('t2',shift = True)
	newt2 = r'F2 / $Hz$'
	data.rename('t2',newt2)
	v = bruker_load_acqu(r'%s%d/'%(path,j))
	dw = 1/v['SW_h']
	dwov = dw/v['DECIM']
	rg = v['RG']
	aq = v['TD']*dw
	if rg>rgmin:
		plotdata = abs(data)
		plotdata.data **= 2
		johnson_factor = 4.0*k_B*T
		plotdata.data /= (aq*johnson_factor)
		t = data.getaxis(newt2)
		mask = logical_and(t>mask_start,
			t<mask_stop)
		avg = plotdata.data[mask].mean() 
		retval = []
		if both or not smoothing:
			plot(plotdata,'-',alpha=0.5,plottype = plottype)
			retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, mean %0.1f'%(dwov*1e6,rg,avg)]
			axis('tight')
		if smoothing:
			# begin convolution
			originalt = plotdata.getaxis(newt2).copy()
			plotdata.ft(newt2,shift = True)
			sigma = smoothing
			siginv = 0.5*sigma**2 # here, sigma is given in the original units (i.e. what we're convolving)
			t = plotdata.getaxis(newt2)
			g = exp(-siginv*t.copy()**2) # we use unnormalized kernel (1 at 0), which is not what I thought!
			plotdata.data *= g
			plotdata.ift(newt2,shift = True)
			t = plotdata.getaxis(newt2).copy()
			t[:] = originalt
			# end convolution
			plot(plotdata,'-',alpha=0.5,plottype = plottype)
			retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, mean %0.1f'%(dwov*1e6,rg,avg)]
			axis('tight')
		return retval
	else:
		return []
#}}}
