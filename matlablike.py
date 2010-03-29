from pylab import *
import textwrap
#rc('image',aspect='auto',interpolation='bilinear')
rc('image',aspect='auto',interpolation='nearest')
rcParams['xtick.direction'] = 'out'
rcParams['xtick.major.size'] = 12
rcParams['xtick.minor.size'] = 6
rcParams['ytick.direction'] = 'out'
rcParams['ytick.major.size'] = 12
rcParams['ytick.minor.size'] = 6
rcParams['legend.fontsize'] = 12
#{{{ errors
class CustomError(Exception):
	def __init__(self, *value):
		if len(value)>1:
			self.value = map(str,value)
		else:
			self.value = value
	def __str__(self):
		retval = ' '.join(self.value)
		return '\n'+'\n'.join(textwrap.wrap(retval,90))
def maprep(*mylist):
	mylist = list(mylist)
	for j in range(0,len(mylist)):
		if type(mylist[j]) is not str:
			mylist[j] = mylist[j].__repr__()
	return ' '.join(mylist)
#}}}
#{{{ indeces to slice
#}}}
#{{{ add slashes for dir's
def dirformat(file):
		#{{{ format strings
		if file[-1]!='/':
			file += '/'
		#}}}
		return file
#}}}
#{{{ old grid and tick
def gridandtick(ax,rotation=(0,0),precision=(2,2),labelstring=('',''),gridcolor=r_[0,0,0],formatonly = False,fixed_y_locator = None,logarithmic = False):
	if not formatonly:
		#{{{x ticks
		# determine the size
		width = abs(diff(ax.get_xlim()))
		if width==0:
			raise CustomError('x axis width is zero')
		widthexp = floor(log(width)/log(10.))-1
		scalefactor = 10**widthexp
		width /= scalefactor
		majorLocator = MultipleLocator(5*scalefactor)
		#majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[0]+'f'+labelstring[0])# labelstring can be used, for instance, for pi
		#ax.xaxis.set_major_formatter(majorFormatter)
		minorLocator   = MultipleLocator(1*scalefactor)
		ax.xaxis.set_major_locator(majorLocator)
		#for the minor ticks, use no labels; default NullFormatter
		ax.xaxis.set_minor_locator(minorLocator)
		#}}}
		#{{{ y ticks
		width = abs(diff(ax.get_ylim()))
		if width==0:
			raise CustomError('y axis width is zero')
		widthexp = floor(log(width)/log(10.))-1
		scalefactor = 10**widthexp
		width /= scalefactor
		if fixed_y_locator == None:
			if logarithmic:
				majorLocator = LogLocator(10)
			else:
				majorLocator   = MultipleLocator(5*scalefactor)
		else:
			majorLocator   = MultipleLocator(fixed_y_locator[4::5])
		#majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[1]+'f'+labelstring[1])# labelstring can be used, for instance, for pi
		#ax.yaxis.set_major_formatter(majorFormatter)
		if fixed_y_locator == None:
			if logarithmic:
				minorLocator = LogLocator(10,subs=r_[0:11])
			else:
				minorLocator   = MultipleLocator(1*scalefactor)
		else:
			minorLocator   = FixedLocator(fixed_y_locator)
		ax.yaxis.set_major_locator(majorLocator)
		#for the minor ticks, use no labels; default NullFormatter
		ax.yaxis.set_minor_locator(minorLocator)
		#}}}
	grid(True,which='major',color=gridcolor,alpha=0.4,linestyle='-')
	grid(True,which='minor',color=gridcolor,alpha=0.2,linestyle='-')
	labels = ax.get_xticklabels()
	setp(labels,rotation=rotation[0],fontsize=10)
	labels = ax.get_yticklabels()
	setp(labels,rotation=rotation[1],fontsize=10)
	fig = gcf()
	fig.autofmt_xdate()
def gridon(gridcolor=r_[0,0,0]):
	grid(True,which='major',color=gridcolor,alpha=0.4,linestyle='-')
	grid(True,which='minor',color=gridcolor,alpha=0.2,linestyle='-')
#}}}
#{{{ a better version?
def othergridandtick(ax,rotation=(0,0),precision=(2,2),labelstring=('',''),gridcolor=r_[0,0,0]):
	#{{{x ticks
	# determine the size
	ax.xaxis.set_major_locator(MaxNLocator(10)) # could use multiplelocator if it keeps try to do multiples of 2
	ax.xaxis.set_minor_locator(MaxNLocator(50))
	#}}}
	#{{{ y ticks
	ax.yaxis.set_major_locator(MaxNLocator(10))
	ax.yaxis.set_minor_locator(MaxNLocator(50))
	#}}}
	grid(True,which='major',color=gridcolor,alpha=0.4,linestyle='-')
	grid(True,which='minor',color=gridcolor,alpha=0.2,linestyle='-')
	labels = ax.get_xticklabels()
	setp(labels,rotation=rotation[0],fontsize=10)
	labels = ax.get_yticklabels()
	setp(labels,rotation=rotation[1],fontsize=10)
#}}}
#{{{ plot wrapper
global OLDplot
OLDplot = plot
global myplotfunc
myplotfunc = OLDplot
def whereblocks(a): # returns contiguous chunks where the condition is true
	parselist = where(a)[0]
	jumps_at = where(diff(parselist)>1)[0]+1
	retlist = []
	lastjump = 0
	for jump in jumps_at:
		retlist += [parselist[lastjump:jump]]
		lastjump = jump
	retlist += [parselist[lastjump:]]
	return retlist
def autolegend(legendstr):
	lg = legend(legendstr,'best')
	lg.get_frame().set_alpha(0.65)
def plot(*args,**kwargs):
	global myplotfunc
	#{{{ all possible properties
	myformat = None 
	myxlabel = None
	myylabel = None
	myx = None
	myy = None
	#}}}
	try:
		#{{{assign all the possible combinations
		if len(args)==1:
			myy = args[0]
		elif (len(args)==2) and (type(args[1]) is str):
			myy = args[0]
			myformat = args[1]
		else:
			myx = args[0]
			myy = args[1]
		if len(args)==3:
			myformat = args[2]
		#if (type(myx) is float64) or (type(myx) is complex128):
		if isscalar(myx):
			myx = array([myx])
		#if (type(myy) is float64) or (type(myy) is complex128):
		if isscalar(myy):
			myy = array([myy])
		#}}}
		#{{{ parse nddata
		if isinstance(myy,nddata):
			if (len(myy.dimlabels)>1):
				mytitle = myy.dimlabels[1]
			if (len(myy.dimlabels)>0):
				myxlabel = myy.dimlabels[0]
			if (myx == None):
				try:
					if (len(myy.axis_coords[0])>0):
						myx = myy.axis_coords[0]
				except:
					raise CustomError("problem, ",shape(myy.axis_coords),"not valid axes for",ndshape(myy))
			myy = myy.data
		#}}}
		#{{{ semilog where appropriate
		if (myx != None) and (len(myx)>1): # by doing this and making myplotfunc global, we preserve the plot style if we want to tack on one point
			b = diff(log10(myx))
			if (size(b)>3) and all(abs((b-b[0])/b[0])<1e-4) and not ('nosemilog' in kwargs.keys()):
				myplotfunc = semilogx
			else:
				myplotfunc = OLDplot
			if ('nosemilog' in kwargs.keys()):
				kwargs.pop('nosemilog')
		if 'plottype' in kwargs.keys():
			if kwargs['plottype'] == 'semilogy':
				myplotfunc = semilogy
			kwargs.pop('plottype')
		#}}}
		#{{{ hsv plots when we have multiple lines
		if len(shape(myy))>1:
			#{{{ hsv plots
			hold(True)
			retval = []
			for j in range(0,myy.shape[1]):
				#{{{ this is the way to assign plot arguments
				plotargs = [myx,myy[:,j],myformat]
				while None in plotargs:
					plotargs.remove(None)
				#}}}
				#{{{ here, i update the kwargs to include the specific color for this line
				newkwargs = kwargs.copy() # kwargs is a dict
				newkwargs.update({'color':cm.hsv(double(j)/double(myy.shape[1]))})
				#}}}
				retval += [myplotfunc(*tuple(plotargs),**newkwargs)]
			#hold(False)
			#}}}
			#}}}
		else:
			plotargs = [myx,myy,myformat]
			while None in plotargs:
				plotargs.remove(None)
			retval = myplotfunc(*tuple(plotargs),**kwargs)
		#{{{ attach labels and such
		if (myxlabel!=None):
			xlabel(myxlabel)
		if (myylabel!=None):
			ylabel(myylabel)
		axis('tight')
		grid(True)
		#}}}
	except:
		message = 'Error plotting,'
		message += ':'
		print r'\begin{verbatim}',CustomError(message,'myformat',myformat,
				'myxlabel',type(myxlabel),shape(myxlabel),
				'myylabel',type(myylabel),shape(myylabel),
				'myy',type(myy),shape(myy),
				'myx',type(myx),shape(myx)).__str__(),r'\end{verbatim}'
		raise
	return retval
#}}}

#{{{nddata
#{{{ shaping and allocating
class ndshape ():
	def __init__(self,*args):
		if len(args) == 2:
			self.shape = list(args[0])
			self.dimlabels = args[1]
		if len(args) == 1: #assum that it's an nddata object
			self.shape = list(args[0].data.shape)
			self.dimlabels = list(args[0].dimlabels)
		return
	def __setitem__(self,reference,setto):
		self.shape[self.dimlabels.index(reference)] = setto
		return
	def copy(self,*args):
		return ndshape(list(self.shape),list(self.dimlabels))
	def __add__(self,arg):
		shape = arg[0]
		dimlabels = arg[1]
		self.shape = shape + self.shape
		self.dimlabels = dimlabels + self.dimlabels
		return self
	def __repr__(self): #how it responds to print
		return zip(self.shape,self.dimlabels).__repr__()
	def __getitem__(self,args):
		mydict = dict(zip(self.dimlabels,self.shape))
		return mydict[args]
	def pop(self,label):
		thisindex = self.dimlabels.index(label)
		self.dimlabels.pop(thisindex)
		self.shape.pop(thisindex)
	def alloc(self,dtype='complex128'):
		return nddata(zeros(tuple(self.shape),dtype=dtype),self.shape,self.dimlabels)
#}}}
#{{{ format out to a certain decimal place
def dp(number,decimalplaces,scientific=False):
	if scientific:
		tenlog = floor(log(number)/log(10.))
		number /= 10**tenlog
		fstring = '%0.'+'%d'%decimalplaces+r'f\times 10^{%d}'%tenlog
	else:
		fstring = '%0.'+'%d'%decimalplaces+'f'
	return fstring%number
#}}}
#{{{ concatenate datalist along dimname
def concat(datalist,dimname):
	#{{{ allocate a new datalist structure	
	t1size = 0
	shapes = map(ndshape,datalist)
	for j in range(0,len(datalist)):
		if dimname in shapes[j].dimlabels:
			t1size += shapes[j][dimname]
			shapetocheck = list(shapes[j].shape)
			shapetocheck.pop(shapes[j].dimlabels.index(dimname))
		else:
			t1size += 1
			shapetocheck = list(shapes[j].shape)
		if j is 0:
			shapetocheckagainst = shapetocheck
		else:
			if any(~(array(shapetocheck) == array(shapetocheckagainst))):
				raise CustomError('For item ',j,'in concat, ',shapetocheck,'!=',shapetocheckagainst,'where all the shapes of the things you\'re trying to concat are:',shapes)
	newdatalist = ndshape(datalist[-1])
	if dimname in newdatalist.dimlabels:
		newdatalist[dimname] = t1size
	else:
		newdatalist += ([t1size],[dimname])
	newdatalist = newdatalist.alloc()
	#}}}
	#{{{ actually contact the datalist
	t1size = 0 # now use it to track to position
	for j in range(0,len(datalist)):
		if dimname in shapes[j].dimlabels:
			newdatalist[dimname,t1size:t1size+shapes[j][dimname]] = datalist[j]
			t1size += shapes[j][dimname]
		else:
			newdatalist[dimname,t1size:t1size+1] = datalist[j]
			t1size += 1
	#}}}
	#{{{ pull the axis labels from the last item in the list
	if len(datalist[-1].axis_coords)>0:
		dimlabels = list(datalist[-1].dimlabels)
		axis_coords = list(datalist[-1].axis_coords)
		if dimname in dimlabels:
			thisindex = dimlabels.index(dimname)
			dimlabels.pop(thisindex)
			axis_coords.pop(thisindex)
		dimlabels += [dimname]
		axis_coords += [r_[0:t1size]]
		newdatalist.labels(dimlabels,axis_coords)
	#}}}
	return newdatalist
#}}}
class nddata ():
	def __init__(self,data,sizes,dimlabels,axis_coords=[],ft_start_time = 0.):
		if not (type(data) is ndarray):
			#if (type(data) is float64) or (type(data) is complex128) or (type(data) is list):
			if isscalar(data) or (type(data) is list):
				data = array(data)
			else:
				raise CustomError('data is not an array, it\'s',type(data),'!')
		if not (type(dimlabels) is list):
			raise CustomError('labels are not a list')
		self.data = reshape(data,sizes)
		self.dimlabels = dimlabels
		self.axis_coords = axis_coords
		self.ft_start_time = ft_start_time
		return
	#{{{ add any dimensions to self that are not present in other
	def matchdims(self,other):
		#print 'diagnose: matching',ndshape(self),'to',ndshape(other)
		addeddims =  list(set(self.dimlabels)^set(other.dimlabels))
		newdims = addeddims + self.dimlabels
		newshape = [1]*len(addeddims) + list(self.data.shape)
		#print 'diagnose: newshape',newshape,'newdims',newdims
		#{{{ reshape to the new dimensions	
		new_axis_coords = [r_[1]]*len(addeddims) + self.axis_coords
		self.data = self.data.reshape(newshape)
		self.dimlabels = newdims
		if len(self.axis_coords)>0:
			self.axis_coords = new_axis_coords
		#}}}
		#{{{ if we are adding dimensions, we will need to reorder to match the order of the other	
		if len(addeddims)>0:
			self.reorder(other.dimlabels)
		#}}}
		return self
	#}}}
	def rename(self,previous,new):
		self.dimlabels[self.dimlabels.index(previous)] = new
		return self
	def __sub__(self,arg):
		a = self.copy()
		b = arg.copy()
		if (shape(arg.data)!=shape(a.data)):
			#print 'diagnose: ndshape(a)',ndshape(a)
			#print '\n\n'
			#print 'diagnose: ndshape(b)',ndshape(b)
			#print '\n\n'
			#print r'diagnose: axis\_coords',map(len,a.axis_coords)
			#print '\n\n'
			a.matchdims(b)
			b.matchdims(a)
			#print 'diagnose:',ndshape(a)
			#print 'diagnose:',ndshape(b)
		try:
			a.data -= b.data
		except:
			print '|-ERROR SUBTRACTING NDDATA-----------'
			print '|','a.data.shape:',a.data.shape
			print '|','arg.data.shape:',arg.data.shape
			print '|------------------------------------'
			raise
		return a
	def __add__(self,arg):
		a = self.copy()
		a.data += arg.data
		return a
	def __mul__(self,arg):
		#{{{ do scalar multiplication
		if isscalar(arg):
			newdata = self.copy()
			newdata.data *= arg
			return newdata
		#}}}
		#{{{ collect the list of dimensions
		newdimlabels = list(self.dimlabels)
		temp = list(arg.dimlabels)
		for item in temp:
			if not (item in newdimlabels):
				newdimlabels += [item]
		#}}}
		#{{{collect the new orderings
		leftorder = []
		leftsize = ones(len(newdimlabels))
		j=0
		for item in newdimlabels:
			if item in self.dimlabels:
				thisindex = self.dimlabels.index(item)
				leftorder += [thisindex]
				leftsize[j] = self.data.shape[thisindex]
			j+=1
		#}}}
		rightorder = []
		rightsize = ones(len(newdimlabels))
		j=0
		for item in newdimlabels:
			if item in arg.dimlabels:
				thisindex = arg.dimlabels.index(item)
				rightorder += [thisindex]
				rightsize[j] = arg.data.shape[thisindex]
			j+=1
		result = self.data.copy().transpose(leftorder).reshape(leftsize) * arg.data.transpose(rightorder).reshape(rightsize)
		if len(self.axis_coords)>0:
			return nddata(result,list(result.shape),newdimlabels,axis_coords=list(self.axis_coords))
		else:
			return nddata(result,list(result.shape),newdimlabels)
	def __div__(self,arg):
		#{{{ collect the list of dimensions
		newdimlabels = list(self.dimlabels)
		temp = list(arg.dimlabels)
		for item in temp:
			if not (item in newdimlabels):
				newdimlabels += [item]
		#}}}
		#{{{collect the new orderings
		leftorder = []
		leftsize = ones(len(newdimlabels))
		j=0
		for item in newdimlabels:
			if item in self.dimlabels:
				thisindex = self.dimlabels.index(item)
				leftorder += [thisindex]
				leftsize[j] = self.data.shape[thisindex]
			j+=1
		rightorder = []
		rightsize = ones(len(newdimlabels))
		j=0
		for item in newdimlabels:
			if item in arg.dimlabels:
				thisindex = arg.dimlabels.index(item)
				rightorder += [thisindex]
				rightsize[j] = arg.data.shape[thisindex]
			j+=1
		result = self.data.copy().transpose(leftorder).reshape(leftsize) / arg.data.transpose(rightorder).reshape(rightsize)
		if len(self.axis_coords)>0:
			return nddata(result,list(result.shape),newdimlabels,axis_coords=list(self.axis_coords))
		else:
			return nddata(result,list(result.shape),newdimlabels)
	def integrate(self,thisaxis):
		self.run_nopop(cumsum,thisaxis)
		if len(self.axis_coords)>0:
			t = self.getaxis(thisaxis)
			dt = t[1]-t[0]
			self.data *= dt
		return self
	def sum(self,axes):
		if (type(axes) is str):
			axes = [axes]
		for j in range(0,len(axes)):
			try:
				thisindex = self.dimlabels.index(axes[j])
			except:
				print '|-ERROR FINDING DIMENSION-----'
				print '| dimlabels is: ',self.dimlabels
				print "| doesn't contain: ",axes[j]
				print '|-----------------------------'
				raise
			self.data = sum(self.data,
					axis=thisindex)
			self.dimlabels.pop(thisindex)
			if self.axis_coords!=[]:
				self.axis_coords.pop(thisindex)
		return self
	def argmax(self,axes):
		if (type(axes) is str):
			axes = [axes]
		for j in range(0,len(axes)):
			try:
				thisindex = self.dimlabels.index(axes[j])
			except:
				print 'error, dimlabels is: ',self.dimlabels
				print "doesn't contain: ",axes[j]
				raise
			self.data = argmax(self.data,
					axis=thisindex)
			self.dimlabels.pop(thisindex)
			if self.axis_coords!=[]:
				self.axis_coords.pop(thisindex)
		return self
	def mean(self,axes):
		if (type(axes) is str):
			axes = [axes]
		for j in range(0,len(axes)):
			try:
				thisindex = self.dimlabels.index(axes[j])
			except:
				print 'error, dimlabels is: ',self.dimlabels
				print "doesn't contain: ",axes[j]
				raise
			self.data = mean(self.data,
					axis=thisindex)
			self.dimlabels.pop(thisindex)
			if self.axis_coords!=[]:
				self.axis_coords.pop(thisindex)
		return self
	def mean_nopop(self,axis):
		self = self.run_nopop(mean,axis=axis)
		return self
	def sum_nopop(self,axes):
		if (type(axes) is str):
			axes = [axes]
		for j in range(0,len(axes)):
			try:
				thisindex = self.dimlabels.index(axes[j])
			except:
				print 'error, dimlabels is: ',self.dimlabels
				print "doesn't contain: ",axes[j]
				raise
			temp = list(self.data.shape)
			temp[thisindex] = 1
			self.data = sum(self.data,
					axis=thisindex)
			self.data = self.data.reshape(temp)
		return self
	def popdim(self,dimname):
		thisaxis = self.dimlabels.index(dimname)
		thisshape = list(self.data.shape)
		if thisshape[thisaxis]!=1:
			raise CustomError("trying to pop a dim that's not length 1")
		thisshape.pop(thisaxis)
		self.data = self.data.reshape(thisshape)
		self.dimlabels.pop(thisaxis)
		self.axis_coords.pop(thisaxis)
		return self
	def ft(self,axes,shiftornot=False,shift=None):
		if shift != None:
			shiftornot = shift
		if (type(axes) is str):
			axes = [axes]
		if not (type(shiftornot) is list):
			shiftornot = [bool(shiftornot)]
		for j in range(0,len(axes)):
			try:
				thisaxis = self.dimlabels.index(axes[j])
			except:
				raise CustomError('error, dimlabels is: ',self.dimlabels)
			self.data = fft(self.data,axis=thisaxis)
			if bool(shiftornot[j]):
				self.data = fftshift(self.data,axes=[thisaxis])
			if len(self.axis_coords)>0:
				t = self.axis_coords[thisaxis]
				dt = t[1]-t[0]
				self.ft_start_time = t[0]
				self.data *= dt
				if bool(shiftornot[j]):
					self.axis_coords[thisaxis] = linspace(-0.5/dt,0.5/dt,size(t))
				else:
					self.axis_coords[thisaxis] = linspace(0,1./dt,size(t))
		return self
	def ift(self,axes,shiftornot=False,shift=None):
		if shift != None:
			shiftornot = shift
		if (type(axes) is str):
			axes = [axes]
		if not (type(shiftornot) is list):
			shiftornot = [bool(shiftornot)]
		for j in range(0,len(axes)):
			try:
				thisaxis = self.dimlabels.index(axes[j])
			except:
				raise CustomError('error, dimlabels is: ',self.dimlabels)
			if bool(shiftornot[j]):
				self.data = ifftshift(self.data,axes=[thisaxis])
			self.data = ifft(self.data,axis=thisaxis)
			if len(self.axis_coords)>0:
				t = self.axis_coords[thisaxis]
				dt = t[1]-t[0]
				self.data *= size(t) * dt # here, the algorithm divides by N, so for integration, we need to not do that
				#{{{ shiftornot specifies the shifting of the initial ft, not this result, so we always return a 0->1 time axis
				self.axis_coords[thisaxis] = linspace(0,1./dt,size(t)) + self.ft_start_time # note that I offset by ft_start_time, which I pull from when I ft'd
				#}}}
		return self
	def __abs__(self):
		out = self.copy()
		out.data = abs(out.data)
		return out
	def circshift(self,axis,amount):
		if amount!=0:
			newdata = ndshape(self).alloc(dtype=self.data.dtype)
			newdata[axis,:-amount] = self[axis,amount:]
			newdata[axis,-amount:] = self[axis,:amount]
			self.data = newdata.data
		return self
	def labels(self,listofstrings,listofaxes):
		self.axis_coords = [[]]*len(self.dimlabels)
		for j in range(0,len(listofstrings)):
			try:
				self.axis_coords[self.dimlabels.index(listofstrings[j])] = listofaxes[j]
			except:
				try:
					raise CustomError("Can't assign the coordinates to"+listofstrings[j]+"as"+listofaxes[j].__repr__())
				except:
					raise CustomError("listofaxes (",len(listofaxes),") isn't same length as ",listofstrings)
		return self
	def sort(self,axisname):
		whichaxis = self.dimlabels.index(axisname)
		order = argsort(self.axis_coords[whichaxis])
		self[axisname,:] = self[axisname,order]
		self.axis_coords[whichaxis] = self.axis_coords[whichaxis][order]
		return self
	def copyaxes(self,other):
		# in the case that the dimensions match, and we want to copy the labels
		self.axis_coords = other.axis_coords
		return self
	def getaxis(self,axisname):
		return self.axis_coords[self.dimlabels.index(axisname)]
	def run(self,func,axis):
		thisaxis = self.dimlabels.index(axis)
		self.data = func(self.data,axis=thisaxis)
		self.dimlabels.pop(thisaxis)
		if self.axis_coords!=[]:
			self.axis_coords.pop(thisaxis)
		return self
	def run_nopop(self,func,axis):
		thisaxis = self.dimlabels.index(axis)
		temp = list(self.data.shape)
		temp[thisaxis] = 1
		self.data = func(self.data,axis=thisaxis)
		#{{{ if the function doesn't rip out the dim, make sure we don't change the dims
		if len(self.data.shape)==len(temp):
			temp[thisaxis] = self.data.shape[thisaxis]
		#}}}
		self.data = self.data.reshape(temp)
		return self
	def smash(self,listofhowmany):
		names = ['error'] * len(listofhowmany)
		newsize = []
		laststop = 0 
		j = 0
		for ndims in listofhowmany:
			thisstop = laststop + ndims
			newsize += [prod(self.data.shape[laststop:thisstop])]
			names[j] = ' x '.join(self.dimlabels[laststop:thisstop])
			laststop = thisstop
			j+=1
		self.data = self.data.reshape(newsize)
		self.dimlabels = names
		return self
	def smashorder(self,listoflists):
		self.reorder(concatenate(listoflists))
		self.smash(map(len,listoflists))
		return self
	def reorder(self,axes):
		neworder = map(self.dimlabels.index,axes)
		self.dimlabels = map(self.dimlabels.__getitem__,neworder)
		if len(self.axis_coords)>0:
			try:
				self.axis_coords = map(self.axis_coords.__getitem__,neworder)
			except:
				raise CustomError('problem mapping',map(len,self.axis_coords),'onto',neworder)
		try:
			self.data = self.data.transpose(neworder)
		except ValueError:
			raise "you didn't supply all the dimensions for the reorder"
		return self
	def __getslice__(self,*args):
		print 'getslice! ',args
	def __setitem__(self,*args):
		if isinstance(args[1],nddata):
			rightdata = args[1].data
			#print 'check size:',self.data.shape,args[1].data.shape
			rightlabels = args[1].dimlabels
		else: # assume it's an ndarray
			rightdata = args[1]
		#{{{ build up the left index list
		leftindex = [slice(None,None,None)]*len(self.dimlabels) # by default, just pass everything, to cut the :,:,: notation
		for j in range(0,len(args[0]),2):
			try:
				thisindex = self.dimlabels.index(args[0][j])
			except:
				print 'ERROR------------------------'
				print args[0][j],' not found'
				print '-----------------------------'
				raise
			leftindex[thisindex] = args[0][j+1]
		#}}}
		try:
			self.data[tuple(leftindex)] = rightdata.squeeze() # assign the data
		except:
			raise CustomError('ERROR ASSIGNING NDDATA:','rightdata.shape:',rightdata.shape,'self.data.shape:',self.data.shape,' dimlabels:',self.dimlabels,' leftindex:',tuple(leftindex),'--> shape of left slice: ',self.data[tuple(leftindex)].shape)
		#print 'STOP SETITEM --------------------------------'
	def copy(self):
		return nddata(self.data.copy(),list(self.data.shape),list(self.dimlabels),axis_coords=map(copy,self.axis_coords),ft_start_time = self.ft_start_time)
	def __getitem__(self,args):
		#print 'START GETITEM -------------------------------'
		#print 'getitem ',args
		#print self.dimlabels
		#print 'args=',args
		if type(args[0]) is str:
			indexlist = [slice(None,None,None)]*len(self.dimlabels) #initialize to all none
			newlabels = list(self.dimlabels)
			if len(self.axis_coords)>0:
				new_axis_coords = list(self.axis_coords)
			#{{{  list of where indeces are on both sides, and what we want to index with
			for j in range(0,len(args),2):
				thisindex = self.dimlabels.index(args[j])
				newlabelindex = newlabels.index(args[j])
				#print 'indexing dimension ',self.dimlabels[thisindex],' with ',args[j+1]
				indexlist[thisindex] = args[j+1]
				if not (type(args[j+1]) is slice):
					newlabels.pop(newlabelindex)
			#}}}
			indexlist = tuple(indexlist)
			if len(self.axis_coords)>0:
				for j in range(0,len(self.axis_coords)):
					try:
						new_axis_coords[j] = new_axis_coords[j][indexlist[j]]
					except:
						raise CustomError('problem reassigning new_axis_coords',new_axis_coords,'indexlist',indexlist)
				k = 0 # we need this, since the index will change as we pop
				for j in range(0,len(self.axis_coords)):
					if (type(new_axis_coords[k]) is int32) or (type(new_axis_coords[k]) is double):
						new_axis_coords.pop(k)
					else:
						k += 1
			try:
				if len(self.axis_coords)>0:
					return nddata(self.data[indexlist],
							self.data[indexlist].shape,
							newlabels,
							axis_coords = new_axis_coords)
				else:
					return nddata(self.data[indexlist],self.data[indexlist].shape,newlabels)
			except:
				print '|-ERROR RETURNING INDEXED NDDATA-------'
				print '| ','self.data.shape:',self.data.shape
				print '| ','indexlist:',indexlist
				print '| ','newlabels:',newlabels
				print '|--------------------------------------'
				raise
		else:
			print 'label your freaking dimensions!'
			#print 'STOP GETITEM --------------------------------'
			raise
#}}}

#{{{subplot_dim
class subplot_dim():
	def __init__(self,firstdim,seconddim):
		self.num = r_[firstdim,seconddim,0]
	def set(self,args,x='',g=True,y='',t='',a=''):
		if type(args) is int:
			number = args
			ax = subplot(*tuple(self.num+r_[0,0,number]))
			xlabel(x)
			ylabel(y)
			title(t)
			grid(g)
		elif (type(args) is tuple) and (len(args) is 3):
			# the second value passed is 
			whichsmall = args[2]
			break_into = args[1]
			number = args[0]
			mydims = self.num*r_[1,break_into,1]+r_[
					0,0,break_into*(number-1)+whichsmall]
			try:
				ax = subplot(*tuple(mydims))
			except:
				print 'failed trying subplots: ', mydims
				raise
			xlabel(x)
			ylabel(y)
			title(t)
			grid(g)
		else:
			print "problem, need to pass either 1 or 3 arguments to set"
			print 'type of args: ',type(args)
		return ax
#}}}
def fa(input,dtype='complex128'):# make a fortran array
	return array(input,order='F',dtype=dtype) # will need transpose reverses the dimensions, since the bracketing still works in C order (inner is last index), but F tells it to store it appropriately in memory
def ndgrid(*input):
	thissize = list([1])
	thissize = thissize * len(input)
	output = list()
	for j in range(0,len(input)):
		tempsize = copy(thissize)
		tempsize[j] = input[j].size
		output.append(input[j].reshape(tempsize))
	return output
def pinvr(C,alpha):
	print 'pinvr called'
	U,S,V = svd(C,full_matrices=0)
	#print 'U S V shapes:'
	#print U.shape
	#print S.shape
	#print V.shape
	if any(~isfinite(U)):
		raise CustomError('pinvr error, U is not finite')
	if any(~isfinite(V)):
		raise CustomError('pinvr error, V is not finite')
	if any(~isfinite(S)):
		raise CustomError('pinvr error, S is not finite')
	S = diag(S / (S**2 + alpha**2))
	if any(~isfinite(S)):
		raise CustomError('pinvr error, problem with S/(S^2+alpha^2) --> set your regularization higher')
	return dot(conj(transpose(V)),
			dot(S,conj(transpose(U))))
def sech(x):
	return 1./cosh(x)
def spectrogram(waveform,f_start,f_stop,npoints_fdom=40,tdom_div=2):
	#npoints_tdom = int(round(double(waveform.len)/double(npoints_fdom)))*npoints_tdom_mult
	npoints_tdom = waveform.len/tdom_div # this seems to be more legible than above 
	resolution = diff(waveform.x[0:2])

	sigma = abs(f_start-f_stop)/double(npoints_fdom)
	#print "sigma = %f resolution = %f"%(sigma,resolution)
	if sigma<4*resolution:
		sigma = 4*resolution

	waveform.def_filter(sigma,npoints_tdom)# define the filter and number of points for the spectrogram windowing (define the filter such that the points are spaced sigma apart)

	# go through and apply the filter for some range of points

	f_axis = linspace(f_start,f_stop,npoints_fdom)

	specgram = zeros((npoints_fdom,npoints_tdom),dtype="complex128")

	for j in range(0,npoints_fdom):

		t_axis, specgram[j,:] = waveform.do_filter(f_axis[j])
		#plot(t_axis,abs(specgram[j,:])) # leave this in for testing what it does in the fdom
	#image(specgram,y=f_axis/1e6,x=t_axis*1e6) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
	imshow(abs(specgram),extent=(t_axis[0]*1e6,t_axis[-1]*1e6,f_axis[-1]/1e6,f_axis[0]/1e6)) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
	return gca()
def image(A,x=[],y=[],**kwargs):
	setlabels = False
	if isinstance(A,nddata):
		setlabels = True
		templabels = list(A.dimlabels)
		x_label = templabels[-1]
		templabels.pop(-1)
		y_label = ''
		while len(templabels)>0:
			y_label += templabels.pop(0)
			if len(templabels)>0:
				y_label += r' $\times$ '
		A = A.data
	if type(x) is list:
		x = array(x)
	if type(y) is list:
		y = array(y)
	if len(x)==0:
		x = [1,A.shape[1]]
	else:
		x = x.flatten()
	if len(y)==0:
		y = [1,A.shape[0]]
	else:
		y = y.flatten()
	myext = (x[0],x[-1],y[-1],y[0])
	extralines = 0
	while A.ndim>2:# to substitude for imagehsvm, etc., so that we just need a ersion of ft
		# order according to how it's ordered in the memory
		# the innermost two will form the image -- first add a line to the end of the images we're going to join up
		tempsize = array(A.shape) # make a tuple the right shape
		tempsize[-2] = 1 # all dims are the same except the image row, of which there is only one
		A = concatenate((A,nan*zeros(tempsize)),axis=(A.ndim-2)) # concatenate along the rows
		tempsize = r_[A.shape[0:-3],A.shape[-2:]]
		tempsize[-2] *= A.shape[-3]
		A = A.reshape(tempsize) # now join them up
		++extralines # keep track of the extra lines at the end
	A = A[:A.shape[0]-extralines,:]
	line_mask = isnan(A)
	#A[line_mask] = A[logical_not(line_mask)].max()
	A[line_mask] = 0
	if iscomplex(A).any():
		A = imagehsv(A)
		imshow(A,extent=myext,**kwargs)
	else:
		imshow(A,extent=myext,**kwargs)
		colorbar()
	if setlabels:
		xlabel(x_label)
		#print y_label
		ylabel(y_label)
	return
def colormap(points,colors,n=256):
	r = interp(linspace(0,1,n),points,colors[:,0].flatten())
	g = interp(linspace(0,1,n),points,colors[:,1].flatten())
	b = interp(linspace(0,1,n),points,colors[:,2].flatten())
	return reshape(r_[r,g,b],(3,n)).T
def imagehsv(A):
	n = 256
	theta = (n-1.)*mod(angle(A)/pi/2.0,1)# angle in 255*cycles
	hsv = colormap(r_[0.,1./3.,2./3.,1.],double(array([
		[1,0,0],
		[0,1,0],
		[0,0,1],
		[1,0,0]])),n=n)
	hsv_norm = sqrt(sum(hsv*hsv,axis=1))
	hsv_norm = reshape(hsv_norm,(hsv_norm.size,1))
	hsv = hsv/hsv_norm
	colors = hsv[ix_(int32(theta.flatten().round()),[0,1,2])]
	colors = reshape(colors,(A.shape[0],A.shape[1],3))
	colors *= abs(A).reshape(A.shape[0],A.shape[1],1)
	colors /= abs(A).max()
	return colors
def myfilter(x,center = 250e3,sigma = 100e3):
	x = (x-center)**2
	x /= sigma**2
	return exp(-x)
class ft():
	def __init__(self,x,data):
		# though this is set up for one dimension, I should really change this so that it can deal with multiple dimensions, where x is a list of the x for the various dimensions of data and we also have an additional argument, which gives text labels to the various axes
		self.x = x[:]
		self.dx = self.x[1]-self.x[0]
		self.data = data[:]
		self.filter_defined = False
	def ft(self):
		self.data = fft(self.data)*self.dx
		self.len = len(self.data)
		self.x = r_[0:ceil(double(self.len)/2.)]
		self.x = r_[self.x,-1.*floor(double(self.len)/2.):0]
		self.x /= self.len * self.dx
		self.dx = self.x[1]-self.x[0]
		#print 'size of f ',len(self.f),' size of fdata ',len(self.fdata)
	def ift(self):
		self.len = len(self.data)
		self.dx = self.x[1]-self.x[0]
		self.data = ifft(self.data)*double(len(self.data))*self.dx
		self.x = linspace(0.,1.,len(self.data))/double(self.dx)
		self.dx = self.x[1]-self.x[0]
	def def_filter(self,sigma,filterpoints):
		''' define the filter for the filter function -- only to be used after ft!!!'''
		self.filter_defined = True
		self.filter=r_[r_[0:ceil(double(filterpoints)/2)+1],r_[-floor(double(filterpoints)/2)+1:0]]
		self.filter *= self.dx
		self.filter = myfilter(self.filter,center=0,sigma=sigma)
		self.filterpoints = filterpoints
	def do_filter(self,f_center):
		''' assuming that we have already ft'd, do an ift about a specific frequency
		do NOT mess with data -- just return the filtered data'''
		fdomtest = False
		self.len = len(self.data)
		f_ind = int(round(double(f_center)/self.dx))
		after = ceil(double(self.filterpoints)/2)
		before = floor(double(self.filterpoints)/2)
		if before>f_ind:
			before = f_ind
		if after+f_ind>self.len:
			after = self.len-f_ind
		thisdata = zeros(self.filterpoints,dtype='complex128')
		if fdomtest:
			smallf_axis = zeros(self.filterpoints,dtype='double') # for testing what it does in the fdom
		# set up the data w/ the center frequency at 0
		if before>0:
			thisdata[-before:] = self.data[f_ind-before:f_ind].copy()
			if fdomtest:
				smallf_axis[-before:] = r_[-before:0] # for testing what it does in the fdom
		if after>0:
			thisdata[0:after] = self.data[f_ind:f_ind+after].copy()
			if fdomtest:
				smallf_axis[0:after] = r_[0:after]
		#print 'filterpoints %d size(filter) %d size(thisdata) %d'%(self.filterpoints,size(self.filter),size(thisdata))
		thisdata *= self.filter
		if fdomtest:
			return (smallf_axis*double(self.dx)+f_center, #this is for testing what it does in the fdom
				(thisdata)*double(self.filterpoints)*self.dx)
		else:
			return (linspace(0.,1.,len(thisdata))/double(self.dx),
				ifft(thisdata)*double(self.filterpoints)*self.dx)
