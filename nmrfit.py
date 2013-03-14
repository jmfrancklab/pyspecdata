from matlablike import *
import sympy

#{{{ different types of fit classes
class t2curve(fitdata):
    def guess(self):
        'this overrides the pinvr guess'
        return [1.0,1.0,-2.0,0.1]
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        for j in range(0,self.multiplicity):
            if j == 0:
                retval = p[2*j]*exp(-x/p[2*j+1])
            else:
                retval = retval + p[2*j]*exp(-x/p[2*j+1])
        return retval
    def fitfunc_raw_symb(self,p,x):
        '''if I'm using a named function, I have to define separately in terms of sympy rather than numpy functions'''
        for j in range(0,self.multiplicity):
            if j == 0:
                retval = p[2*j]*sympy.exp(-x/p[2*j+1])
            else:
                retval = retval + p[2*j]*sympy.exp(-x/p[2*j+1])
        return retval
    def __init__(self,*args,**kwargs):
        if 'multiplicity' in kwargs:
            multiplicity = kwargs.pop('multiplicity')
        else:
            multiplicity = 1
        fitdata.__init__(self,*args,**kwargs)
        alphaletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
        self.symbol_list = [k for j in range(0,multiplicity)
                for k in r'M_%s(0)'%alphaletter[j],r'T_2%s'%alphaletter[j]]
        #self.starting_guesses = map(double,[array([1,0.1]*multiplicity)])
        self.starting_guesses = map(double,[array([1,2,-1,0.1])])
        self.guess_lb = array([0.,-inf]*multiplicity)
        self.guess_ub = array([+inf,+inf]*multiplicity)
        print 'symbol list is:',self.symbol_list
        print 'starting guesses is:',self.starting_guesses
        print 'guess lb is:',self.guess_lb
        print 'guess ub is:',self.guess_ub
        self.multiplicity = multiplicity
        self.gen_symbolic(r'M(t)')
        return
class t1curve(fitdata):
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        return p[0]+(p[1]-p[0])*exp(-x/p[2])
    def fitfunc_raw_symb(self,p,x):
        '''if I'm using a named function, I have to define separately in terms of sympy rather than numpy functions'''
        return p[0]+(p[1]-p[0])*sympy.exp(-x/p[2])
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of T1, we want to return ln(y-M(\infty)) = ln(M(0)-M(\infty)) - t/T_1
        '''
        #print 'DEBUG: y is',y
        #print 'DEBUG: M(\infty) is',self.output(r'M(\infty)')
        temp = self.output(r'M(\infty)')-y # the argument for log
        #print 'DEBUG: temp is',temp
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        rety = log(temp)
        if yerr != None:
            reterr = yerr/abs(temp)
        mask = isfinite(rety)
        retx = x # for instance, in emax, this is not just x
        xname = self.fit_axis # same as the fit axis
        yname = r'$ln(M(\infty)-M(t))$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([self.fit_axis],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
    def linerror(self,x,y):
        '''propagate the error for linfunc
        '''
        rety = log(y-self.output(r'M(\infty)'))
        mask = isfinite(rety)
        x_axis_of_linear_plot = x # for instance, in emax, this is not just x
        retval = nddata(rety,
                [size(rety),1],
                [self.fit_axis,r'$ln(M(t)-M(\infty))$'])
        retval.labels([self.fit_axis],
                [x_axis_of_linear_plot.copy()])
        return retval
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        self.symbol_list = [r'M(\infty)',r'M(0)',r'T_1'] # note that it must notbe possible to find part of one of the later strings by searching for one of the earlier strings
        self.starting_guesses = map(double,[r_[3,-3,0.1],r_[1,1,1],r_[0,0,1],r_[-100,100,0.03],r_[0.001,0.001,0.001],r_[1,-1,4.0]])
        self.guess_lb = r_[-inf,-inf,1e-4]
        self.guess_ub = r_[+inf,+inf,20.]
        self.gen_symbolic(r'M(t)')
        return
class emax_legacy(fitdata):
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        newdata = self.copy()
        newdata.sort(self.fit_axis)
        power = newdata.getaxis(self.fit_axis)
        integral = newdata.data
        largest = len(power)
        initial_slope = (integral[largest/4]-integral[0])/(power[largest/4]-power[0])
        approx_emax = integral[-1]
        #print 'DEBUG: guessed initial slope',initial_slope,'approx emax',approx_emax
        return [1,-initial_slope,-initial_slope/(1-approx_emax)]
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        return (p[0]-(p[1]*x/(1.+p[2]*x)))
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of E_max, we want 1/(1-E) = 
        '''
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        #print "linfunc passed x=",x,"and y=",y
        rety = 1./(1.-y)
        if yerr != None:
            reterr = yerr/((1.-y)**2)
        mask = isfinite(rety)
        retx = 1./x # for instance, in emax, this is not just x
        xname = r'1 / '+self.fit_axis # same as the fit axis
        yname = r'$\frac{1}{1-E(p)}$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([xname],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        self.function_string = r'$E(p)=c_0-Ap/(1+Bp)$'
        self.symbol_list = [r'c_0',r'A',r'B'] # note that it must notbe possible to find part of one of the later strings by searching for one of the earlier strings
        return
class emax(fitdata):
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        newdata = self.copy()
        newdata.sort(self.fit_axis)
        power = newdata.getaxis(self.fit_axis)
        integral = newdata.data
        largest = len(power)
        initial_slope = (integral[largest/4]-integral[0])/(power[largest/4]-power[0])
        approx_emax = integral[-1]/integral[0]
        return [approx_emax,integral[0],-initial_slope] #guess [r'E_{max}',r'v',r'A']
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        #self.function_string = r'$E(p)=v-Apv(1-E_{max})/(1-E_{max}+Ap)$'
        #self.symbol_list = [r'E_{max}',r'v',r'A'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        return (p[1]-(p[2]*x*p[1]*(1.-p[0])/(1.-p[0]+p[2]*x)))
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of E_max, we want 1/(1-E) = 
        '''
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        #print "linfunc passed x=",x,"and y=",y
        rety = 1./(1.-y)
        if yerr != None:
            reterr = yerr/((1.-y)**2)
        mask = isfinite(rety)
        retx = 1./x # for instance, in emax, this is not just x
        xname = r'1 / '+self.fit_axis # same as the fit axis
        yname = r'$\frac{1}{1-E(p)}$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([xname],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        #self.function_string = r'$E(p)=v-Apv(1-E_{max})/(1-E_{max}+Ap)$'
        self.symbol_list = [r'E_{max}',r'v',r'A'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        self.gen_symbolic(r'E(p)')
        return
class one_minus_emax(fitdata):
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        self.symbol_list = [r'E_{max}',r'A'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        self.gen_symbolic(r'1-E(p)')
        return
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        return ((x*p[1]*(1.-p[0])/(1.-p[0]+p[1]*x)))
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        newdata = self.copy()
        newdata.sort(self.fit_axis)
        power = newdata.getaxis(self.fit_axis)
        one_minus_E = newdata.data
        largest = len(power)
        initial_slope = (one_minus_E[largest/4]-one_minus_E[0])/(one_minus_E[largest/4]-one_minus_E[0])
        approx_emax = (1.0-one_minus_E[-1])
        return [approx_emax,initial_slope] #guess [r'E_{max}',r'A']
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of E_max, we want 1/(1-E) = 
        '''
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        #print "linfunc passed x=",x,"and y=",y
        rety = 1./(y)
        if yerr != None:
            reterr = yerr/(y**2) # check this later
        mask = isfinite(rety)
        retx = 1./x # for instance, in emax, this is not just x
        xname = r'1 / '+self.fit_axis # same as the fit axis
        yname = r'$\frac{1}{1-E(p)}$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([xname],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
class xismax(fitdata):
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        self.function_string = r'$\xi s(p)=B p \xi s_{max}/(\xi s_{max} + B p$'
        self.symbol_list = [r'B',r'\xi s_{max}'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        return
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        #self.function_string = r'$\xi s(p)=B p \xi s_{max}/(\xi s_{max} + B p$'
        #self.symbol_list = [r'B',r'\xi s_{max}'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        return (p[0]*x*p[1])/(p[1]+p[0]*x)
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        newdata = self.copy()
        newdata.sort(self.fit_axis)
        power = newdata.getaxis(self.fit_axis)
        integral = newdata.data
        largest = len(power)
        initial_slope = (integral[largest/4]-integral[0])/(power[largest/4]-power[0])
        #{{{ use the indexing functions to more legibly set the return values
        retval = zeros(len(self.symbol_list))
        retval[self._pn(r'\xi s_{max}')] = integral[-1]
        retval[self._pn(r'B')] = initial_slope
        return retval
        #}}}
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of E_max, we want 1/(1-E) = 
        '''
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        #print "linfunc passed x=",x,"and y=",y
        rety = 1./(1.-y)
        if yerr != None:
            reterr = yerr/((1.-y)**2)
        mask = isfinite(rety)
        retx = 1./x # for instance, in emax, this is not just x
        xname = r'1 / '+self.fit_axis # same as the fit axis
        yname = r'$\frac{1}{1-E(p)}$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([xname],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
class smax(fitdata):
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        self.symbol_list = [r'x',r'\xi'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        self.gen_symbolic(r'xismax(C)')
        return
    def fitfunc_raw(self,p,C):
        '''just the actual fit function to return the array y as a function of p and x'''
        x = p[0]*C
        xi = p[1]
        snrmax = x/(x+3)
        return xi*(1./3.+2./3.*snrmax)
    def set_guess(self,**kwargs):
        newdata = self.copy()
        maxC = lambda x: x == self.getaxis(self.fit_axis).max()
        self.guesses = {r'xi':newdata[self.fit_axis,maxC].data[-1],
            'x':1.}
        self.guesses.update(kwargs)
        return
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        if not hasattr(self,'guesses'):
            self.set_guess()
        return [self.guesses[x] for x in ['x','xi']] #guess [r'E_{max}',r'A']
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of E_max, we want 1/(1-E) = 
        '''
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        #print "linfunc passed x=",x,"and y=",y
        rety = 1./(y)
        if yerr != None:
            reterr = yerr/(y**2) # check this later
        mask = isfinite(rety)
        retx = 1./x # for instance, in emax, this is not just x
        xname = r'1 / '+self.fit_axis # same as the fit axis
        yname = r'$\frac{1}{1-E(p)}$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([xname],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
#{{{ k_sigma s(p)
class ksp (fitdata):
    def __init__(self,*args,**kwargs):
        fitdata.__init__(self,*args,**kwargs)
        #self.symbol_list = [r'ks_{max}',r'p_{half}']
        self.symbol_list = [r'ksmax',r'phalf']
        self.starting_guesses = map(double,[r_[60,0.5],r_[20,0.05],r_[1.,1.]])
        self.guess_lb = r_[0.1,1e-4]
        self.guess_ub = r_[500.,30.]
        self.zero_guess = r_[True,True]
        self.gen_symbolic(r'k_{\sigma} s(p)')
        return
    def fitfunc_raw(self,p,x):
        return p[0]/(p[1]+x)*x
    def fitfunc_raw_symb(self,p,x):
        try:
            return p[0]/(p[1]+x)*x
        except:
            raise CustomError('types of p',p,map(type,p),'type of x',type(x),'x=',x)
#}}}
#}}}
