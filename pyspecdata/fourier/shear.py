from ..general_functions import *
if not inside_sphinx():
    from pylab import r_

def extend_for_shear(self,altered_axis,propto_axis,skew_amount,verbose = False):
    "this is propto_axis helper function for `.fourier.shear`"
    #{{{ in the time domain, altered_axis is the one that's altered (and
    #       needs to be extended), while the shearing is proportional to
    #       -by_amount*propto_axis
    if verbose: print("extending to account for the shear along ",altered_axis,"by",skew_amount,"which gives lesser and greater expansion amounts of", end=' ')
    try:
        shear_displacement = skew_amount * self.getaxis(propto_axis
                )[r_[0,-1]]
    except Exception as e:
        if self.getaxis(propto_axis) is None:
            raise RuntimeError("propto_axis ("+propto_axis+") is not set")
        else:
            raise e
    shear_displacement = sort(shear_displacement) # this gives the lesser
    #       and greater, respectively (i.e. le, gt -- not smaller/bigger),
    #       of the two shear displacements, so that I know how to extend.
    #       I need to manually sort, because I don't know if skew_amount is
    #       negative or positive.
    #{{{ actually extend: leave alone if zero.
    if verbose: print(" and ".join(map(str,shear_displacement)))
    for j in [0,-1]:
        if shear_displacement[j] != 0.:
            if verbose: print(' '.join(map(str,("preparing to extend altered_axis (",
                    altered_axis,")",
                    self.getaxis(altered_axis)[r_[0,-1]], "to",
                    self.getaxis(altered_axis)[j] + shear_displacement[j],
                    "along the", ['lesser','greater'][j],"side of",
                    "altered_axis (",self.getaxis(altered_axis)[j],") by adding",
                    shear_displacement[j],"to it"))))
            if ((self.getaxis(altered_axis)[j]>0) ^
                    (shear_displacement[j]>0)):
                if verbose: print("skipping this extension, since the shear seems to be trying to push the axis back towards zero")
            else:
                self.extend(altered_axis, self.getaxis(altered_axis)[j] +
                        shear_displacement[j])# match greater with greater and
                #       lesser with lesser (doesn't matter to me which side of
                #       propto_axis that they came from)
    #}}}
    #}}}
    return self

def shear(self,altered_axis,propto_axis,by_amount,zero_fill = False,start_in_conj = False):
    'the fourier shear method -- see .shear() documentation'
    #{{{ see if it's in the frequency or time domain
    if self.get_ft_prop(altered_axis) and self.get_ft_prop(propto_axis):
        frequency_domain = True
    elif not self.get_ft_prop(altered_axis) and not self.get_ft_prop(propto_axis):
        frequency_domain = False
    else:
        raise ValueError("In order to shear, both dimensions must be in the same (time vs. frequency) domain.  Currently, they are {:s} for {:s} and {:s} for {:s}.".format(
            self.get_ft_prop(altered_axis),altered_axis,
            self.get_ft_prop(propto_axis),propto_axis
            ))
    if start_in_conj:# then I want to actually specify the shear in the other
        #               domain than the one I start in, and just skip the
        #               beginning
        frequency_domain = not frequency_domain
    #}}}
    if frequency_domain:
        print("entering time domain for",altered_axis)
        if not start_in_conj:
            if zero_fill:
                self.extend_for_shear(altered_axis,propto_axis,by_amount)
            self.ift(altered_axis)
            self.ift(propto_axis) # before expansion
        else:
            if zero_fill:
                raise ValueError("I can't zero fill  because you chose to start in the conjugate dimension")
        print("conjugate domain extension:")
        self.extend_for_shear(propto_axis,altered_axis,-by_amount) # in
        #       the time domain, propto_axis is the one that's altered
        #       (and needs to be extended), while the shearing is
        #       proportional to -by_amount*altered_axis
        self.ft(propto_axis) # after expansion
        print("applying phase shift")
        phaseshift = self.fromaxis([altered_axis,propto_axis],
                lambda x,y: exp(2j*pi*by_amount*x*y))
        self.data *= phaseshift.data
        print("back to frequency domain")
        self.ft(altered_axis)
    else:
        print("entering time domain for",altered_axis)
        if not start_in_conj:
            if zero_fill:
                self.extend_for_shear(altered_axis,propto_axis,by_amount)
            self.ft(altered_axis)
            self.ft(propto_axis) # before expansion
        else:
            if zero_fill:
                raise ValueError("I can't zero fill  because you chose to start in the conjugate dimension")
        print("conjugate domain extension:")
        self.extend_for_shear(propto_axis,altered_axis,-by_amount) # in
        #       the time domain, propto_axis is the one that's altered
        #       (and needs to be extended), while the shearing is
        #       proportional to -by_amount*altered_axis
        self.ift(propto_axis) # after expansion
        print("applying phase shift")
        phaseshift = self.fromaxis([altered_axis,propto_axis],
                lambda x,y: exp(-2j*pi*by_amount*x*y))
        self.data *= phaseshift.data
        print("back to frequency domain")
        self.ift(altered_axis)
    return self
