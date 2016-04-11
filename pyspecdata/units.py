from numpy import *
import re
from itertools import chain

class units (object):
    r"""Each instance of this object stores a numerical representation of a single set of units, and there are routines to set units by
    (*i.e.* :func:`parsing <pyspecdata.units.parse>`) strings to units
    and routines to convert units to an
    :func:`str <pyspecdata.units.part>` representation.

    At its core, the units are represented by three numpy structured arrays per instance:
    one for the coefficients in terms of base units,
    one for the order of magnitude used to determine any prefixes added to the base units,
    and one for any scaling factors needed to convert to base units.

    An array of these structured arrays can be converted into a row
    vector with ``.view((float16,len(base_dtype.names)))``.
    "Base Units" here are the same as SI base units **except** that it uses g instead of kg (so we can do  the prefixes correctly), and we have added rad.
    """
    def __init__(self):
        self.base_dtype = dtype(map(lambda x: (x, float16), ['m', 'g', 's', 'A',
            'K', 'mol', 'cd', 'rad']))
        self.oom_names =   ['T' , 'G' , 'M' , 'k' , 'c' , 'm' , u'\u03bc' , 'n' , 'p']
        oom_values = r_[12 , 9   , 6   , 3   , -2  , -3  , -6     , -9  , -12]
        self.oom_dict = dict(zip(self.oom_names,oom_values))
    def parse(self,in_str):
        u"""Take `in_str` and parse it as a unit or series of units, and set the units associated with the current instance to the result.

        Addition, subtraction, and parentheses are not allowed, and we define a
        non-standard order of operations, as follows:

        #. ``\\\\mu`` (including a trailing space) and ``u`` are converted to the utf-8 equivalent (\u03bc)
        #. ``...^-0.125`` any number is assumed to be part of the exponent, and only numbers are allowed.
        #. ``*`` multiplication
        #. a space also represents multiplication
        #. ``.../...`` comes **after** all other operations, as is typical for single-line text
        #. ``sqrt(...)`` comes "last" in the sense that we take care of everything both inside and outside the sqrt first, before applying the sqrt.
        
        At this point, I use split to break up according to the order of operations, assign powers to each, and determine the prefix.
        However, I'm currently only using base units, and I will eventually want to use derived units as well.
        """
        in_str = unicode(in_str)
        mu_re = re.compile(r'\bu|\\mu ')
        in_str = mu_re.sub(u'\u03bc',in_str)
        working_list = [in_str]
        unit_powers = [1]
        for op_num,this_regexp in enumerate([r'\bsqrt\((.*)\)',
            r'(?:\*| *)',
            r'/',
            r'(?:\^|\*\*)([\-\.0-9]+)', ]):# the first in the order of operations are actually inside
            last_len = 0
            while len(working_list) > last_len:
                last_len = len(working_list)
                working_list = map(lambda x: re.split(this_regexp,x,maxsplit = 1),
                        working_list)
                print "broken according to",this_regexp
                new_unit_powers = [[unit_powers[j]]*len(working_list[j]) for j
                        in range(len(working_list))]
                for j,v in enumerate(new_unit_powers):
                    if len(v) > 1:# then, it split this element
                        if op_num == 0:
                            new_unit_powers[j][1] = float(new_unit_powers[j][1]) / 2.0 # middle element is inside parens, and goes to half power
                        elif op_num == 1:
                            pass # for multiplication, we don't need to do anything
                        elif op_num == 2:
                            print "unit powers before messing:",new_unit_powers
                            new_unit_powers[j][1] *= -1 # second half of this expression becomes negative
                            # {{{ all following elements are also negative
                            if j+1 < len(new_unit_powers):
                                print "setting all subsequent to -1"
                                for k in range(j+1,len(new_unit_powers)):
                                    new_unit_powers[k] = map(lambda x: -1*x,new_unit_powers[k])
                            # }}}
                        elif op_num == 3:
                            print "unit powers before messing:",new_unit_powers
                            print "working list before messing",working_list
                            new_unit_powers[j][0] *= float(working_list[j][1])
                            working_list[j].pop(1)
                            new_unit_powers[j].pop(1)
                            print "unit power after messing:",new_unit_powers
                            print "working list after messing:",working_list
                        working_list[j] = filter(lambda x: len(x)>0,working_list[j])
                        new_unit_powers[j] = [new_unit_powers[j][k] for k in
                                range(len(working_list[j])) if
                                len(working_list[j][k]) > 0]
                unit_powers = new_unit_powers
                print working_list
                print unit_powers
                working_list = list(chain(*working_list))
                unit_powers = list(chain(*unit_powers))
        print "Final result:"
        print unit_powers
        oom_regexp = r'\b('+'|'.join(self.oom_names)+'){0,1}'+'('+'|'.join(self.base_dtype.names)+r')\b'
        print "oom_regexp",oom_regexp.encode('ascii','replace')
        oom_regexp = re.compile(oom_regexp)
        print working_list
        print "breaks into:"
        for n,j in enumerate(working_list):
            m = oom_regexp.match(j)
            if m:
                prefix,unit = m.groups()
                un_np = zeros(1,dtype = self.base_dtype)
                oom_np = zeros(1,dtype = self.base_dtype)
                un_np[unit] = unit_powers[n]
                if prefix is not None:
                    oom_np[unit] = self.oom_dict[prefix]
                print r_[un_np,oom_np] 
            else:
                raise ValueError(' '.join(["couldn't match:",j,"to a unit"]))
        return self
    def str(self,number):
        r"""Give a string that prints `number`, which has the units
        given by the current instance of the class.
        Choose the simplest possible expression for the units.

        When printing, we have a matrix that give all our "representation" units,
        and we use a pseudoinverse to give us the simplest possible expression of our units.
        (This is assuming that all derived units are defined in terms of powers
        greater than or equal to 1 of the base units, because powers of
        magnitude less than 1 would complicate things by allowing us to reduce
        the norm by spreading across several derived units -- in that case, we
        might want to add a threshold before taking
        the pinv.)
        """
        return

