r'''The :class:`ndshape` class allows you to allocate arrays and determine the shape of existing arrays.'''
import numpy as np
from .general_functions import *

class ndshape_base ():
    r'''The base ndshape class, which doesn't include an allocation method.'''
    def __init__(self,*args):
        """Create an nddata in one of various ways.  Values are
        always copies of the values they were created from
        (changing the shape of the initialization parameters will
        not change/update the new nddata that's created).

        >>> nddata_instance = ndshape(shapes,dimlabels)

        >>> nddata_instance = ndshape(list_of_pairs)

        or 

        >>> nddata_instance = ndshape(nddata_instance)

        Parameters
        ==========
        
        shapes : list of int

            the sizes of the dimensions, in order

        dimlabels : list of str
            
            the names of the dimensions, in order

        list_of_pairs : list of tuples

            zip(dimlabels,shapes)

        nddata_instance : nddata
            
            determine the shape of the nddata,
            and return it
        """
        self.zero_dimensional = False
        if len(args) == 2:
            self.shape = list(args[0])
            self.dimlabels = args[1]
        if len(args) == 1: #assum that it's an nddata object
            if hasattr(args[0],'dimlabels') and hasattr(args[0],'data'):# test that it's nddata without using any of the functions in core.py
                self.shape = list(args[0].data.shape)
                self.dimlabels = list(args[0].dimlabels)
                if len(self.shape) == 0 and len(self.dimlabels) == 0:
                    self.zero_dimensional = True
                    return
            elif isinstance(args[0], list):
                self.dimlabels, self.shape = list(map(list,list(zip(*args[0]))))
            else:
                raise ValueError('If you pass a single argument, it must be an nddata')
        return
    def __setitem__(self,reference,setto):
        self.shape[self.axn(reference)] = setto
        return
    def axn(self,axis):
        r'return the number for the axis with the name "axis"'
        try:
            return self.dimlabels.index(axis)
        except:
            raise ValueError(' '.join(map(repr,['there is no axis named',axis,'all axes are named',self.dimlabels])))
    def copy(self):
        try:
            return deepcopy(self)
        except:
            raise RuntimeError('Some type of error trying to run deepcopy on'+repr(self))
    def matchdims(self,arg):
        r'returns shape with [not in self, len 1] + [overlapping dims between arg + self] + [not in arg] --> this is better accomplished by using sets as I do in the matchdims below'
        for k in set(self.dimlabels) & set(arg.dimlabels):
            a = arg.shape[arg.axn(k)]
            b = self.shape[self.axn(k)]
            if a != b:
                raise CustomError('the',k,'dimension is not the same for self',self,'and arg',arg)
        if hasattr(arg,'dimlabels') and hasattr(arg,'data'):# test that it's nddata without using any of the functions in core.py
            arg = ndshape(arg)
        #{{{ add extra 1-len dims
        addeddims = set(self.dimlabels) ^ set(arg.dimlabels) & set(arg.dimlabels)
        self.dimlabels = list(addeddims) + self.dimlabels
        self.shape = [1] * len(addeddims) + list(self.shape)
        #}}}
        return self
    def __radd__(self,arg):
        'take list of shape,dimlabels'
        shape = arg[0]
        dimlabels = arg[1]
        if isinstance(shape, str):
            shape,dimlabels = dimlabels,shape
        if np.isscalar(self.shape):
            self.shape = [self.shape]
        if np.isscalar(self.dimlabels):
            self.dimlabels = [self.dimlabels]
        if np.isscalar(shape):
            shape = [shape]
        if np.isscalar(dimlabels):
            dimlabels = [dimlabels]
        self.shape = shape + self.shape
        self.dimlabels = dimlabels + self.dimlabels
        return self
    def __add__(self,arg):
        '''take list of shape,dimlabels
        this is the correct function, until I can fix my back-references for add, which does it backwards'''
        shape = arg[0]
        dimlabels = arg[1]
        if isinstance(shape, str):
            shape,dimlabels = dimlabels,shape
        if np.isscalar(self.shape):
            self.shape = [self.shape]
        if np.isscalar(self.dimlabels):
            self.dimlabels = [self.dimlabels]
        if np.isscalar(shape):
            shape = [shape]
        if np.isscalar(dimlabels):
            dimlabels = [dimlabels]
        self.shape = self.shape + shape
        self.dimlabels = self.dimlabels + dimlabels
        return self
    def __repr__(self): #how it responds to print
        return list(zip(self.shape,self.dimlabels)).__repr__()
    def __iter__(self):
        self._index = 0
        return self
    def max(self):
        idx = np.argmax(self.shape)
        return self.dimlabels[idx]
    def min(self):
        idx = np.argmin(self.shape)
        return self.dimlabels[idx]
    def __next__(self):
        if self._index < len(self.shape):
            k,v = self.dimlabels[self._index],self.shape[self._index]
            self._index += 1
            return k,v
        else:
            raise StopIteration
    def __getitem__(self,args):
        try:
            mydict = dict(list(zip(self.dimlabels,self.shape)))
        except Exception as e:
            raise ValueError(strm("either dimlabels=",self.dimlabels,"or shape",
                self.shape,"not in the correct format") + explain_error(e))
        try:
            return mydict[args]
        except:
            raise ValueError(strm("one or more of the dimensions named",args,"do not exist in",self.dimlabels))
    def rename(self,before,after):
        r'rename a dimension'
        thisindex = self.axn(before)
        self.dimlabels[thisindex] = after
        return self
    def pop(self,label):
        r'remove a dimension'
        thisindex = self.axn(label)
        self.dimlabels.pop(thisindex)
        self.shape.pop(thisindex)
        return self
