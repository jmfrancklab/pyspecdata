from .core import nddata

# from here: https://stackoverflow.com/questions/7388258/replace-property-for-perfomance-gain?noredirect=1&lq=1
class CachedAttribute(object):    
    '''Computes attribute value and caches it in the instance.
    From the Python Cookbook (Denis Otkidach)
    This decorator allows you to create a property which can be computed once and
    accessed many times. Sort of like memoization.
    '''
    def __init__(self, method, name=None):
        # record the unbound-method and the name
        self.method = method
        self.name = name or method.__name__
        self.__doc__ = method.__doc__
    def __get__(self, inst, cls):
        # self: <__main__.cache object at 0xb781340c>
        # inst: <__main__.Foo object at 0xb781348c>
        # cls: <class '__main__.Foo'>       
        if inst is None:
            # instance attribute accessed on class, return self
            # You get here if you write `Foo.bar`
            return self
        # compute, cache and return the instance's attribute value
        result = self.method(inst)
        # setattr redefines the instance's attribute so this doesn't get called again
        setattr(inst, self.name, result)
        return result

class nddata_placeholder(nddata):
    r"""When you open a file with many datasets, and only want to load one, this is used to store a placeholder to the different datasets.
    Rather than passing the initialization any type of data, you simply pass it one argument: the "data hook":
    a function that accepts a single argument (self), alters the dimlabels and axes as needed, and returns an ndarray containing the data.

    You can attach only the necessary information to the placeholder, and then
    load it into a dictionary that can be explored by the user as they look for
    the right dataset.
    """
    def __init__(self,_data_hook):
        self._data_hook = _data_hook
    @CachedAttribute
    def data(self):
        return self._data_hook()
