"""provides object fl_dummy that can take the place of a figure list, but doesn't do anything.

`fl=fl_dummy` is a better default than `fl=None`, because it doesn't require any conditionals
"""
class fl_dummy_class (object):
    def plot(*args):
        pass
    def next(*args):
        pass
    def push_marker(*args):
        pass
    def pop_marker(*args):
        pass
fl_dummy = fl_dummy_class

