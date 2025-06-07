import numpy as np
from conftest import load_module

# load dependencies in order
load_module("general_functions")
ns = load_module("ndshape")


def test_basic_init_and_properties():
    s = ns.ndshape_base([3, 4], ["x", "y"])
    assert s.shape == [3, 4]
    assert s.dimlabels == ["x", "y"]
    assert s.axn("y") == 1
    assert s["x"] == 3


def test_list_of_pairs_init():
    s = ns.ndshape_base([("x", 2), ("y", 5)])
    assert s.shape == [2, 5]
    assert s.dimlabels == ["x", "y"]


def test_nddata_like_init_and_methods():
    class Dummy:
        def __init__(self):
            self.data = np.zeros((2, 3))
            self.dimlabels = ["a", "b"]

    s = ns.ndshape_base(Dummy())
    assert s.shape == [2, 3]
    assert s.dimlabels == ["a", "b"]

    s.rename("b", "c")
    assert s.dimlabels == ["a", "c"]
    s["a"] = 5
    assert s.shape[0] == 5
    s.pop("c")
    assert s.shape == [5]
    assert s.dimlabels == ["a"]


def test_iteration():
    s = ns.ndshape_base([3, 4], ["x", "y"])
    items = list(s)
    assert items == [("x", 3), ("y", 4)]
