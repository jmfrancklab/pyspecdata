import numpy as np
import pytest
from conftest import load_module

gf = load_module("general_functions")


def test_dp_formatting():
    assert gf.dp(3.14159, 2) == "3.14"
    assert gf.dp(1234.567, 2) == "1.23\\times 10^{3}"
    assert gf.dp(1234.567, 2, scientific=True) == "1.23\\times 10^{3}"


def test_emptytest():
    assert gf.emptytest([])
    assert gf.emptytest(np.array([]))
    assert gf.emptytest(None)
    assert not gf.emptytest([1])


def test_autostringconvert():
    assert gf.autostringconvert("abc") == "abc"
    obj = 5
    assert gf.autostringconvert(obj) is obj


def test_process_kwargs_defaults():
    kwargs = {}
    a, b = gf.process_kwargs([("a", 1), ("b", 2)], kwargs)
    assert (a, b) == (1, 2)
    assert kwargs == {}


def test_process_kwargs_override():
    kwargs = {"a": 10}
    a, b = gf.process_kwargs([("a", 1), ("b", 2)], kwargs)
    assert (a, b) == (10, 2)
    assert kwargs == {}


def test_process_kwargs_passthrough():
    kwargs = {"c": 3}
    result = gf.process_kwargs([("a", 1)], kwargs, pass_through=True)
    assert result == 1
    assert kwargs == {"c": 3}


def test_process_kwargs_unknown_error():
    kwargs = {"c": 3}
    with pytest.raises(ValueError):
        gf.process_kwargs([("a", 1)], kwargs)
