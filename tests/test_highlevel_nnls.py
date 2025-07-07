import numpy as np
from numpy import linspace, r_, exp, sqrt
from numpy.random import seed
from conftest import load_module
import importlib.util
import sys

# ensure submodules are reloaded fresh for this test
for name in ["pyspecdata.nnls", "pyspecdata.matrix_math.nnls", "pyspecdata.core"]:
    sys.modules.pop(name, None)

load_module("nnls")
load_module("matrix_math.nnls")

load_module("general_functions")
core = load_module("core")
nddata = core.nddata
init_logging = load_module("general_functions").init_logging


def test_highlevel_nnls():
    seed(1234)
    init_logging("debug")
    vd_list = nddata(linspace(5e-4, 10, 25), "vd")
    t1_name = r"$\\log(T_1)$"
    logT1 = nddata(r_[-4:2:100j], t1_name)

    def Gaussian_1d(axis, mu1, sigma1):
        this_G = exp(-((axis - mu1) ** 2) / 2 / sigma1**2)
        return this_G

    true_F = Gaussian_1d(logT1.C.run(lambda x: 10 ** (x)), 6, 0.3)

    K = 1.0 - 2 * exp(-vd_list / 10 ** (logT1))
    K.reorder("vd")

    M = K @ true_F
    M.add_noise(0.2)
    M /= 0.2

    solution = M.C.nnls(
        "vd", logT1, lambda x, y: 1 - 2 * exp(-x / 10 ** (y)), l="BRD"
    )
    solution_confirm = M.C.nnls(
        "vd",
        logT1,
        lambda x, y: 1 - 2 * exp(-x / 10 ** (y)),
        l=sqrt(solution.get_prop("opt_alpha")),
    )

    diff = np.linalg.norm(solution.data - solution_confirm.data)
    assert diff < 0.2

    axis_T1 = 10 ** logT1.data
    max_T1 = axis_T1[np.argmax(solution.data)]
    avg_T1 = np.average(axis_T1, weights=solution.data)
    assert abs(max_T1 - 6) < 1
    assert abs(avg_T1 - 6) < 5
