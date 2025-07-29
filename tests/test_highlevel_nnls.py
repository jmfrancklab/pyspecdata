import numpy as np
from numpy import linspace, r_, exp, sqrt
from numpy.random import seed
from conftest import load_module
from pyspecdata.matrix_math import venk_nnls
core = load_module("core")
nddata = core.nddata
init_logging = load_module("general_functions").init_logging

# ensure submodules are reloaded fresh for this test
mu1 = 0.5
sigma1 = 0.3


def test_highlevel_nnls():
    seed(1234)
    init_logging("debug")
    vd_list = nddata(linspace(5e-4, 10, 25), "vd")
    t1_name = r"$\\log(T_1)$"
    logT1 = nddata(r_[-4:2:100j], t1_name)

    true_F = (
        1
        / np.sqrt(2 * np.pi * sigma1**2)
        * exp(-((logT1 - mu1) ** 2) / 2 / sigma1**2)
    )

    K = 1.0 - 2 * exp(-vd_list / 10 ** (logT1))
    K.reorder("vd")

    M = K @ true_F
    M.add_noise(0.2)
    M /= 0.2

    solution = M.C.nnls(
        "vd", logT1, lambda x, y: 1 - 2 * exp(-x / 10 ** (y)), l="BRD"
    )
    solution_stackcalc = M.C.nnls(
        "vd",
        logT1,
        lambda x, y: 1 - 2 * exp(-x / 10 ** (y)),
        l=sqrt(solution.get_prop("opt_alpha")),
    )
    diff = np.linalg.norm(solution.data - solution_stackcalc.data)
    assert diff < 0.01 * np.linalg.norm(solution.data)

    solution_venk = M.C.nnls(
        "vd",
        logT1,
        lambda x, y: 1 - 2 * exp(-x / 10 ** (y)),
        l=sqrt(solution.get_prop("opt_alpha")),
        method=venk_nnls,
    )
    diff = np.linalg.norm(solution_stackcalc.data - solution_venk.data)
    assert diff < 0.054 * np.linalg.norm(
        solution_stackcalc.data
    ), f"venk diff is {diff/np.linalg.norm(solution_stackcalc.data)}"

    max_log_T1 = logT1.data[np.argmax(solution.data)]
    avg_log_T1 = np.average(logT1.data, weights=solution.data)
    betteravg_log_T1 = np.average(
        logT1[t1_name:(-2, None)].data,
        weights=solution[t1_name:(-2, None)].data,
    )
    assert (
        abs(max_log_T1 - mu1) < 0.1
    ), f"max_log_T1 is {max_log_T1} while mu1 is {mu1}"
    assert (
        abs(avg_log_T1 - mu1) < 0.2
    ), f"avg_log_T1 is {avg_log_T1} while mu1 is {mu1}"
    assert (
        abs(betteravg_log_T1 - mu1) < 0.2
    ), f"betteravg_log_T1 is {betteravg_log_T1} while mu1 is {mu1}"
