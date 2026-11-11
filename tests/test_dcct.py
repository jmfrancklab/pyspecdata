from collections import OrderedDict
import sys
import types
import matplotlib.pyplot as plt
import numpy as np
import pytest
import sympy as s
from conftest import load_module

plt.switch_backend("Agg")

sys.modules.setdefault("_nnls", types.ModuleType("_nnls"))

core = load_module("core", use_real_pint=True)
pkg = sys.modules["pyspecdata"]
pkg.nddata = core.nddata
pkg.ndshape = core.ndshape
dcct_mod = load_module("plot_funcs.DCCT_function", use_real_pint=True)
fake_data_mod = load_module("generate_fake_data", use_real_pint=True)
nddata = core.nddata


def test_getitem_preserves_none_axis_coords_after_chunk():
    idx = nddata(np.r_[0:8], [-1], ["smooshed"])
    idx.chunk("smooshed", ["ph1", "ph2"], [2, 4])
    sliced = idx["ph1", 0]
    assert sliced.dimlabels == ["ph2"]
    assert sliced.getaxis("ph2") is None
    np.testing.assert_array_equal(sliced.data, np.r_[0:4])


def phase_cycled_data():
    data = np.arange(2 * 4 * 3 * 16, dtype=float).reshape(2, 4, 3, 16)
    data = data + 1j * data[::-1]
    d = nddata(data, [2, 4, 3, 16], ["ph1", "ph2", "vd", "t2"])
    d.set_axis("ph1", np.r_[0, 2] / 4)
    d.set_axis("ph2", np.r_[0:4] / 4)
    d.set_axis("vd", np.linspace(0.0, 1.0, 3))
    d.set_axis("t2", np.linspace(0.0, 15e-3, 16))
    d.set_units("t2", "s")
    return d.ft("t2", shift=True).ft(["ph1", "ph2"])


def test_dcct_smoke_with_phase_dimensions():
    fig = plt.figure()
    ax_list, _ = dcct_mod.DCCT(phase_cycled_data(), fig=fig)
    assert len(ax_list) == 8
    # the brackets label the (aliased) coherence changes
    all_text = [j.get_text() for j in fig.texts]
    assert all_text.count("-1, +1") == 1  # ph1
    assert all_text.count("-2, +2") == 2  # ph2, once for each ph1
    # each image is centered on its coordinates, extending half a pixel
    # past the first and last
    left, right, bottom, top = ax_list[0].images[0].get_extent()
    np.testing.assert_allclose([bottom, top], [-0.25, 1.25])
    plt.close("all")


def test_dcct_inner_phase_dimension_ticks():
    d = phase_cycled_data().reorder(["vd", "ph1", "ph2", "t2"])
    fig = plt.figure()
    ax_list, _ = dcct_mod.DCCT(d, fig=fig)
    assert len(ax_list) == 6
    assert [j.get_text() for j in ax_list[0].get_yticklabels()] == [
        "0",
        "+1",
        "-2, +2",
        "-1",
    ]
    plt.close("all")


def test_dcct_arrow_width_and_cmap_check():
    d = phase_cycled_data()
    fig = plt.figure()
    dcct_mod.DCCT(d, fig=fig, arrow_width_px=8)
    arrows = [j for j in fig.artists if hasattr(j, "get_mutation_scale")]
    assert len(arrows) == 3  # an arrow for ph1, ph2, and vd
    # the shaft of the arrow is 0.2 times the mutation scale (in points)
    np.testing.assert_allclose(
        [j.get_mutation_scale() for j in arrows], 8 * 72 / fig.dpi / 0.2
    )
    with pytest.raises(ValueError):
        dcct_mod.DCCT(d, fig=fig, cmap="viridis")
    for removed_option in [
        "pass_frq_slice",
        "diagnostic",
        "ax",
        "spacing",
        "x_first",
        "renumber",
    ]:
        with pytest.raises(ValueError, match=removed_option + ".*deprecated"):
            dcct_mod.DCCT(d, fig=fig, **{removed_option: True})
    plt.close("all")


def test_dcct_alias_labels_and_incomplete_phase_cycle():
    d = phase_cycled_data()
    fig = plt.figure()
    dcct_mod.DCCT(d, fig=fig, max_coh_jump={"ph1": 1, "ph2": 1})
    # with max_coh_jump of 1, nothing aliases onto the ph2 element at 2
    assert [j.get_text() for j in fig.texts].count("X") == 2
    # 2 steps of 1/4 cycle don't go all the way around the circle, whether
    # we are in the phase cycling domain or the coherence domain
    d.ift("ph1")
    d.set_axis("ph1", np.r_[0:2] / 4)
    with pytest.raises(ValueError, match="once around the circle"):
        dcct_mod.DCCT(d, fig=fig)
    d.ft("ph1")
    with pytest.raises(ValueError, match="once around the circle"):
        dcct_mod.DCCT(d, fig=fig)
    plt.close("all")


def test_dcct_coherence_mask_and_shareaxis():
    def unmask_signal(coh_array):
        assert coh_array.dimlabels == ["ph1", "ph2"]
        assert coh_array.get_prop("coherence_pathway") == {"ph1": 0, "ph2": 1}
        coh_array["ph1", 0]["ph2", 1] = 1
        return coh_array

    d = phase_cycled_data()
    d.set_prop("coherence_pathway", {"ph1": 0, "ph2": 1})
    fig = plt.figure()
    ax_list, _ = dcct_mod.DCCT(d, fig=fig, coherence_unmask_fn=unmask_signal)
    # axes are numbered with ph2 varying fastest, so only number 1 is not
    # hatched
    assert [len(j.patches) for j in ax_list] == [1, 0, 1, 1, 1, 1, 1, 1]
    ax_list[3].set_xlim(-10, 10)
    assert ax_list[0].get_xlim() == (-10, 10)
    fig = plt.figure()
    ax_list, _ = dcct_mod.DCCT(d, fig=fig, shareaxis=False)
    ax_list[3].set_xlim(-10, 10)
    assert ax_list[0].get_xlim() != (-10, 10)
    plt.close("all")


def test_fake_data_uses_scan_based_temp_units():
    t2, vd = s.symbols("t2 vd")
    data = fake_data_mod.fake_data(
        s.exp(-vd / 0.2) * s.exp(1j * 2 * s.pi * 10 * t2),
        OrderedDict(
            [
                ("vd", nddata(np.linspace(0.0, 1.0, 4), "vd")),
                ("t2", nddata(np.linspace(0.0, 10e-3, 16), "t2")),
            ]
        ),
        {},
    )
    assert data.get_units("t2") == "s"
    assert data.dimlabels == ["vd", "t2"]
