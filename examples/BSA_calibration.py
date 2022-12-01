from pylab import *
from pyspecdata import *
import numpy as np
import matplotlib.pyplot as plt

# HERE we get dictionaries of nddata
dataWKBSA = find_file(
    "221110_BSAexerciseWK_0p07-0percentBSAcalibration.BSW", exp_type="UV_Vis/BSA_calib"
)
# print("the experiments present in this file are:", data.keys())
wv = "$\\lambda$"
# HERE make a dictionary with keys equal to labels (using concentration,
# and with values equal to the actual nddata objects for the spectra
list_of_runs = [
    {
        "Water": dataWKBSA["UP_H2O"],
        "0 %": dataWKBSA["K-PBSbuffer"],
        "0.0078 %": dataWKBSA["J-0p0078%BSA"],
        "0.0156 %": dataWKBSA["H-0p0156%BSA"],
        "0.0234 %": dataWKBSA["G-0p0234%BSA_actual"],
        "0.0311 %": dataWKBSA["F-0p0311%BSA"],
        "0.0389 %": dataWKBSA["E-0p0389%BSA"],
        "0.0466 %": dataWKBSA["D-0p0466%BSA"],
        "0.0544 %": dataWKBSA["C-0p0544%BSA"],
        "0.0622 %": dataWKBSA["B-0p0622%BSA"],
        "0.0700 %": dataWKBSA["A-0p0700%BSA"],
    }
]
# Dictionary with labels of each run in the UV file
# HERE add a dilution factor -- we'll multiply in the final plot
dilution_factor = [1]
# HERE give labels of equal length
label_of_run = ["WK BSA run"]
finalcurve = [{} for j in range(len(list_of_runs))]  # two dicts, one for each
A280 = [{} for j in range(len(list_of_runs))]  # two dicts, one for each
print(range(len(list_of_runs)))
# }}}
whichunits = "%"
# prop_cycle = plt.rcParams["axes.prop_cycle"]
# prop_cycle = prop_cycle * plt.cycler(alpha=[0.5])
# print(plt.rcParams["ytick.major.size"])
# plt.rcParams["ytick.major.size"] = 3
# print(plt.rcParams["ytick.major.size"])
# plt.rcParams["font.size"] = 8
# plt.rcParams["legend.fontsize"] = 8
for j, thisset in enumerate(list_of_runs):
    # thecycler = prop_cycle()
    for k, v in thisset.items():
        v.name(label_of_run[j])
        # kw = next(thecycler)
        thisd = v[wv:(280, None)]
        if k not in ["Water"]:
            A280[j][k] = v[wv:(276, 281)].C.mean(wv).item()
        if k not in [f"0 %", "Water"]:
            BSA_diff = thisd[wv:(280, None)] - thisset[f"0 %"][wv:(280, None)]
            finalcurve[j][k] = BSA_diff[wv:(276, 281)].C.mean(wv).item()
figure()
for j in range(len(finalcurve)):
    thisdil = dilution_factor[j]
    x, y = zip(*finalcurve[j].items())
    x = [float(j.replace(f" {whichunits}", "")) for j in x]
    y = list(y)
    plot(
        x, array(y) * thisdil, "o", alpha=0.5, label=f"{label_of_run[j]} bsln_adj_A280"
    )
    plt.plot(x, y)
    # }}}
    # {{{ A280 standin
    x, y = zip(*A280[j].items())
    x = [float(j.replace(f" {whichunits}", "")) for j in x]
    y = list(y)
    plot(x, array(y) * thisdil, "x", alpha=0.5, label=f"{label_of_run[j]} A280")
    plt.plot(x, y)
    # }}}
plt.legend(**dict(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0))
plt.tight_layout()
title("final calibration curve")
print(A280)
print(finalcurve)
# }}}
show()
