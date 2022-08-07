import matplotlib as mplt
import load_scan_data as lsd
import make_extended_data as med

import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib import rc
import matplotlib.transforms as transforms

import os
import re
import pprint
default_linewidth = 1.0;
default_ticksize = 10.0;

mplt.rcParams['lines.linewidth'] =   default_linewidth;
mplt.rcParams['axes.linewidth'] =    default_linewidth;
mplt.rcParams['xtick.major.size'] =  default_ticksize;
mplt.rcParams['xtick.major.width'] = default_linewidth;
mplt.rcParams['ytick.major.size'] =  default_ticksize;
mplt.rcParams['ytick.major.width'] = default_linewidth;

rc('font', **{'size': 15.0});
rc('axes', **{'labelsize': 15.0});
rc('mathtext', **{'fontset':'stixsans'});

import matplotlib.pyplot as plt

import os,sys, argparse
from netCDF4 import Dataset
import numpy as np


def repNaN2None(arr):
    for i, elm in enumerate(arr):
        if np.isnan(elm):
            arr[i] = None

varying_axis = {"gamma" : "xi", "xi" : "gamma"}

param_map = {"gamma": "Q", "xi" : "xi"}
param_rng = {
    "Q"     : [-0.01, 2.0],
    "xi"    : [-2, 2],
}

detailed = False

if detailed == True:
    # The following is used to measure the bifurcation regime.
    selected_data = {
        'gamma'    : np.array([5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200], dtype=float) / 1000,  # Sv
        'xi'       : np.array([-1.50, -1.00, -0.50, 0, 0.50], dtype=float),
    }
else: 
    # The following is used to output figure in the paper
    selected_data = {
        #'gamma'    : list(np.array([0, 100, 200, 300, 400, 500], dtype=float) / 1000),  # Sv
        'gamma'    : list(np.array([750 + i*25 for i in range(11)], dtype=float) / 1000),  # Sv
        #'xi'       : list(np.array([-1.50, -1.00, -0.9, -0.8, -0.7, -0.6, -0.50, 0, 0.50], dtype=float)),
        'xi'       : list(np.linspace(-0.5, 0.5, 11)),
    }



for k, v in selected_data.items():
    selected_data[k] = selected_data[k][::-1]

residue_threshold = 3e-13

folders = {}

for axis, pts in selected_data.items():

    folders[axis] = []

    for i, pt in enumerate(pts):



        if axis == "xi":
            folders[axis].append("data/continuation_data_H295m_073122/batch_H295m_fixed_xi/output_redo_balanced_tanh/CM_xi%08d_pos/lb8" % (round(pt * 100), ))
        elif axis == "gamma":
            folders[axis].append("data/continuation_data_H295m_073122/batch_H295m_fixed_gamma/output_redo_balanced_tanh/CM_gamma%08d_neg/lb8" % (round(pt * 1000), ))
    


#target_vars = ["psi1000", "db_ew", "s1000_hlat", "chi1000", "d_cvt", "chi_dbdz"] #"cvt_e", "cvt_w"]

# Data contains a list of folders (cases)
# and each case has multiple data files
data_to_delete = {}
data = {}
coor = None
loaded_varnames = ["Psib", "chi", "be", "bw", "res", "stable"]

for axis in ["gamma", "xi"]:

    data[axis] = []
    data_to_delete[axis] = []

    
    actual_loaded_varnames = loaded_varnames.copy()
    actual_loaded_varnames.append(param_map[varying_axis[axis]])


    for i, folder in enumerate(folders[axis]):

        print("Loading the folder: %s" % (folder,))

        #_data, _coor = lsd.loadScanData(folder, actual_loaded_varnames, load_coor=(coor is None))

        try:
            _data, _coor = lsd.loadScanData(folder, actual_loaded_varnames, load_coor=(coor is None))
        except Exception as e:
            print("Error occurs. Skip this one.")
            data_to_delete[axis].append(i)
            data[axis].append(None)
            continue
     
        print("Number of records: %d" % (_data[actual_loaded_varnames[0]].shape[0])) 
        if coor is None:
            coor = _coor

        med.makeExtendedData(_data, coor)

        data[axis].append(_data)

for axis in ["gamma", "xi"]:
    data_to_delete[axis].reverse() # important. delete from last to first to avoid reordering

for axis in ["gamma", "xi"]:
    print("[%s] Delete index: " % (axis,) , data_to_delete[axis])
    for rm_idx in data_to_delete[axis]:
        del folders[axis][rm_idx]
        del data[axis][rm_idx]
        del selected_data[axis][rm_idx]


for axis in ["gamma", "xi"]:
    for d in data[axis]:
        if "Q" in d:
            d["Q"] /= 1e6
    
        d["chi1000"] /= 1e6 * 1e-6
        d["psi1000"] /= 1e6 
        d["Psib"] /= 1e6
        d["db_ns"] *= 1e3
        d["db_ew"] *= 1e3
        d["s1000"] *= 1e5
        d["s1000_hlat"] *= 1e5

# remove data points that is not converging

for axis in ["gamma", "xi"]:

    param = param_map[varying_axis[axis]]

    for d in data[axis]:
        nan_idx = d["res"] > residue_threshold
        d[param][nan_idx] = np.nan
        d[param][d[param] < param_rng[param][0]] = np.nan
        d[param][d[param] > param_rng[param][1]] = np.nan
        for v in ["psi1000"]:
            d[v][nan_idx] = np.nan
     


print("Data loaded. Plotting now...")

fig, ax = plt.subplots(1, 2, figsize=(12, 6), squeeze=True, constrained_layout=True)

axis_i_map = {
    "xi" : 0,
    "gamma" : 1,
}

other_axis_i_map = {
    "xi" : 1,
    "gamma" : 0,
}


for axis in ["xi", "gamma"]:

    Ncases = len(selected_data[axis])

    cmap = plt.cm.get_cmap("gist_rainbow", Ncases)
    colors = [ cmap(i) for i in range(Ncases)]


    axis_i       = axis_i_map[axis] 
    other_axis_i = other_axis_i_map[axis] 
    
    _ax = ax[axis_i]
    _other_ax = ax[other_axis_i]

    param = param_map[varying_axis[axis]]


   
    # Set title
    #_ax.text(0.1, 0.9, "(%s)" % ("abc"[axis_i]), transform=_ax.transAxes, va="top", ha="left", size=25)
 
    # Plot the scanned line
    trans = transforms.blended_transform_factory(_other_ax.transData, _other_ax.transAxes)

    for fixed_pt in selected_data[axis]:
        _other_ax.plot([fixed_pt] * 2, [0, 1], ls="--", color="#aaaaaa", transform=trans, zorder=1)
        

    for i, d in enumerate(data[axis]):
    
        print("Plotting %s (%s)" % (axis, folders[axis][i]))

        vals, rngs = lsd.detectRanges(d["stable"])

        for j, val in enumerate(vals):
            
            rng = rngs[j]

            x   = d[param][rng]
            res = d["res"][rng]


            for k, var in enumerate(["psi1000",]):
        
                y   = d[var][rng]

                linestyle = "solid" if val == 1.0 else "dotted"
     
                kw_label = {
                    "color" : colors[i],
                }

                if j == 0 and k == 0:
                    fixed_value = selected_data[axis][i]
                    if axis == "gamma":
                        kw_label['label'] = "$ \\gamma = %.3f $ Sv" % (fixed_value,)
                    elif axis == "xi":
                        kw_label['label'] = "$ \\xi = %.2f $" % (fixed_value,)

                _ax.plot(x, y, zorder=10, linestyle=linestyle, **kw_label)

    _ax.set_title("(%s)" % ("abcdefg"[axis_i] ))
    _ax.set_xlabel(["$\\gamma$ [Sv]",  "$\\xi$" ][axis_i])

    _ax.set_ylabel("$\\left\\langle\\psi\\right\\rangle$ [Sv]")
    _ax.set_xlim([[0.0, 1.0], [-1.8, 0.6]][axis_i])
    _ax.set_ylim([-1, 6.0])

    _ax.legend(loc="upper left", fontsize=12, handlelength=1.0)

fig.savefig("figures/ZATOM_bifur_gamma_xi.png", dpi=300)

plt.show()

