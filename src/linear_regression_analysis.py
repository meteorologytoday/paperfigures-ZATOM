import matplotlib as mplt
import load_scan_data as lsd
import make_extended_data as med

import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib import rc

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

parser = argparse.ArgumentParser()
parser.add_argument('--folder', nargs='+')
parser.add_argument('--legend', nargs='*')
parser.add_argument('--colors', nargs='*')
parser.add_argument('--auto-color', action="store_true")
parser.add_argument('--output-bifur', default="")
parser.add_argument('--gamma-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--residue-threshold', type=float, default=1e-12)
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)

repNaN2None(args.gamma_rng)

data = []
coor = {}

target_vars = ["psi1000", "db_ew", "s1000", "chi1000", "d_cvt"] #"cvt_e", "cvt_w"]

folders = args.folder
legends   = args.legend

if (legends is None) or (legends is not None and len(legends) < len(folders)):
    print("Legends do not have the same length of folders. Use numbering instead")
    legends = ["%d" % (i,) for i in range(len(folders))]

print("Legends are")
print(legends)



# Data contains a list of folders (cases)
# and each case has multiple data files
data_to_delete = []
data = []
coor = None
loaded_varnames = ["Psib", "chi", "Q", "be", "bw", "res", "stable"]
for i, folder in enumerate(folders):

    print("Loading the folder: %s" % (folder,))

    _data, _coor = lsd.loadScanData(folder, loaded_varnames, load_coor=(coor is None))

    try:
        _data, _coor = lsd.loadScanData(folder, loaded_varnames, load_coor=(coor is None))
    except Exception as e:
        print("Error occurs. Skip this one.")
        data_to_delete.append(i)
        data.append(None)
        continue
 
    print("Number of records: %d" % (len(_data["Q"]))) 
    if coor is None:
        coor = _coor

    med.makeExtendedData(_data, coor)

    data.append(_data)

data_to_delete.reverse() # important. delete from last to first to avoid reordering
print("Delete index: ", data_to_delete)

for rm_idx in data_to_delete:
    del folders[rm_idx]
    del data[rm_idx]
    del legends[rm_idx]


for d in data:
    
    d["Q"] /= 1e6
    d["chi1000"] /= 1e6 * 1e-6
    d["psi1000"] /= 1e6 
    d["Psib"] /= 1e6
    d["db_ns"] *= 1e3
    d["db_ew"] *= 1e3
    d["s1000"] *= 1e5
    d["s1000_hlat"] *= 1e5

    y_ind = np.argmin(np.abs(coor["y_T"] - 60.0))
    z_ind = np.argmin(np.abs(coor["z_W"] - 1000.0))

    print("z_ind = %d, y_ind = %d" % (z_ind, y_ind)) 
    d["psi_fixed"] = d["Psib"][:, y_ind, z_ind] 

# remove data points that is not converging
for d in data:
    nan_idx = d["res"] > args.residue_threshold
    d["Q"][nan_idx] = np.nan
    d["Q"][d["Q"] > args.gamma_rng[1]] = np.nan
    d["Q"][d["Q"] < args.gamma_rng[0]] = np.nan
    for v in target_vars:
        d[v][nan_idx] = np.nan
 
if args.auto_color:

    print("Option --auto-color is on. Use self-generated colors.")
    cmap = plt.cm.get_cmap("Set1", len(legends))
    colors = [ cmap(i) for i in range(len(legends)) ]
    
else:

    if args.colors is None:
        colors = [
            'red',
            'darkorange',
            'darkgreen', 
            'blue',
            'violet',
            'black',
            "grey",
            "dodgerblue",
            "brown",
            "purple",
            "pink",
        ]
    else:

        colors = args.colors

        if len(colors) < len(legends):
            raise Exception("Colors provided have to be more than number of legends")
                   
print("Data loaded. Plotting now...")

# Multi linear regression
from sklearn.linear_model import LinearRegression
   
fig, ax = plt.subplots(1, 4, figsize=(18, 4), constrained_layout=True)

marker = ["o", "^", "+", "x", "s"]

for i in range(len(legends)):
    
    print("Doing linear regression of %s (%s)" % (legends[i], folders[i]))

    d = data[i]

    Q  = d["Q"]
    chi_dbdz = d["chi1000"] * d["s1000"]
    dq = d["cvt_w"] - d["cvt_e"]
    db_ew = d["db_ew"]
    psi = d["psi1000"]

    used_idx = np.isfinite(Q)
    Q        = Q[used_idx]
    chi_dbdz = chi_dbdz[used_idx]
    dq       = dq[used_idx]
    db_ew    = db_ew[used_idx]
    psi      = psi[used_idx]

    if len(Q) == 0:
        print("This case has no valid point")
        continue
    """
    X = np.zeros((len(Q), 3))
    X[:, 0] = chi_dbdz[:]
    X[:, 1] = dq[:]
    X[:, 2] = Q[:]

    X = np.zeros((len(Q), 2))
    X[:, 0] = chi_dbdz[:]
    X[:, 1] = dq[:]
    """

    # ======================
    X = np.zeros((len(Q), 3))
    X[:, 0] = chi_dbdz[:]
    X[:, 1] = dq[:]
    X[:, 2] = Q[:]
    
    y = db_ew
    reg = LinearRegression(normalize=True).fit(X, y)
    y_predict = reg.predict(X) 
    print(reg.score(X, y))
    print(reg.coef_)

    ax[0].scatter(y, y_predict, s=10, marker=marker[i], label=legends[i])


    # ======================
    X = np.zeros((len(Q), 1))
    X[:, 0] = chi_dbdz[:]
    
    y = db_ew
    reg = LinearRegression(normalize=True).fit(X, y)
    y_predict = reg.predict(X) 
    print(reg.score(X, y))
    print(reg.coef_)
    ax[1].scatter(y, y_predict, s=10, marker=marker[i], label=legends[i])


    # ======================
    X = np.zeros((len(Q), 1))
    X[:, 0] = dq[:]
    
    y = db_ew
    reg = LinearRegression(normalize=True).fit(X, y)
    y_predict = reg.predict(X) 
    print(reg.score(X, y))
    print(reg.coef_)
    ax[2].scatter(y, y_predict, s=10, marker=marker[i], label=legends[i])



    # ======================
    X = np.zeros((len(Q), 1))
    X[:, 0] = Q[:]
    
    y = db_ew
    reg = LinearRegression(normalize=True).fit(X, y)
    y_predict = reg.predict(X) 
    print(reg.score(X, y))
    print(reg.coef_)
    ax[3].scatter(y, y_predict, s=10, marker=marker[i], label=legends[i])

for i in range(4):
    ax[i].plot([-3, 3], [-3, 3], ls="dashed", color="gray")


ax[0].set_title("All comp")
ax[1].set_title("chi_dbdz")
ax[2].set_title("dq")
ax[3].set_title("Q")


ax[0].legend()

plt.show()

