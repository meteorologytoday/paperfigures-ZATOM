import load_scan_data as lsd
import make_extended_data_v3 as med

import os
import re
import pprint
default_linewidth = 1.0;
default_ticksize = 10.0;


import os,sys, argparse
from netCDF4 import Dataset
import numpy as np


def repNaN2None(arr):
    for i, elm in enumerate(arr):
        if np.isnan(elm):
            arr[i] = None

parser = argparse.ArgumentParser()
parser.add_argument('--folder', type=str, nargs=2)
parser.add_argument('--diff-idx', type=int, nargs=2)
parser.add_argument('--titles', type=str, nargs=2, default=["State 0", "State 1 - State 0"])
parser.add_argument('--output', type=str, default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)

data = []
coor = {}

folders = args.folder


# Data contains a list of folders (cases)
# and each case has multiple data files
data_to_delete = []
data = []
coor = None
loaded_varnames = ["Q", "Psib", "chi", "be", "bw", "qw", "qe", "res", "stable"] + med.necessary_variables

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
 
    print("Number of records: %d" % (len(_data[list(_data.keys())[0]]))) 
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


# find diff
print("Taking difference")
print("diff_idx[0] = ", args.diff_idx[0])
print("diff_idx[1] = ", args.diff_idx[1])
d_diff = {
    varname : data[1][varname][args.diff_idx[1]] - data[0][varname][args.diff_idx[0]] for varname in data[1].keys()
}


d_0 = {
    varname : data[0][varname][args.diff_idx[0]] for varname in data[0].keys()
}
d_1 = {
    varname : data[1][varname][args.diff_idx[1]] for varname in data[1].keys()
}


print("Max Psi in d_0: %f Sv " % (np.amax(d_0["Psib"])/1e6,))
print("Max Psi in d_1: %f Sv " % (np.amax(d_1["Psib"])/1e6,))
print("Max dPsi in d_diff: %f Sv " % (np.amax(d_diff["Psib"]) / 1e6,))
print("Min dPsi in d_diff: %f Sv " % (np.amin(d_diff["Psib"]) / 1e6,))

print("Loading matplotlib...")
import matplotlib as mplt

if args.no_display:
    print("`--no-display` is set.")
    mplt.use("Agg")

else:
    print("`--no-display` is not set. Will show figures...")
    mplt.use("TkAgg")

import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib import rc

import tool_fig_config
import cmocean
print("Done")
mplt.rcParams['lines.linewidth'] =   default_linewidth;
mplt.rcParams['axes.linewidth'] =    default_linewidth;
mplt.rcParams['xtick.major.size'] =  default_ticksize;
mplt.rcParams['xtick.major.width'] = default_linewidth;
mplt.rcParams['ytick.major.size'] =  default_ticksize;
mplt.rcParams['ytick.major.width'] = default_linewidth;

rc('font', **{'size': 15.0});
rc('axes', **{'labelsize': 15.0});
rc('mathtext', **{'fontset':'stixsans'});



plot_z_W = - coor["z_W"] 
plot_z_T = - coor["z_T"] 

ncol = 2
nrow = 1

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 6,
    h = 4,
    wspace = 2.5,
    hspace = 1.5,
    w_left = 1.0,
    w_right = 1.5,
    h_bottom = 1.0,
    h_top = 1.0,
    ncol = ncol,
    nrow = nrow,
)

fig, ax = plt.subplots(
    nrow, ncol,
    figsize=figsize,
    subplot_kw=dict(aspect="auto"),
    gridspec_kw=gridspec_kw,
    constrained_layout=False,
    squeeze=False,
    sharex=False,
)


bdiff_scale_power = -3
dbdiff_scale_power = -3

bdiff = (d_0["be"] - d_0["bw"]) / 10**bdiff_scale_power
dbdiff = (d_diff["be"] - d_diff["bw"]) / 10**dbdiff_scale_power

bdiff_ticks = [-2, -1, 0, 1, 2]
dbdiff_ticks = [-2, -1, 0, 1, 2]

levels_dbdiff = np.linspace(-1, 1, 11) * 2 
levels_bdiff = np.linspace(-1, 1, 11) * 2

levels_ui = np.arange(-100, 100, 0.5)
levels_we = np.arange(-100, 100, .1)
levels_ww = np.arange(-100, 100, .5)
levels_dpsi = np.arange(-20, 20, 0.5)
levels_psi = np.arange(-20, 20, 1)
dwf_color = (0.6, 0.75, 1.0)
b_cmap = "cmo.balance"


# ======================

_ax = ax.flatten()[0]

CS = _ax.contour(coor["y_V"], plot_z_W, d_0["Psib"].T / 1e6, levels_psi, colors="black")
_ax.clabel(CS, CS.levels, inline=True, fmt="%.1f", fontsize=10, inline_spacing=4)

b_mappable_d_0 = _ax.contourf(coor["y_T"], plot_z_T, bdiff.T, levels_bdiff, cmap=b_cmap, extend="both")

cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
plt.colorbar(mappable=b_mappable_d_0, cax=cax, cmap=b_cmap, orientation="vertical", label="$b_e- b_w$ [$ \\times 10^{%d} \\mathrm{m}^2 / \\mathrm{s}$]" % (bdiff_scale_power), ticks=bdiff_ticks)

#dwf_west = np.zeros_like(d_diff['dwf_west'].copy())
#ddq = np.zeros_like(d_diff["dq"])
#ddq[d_diff["dq"] < 0] = 1.0

#cs_dwf_west = _ax.contourf(coor['y_T'], plot_z_T, ddq, [0, 0.5, 1.5], colors="none", hatches=[None, ".."])

_ax.set_title("(a) %s" % (args.titles[0]),)


# ======================

_ax = ax.flatten()[1]

CS = _ax.contour(coor["y_V"], plot_z_W, d_diff["Psib"].T / 1e6, levels_dpsi, colors="black")
_ax.clabel(CS, CS.levels, inline=True, fmt="%.1f", fontsize=10, inline_spacing=4)

b_mappable_d_diff = _ax.contourf(coor["y_T"], plot_z_T, dbdiff.T, levels_dbdiff, cmap=b_cmap, extend="both")

cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
plt.colorbar(mappable=b_mappable_d_diff, cax=cax, cmap=b_cmap, orientation="vertical", label="$\\Delta \\left( b_e- b_w \\right)$ [$ \\times 10^{%d} \\mathrm{m}^2 / \\mathrm{s}$]" % (dbdiff_scale_power), ticks=dbdiff_ticks)


#dwf_west = np.zeros_like(d_diff['dwf_west'].copy())
#ddq = np.zeros_like(d_diff["dq"])
#ddq[d_diff["dq"] < 0] = 1.0

#cs_dwf_west = _ax.contourf(coor['y_T'], plot_z_T, ddq, [0, 0.5, 1.5], colors="none", hatches=[None, ".."])
_ax.set_title("(b) %s" % (args.titles[1]),)


for _ax in ax.flatten():
    _ax.set_xlabel("Latitude [ deg ]")
    _ax.set_ylabel("Depth [ km ]")
    _ax.invert_yaxis()
    _ax.grid()



#be_plot = ( d["be"][s, :, :].transpose() - b_neworigin ) * b_factor 
#bw_plot = ( d["bw_bnd"][s, :, :].transpose() - b_neworigin ) * b_factor 

#print("Max and min of be: %.2e ; %.2e" % (np.amax(be_plot), np.amin(be_plot)))
#print("Max and min of bw: %.2e ; %.2e" % (np.amax(bw_plot), np.amin(bw_plot)))

#print("Max and min of bw: %.2e ; %.2e" % (np.amax(d["bw"]), np.amin(d["bw"])))
#print("Max and min of ui: %.2e ; %.2e" % (np.amax(d["ui"]), np.amin(d["ui"])))
#print("Max and min of we: %.2e ; %.2e" % (np.amax(d["we"]), np.amin(d["we"])))
#print("Max and min of ww: %.2e ; %.2e" % (np.amax(d["ww"]), np.amin(d["ww"])))

    #_ax = ax2[:, k]
   


#b_mappable = _ax[1].contourf(coor["y_T"], plot_z_T, bw_plot, levels_b, cmap=b_cmap, extend='both')

#CS = _ax[1].contour(coor["y_V"], plot_z_W, d["Psib"][s, :, :].transpose(), levels_psi, colors="black")
#_ax[1].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=1)

#b_mappable = _ax[2].contourf(coor["y_T"], plot_z_T, be_plot, levels_b, cmap=b_cmap, extend="both")

#CS = _ax[2].contour(coor["y_T"], plot_z_W, d["chi"][s, :, :].transpose(), levels_psi, colors="black")
#_ax[2].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=1)

#cs_dwf_west = _ax[1].contourf(coor['y_T'], plot_z_T, d['dwf_west'][s, :, :].transpose(), [0, 0.5,1.5], colors="none", hatches=[None, ".."])
#cs_dwf_east = _ax[2].contourf(coor['y_T'], plot_z_T, d['dwf_east'][s, :, :].transpose(), [0, 0.5,1.5], colors="none", hatches=[None, ".."])


#for _cs in [cs_dwf_east, cs_dwf_west]:
#    for _, collection in enumerate(_cs.collections):
#        collection.set_edgecolor((1.0, 0.5, 0.5))
#        collection.set_linewidth(0.0)

"""

        for __ax in _ax:
            __ax.set_ylabel("Depth [km]")
    
    _ax[-1].set_xlabel("Latitude [${}^\\circ\\mathrm{N}$]")

    for __ax in _ax:
        __ax.set_xticks([20,30,40,50,60])

    _ax[0].set_yticks([0, 1, 2, 3, 4])
    _ax[1].set_yticks([0, 0.5, 1])
    _ax[2].set_yticks([0, 0.5, 1])
    _ax[1].set_yticklabels(["0", ".5", "1"])
    _ax[2].set_yticklabels(["0", ".5", "1"])

    if k != 0:
        __ax.set_yticklabels([""] * len(__ax.get_yticks()))

    _ax[0].set_ylim([4.5, 0])
    _ax[1].set_ylim([1, 0])
    _ax[2].set_ylim([1, 0])

""" 
#cax = fig2.add_subplot(gs0[:, -1])
if args.output != "":
    print("Output file: ", args.output)
    fig.savefig(args.output, dpi=300)

if not args.no_display:
    print("Showing figures...")
    plt.show()

