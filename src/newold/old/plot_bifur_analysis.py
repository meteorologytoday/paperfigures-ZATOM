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
parser.add_argument('--output-marks', default="")
parser.add_argument('--offset-marks', type=float, default=0)
parser.add_argument('--gamma-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--psi-fixed-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--s1000-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--db_ns-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--db_ew-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--cvt_e-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--cvt_w-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--marks', nargs='*', type=float, default=[])
parser.add_argument('--marks-pos', nargs='*', type=str, default=[])
parser.add_argument('--residue-threshold', type=float, default=1e-12)
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)

repNaN2None(args.gamma_rng)
repNaN2None(args.psi_fixed_rng)
repNaN2None(args.s1000_rng)
repNaN2None(args.db_ns_rng)
repNaN2None(args.db_ew_rng)
repNaN2None(args.cvt_w_rng)
repNaN2None(args.cvt_e_rng)

data = []
coor = {}

yvars = ["Psib", "s1000", "db_ns", "db_ew"]

folders = args.folder
legends   = args.legend

if (legends is None) or (legends is not None and len(legends) < len(folders)):
    print("Legends do not have the same length of folders. Use numbering instead")
    legends = ["%d" % (i,) for i in range(len(folders))]

print("Legends are")
print(legends)

if len(args.marks) % 2 != 0:
    raise Exception("Length of mark is not multiple of 2. It has to be a list of (x, y) pairs. ")


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

# remove data points that is not converging
for d in data:
    for v in yvars:
        d[v][d["res"] > args.residue_threshold] = np.nan

for d in data:
    
    d["Q"] /= 1e6
    d["Psib"] /= 1e6
    d["db_ns"] *= 1e3
    d["db_ew"] *= 1e3
    d["s1000"] *= 1e5
    d["s1000_hlat"] *= 1e5

    y_ind = np.argmin(np.abs(coor["y_T"] - 60.0))
    z_ind = np.argmin(np.abs(coor["z_W"] - 1000.0))

    print("z_ind = %d, y_ind = %d" % (z_ind, y_ind)) 
    d["psi_fixed"] = d["Psib"][:, y_ind, z_ind] 

# Pick out the marks
def dist(x1,y1,x2,y2):
    return ((x1-x2)**2 + (y1-y2)**2)**0.5

nmarkpairs = int( len(args.marks) / 2 )
shortest_dist = np.zeros((nmarkpairs,)) + 1e5
mark_index = np.zeros((nmarkpairs, 2), dtype=int) - 1
        
for k in range(nmarkpairs):
    print("Search for the marker case closest to : (gamma, psi) = (%f, %f)" % (args.marks[k*2], args.marks[k*2+1]))
    for i, d in enumerate(data):
        for s in range(len(d["Q"])):
            _dist = dist(d["Q"][s]*10, d["psi_fixed"][s], args.marks[k*2]*10, args.marks[k*2+1])
            if _dist < shortest_dist[k]:
                shortest_dist[k] = _dist
                mark_index[k, 0] = i
                mark_index[k, 1] = s
    
    print("Found: (%d, %d)" % (mark_index[k, 0], mark_index[k, 1]))



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

plot_z_W = coor["z_W"] / 1e3
plot_z_T = coor["z_T"] / 1e3

fig, ax = plt.subplots(2, 3, figsize=(12, 8), squeeze=False, constrained_layout=True)

ax_flat = ax.flatten(order='F')
       
target_vars = ["psi_fixed", "db_ew", "s1000", "db_ns", "cvt_e", "cvt_w"]

for i in range(len(legends)):
    
    print("Plotting %s (%s)" % (legends[i], folders[i]))

    d = data[i]
    offset = 0

    vals, rngs = lsd.detectRanges(d["stable"])

    for j, val in enumerate(vals):
        
        rng = rngs[j]

        x   = d["Q"][rng]
        res = d["res"][rng]


        for k, var in enumerate(target_vars):
    
            y   = d[var][rng]

            linestyle = "solid" if val == 1.0 else "dashed"
 
            kw_label = {
                "color" : colors[i],
            }

            if j == 0 and k == 0:
                kw_label['label'] = legends[i]

            ax_flat[k].plot(x, y, zorder=1, linestyle=linestyle, **kw_label)



for k in range(nmarkpairs):
    print("Marking the %d-th marker." % k)
    print(mark_index[k])
    d = data[mark_index[k,0]]
    s = mark_index[k, 1]
    
    for l, var in enumerate(target_vars):
        #ax_flat[l].scatter(d["Q"][s], d[var][s], s=10, marker="o", c="yellow", edgecolor="k", zorder=50)
        ax_flat[l].scatter(d["Q"][s], d[var][s], s=10, marker="o", c="k", edgecolor="k", zorder=50)
        ax_flat[l].text(d["Q"][s] + args.offset_marks, d[var][s], "123456789"[k], size=12, va="center", ha="left")

ax[0, 0].legend(fontsize=12, handlelength=1.0, labelspacing=0.25)


labels = {
    "psi_fixed"  : r"$\psi$",
    "s1000"  : r"$s$",
    "s1000_hlat"  : r"$s_{\mathrm{hlat}}$",
    "db_ns"  : r"$\delta b$",
    "db_ew"  : r"$b_e^* - b_w^*$",
    "cvt_e"  : r"$\tilde{q}_e$",
    "cvt_w"  : r"$\tilde{q}_w$",
}

units = {
    "psi_fixed" : "[Sv]",
    "s1000"  : r"[ $ \times 10^{-5} \mathrm{s}^{-2} $]",
    "s1000_hlat"  : r"[ $ \times 10^{-5} \mathrm{s}^{-2} $]",
    "db_ns"  : r"[ $ \times 10^{-3} \mathrm{m} / \mathrm{s}^{2} $]",
    "db_ew"  : r"[ $ \times 10^{-3} \mathrm{m} / \mathrm{s}^{2} $]",
    "cvt_e"  : r"",
    "cvt_w"  : r"",
}

for l, var in enumerate(target_vars):
    ax_flat[l].set_ylabel("%s" % (units[var],))
    
    ylim_attr = "%s_rng" % (var,)
    if hasattr(args, ylim_attr):
        ax_flat[l].set_ylim(getattr(args, ylim_attr))

    ax_flat[l].set_xlabel("$\\gamma$ [Sv]")
    ax_flat[l].set_xlim(args.gamma_rng)
    ax_flat[l].set_title("(%s) %s" % ("abcdefg"[l], labels[var],))

if args.output_bifur != "":
    fig.savefig(args.output_bifur, dpi=300)


if nmarkpairs != 0:

    # finding a proper range for buoyancy
    b_factor = 1e3
    b_rng = [ np.nan, np.nan ]
    b_concat = np.array([], dtype=float)
    for k in range(nmarkpairs):
        d = data[mark_index[k,0]]
        s = mark_index[k, 1]
        b_concat = np.concatenate((b_concat, d["bw_bnd"][s, :, :].flatten()))
        b_concat = np.concatenate((b_concat, d["be"][s, :, :].flatten()))

    b_mean = np.mean(b_concat)
    b_std  = np.std(b_concat)

    nstd = 3
    b_neworigin = b_mean - nstd * b_std
    b_rng = np.array([0.0, 2 * nstd * b_std])

    # need b_factor because later on all the buoyancy data will be multiplied by b_factor
    b_rng[0] = np.floor(b_rng[0] * b_factor / 2) * 2
    b_rng[1] = np.floor(b_rng[1] * b_factor / 2) * 2

    levels_b = np.arange(b_rng[0], b_rng[1], 2)
    #levels_b = np.arange(.24, .28, 0.005)
    levels_ui = np.arange(-100, 100, 0.5)
    levels_we = np.arange(-100, 100, .1)
    levels_ww = np.arange(-100, 100, .5)
    levels_psi = np.arange(-20, 20, 1)
    dwf_color = (0.6, 0.75, 1.0)

    b_cmap = "bwr"
    
    fig2 = plt.figure(figsize=(4*nmarkpairs+1, 4*3)) 
    gs0 = gridspec.GridSpec(3, nmarkpairs+1, figure=fig, width_ratios=nmarkpairs * [1] + [0.05], wspace=0.3)
    #fig2, ax2 = plt.subplots(3, nmarkpairs, sharey=False, sharex=True, figsize=(4*nmarkpairs, 4*3), gridspec_kw={'bottom':0.2}, squeeze=False) 
 
    for k in range(nmarkpairs):
 
        print("Plotting the %d-th marker streamfunction" % k)
        d = data[mark_index[k,0]]
        s = mark_index[k, 1]
    
        _ax = [ fig2.add_subplot(gs0[i, k]) for i in range(3) ]
      
        be_plot = ( d["be"][s, :, :].transpose() - b_neworigin ) * b_factor 
        bw_plot = ( d["bw_bnd"][s, :, :].transpose() - b_neworigin ) * b_factor 

        #print("Max and min of be: %.2e ; %.2e" % (np.amax(be_plot), np.amin(be_plot)))
        #print("Max and min of bw: %.2e ; %.2e" % (np.amax(bw_plot), np.amin(bw_plot)))

        #print("Max and min of bw: %.2e ; %.2e" % (np.amax(d["bw"]), np.amin(d["bw"])))
        #print("Max and min of ui: %.2e ; %.2e" % (np.amax(d["ui"]), np.amin(d["ui"])))
        #print("Max and min of we: %.2e ; %.2e" % (np.amax(d["we"]), np.amin(d["we"])))
        #print("Max and min of ww: %.2e ; %.2e" % (np.amax(d["ww"]), np.amin(d["ww"])))

        #_ax = ax2[:, k]
       

        CS = _ax[0].contour(coor["y_V"], plot_z_W, d["Psib"][s, :, :].transpose(), levels_psi, colors="black")
        _ax[0].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=4)
        
        b_mappable = _ax[1].contourf(coor["y_T"], plot_z_T, bw_plot, levels_b, cmap=b_cmap, extend='both')
        #CS = _ax[1].contour(coor["y_V"], plot_z_W, d["Psib"][s, :, :].transpose(), levels_psi, colors="black")
        #_ax[1].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=1)

        b_mappable = _ax[2].contourf(coor["y_T"], plot_z_T, be_plot, levels_b, cmap=b_cmap, extend="both")
        #CS = _ax[2].contour(coor["y_T"], plot_z_W, d["chi"][s, :, :].transpose(), levels_psi, colors="black")
        #_ax[2].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=1)
        
        cs_dwf_west = _ax[1].contourf(coor['y_T'], plot_z_T, d['dwf_west'][s, :, :].transpose(), [0, 0.5,1.5], colors="none", hatches=[None, ".."])
        cs_dwf_east = _ax[2].contourf(coor['y_T'], plot_z_T, d['dwf_east'][s, :, :].transpose(), [0, 0.5,1.5], colors="none", hatches=[None, ".."])


        #for _cs in [cs_dwf_east, cs_dwf_west]:
        #    for _, collection in enumerate(_cs.collections):
        #        collection.set_edgecolor((1.0, 0.5, 0.5))
        #        collection.set_linewidth(0.0)
 


        _ax[0].set_title("State %s ($%.2f \\, \\mathrm{Sv}$)" % ("123456789"[k], d["Q"][s]))

        if k == 0:
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

        
    cax = fig2.add_subplot(gs0[:, -1])
    plt.colorbar(mappable=b_mappable, cax=cax, cmap=b_cmap, orientation="vertical", ticks=[], label="$b_w^*$ and $b_e^*$ [$\\mathrm{m}^2 / \\mathrm{s}$]")
    
    if args.output_marks != "":
        fig2.savefig(args.output_marks, dpi=300)



plt.show()

