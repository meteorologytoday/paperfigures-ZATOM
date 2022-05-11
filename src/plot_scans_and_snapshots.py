import matplotlib as mplt
import load_scan_data as lsd
#mplt.use("TkAgg")

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
parser.add_argument('--output-bifur', default="")
parser.add_argument('--output-marks', default="")
parser.add_argument('--F-rng', nargs=2, type=float)
parser.add_argument('--psi-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--numbering-jmp', type=int, default=10)
parser.add_argument('--show-numbering', action="store_true")
parser.add_argument('--show-filename', action="store_true")
parser.add_argument('--metric', type=str, default="max", choices=["max", "60N", "40N", "20N", "mean"])
parser.add_argument('--marks', nargs='*', type=float, default=[])
parser.add_argument('--legend-coor', nargs='*', type=float, default=[])
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)

repNaN2None(args.psi_rng)
data = []
coor = {}


folders = args.folder
legends   = args.legend

if (legends is None) or (legends is not None and len(legends) < len(folders)):
    print("Legends do not have the same length of folders. Use numbering instead")
    legends = ["%d" % (i,) for i in range(len(folders))]

print(legends)

if len(args.marks) % 2 != 0:
    raise Exception("Length of mark is not multiple of 2. It has to be a list of (x, y) pairs. ")


# Data contains a list of folders (cases)
# and each case has multiple data files

data = []
coor = None
loaded_varnames = ["Psib", "Q", "be", "bw", "we", "res", "stable", "ui", "ww"]
for i, folder in enumerate(folders):

    print("Loading the folder: %s" % (folder,))

    _data, _coor = lsd.loadScanData(folder, loaded_varnames, load_coor=(coor is None))
    data.append(_data)

    print("Number of records: %d" % (len(_data["Q"]))) 
    if coor is None:
        coor = _coor


    _data["chi_max"] = np.amax(- _data["we"] * np.cos(coor["y_T"] * np.pi/180.0)[None, :, None] * (6.4e6 * 15*np.pi/180.0) * (6.4e6 * 1*np.pi/180.0), axis=2) / 1e6
    _data["chi_min"] = np.amin(- _data["we"] * np.cos(coor["y_T"] * np.pi/180.0)[None, :, None] * (6.4e6 * 15*np.pi/180.0) * (6.4e6 * 1*np.pi/180.0), axis=2) / 1e6
           
    print("max chi: ", np.amax(_data["chi_max"])) 
    print("min chi: ", np.amax(_data["chi_min"])) 

res_threshold = 1e-12#1.0 * (9.81*2e-4) / (86400*360*10)
print("res_threshold = %.3e" % (res_threshold,))

y_ind = 0
if args.metric == "20N":
    y_ind = np.argmin(np.abs(coor["y_T"] - 20.0))
elif args.metric == "40N":
    y_ind = np.argmin(np.abs(coor["y_T"] - 40.0))
elif args.metric == "60N":
    y_ind = np.argmin(np.abs(coor["y_T"] - 60.0))


# look for proper z coordinate
z_ind = np.argmin(np.abs(coor["z_W"] - 1000.0))

print("z_ind = %d, y_ind = %d" % (z_ind, y_ind)) 

for d in data:
    d["Q"] /= 1e6
    d["Psib"] /= 1e6

    if args.metric == "max":
       d["psi_measure"] = np.amax(d["Psib"][:, :, :], axis=(1,2))
    elif args.metric in ["20N", "40N", "60N"]:
        #psi = np.amax(d["Psib"][:, :, y_ind], axis=1)
       d["psi_measure"] = d["Psib"][:, y_ind, z_ind]
    elif args.metric == "mean":
       d["psi_measure"] = np.mean(d["Psib"], axis=(1,2))

# Compute the deep water formation grids

def computeDeepWaterFormationGrids(b):

    flags = b * 0
    b_lower = np.roll(b, -1, axis=1)
    
    flags[b_lower > b] = 1.0
    flags[:, -1] = flags[:, -2]
    
    return flags

for d in data:

    Q = d["Q"]
    be = d["be"]
    bw = d["bw"]
   
    dwf_east = be * 0.0
    dwf_west = be * 0.0
    for j in range(len(Q)-1):

        dwf_east[j, :, :] = computeDeepWaterFormationGrids(be[j, :, :])
        dwf_west[j, :, :] = computeDeepWaterFormationGrids(bw[j, :, :])

    d["dwf_east"] = dwf_east
    d["dwf_west"] = dwf_west

    #print(np.sum(dwf_east))

# Pick out the marks
def dist(x1,y1,x2,y2):
    return ((x1-x2)**2 + (y1-y2)**2)**0.5

nmarkpairs = int( len(args.marks) / 2 )
shortest_dist = np.zeros((nmarkpairs,)) + 1e5
mark_index = np.zeros((nmarkpairs, 2), dtype=int) - 1
        
for k in range(nmarkpairs):
    print("Search for the marker case closest to : (gamma, psi) = (%f, %f)" % (args.marks[k*2], args.marks[k*2+1]))
    for i, d in enumerate(data):
        print(len(d["Q"]))
        for s in range(len(d["Q"])):
            _dist = dist(d["Q"][s], d["psi_measure"][s], args.marks[k*2], args.marks[k*2+1])
            if _dist < shortest_dist[k]:
                shortest_dist[k] = _dist
                mark_index[k, 0] = i
                mark_index[k, 1] = s
    print("Found: (%d, %d)" % (mark_index[k, 0], mark_index[k, 1]))



# Decide color and legend
legend_coor = False
if len(args.legend_coor) != 0:
    legend_coor = True
    print("Legend will be plotted using text with the specified coordinate in --legend-coor. Color will be set to black.")

                
print("Data loaded. Plotting now...")

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

plot_z_W = coor["z_W"] / 1e3
plot_z_T = coor["z_T"] / 1e3

fig, ax = plt.subplots(1, 1, figsize=(6, 6), squeeze=False, constrained_layout=True)

for i in range(len(legends)):
    
    print("Plotting %s (%s)" % (legends[i], folders[i]))

    d = data[i]
    offset = 0

    vals, rngs = lsd.detectRanges(d["stable"])

    for j, val in enumerate(vals):
        
        rng = rngs[j]

        Q   = d["Q"][rng]
        res = d["res"][rng]
        psi_measure = d["psi_measure"][rng]

        linestyle = "solid" if val == 1.0 else "dashed"
 
        kw_label = {
            "color" : colors[i],
        }

        if j == 0:
            kw_label['label'] = legends[i]


        if legend_coor == True:
        
            _x = args.legend_coor[i*2]
            _y = args.legend_coor[i*2+1]
            ax[0, 0].text(_x, _y, legends[i], color=colors[i], size=15, ha="center", va="center")
   

        ax[0, 0].plot(Q, psi_measure, zorder=1, linestyle=linestyle, **kw_label)




for k in range(nmarkpairs):
    print("Marking the %d-th marker." % k)
    print(mark_index[k])
    d = data[mark_index[k,0]]
    s = mark_index[k, 1]
    ax[0,0].scatter(d["Q"][s], d["psi_measure"][s], s=80, marker="*", c="yellow", edgecolor="k", zorder=50)
    ax[0,0].text(d["Q"][s], d["psi_measure"][s] + 0.08, "123456789"[k], size=15, va="bottom", ha="center")

# if True then it is labeled directly on the specified coordinate
if legend_coor == False:
    ax[0, 0].legend()

if args.metric == "30N":
    ax[0, 0].set_ylabel(r"$\Psi$ at $30^{\circ}\mathrm{N}$ 1000m depth [Sv]")

elif args.metric == "max":
    ax[0, 0].set_ylabel(r"$\Psi$ maximum [Sv]")

ax[0, 0].set_ylim(args.psi_rng)
    
ax[0, 0].set_xlabel("$\\gamma$ [Sv]")
ax[0, 0].set_xlim(args.F_rng)

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
        b_concat = np.concatenate((b_concat, d["bw"][s, :, :].flatten()))
        b_concat = np.concatenate((b_concat, d["be"][s, :, :].flatten()))

    b_mean = np.mean(b_concat)
    b_std  = np.std(b_concat)

    nstd = 2
    b_neworigin = b_mean - nstd * b_std
    b_rng = np.array([0.0, 2 * nstd * b_std])

    print("b_neworigin: %f" % (b_neworigin,))

    # need b_factor because later on all the buoyancy data will be multiplied by b_factor
    b_rng[0] = np.floor(b_rng[0] * b_factor / 2) * 2
    b_rng[1] = np.floor(b_rng[1] * b_factor / 2) * 2


    fig2, ax2 = plt.subplots(3, nmarkpairs, sharey=False, sharex=True, figsize=(4*nmarkpairs, 4*3), gridspec_kw={'bottom':0.2}, squeeze=False)
   
    levels_b = np.arange(b_rng[0], b_rng[1], 2)
    #levels_b = np.arange(.24, .28, 0.005)
    levels_ui = np.arange(-100, 100, 0.5)
    levels_we = np.arange(-100, 100, .1)
    levels_ww = np.arange(-100, 100, .5)
    levels_psi = np.arange(-20, 20, 1)
    dwf_color = (0.6, 0.75, 1.0)
 
    for k in range(nmarkpairs):
      
        be_plot = ( d["be"][s, :, :].transpose() - b_neworigin ) * b_factor 
        bw_plot = ( d["bw"][s, :, :].transpose() - b_neworigin ) * b_factor 

        print("Max and min of be: %.2e ; %.2e" % (np.amax(be_plot), np.amin(be_plot)))
        print("Max and min of bw: %.2e ; %.2e" % (np.amax(bw_plot), np.amin(bw_plot)))

        print("Max and min of bw: %.2e ; %.2e" % (np.amax(d["bw"]), np.amin(d["bw"])))
        #print("Max and min of ui: %.2e ; %.2e" % (np.amax(d["ui"]), np.amin(d["ui"])))
        #print("Max and min of we: %.2e ; %.2e" % (np.amax(d["we"]), np.amin(d["we"])))
        #print("Max and min of ww: %.2e ; %.2e" % (np.amax(d["ww"]), np.amin(d["ww"])))

        _ax = ax2[:, k]
        print("Plotting the %d-th marker streamfunction" % k)
        d = data[mark_index[k,0]]
        s = mark_index[k, 1]
        
        CS = _ax[0].contour(coor["y_V"], plot_z_W, d["Psib"][s, :, :].transpose(), levels_psi, colors="black")
        _ax[0].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=4)
        
           
        #CS = _ax[1].contour(coor["y_T"], plot_z_W, d["ww"][s, :, :].transpose() * 1e5, levels_ww, colors="red")
        #_ax[1].clabel(CS, CS.levels, inline=True, fmt="%.1f", fontsize=10, inline_spacing=2)

        CS = _ax[1].contour(coor["y_T"], plot_z_T, bw_plot, levels_b, colors="black")
        _ax[1].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=2)
        cs_dwf_west = _ax[1].contourf(coor['y_T'], plot_z_T, d['dwf_west'][s, :, :].transpose(), [0, 0.5,1.5], colors=['none', dwf_color])


        #CS = _ax[2].contour(coor["y_T"], plot_z_W, d["we"][s, :, :].transpose() * 1e5, levels_we, colors="red")
        #_ax[2].clabel(CS, CS.levels, inline=True, fmt="%.1f", fontsize=10, inline_spacing=2)
        
        CS = _ax[2].contour(coor["y_T"], plot_z_T, be_plot, levels_b, colors="black")
        _ax[2].clabel(CS, CS.levels, inline=True, fmt="%d", fontsize=10, inline_spacing=2)
        cs_dwf_east = _ax[2].contourf(coor['y_T'], plot_z_T, d['dwf_east'][s, :, :].transpose(), [0, 0.5,1.5], colors=['none', dwf_color])


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
            __ax.set_yticks([0, 1, 2, 3, 4])

            if k != 0:
                __ax.set_yticklabels([""] * len(__ax.get_yticks()))

        _ax[0].set_ylim([4.5, 0])
        _ax[1].set_ylim([1.0, 0])
        _ax[2].set_ylim([1.0, 0])

    if args.output_marks != "":
        fig2.savefig(args.output_marks, dpi=300)



plt.show()









