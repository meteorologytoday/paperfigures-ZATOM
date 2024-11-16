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
parser.add_argument('--folder', nargs='+')
parser.add_argument('--legend', nargs='*')
parser.add_argument('--title', type=str, nargs='*', default=None)
parser.add_argument('--colors', nargs='*')
parser.add_argument('--auto-color', action="store_true")
parser.add_argument('--no-legend', action="store_true")
parser.add_argument('--output-bifur', default="")
parser.add_argument('--output-marks', default="")
parser.add_argument('--offset-marks', type=float, default=0)
parser.add_argument('--legend-loc', type=str, default="best")
parser.add_argument('--varnames', type=str, nargs="+", required=True)
parser.add_argument('--param', type=str, choices=["gamma", "xi"])
parser.add_argument('--param-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--psi-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--mode-psi-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--s1000-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--db_ns-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--db_ew-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--cvt_e-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--cvt_w-rng', nargs=2, type=float, default=[np.nan, np.nan])
parser.add_argument('--coe-gamma-xi', type=float, default=9.57e-19)
parser.add_argument('--marks', nargs='*', type=float, default=[])
parser.add_argument('--marks-pos', nargs='*', type=str, default=[])
parser.add_argument('--mark-labels', nargs='*', type=str, default=[])
parser.add_argument('--mark-sides', nargs='*', type=str, default=[], choices=["R", "L"])
parser.add_argument('--residue-threshold', type=float, default=1e-12)
parser.add_argument('--no-display', action="store_true")
parser.add_argument('--ncol', type=int, default=1)
parser.add_argument('--thumbnail-skip', type=int, default=0)
parser.add_argument('--text', nargs='*', type=str, default=[])
parser.add_argument('--text-pos', nargs='*', type=float, default=[])
parser.add_argument('--text-ax', nargs='*', type=int, help="Specify which axes to put text", default=[])
parser.add_argument('--put-var-on-yaxis', action="store_true")
parser.add_argument('--mode', type=int, help="Mode to compute.", default=1)


args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)

repNaN2None(args.param_rng)
repNaN2None(args.psi_rng)
repNaN2None(args.s1000_rng)
repNaN2None(args.db_ns_rng)
repNaN2None(args.db_ew_rng)
repNaN2None(args.cvt_w_rng)
repNaN2None(args.cvt_e_rng)

data = []
coor = {}

#target_vars = ["mode_psi", "mode_db_ew", "mode_s", "mode_chi", "mode_dq", "mode_chi_dbdz"]
#target_vars = ["mode_psi", "mode_chi_dbdz", "mode_chi_dbdz_product", "mode_chi", "mode_s_eff", "mode_dq"]
target_vars = args.varnames 

#["mode_psi", "mode_chi_dbdz", "mode_dq",]# "mode_chi_dbdz_product", "mode_chi", "mode_s_eff",]

folders = args.folder
legends   = args.legend

if args.no_legend:
    print("The flag `--no-legend ` is set.")
else:
    if (legends is None) or (legends is not None and len(legends) < len(folders)):
        print("Legends do not have the same length of folders. Use numbering instead")
        legends = ["%d" % (i,) for i in range(len(folders))]

    print("Legends are")
    print(legends)

if len(args.marks) % 2 != 0:
    raise Exception("Length of mark is not multiple of 2. It has to be a list of (x, y) pairs. ")

if len(args.text_pos) != 2*len(args.text):
    raise Exception("Length of `--text-pos` does not match length of `--text` * 2. ")

if len(args.text_ax) != len(args.text):
    raise Exception("Length of `--text-ax` does not match length of `--text`. ")


#
if args.param == "gamma":
    print("Using freshwater forcing as bifurcation parameter.")
    param = "Q"

elif args.param == "xi":
    print("Using zonal asymmetry as bifurcation parameter.")
    param = "xi"


# Data contains a list of folders (cases)
# and each case has multiple data files
data_to_delete = []
data = []
coor = None
loaded_varnames = ["Psib", "chi", param, "be", "bw", "qw", "qe", "res", "stable", "ui"] + med.necessary_variables
for i, folder in enumerate(folders):

    print("Loading the folder: %s" % (folder,))

    try:
        _data, _coor = lsd.loadScanData(
            folder,
            loaded_varnames,
            load_coor=(coor is None), 
            residue_threshold = args.residue_threshold,
        )
    except Exception as e:
        print("Error occurs. Skip this one.")
        data_to_delete.append(i)
        data.append(None)
        continue
 
    print("Number of records: %d" % (len(_data[param]))) 
    if coor is None:
        coor = _coor

    med.makeExtendedData(_data, coor, lat_s=40.0, lat_n=70.0, merge=True, verbose=True, mode=args.mode)

    data.append(_data)
    #print("mode_dq = ", _data["mode_dq"])
    
data_to_delete.reverse() # important. delete from last to first to avoid reordering
print("Delete index: ", data_to_delete)

for rm_idx in data_to_delete:
    del folders[rm_idx]
    del data[rm_idx]
    del legends[rm_idx]

    
for d in data:
   
    if param == "Q": 
        d[param] /= 1e6

    d["mode_psi"] /= 1e6 

    y_ind = np.argmin(np.abs(coor["y_T"] - 60.0))
    z_ind = np.argmin(np.abs(coor["z_W"] - 1000.0))


# Pick out the marks
def dist(x1,y1,x2,y2):
    return ((x1-x2)**2 + (y1-y2)**2)**0.5

nmarkpairs = int( len(args.marks) / 2 )
shortest_dist = np.zeros((nmarkpairs,)) + 1e5
mark_index = np.zeros((nmarkpairs, 2), dtype=int) - 1
        
for k in range(nmarkpairs):
    print("Search for the marker case closest to : (%s, psi) = (%f, %f)" % (param, args.marks[k*2], args.marks[k*2+1]))
    for i, d in enumerate(data):
        for s in range(len(d[param])):
            _dist = dist(d[param][s]*10, d["mode_psi"][s], args.marks[k*2]*10, args.marks[k*2+1])
            if _dist < shortest_dist[k]:
                shortest_dist[k] = _dist
                mark_index[k, 0] = i
                mark_index[k, 1] = s
    
    print("Found: (%d, %d)" % (mark_index[k, 0], mark_index[k, 1]))



          
print("Data loaded. Plotting now...")

import tool_fig_config

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

import colorblind
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

if args.auto_color:

    print("Option --auto-color is on. Use self-generated colors.")
    #cmap = plt.cm.get_cmap("cmo.phase", len(legends))
    cmap = mplt.colors.LinearSegmentedColormap.from_list('my_cmap', [
        '#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442'
#        (255, 0, 0),
#        colorblind.BW8color["vermillion"],
#        colorblind.BW8color["orange"],
#        colorblind.BW8color["reddishpurple"],
#        colorblind.BW8color["skyblue"],
#        colorblind.BW8color["blue"],
    ])
    
    colors = [ cmap(i/len(args.folder)) for i in range(len(args.folder)) ]
    
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
         

plot_z_W = - coor["z_W"] / 1e3
plot_z_T = - coor["z_T"] / 1e3

ncol = args.ncol
nrow = int(np.ceil(len(target_vars) / ncol)) 

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 4,
    h = 4,
    wspace = 1.0,
    hspace = 1.5,
    w_left = 1.0,
    w_right = 1.0,
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

#fig, ax = plt.subplots(2, 3, figsize=(12, 8), squeeze=False, constrained_layout=True)

ax_flat = ax.flatten(order='C')
       


for i in range(len(args.folder)):
    
    print("Plotting %s" % (folders[i]))

    d = data[i]
    offset = 0

    vals, rngs = lsd.detectRanges(d["stable"])

    for j, val in enumerate(vals):
        
        rng = rngs[j]

        x   = d[param][rng]
        res = d["res"][rng]


        for k, var in enumerate(target_vars):
           
            #print(var, "=", d[var]) 
            y   = d[var][rng]

            linestyle = "solid" if val == 1.0 else "dashed"
 
            kw_label = {
                "color" : colors[i],
            }

            if (not args.no_legend) and j == 0 and k == 0:
                kw_label['label'] = legends[i]

            ax_flat[k].plot(x, y, zorder=1, linestyle=linestyle, **kw_label)


for k in range(len(args.text)):

    print("Marking the %d-th text." % k)
    _ax = ax_flat[args.text_ax[k]]
    x = args.text_pos[2*k+0]
    y = args.text_pos[2*k+1]
    s = args.text[k]
    _ax.text(x, y, s, size=12, va="center", ha="center")
        



for k in range(nmarkpairs):
    print("Marking the %d-th marker." % k)
    print(mark_index[k])
    d = data[mark_index[k,0]]
    s = mark_index[k, 1]

    if args.mark_sides[k] == "R":
        offset_x = args.offset_marks
        ha = "left"
    elif args.mark_sides[k] == "L":
        offset_x = - args.offset_marks
        ha = "right"

    for l, var in enumerate(target_vars):


        #ax_flat[l].scatter(d["Q"][s], d[var][s], s=10, marker="o", c="yellow", edgecolor="k", zorder=50)
        ax_flat[l].scatter(d[param][s], d[var][s], s=10, marker="o", c="k", edgecolor="k", zorder=50)
        ax_flat[l].text(d[param][s] + offset_x, d[var][s], args.mark_labels[k], size=12, va="center", ha=ha)
        

if not args.no_legend:
    ax[0, 0].legend(fontsize=12, handlelength=1.0, labelspacing=0.25, borderpad=0.2, loc=args.legend_loc)


labels = {
    "psi_fixed" : r"$\overline{\left\langle\psi\right\rangle}$",
    "mode_psi" : r"$ \left \langle \tilde{\psi}^{1} \right \rangle $",
    "mode_chi"   : r"$\left\langle\chi\right\rangle$",
    "mode_s_eff"  : r"$\left\langle\partial_z \overline{b}^{\mathrm{eff}} \right\rangle$",
    "s1000_hlat"  : r"$s_{\mathrm{hlat}}$",
    "db_ns"  : r"$\left\langle\partial_y b_e^* \right\rangle$",
    "mode_db_ew"  : r"$\left\langle b_e^* - b_w^* \right\rangle$",
    "cvt_e"  : r"$\tilde{q}_e$",
    "cvt_w"  : r"$\tilde{q}_w$",
    "d_cvt"  : r"$\Delta \tilde{q}$",
    "mode_dq"     : r"$\left\langle \Delta q \right\rangle$",
    "mode_chi_dbdz"  : r"$\left\langle\chi \partial_z \overline{b} \right\rangle$",
    "mode_chi_dbdz_product"  : r"$\left\langle\chi\right\rangle \left\langle\partial_z \overline{b} \right\rangle$",
    "mode_ui_adv"  : r"$\left\langle u_i \frac{b_e - b_w}{L_b} \right\rangle$",
    "mode_ZOC"     : r"$\left\langle\partial_z \overline{b}^{\mathrm{eff}} \right\rangle + \left\langle u_i \frac{b_e - b_w}{L_b} \right\rangle$",
}

units = {
    "psi_fixed"   : "[Sv]",
    "mode_chi"     : r"[$ \mathrm{Sv} / \left( 1000 \mathrm{km}\right)$]",
    "mode_psi"     : "[Sv]",
    "mode_ui_adv"  : r"[Sv]",
    "mode_ZOC"     : r"[Sv]",

    "mode_s_eff"   : r"[ $ \times 10^{-5} \mathrm{s}^{-2} $]",
    "s1000_hlat"  : r"[ $ \times 10^{-5} \mathrm{s}^{-2} $]",
    "db_ns"       : r"[ $ \times 10^{-3} \mathrm{m} / \mathrm{s}^{2} $]",
    "mode_db_ew"       : r"[ $ \times 10^{-3} \mathrm{m} / \mathrm{s}^{2} $]",
    "cvt_e"       : r"",
    "cvt_w"       : r"",
    "d_cvt"       : r"",
    "mode_dq"                : r"[ $ \times 10^{-12} \mathrm{m} / \mathrm{s}^{2} $]",
    "mode_chi_dbdz"          : r"[ $ \times 10^{-12} \mathrm{m}^2 / \mathrm{s}^{3} $]",
    "mode_chi_dbdz_product"  : r"[ $ \times 10^{-12} \mathrm{m}^2 / \mathrm{s}^{3} $]",

}

for l, var in enumerate(target_vars):
    
    if args.put_var_on_yaxis:
        ax_flat[l].set_ylabel("%s %s" % (labels[var], units[var],))
    else:
        ax_flat[l].set_ylabel("%s" % (units[var],))
    
    ylim_attr = "%s_rng" % (var,)
    if hasattr(args, ylim_attr):
        ax_flat[l].set_ylim(getattr(args, ylim_attr))

    if args.param == "gamma":
        ax_flat[l].set_xlabel("$\\gamma$ [Sv]")
    elif args.param == "xi":
        ax_flat[l].set_xlabel("$\\xi$")


    ax_flat[l].set_xlim(args.param_rng)

    if args.title is not None:
        _title = args.title[l]
    else:
        _title = " %s" % (labels[var],)

    ax_flat[l].set_title("(%s)%s" % ("abcdefg"[l+args.thumbnail_skip], _title,))

    ax_flat[l].grid()

if args.output_bifur != "":
    print("Saving output: ", args.output_bifur)
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

    ncol = nmarkpairs
    nrow = 3 
    figsize, gridspec_kw = tool_fig_config.calFigParams(
        w = 4,
        h = 4,
        wspace = 1.0,
        hspace = 0.7,
        w_left = 1.0,
        w_right = 1.0,
        h_bottom = 1.0,
        h_top = 1.0,
        ncol = ncol,
        nrow = nrow,
    )

    fig2, ax = plt.subplots(
        nrow, ncol,
        figsize=figsize,
        subplot_kw=dict(aspect="auto"),
        gridspec_kw=gridspec_kw,
        constrained_layout=False,
        squeeze=False,
        sharex=False,
    )
   
    #fig2 = plt.figure(figsize=(4*nmarkpairs+1, 4*3)) 
    #gs0 = gridspec.GridSpec(3, nmarkpairs+1, figure=fig, width_ratios=nmarkpairs * [1] + [0.05], wspace=0.3)
    #fig2, ax2 = plt.subplots(3, nmarkpairs, sharey=False, sharex=True, figsize=(4*nmarkpairs, 4*3), gridspec_kw={'bottom':0.2}, squeeze=False) 
 
    for k in range(nmarkpairs):
 
        print("Plotting the %d-th marker streamfunction" % k)
        d = data[mark_index[k,0]]
        s = mark_index[k, 1]
    
        #_ax = [ fig2.add_subplot(gs0[i, k]) for i in range(3) ]
        _ax = ax[:, k]#[ fig2.add_subplot(gs0[i, k]) for i in range(3) ]
      
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
 


        if param == "Q":
            _ax[0].set_title("State %s ($%.2f \\, \\mathrm{Sv}$)" % ("123456789"[k], d[param][s]))
        elif param == "xi":
            _ax[0].set_title("State %s ($%.2f$)" % ("123456789"[k], d[param][s]))

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

        #_ax[0].set_ylim([4.5, 0])
        #_ax[1].set_ylim([1, 0])
        #_ax[2].set_ylim([1, 0])

        
    #cax = fig2.add_subplot(gs0[:, -1])
    #plt.colorbar(mappable=b_mappable, cax=cax, cmap=b_cmap, orientation="vertical", ticks=[], label="$b_w^*$ and $b_e^*$ [$\\mathrm{m}^2 / \\mathrm{s}$]")
    
    if args.output_marks != "":
        fig2.savefig(args.output_marks, dpi=300)

if not args.no_display:
    print("Showing figures...")
    plt.show()

