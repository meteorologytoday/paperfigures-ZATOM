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


def findSamplingIdxBySpacing(s, spacing):
 
    if not np.all(np.isfinite(s)):
        raise Exception("The array `s` is not all finite")
  
    ds = s[1:] - s[:-1]

    if np.any(ds <= 0):
        raise Exception("The array `s` is not monotonically increasing")

 
    idx = []
    anchored_idx = 0

    for i in range(len(s)):
        
        ds = s[i] - s[anchored_idx]
        
        if ds >= spacing:

            idx.append(i)
            anchored_idx = i

    return np.array(idx)        
    
    


def repNaN2None(arr):
    for i, elm in enumerate(arr):
        if np.isnan(elm):
            arr[i] = None

parser = argparse.ArgumentParser()
parser.add_argument('--folder', required=True)
parser.add_argument('--title', type=str, default="")
parser.add_argument('--output', default="")
parser.add_argument('--residue-threshold', type=float, default=1e-12)

parser.add_argument('--sampling-spacing', type=float, help="Sampling spacing after scaled", default=1.0)

parser.add_argument('--param', type=str, choices=["gamma", "xi"])
parser.add_argument('--param-scale', type=float, help="Scale of gamma in Sv", required=True)
parser.add_argument('--param-rng', nargs=2, type=float, default=None)

parser.add_argument('--plot-param-rng', nargs=2, type=float, default=[None, None])
parser.add_argument('--plot-psi-rng', nargs=2, type=float, default=[None, None])
parser.add_argument('--plot-chi-rng', nargs=2, type=float, default=[None, None])

parser.add_argument('--psi-scale', type=float, help="Scale of psi in Sv", required=True)
parser.add_argument('--mode', type=int, help="Mode to compute.", default=1)

parser.add_argument('--label-offset', type=float, help="Offset in Sv", default=0.05)
parser.add_argument('--label-sides', nargs='*', type=str, default=[], choices=["R", "L"])
parser.add_argument('--label-idx', type=int, nargs="*", help="Position to label", default=[])
parser.add_argument('--labels', type=str, nargs="*", help="Labels")

parser.add_argument('--suptitle', type=str, help="Suptitle", default="")

parser.add_argument('--dont-plot-psi', action="store_true", )
parser.add_argument('--dont-plot-chi', action="store_true", )
parser.add_argument('--dont-plot-dpsidt-decomp', action="store_true", )
parser.add_argument('--dont-plot-dchidt-decomp', action="store_true", )
parser.add_argument('--thumbnail-numbering', type=str, default="abcdefghijklmn")

parser.add_argument('--nrow', type=int, default=1)


parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)


coor = {}

folder = args.folder

if args.param == "gamma":
    print("Using freshwater forcing as bifurcation parameter.")
    param = "Q"
    param_factor = 1e-6

elif args.param == "xi":
    print("Using zonal asymmetry as bifurcation parameter.")
    param = "xi"
    param_factor = 1.0


# Data contains a list of folders (cases)
# and each case has multiple data files
loaded_varnames = ["Psib", "we", "ww", "vw", "chi", "be", "bw", "qw", "qe", "res", "stable", "ui", param,] + med.necessary_variables

print("Loading the folder: %s" % (folder,))
data, coor = lsd.loadScanData(
    folder,
    loaded_varnames,
    load_coor=True,
    residue_threshold = args.residue_threshold,
)
    
print("Number of records: %d" % (len(data[param]))) 
med.makeExtendedData(data, coor, lat_s=40.0, lat_n=70.0, merge=True, verbose=True, mode=args.mode)

"""
print("Finding the indexes to to subsampling")
param_tmp = data[param] / (1e6*args.param_scale)
psi_tmp = data["mode_psi"] / (1e6*args.psi_scale)

ds = np.zeros_like(param_tmp)
ds[1:] = ((param_tmp[1:] - param_tmp[:-1])**2 + (psi_tmp[1:] - psi_tmp[:-1])**2)**0.5
s = np.cumsum(ds)

idx = findSamplingIdxBySpacing(s, args.sampling_spacing) 

param_sub = data["Q"][idx]
#idx = idx[ (Q_sub >= args.gamma_rng[0]) & (Q_sub <= args.gamma_rng[1]) ]

data["subsampling_idx"] = idx
"""

# Remove data out of range of interest
if args.param_rng is not None:
    new_data = dict()
    sub_idx = (data[param]*param_factor >= args.param_rng[0]) & (data[param]*param_factor <= args.param_rng[1])
    for k, arr in data.items():
        new_data[k] = arr[sub_idx]

    data = new_data 
        

vals, rngs = lsd.detectRanges(data["stable"])



plot_infos = dict(

    mode_dbdiffdt_dueto_ZOC_T = dict(
        color = "",
        label = r"adv T ZOC",
    ),
    mode_dbdiffdt_dueto_MOC_T = dict(
        color = "",
        label = r"adv T MOC",
    ),

    mode_dbdiffdt_dueto_ZOC_S = dict(
        color = "",
        label = r"adv S ZOC",
    ),
    mode_dbdiffdt_dueto_MOC_S = dict(
        color = "",
        label = r"adv S MOC",
    ),

    mode_chi = dict(
        color = "",
        label = r"$\left\langle \chi \right\rangle$",
    ),
    
    mode_psi = dict(
        color = "",
        label = r"$\left\langle \Psi \right\rangle$",
    ),
    
    mode_dbdiffdt_dueto_ZOCMOC = dict(
        color = "",
        label = r"ZOCMOC",
    ),
    
    mode_dbdiffdt_dueto_ADV_ZOC = dict(
        color = "",
        label = r"$\dot{{\tilde{{\psi}}^{mode:d}_{{\mathrm{{ZOCA}}}}$",
    ),
    mode_dbdiffdt_dueto_ADV_MOC = dict(
        color = "",
        label = r"$\dot{{\tilde{{\psi}}^{mode:d}_{{\mathrm{{MOCA}}}}$",
    ),

    mode_dbdiffdt_dueto_SS = dict(
        color = "",
        label = r"$\left\langle \dot{\psi} \right\rangle_{\mathrm{FWF}}$",
    ),

    mode_dbdiffdt_dueto_FRC = dict(
        color = "",
        label = r"$\left\langle \dot{\psi} \right\rangle_{\mathrm{TFRC}}$",
    ),

    mode_dbdiffdt_dueto_CVA = dict(
        color = "",
        label = r"$\left\langle \dot{\psi} \right\rangle_{\mathrm{CVA}}$",
    ),

    mode_dbdiffdt_dueto_BGDIFU = dict(
        color = "",
        label = r"$\left\langle \dot{\psi} \right\rangle_{\mathrm{BGDF}}$",
    ),
 
    mode_dbdiffdt_dueto_ZNLHDIFU = dict(
        color = "",
        label = r"$\left\langle \dot{\psi} \right\rangle_{\mathrm{BDYDF}}$",
    ),
 
    mode_dbdiffdt_sum = dict(
        color = "",
        label = r"SUM",
    ),

    mode_dchidt_dueto_ADV_ZOC = dict(
        color = "",
        label = r"$\left\langle \dot{\chi} \right\rangle_{\mathrm{ZOCA}}$",
    ),
    mode_dchidt_dueto_ADV_MOC = dict(
        color = "",
        label = r"$\left\langle \dot{\chi} \right\rangle_{\mathrm{MOCA}}$",
    ),

    mode_dchidt_dueto_SS = dict(
        color = "",
        label = r"$\left\langle \dot{\chi} \right\rangle_{\mathrm{FWF}}$",
    ),

    mode_dchidt_dueto_FRC = dict(
        color = "",
        label = r"$\left\langle \dot{\chi} \right\rangle_{\mathrm{TFRC}}$",
    ),

    mode_dchidt_dueto_CVA = dict(
        color = "",
        label = r"$\left\langle \dot{\chi} \right\rangle_{\mathrm{CVA}}$",
    ),

    mode_dchidt_dueto_BGDIFU = dict(
        color = "",
        label = r"$\left\langle \dot{\chi} \right\rangle_{\mathrm{BGDF}}$",
    ),
 
    mode_dchidt_dueto_ZNLHDIFU = dict(
        color = "",
        label = r"$\left\langle \dot{\chi} \right\rangle_{\mathrm{BDYDF}}$",
    ),
 
    mode_dchidt_sum = dict(
        color = "",
        label = r"SUM",
    ),


)


print("Data loaded. Plotting now...")

print("Loading matplotlib...")
import tool_fig_config
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
mplt.rcParams['lines.linewidth'] =   default_linewidth;
mplt.rcParams['axes.linewidth'] =    default_linewidth;
mplt.rcParams['xtick.major.size'] =  default_ticksize;
mplt.rcParams['xtick.major.width'] = default_linewidth;
mplt.rcParams['ytick.major.size'] =  default_ticksize;
mplt.rcParams['ytick.major.width'] = default_linewidth;

rc('font', **{'size': 15.0});
rc('axes', **{'labelsize': 15.0});
rc('mathtext', **{'fontset':'stixsans'});

import colorblind

print("Done")

colors = [
    'black',
    'orange',
    'skyblue',
    "bluishgreen",
    'yellow',
    "blue",
    'vermillion',
    'reddishpurple',
]



ncol = 4
nrow = 1

if args.dont_plot_psi:
    ncol -= 1

if args.dont_plot_chi:
    ncol -= 1

if args.dont_plot_dpsidt_decomp:
    ncol -= 1

if args.dont_plot_dchidt_decomp:
    ncol -= 1


nrow = args.nrow
ncol = int(np.ceil(ncol / args.nrow))



d = data  # shorthand

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 4,
    h = 4,
    wspace = 1.2,
    hspace = 1.5,
    w_left = 2.0,
    w_right = 2.0,
    h_bottom = 2.0,
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

fig.suptitle(args.suptitle)

ax_flat = ax.flatten(order='C')
ax_idx = 0

if not args.dont_plot_psi:

    # Plot bifurcation
    print("Plotting bifurcation diagram")
    _ax = ax_flat[ax_idx] ; ax_idx += 1

    varname = "mode_psi"
    factor = 1e-6
    for j, val in enumerate(vals):
        
        rng = rngs[j]
        x   = d[param][rng] * param_factor
        y   = d[varname][rng] * factor
        res = d["res"][rng]

        linestyle = "solid" if val == 1.0 else "dashed"

        kw_label = {
            "color" : "black",
        }

        _ax.plot(x, y, zorder=1, linestyle=linestyle, **kw_label)

    #subsampling_idx = data["subsampling_idx"]
    #subset_param = data[param][subsampling_idx] * param_factor
    #subset_psi   = data[varname][subsampling_idx] / 1e6 # Sv
    #subset_s     = s[subsampling_idx]
    #_ax.scatter(subset_param, subset_psi, s=10, c="black", marker="o")

    for k, idx in enumerate(args.label_idx):

        subset_param = data[param][idx] * param_factor
        subset_psi   = data[varname][idx] / 1e6 # Sv
        label_side = args.label_sides[k]

        if label_side == "R":
            offset_x = args.label_offset
            ha = "left"
        elif label_side == "L":
            offset_x = - args.label_offset
            ha = "right"

        _ax.scatter(subset_param, subset_psi, s=10, c="black", marker="o")
        _ax.text(subset_param + offset_x, subset_psi, args.labels[k], size=12, va="center", ha=ha)
        
    _ax.set_axisbelow(True)
    _ax.grid()
    _ax.set_ylabel("[Sv]")
    
    _ax.set_title("(%s) $ \\left \\langle \\tilde{\\psi}^{%d} \\right \\rangle$ " % (args.thumbnail_numbering[ax_idx-1], args.mode,))

    _ax.set_ylim(args.plot_psi_rng)

if not args.dont_plot_dpsidt_decomp:

    print("Plotting parametric dpsi/dt diagram")

    factor = 86400*360*1 * 1e-6 # Sv per year

    _ax = ax_flat[ax_idx] ; ax_idx += 1
    for i, (varname_suffix, tendency_name, colorname) in enumerate([
        ("ADV_ZOC", "ZOCA", "black"), 
        ("ADV_MOC", "MOCA", "orange"), 
        ("SS", "FWF", "skyblue"),
        ("FRC", "TFRC", "bluishgreen"),
        ("BGDIFU", "BGDF", "yellow"),
        ("ZNLHDIFU", "BDYDF", "vermillion"),
        ("CVA", "CVA", "reddishpurple"),
    #    "mode_dbdiffdt_sum",
    ]):
        varname = "mode_dbdiffdt_dueto_%s" % (varname_suffix,)
        label = r"$ \dot{\tilde{\psi}^{%d}}_{\mathrm{%s}}$" % (args.mode, tendency_name)
        color = colorblind.BW8color[colorname] if colorname in colorblind.BW8color else colorname
        for j, val in enumerate(vals):
            
            rng = rngs[j]
            #x   = s[rng]
            x   = d[param][rng] * param_factor
            y   = data[varname][rng] * factor

            linestyle = "solid" if val == 1.0 else "dashed"
                
            kw_label = {
                "color" : color,
            }

            if j == 0:
                kw_label["label"] = label

            _ax.plot(x, y, zorder=1, linestyle=linestyle, **kw_label)

        for k, idx in enumerate(args.label_idx):
            subset_param = data[param][idx] * param_factor
            subset_psi   = data[varname][idx] * factor
            _ax.scatter(subset_param, subset_psi, s=10, color=color, marker="o")


    #_trans = mplt.transforms.blended_transform_factory(_ax.transData, _ax.transAxes)
    #for i, _s in enumerate(subset_s):
    #    _ax.plot([_s, _s], [0, 1],color="gray", linestyle="solid", transform=_trans)
    _ax.set_axisbelow(True)
    _ax.grid()
    _ax.legend(bbox_to_anchor=(1.02, 1.0), ncol=1, loc='upper left', fontsize=12)
    _ax.set_ylabel("[$ \\mathrm{Sv} \\cdot \\mathrm{year}^{-1}$]")

    """
    _ax_twinx = _ax.twinx()
    varname = "mode_chi"
    plot_info = plot_infos[varname]
    label = plot_info["label"]
    for j, val in enumerate(vals):
        
        rng = rngs[j]
        x   = s[rng]
        y   = data[varname][rng]

        linestyle = "solid" if val == 1.0 else "dashed"

        kw_label = {
            "color" : "black",
        }

        if j == 0:
            kw_label["label"] = label

        _ax_twinx.plot(x, y, zorder=1, linestyle=linestyle, **kw_label, linewidth=2)
    """

    _ax.set_title("(%s) Decomp of $  \\left \\langle \\partial \\tilde{\\psi}^{%d} / \\partial t \\right \\rangle$" % (args.thumbnail_numbering[ax_idx-1], args.mode, ))


if not args.dont_plot_chi:
    
    print("Plot Chi")

    _ax = ax_flat[ax_idx] ; ax_idx += 1
    for i, (varname, factor) in enumerate([
        ("mode_chi", 1e-6 * 1e6) # Sv / 1000km
    ]):

        plot_info = plot_infos[varname]
        label = plot_info["label"]
        for j, val in enumerate(vals):
            
            rng = rngs[j]
            #x   = s[rng]
            x   = d[param][rng] * param_factor
            y   = data[varname][rng] * factor

            linestyle = "solid" if val == 1.0 else "dashed"
            color = colorblind.BW8color[colorname] if (colorname := colors[i]) in colorblind.BW8color else colorname
            kw_label = {
                "color" : color,
            }

            if j == 0:
                kw_label["label"] = label

            _ax.plot(x, y, zorder=1, linestyle=linestyle, **kw_label)

        for k, idx in enumerate(args.label_idx):
            subset_param = data[param][idx] * param_factor
            subset_psi   = data[varname][idx] * factor
            _ax.scatter(subset_param, subset_psi, s=10, color=color, marker="o")


    _ax.set_ylabel("[$ \\mathrm{Sv} \\cdot \\left( 1000 \\, \\mathrm{km} \\right)^{-1}$]")

    #_trans = mplt.transforms.blended_transform_factory(_ax.transData, _ax.transAxes)
    #for i, _s in enumerate(subset_s):
    #    _ax.plot([_s, _s], [0, 1],color="gray", linestyle="solid", transform=_trans)

    #_ax.legend()
    _ax.set_axisbelow(True)
    _ax.grid()

    _ax.set_title("(%s) $ \\left\\langle \\chi \\right\\rangle$ " % (args.thumbnail_numbering[ax_idx-1],))
    _ax.set_title("(%s) $ \\left \\langle \\tilde{\\chi}^{%d} \\right \\rangle$" % (args.thumbnail_numbering[ax_idx-1], args.mode, ))
    
    _ax.set_ylim(args.plot_chi_rng)




if not args.dont_plot_dchidt_decomp:

    print("Plotting parametric dchi/dt diagram")

    factor = 86400*360*1 * 1e-6 * 1e6 # Sv per year per 1000km

    _ax = ax_flat[ax_idx] ; ax_idx += 1
    for i, (varname_suffix, tendency_name, colorname) in enumerate([
        ("ADV_ZOC", "ZOCA", "black"), 
        ("SS", "FWF", "skyblue"),
        ("FRC", "TFRC", "bluishgreen"),
        ("BGDIFU", "BGDF", "yellow"),
        ("ZNLHDIFU", "BDYDF", "vermillion"),
        ("CVA", "CVA", "reddishpurple"),
    ]):
 
        varname = "mode_dchidt_dueto_%s" % (varname_suffix,)
        label = r"$ \dot{\tilde{\chi}^{%d}}_{\mathrm{%s}}$" % (args.mode, tendency_name)

        color = colorblind.BW8color[colorname] if colorname in colorblind.BW8color else colorname
        for j, val in enumerate(vals):
            
            rng = rngs[j]
            #x   = s[rng]
            x   = d[param][rng] * param_factor
            y   = data[varname][rng] * factor

            linestyle = "solid" if val == 1.0 else "dashed"
                
            kw_label = {
                "color" : color,
            }

            if j == 0:
                kw_label["label"] = label

            _ax.plot(x, y, zorder=1, linestyle=linestyle, **kw_label)

        for k, idx in enumerate(args.label_idx):
            subset_param = data[param][idx] * param_factor
            subset_psi   = data[varname][idx] * factor
            _ax.scatter(subset_param, subset_psi, s=10, color=color, marker="o")


    #_trans = mplt.transforms.blended_transform_factory(_ax.transData, _ax.transAxes)
    #for i, _s in enumerate(subset_s):
    #    _ax.plot([_s, _s], [0, 1],color="gray", linestyle="solid", transform=_trans)
    _ax.set_axisbelow(True)
    _ax.grid()
    _ax.legend(bbox_to_anchor=(1.02, 1.0), ncol=1, loc='upper left', fontsize=12)
    _ax.set_ylabel("[$ \\mathrm{Sv} \\cdot \\mathrm{year}^{-1} \\cdot \\left(1000 \\, \\mathrm{km}\\right)^{-1}$]")
     
    _ax.set_title("(%s) Decomp of $ \\left \\langle \\partial \\tilde{\\chi}^{%d} / \\partial t \\right \\rangle $" % (args.thumbnail_numbering[ax_idx-1], args.mode, ))


for _ax in ax_flat:
    _ax.set_xlabel("$ \\gamma $ [Sv]")
    _ax.set_xlim(args.plot_param_rng)
    _ax.set_xticks(np.arange(0, 0.125, 0.02))


if args.output != "":
    print("Saving output: ", args.output)
    fig.savefig(args.output, dpi=300)

if not args.no_display:
    print("Showing figures...")
    plt.show()

