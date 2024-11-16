from pathlib import Path
import os,sys, argparse
import xarray as xr
import numpy as np
import pprint
import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output',    default="")
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)

coor = {}
data = {}

ds = xr.open_dataset(args.input, decode_timedelta=False)

def normalizeTo01(d):

    M = np.amax(d)
    m = np.amin(d)

    d = (d - m) / (M - m)

    return d


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

import cmocean as cmo
import colorblind
print("Done")

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

ncol = 1
nrow = 1

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

print("Creating canvas...")
fig, ax = plt.subplots(
    nrow, ncol,
    figsize=figsize,
    subplot_kw=dict(aspect="auto"),
    gridspec_kw=gridspec_kw,
    constrained_layout=False,
    squeeze=False,
    sharex=False,
)
print("done.")

_ax = ax[0, 0]

cmap = cmo.cm.ice_r

y = ds.coords["beta"].to_numpy()
x = ds.coords["tau"].to_numpy() / 360

data = ds["ln_posterior"].to_numpy()
data = normalizeTo01(data)

levels = np.linspace(0.999, 1, 101)

mappable = _ax.contourf(y, x, data, levels, cmap=cmap, extend="both")

cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
plt.colorbar(mappable=mappable, cax=cax, orientation="vertical", label="Log Posterior")#, ticks = )

if args.output != "":
    print("Outputting figure: %s" % args.output)
    fig.savefig(args.output, dpi=200)

if not args.no_display:
    print("Showin figure...")
    plt.show()

plt.close(fig)








