import matplotlib as mplt

#mplt.use('Agg')

import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib import rc
from pathlib import Path


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
#from pprint import pprint
import pprint
import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--data-file', required=True)
parser.add_argument('--output',    default="")
parser.add_argument('--idx', type=int, default=-1)
parser.add_argument('--display', action="store_true")
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)

coor = {}
data = {}

with Dataset(args.data_file, "r") as f:

    for varname in ["bw", "be", "Psib", "ui", "vw", "ww", "we"]:
        data[varname] = f.variables[varname][:][args.idx]


    for varname in ["y_V", "z_W", "y_T", "z_T"]:
        coor[varname] = f.variables[varname][:]

    data["vw"] = (data["vw"][1:, :] + data["vw"][:-1, :])/2.0
    data["ww"] = (data["ww"][:, 1:] + data["ww"][:, :-1])/2.0
    data["we"] = (data["we"][:, 1:] + data["we"][:, :-1])/2.0

#Ns = data["Psib"].shape[0]
#print("Number of scans: %d" % Ns)

coor["y_V"] *= 180.0/np.pi
coor["y_T"] *= 180.0/np.pi
coor["z_T"] *= -0.001
coor["z_W"] *= -0.001

ui_levels=np.linspace(-1, 1, 21) * 0.01
psi_levels=np.linspace(-50, 50, 101) * 2
w_levels=np.linspace(-100, 100, 101)

cmap1 = cm.get_cmap("prism")
cmap2 = cm.get_cmap("bwr")

clevels1 = np.linspace(0, 4, 101)
clevels2 = np.linspace(-1, 1, 11) * 0.2

yticks = np.linspace(0, 4, 5)
xticks = np.array([10,40,70,])

fig, ax = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True, sharey=False)

#fig.suptitle("[%s] s = %d" % (args.casename, s+1))

aspect_ratio = 6.4e6 * (70-10)/360 / 4e3
print("Aspect ratio: %.1e" % aspect_ratio)

ui   = data["ui"][:, :].transpose()
ww   = data["ww"][:, :].transpose() * aspect_ratio
we   = data["we"][:, :].transpose() * aspect_ratio
vw   = data["vw"][:, :].transpose()
ve   = vw * 0
bw   = (data["bw"][:, :].transpose() - 0.24) * 1e2
be   = (data["be"][:, :].transpose() - 0.24) * 1e2
Psib = data["Psib"][:, :].transpose() / 1e6

print(bw.shape)

mappable1 = ax[0].contourf(coor["y_T"], coor["z_T"], bw, clevels1, cmap=cmap1, extend="both")
mappable1 = ax[1].contourf(coor["y_T"], coor["z_T"], be, clevels1, cmap=cmap1, extend="both")

skip = slice(None, None, 2)
ax[0].quiver(coor["y_T"][skip], coor["z_T"][skip], vw[skip, skip], ww[skip, skip] * 10, width=0.005)
ax[1].quiver(coor["y_T"][skip], coor["z_T"][skip], ve[skip, skip], we[skip, skip] * 10, width=0.005)
#ax[0].contour(coor["y_T"], coor["z_T"], ui, ui_levels, colors="black")
#ax[1].contour(coor["y_T"], coor["z_T"], ui, ui_levels, colors="black")

mappable2 = ax[2].contourf(coor["y_T"], coor["z_T"], be - bw, clevels2, cmap=cmap2, extend="both")
ax[2].contour(coor["y_V"], coor["z_W"], Psib, psi_levels, colors="black")

cb1 = plt.colorbar(mappable1, ax=ax[1], orientation="vertical")
cb2 = plt.colorbar(mappable2, ax=ax[2], orientation="vertical")

cb1.set_label(r"[ $\times\,10^{-2}\,\mathrm{m}^2 / \mathrm{s}^2$ ]")
cb2.set_label(r"[ $\times\,10^{-2}\,\mathrm{m}^2 / \mathrm{s}^2$ ]")

ax[0].set_title(r"$b_w$")
ax[1].set_title(r"$b_e$")
ax[2].set_title(r"$b_e-b_w$ ; $\Psi_b$ (intvl = $2 \mathrm{Sv}$)")

ax[0].set_ylabel("Depth [km]")

for _ax in ax.flatten():
    _ax.set_xticks(xticks)
    _ax.set_xlabel("Latitude [ deg ]")
    
ax[0].set_yticks(yticks)
ax[1].set_yticks(yticks)
ax[2].set_yticks(yticks)

#ax[0].set_ylim([0, 1])
#ax[1].set_ylim([0, 1])
ax[0].set_ylim([0, coor["z_W"][-1]])
ax[1].set_ylim([0, coor["z_W"][-1]])
ax[2].set_ylim([0, coor["z_W"][-1]])


ax[0].invert_yaxis()
ax[1].invert_yaxis()
ax[2].invert_yaxis()


# clear unwanted Y axis labels
for _ax in ax[1:2].flatten():
    _ax.set_yticklabels(["" for _ in yticks])

if args.output != "":
    print("Outputting figure: %s" % args.output)
    fig.savefig(args.output, dpi=200)

if args.display == True:
    plt.show()

plt.close(fig)








