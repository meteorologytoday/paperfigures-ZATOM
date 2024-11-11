import pprint 
import xarray as xr
import load_scan_data as lsd
import make_extended_data as med

import matplotlib.pyplot as plt

import os,sys, argparse
from netCDF4 import Dataset
import numpy as np

def computeError(data, tau, A, C):
    err = data["mode1_psi"] - (tau / A) * ( data["mode1_chi_dbdz"] + data["mode1_dq"] + C * data["xi"] * data["Q"] )
    return err



parser = argparse.ArgumentParser()
parser.add_argument('--folder', nargs='+')
parser.add_argument('--beta-rng', nargs=2, type=float, default=[0.0, 1.0])
parser.add_argument('--gamma-rng', nargs=2, type=float, required=True)
parser.add_argument('--gamma-scale', type=float, help="Scale of gamma in Sv", required=True)
parser.add_argument('--psi-scale', type=float, help="Scale of psi in Sv", required=True)
parser.add_argument('--sampling-spacing', type=float, help="Sampling spacing after scaled", required=True)
parser.add_argument('--tau-rng', nargs=2, type=float, help="Tau range in days", default=[0, 20])
parser.add_argument('--N-tau', type=int, default=4)
parser.add_argument('--N-beta', type=int, default=4)
parser.add_argument('--const-A', type=float, help="Const A", required=True)
parser.add_argument('--const-C', type=float, help="Const C", required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--no-display', action="store_true")
parser.add_argument('--residue-threshold', type=float, default=1e-12)
args = parser.parse_args()

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(args)


data = []
coor = {}

folders = args.folder
N_data = len(folders)

data = []
coor = None
loaded_varnames = ["Psib", "chi", "be", "bw", "qw", "qe", "res", "stable", "Q", ]
for i, folder in enumerate(folders):

    print("Loading the folder: %s" % (folder,))
    _data, _coor = lsd.loadScanData(
        folder,
        loaded_varnames,
        load_coor=(coor is None),
        residue_threshold = args.residue_threshold,
    )

    print("Number of records: %d" % (len(_data["Q"]))) 
    if coor is None:
        coor = _coor

    med.makeExtendedData(_data, coor)

    data.append(_data)

print("Finding the indexes to to subsampling")
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
    
    
for d in data:
    
    Q_tmp = d["Q"] / (1e6*args.gamma_scale)
    psi_tmp = d["mode1_psi"] / (1e6*args.psi_scale)
  
    ds = np.zeros_like(Q_tmp)
    ds[1:] = ((Q_tmp[1:] - Q_tmp[:-1])**2 + (psi_tmp[1:] - psi_tmp[:-1])**2)**0.5
    s = np.cumsum(ds)
    
    idx = findSamplingIdxBySpacing(s, args.sampling_spacing) 

    Q_sub = d["Q"][idx]/1e6
    idx = idx[ (Q_sub >= args.gamma_rng[0]) & (Q_sub <= args.gamma_rng[1]) ]

    d["subsampling_idx"] = idx


if not args.no_display:

    print("`--no-display` is not set. Will show figures...")

    import tool_fig_config

    print("Loading matplotlib...")
    import matplotlib as mplt

    mplt.use("TkAgg")

    import matplotlib.pyplot as plt

    import matplotlib.gridspec as gridspec
    from matplotlib import cm
    from matplotlib import rc

    import colorblind
    print("Done")

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

    fig, ax = plt.subplots(
        nrow, ncol,
        figsize=figsize,
        subplot_kw=dict(aspect="auto"),
        gridspec_kw=gridspec_kw,
        constrained_layout=False,
        squeeze=False,
        sharex=False,
    )

    _ax = ax[0, 0]
    param = "Q"
    var = "mode1_psi"

    for i in range(len(args.folder)):
        
        print("Plotting %s" % (folders[i]))

        color = colorblind.BW8color[list(colorblind.BW8color.keys())[i]]
        d = data[i]
        vals, rngs = lsd.detectRanges(d["stable"])

        for j, val in enumerate(vals):
           
 
            rng = rngs[j]

            x   = d[param][rng]
            res = d["res"][rng]

            y   = d[var][rng]

            linestyle = "solid" if val == 1.0 else "dashed"
 
            kw_label = {
                "color" : color,
            }


            _ax.plot(x, y, zorder=1, linestyle=linestyle, **kw_label)

        subsampling_idx = d["subsampling_idx"]

        x = d[param][subsampling_idx]
        y = d[var][subsampling_idx]
        _ax.scatter(x, y, s=10, marker="o", c=color, edgecolor=color, zorder=50)

    print("Showing figures...")
    plt.show()



print("Fitting data now...")



# Step 1: Construct parameter space
N_beta = args.N_beta
N_tau = args.N_tau

param_tau  = np.linspace(args.tau_rng[0]*86400, args.tau_rng[1]*86400, N_tau)
param_beta = np.linspace(args.beta_rng[0], args.beta_rng[1], N_beta)
sigma_psi = 0.5e6
ln_post = np.zeros( (N_tau, N_beta) )


# Step 2:
# scan through possible beta and tau
for j, tau in enumerate(param_tau):
    for i, beta in enumerate(param_beta):
        for k, d in enumerate(data):
            
            med.makeExtendedData(d, coor, beta=beta)
            err = computeError(d, tau, args.const_A, args.const_C)
            ln_post[j, i] += - np.sum(err**2 / (2 * sigma_psi**2))

ln_post -= np.amax(ln_post)

post = np.exp(ln_post)

# Step 3: output posterior
print("Output posterior file: ", args.output)

ds = xr.Dataset(
    data_vars = dict(
        ln_posterior = (["tau", "beta"], ln_post),
        posterior = (["tau", "beta"], post),
    ),
    coords = dict(
        tau  = (["tau", ], param_tau / 86400.0, dict(units="days")),
        beta = (["beta", ], param_beta),
    ),
)


ds.to_netcdf(args.output)



