from netCDF4 import Dataset
import os, re
import numpy as np
    
omega = 7.3e-5   # 1/s
alpha_T = 2e-1   # kg/m^3/K
alpha_S = 7e-1   # kg/m^3/PSU
g0 = 10.0        # m/s^2
rho0 = 1e3       # kg/m^3
beta = 0.9

"""
    This function load the scan data.
    It load every xxxx.nc file in the given folder and
    connect them.

    It also load the coordinate data. It selects the first
    xxxx.nc file.
"""
def makeExtendedData(data, coor):
    
    Ns = data["Psib"].shape[0]

    data["bw_bnd"] = data["be"] + 2 * ( data["bw"] - data["be"] )
    psi_W = (data["Psib"][:, :-1, :] + data["Psib"][:, 1:, :]) / 2

    # Next few line is to remove integration of certain latitudes
    lat_s = 40.0
    lat_s_ind = np.argmin(np.abs(coor["y_T"] - lat_s))

    lat_n = 68.0
    lat_n_ind = np.argmin(np.abs(coor["y_T"] - lat_n))

    print("Integrate only the range between lat [%f, %f], with idx founc [%d, %d]" % (
        lat_s, lat_n,
        lat_s_ind, lat_n_ind,
    ))

    coor["cos_lat"][0:lat_s_ind] = 0.0
    coor["cos_lat"][lat_n_ind:] = 0.0
    print( coor["cos_lat"])

    Lw = coor["dx_T"][0, 0]
    Le = coor["dx_T"][1, 0]
    Lambda = Lw / Le

    H = np.sum(coor["dz_T"])

    mode_1_dz_weight_T = - (2 / H) * np.sin( coor["z_T"] * np.pi / H) * coor["dz_T"]
    mode_1_dz_weight_W = - (2 / H) * np.sin( coor["z_W"] * np.pi / H) * coor["dz_W"]
    
    print("Total depth H = ", H)
    print("Beta = ", beta)
    #b_mean = ( data["bw"] * Lw + data["be"] * Le ) / ( Lw + Le )
    b_eff  = (1 - beta) * data["bw"] + data["be"] * Lambda
    s_eff_W = ( b_eff[:, :, :-1] - b_eff[:, :, 1:] ) / ( coor["dz_W"][1:-1][None, None, :] )
    s_eff_vint = np.sum(s_eff_W * mode_1_dz_weight_W[None, None, 1:-1], axis=2)
    data["mode1_s_eff"] = np.average( s_eff_vint, weights=coor["cos_lat"], axis=1)
    

    chi_vint = np.sum(data["chi"] * mode_1_dz_weight_W[None, None, :], axis=2)
    data["mode1_chi"] = np.average( chi_vint, weights=coor["cos_lat"], axis=1)
    
    data["mode1_chi_dbdz_product"] = data["mode1_chi"] * data["mode1_s_eff"] / Lw 

    chi_dbdz = s_eff_W * data["chi"][:, :, 1:-1] / Lw
    chi_dbdz_vint = np.sum( chi_dbdz * mode_1_dz_weight_W[None, None, 1:-1], axis=2)
    data["mode1_chi_dbdz"] = np.average( chi_dbdz_vint[:, :], weights=coor["cos_lat"], axis=1)


    f_co = 2 * omega * coor["sin_lat"]
    
    f_W = (psi_W[0, :, :] * 0 + 1.0) * f_co[:, None]  # borrow the shape from psi_W
    mode1_f = np.average( np.sum( f_W * mode_1_dz_weight_W[None, :], axis=1), weights=coor["cos_lat"], axis=0 )
    print("mode1_f = ", mode1_f)

    fpsi_W = f_co[None, :, None] * psi_W
    fpsi_vint = np.sum( fpsi_W * mode_1_dz_weight_W[None, None, :], axis=2)
    data["mode1_psi"] = np.average( fpsi_vint, weights=coor["cos_lat"], axis=1) / mode1_f
    
    db_ew = data["be"] - data["bw_bnd"]
    db_ew = np.sum(db_ew * mode_1_dz_weight_T[None, None, :], axis=2)
    db_ew = np.average(db_ew, weights=coor["cos_lat"], axis=1) 
    data["mode1_db_ew"] = db_ew

    # Convective mixing and convectivity
  
    be = data["be"]
    bw = data["bw"] 
    dwf_east = np.zeros(be.shape)
    dwf_west = np.zeros(be.shape)

    cvt_e = np.zeros( (Ns,) )
    cvt_w = np.zeros( (Ns,) )

    for j in range(Ns):

        flags_e = computeDeepWaterFormationGrids(be[j, :, :])
        flags_w = computeDeepWaterFormationGrids(bw[j, :, :])

        dwf_east[j, :, :] = flags_e
        dwf_west[j, :, :] = flags_w

        cvt_e[j] = computeConvectivity(flags_e) 
        cvt_w[j] = computeConvectivity(flags_w) 

    data["dwf_east"] = dwf_east
    data["dwf_west"] = dwf_west
    data["cvt_e"] = cvt_e
    data["cvt_w"] = cvt_w
    data["d_cvt"] = cvt_w - cvt_e


    if ("qe" in data) and ("qw" in data):
        print("Compute dq")
        
        # notice that q is separated into T and S components. So z is now the fourth axis
        qw_vint = np.sum( data["qw"][:, :, :, :] * mode_1_dz_weight_T[None, None , None, :], axis=3)
        qw = np.average(qw_vint, weights=coor["cos_lat"], axis=2) 
 
        qe_vint = np.sum( data["qe"][:, :, :, :] * mode_1_dz_weight_T[None, None , None, :], axis=3)
        qe = np.average(qe_vint, weights=coor["cos_lat"], axis=2) 
        
        # Convert T and S into buoyancy
        qw_b = g0 / rho0 * (qw[:, 0] * alpha_T - qw[:, 1] * alpha_S) 
        qe_b = g0 / rho0 * (qe[:, 0] * alpha_T - qe[:, 1] * alpha_S) 

        data["mode1_dq"] = qe_b - qw_b

# Compute the deep water formation grids
def computeDeepWaterFormationGrids(b):

    flags = np.zeros(b.shape)
    b_lower = np.roll(b, -1, axis=1)
    
    flags[b_lower > b] = 1.0
    flags[:, -1] = flags[:, -2]
    
    return flags

# Compute the fraction of convective mixing. 
# It has to extend vertically from the surface downward.
def computeConvectivity(flags):
    
    cnt = 0

    # Search from surface
    for j in range(flags.shape[0]):
        for k in range(flags.shape[1]):
            if flags[j, k]:
                cnt += 1
            else:
                pass    # count all flags
                #break  # only count if it is unstable all the way to the top
    
    return cnt / flags.size


