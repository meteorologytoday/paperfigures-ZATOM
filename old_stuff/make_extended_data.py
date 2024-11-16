from netCDF4 import Dataset
import os, re
import numpy as np

R_earth = 6.4e6  # Radius of earth in meter    
omega = 7.3e-5   # 1/s
alpha_T = 2e-1   # kg/m^3/K
alpha_S = 7e-1   # kg/m^3/PSU
g0 = 10.0        # m/s^2
rho0 = 1e3       # kg/m^3

T2b_factor = alpha_T * g0 / rho0
S2b_factor = - alpha_S * g0 / rho0

def TS2b(T, S):
    return T * T2b_factor + S * S2b_factor

def mode_integration(d, wgt, axis=-1):

    if axis == -1:
        axis = len(d.shape)-1

    return np.sum(d * wgt, axis=axis)

def yavg(d, wgt, axis=-1):

    if axis == -1:
        axis = len(d.shape)-1

    return np.average(d, weights=wgt, axis=axis)

# Variables dimension assumes to be:
# (Sample, X, Y, Z)



"""
    This function load the scan data.
    It load every xxxx.nc file in the given folder and
    connect them.

    It also load the coordinate data. It selects the first
    xxxx.nc file.
"""
def makeExtendedData(data, coor, beta = 0.9, lat_s=40.0, lat_n=68.0):
    
    Ns = data["Psib"].shape[0]

    data["bw_bnd"] = data["be"] + 2 * ( data["bw"] - data["be"] )
    psi_W = (data["Psib"][:, :-1, :] + data["Psib"][:, 1:, :]) / 2

    # Next few line is to remove integration of certain latitudes
    lat_s_ind = np.argmin(np.abs(coor["y_T"] - lat_s))
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

    # ui component
    dbdiffdt_dueto_ui = data["ui"] * (data["be"] - data["bw"]) / Lw
    data["mode1_ui_adv"] = yavg(
        mode_integration(
            dbdiffdt_dueto_ui,
            wgt = mode_1_dz_weight_T,
        ),
        wgt=coor["cos_lat"],
    )

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

    data["mode1_ZOC"] = data["mode1_ui_adv"] + data["mode1_chi_dbdz"]

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
#        print("mode1_dq = ", data["mode1_dq"])
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

# This function computes the exact terms for diagnostic relationship
def makeExtendedData2(data, coor, lat_s = 40.0, lat_n = 68.0, verbose = False, merge = False):
   
    new_data = {}
 
    Ns = data["Psib"].shape[0]
    
    # Compute bw^*
    new_data["bw_bnd"] = data["be"] + 2 * ( data["bw"] - data["be"] )
    psi_W = (data["Psib"][:, :-1, :] + data["Psib"][:, 1:, :]) / 2


    wgt_cos_lat_T = np.copy(coor["cos_lat_T"])
    wgt_cos_lat_V = np.copy(coor["cos_lat_V"])

    # Remove integration of certain latitudes
    lat_s_ind = np.argmin(np.abs(coor["y_T"] - lat_s))
    lat_n_ind = np.argmin(np.abs(coor["y_T"] - lat_n))
    verbose and print("wgt_cos_lat_T: Integrate only the range between lat [%f, %f], with idx founc [%d, %d]" % (
        lat_s, lat_n,
        lat_s_ind, lat_n_ind,
    ))
    wgt_cos_lat_T[0:lat_s_ind] = 0.0
    wgt_cos_lat_T[lat_n_ind:] = 0.0

    lat_s_ind = np.argmin(np.abs(coor["y_V"] - lat_s))
    lat_n_ind = np.argmin(np.abs(coor["y_V"] - lat_n))
    verbose and print("wgt_cos_lat_V: Integrate only the range between lat [%f, %f], with idx founc [%d, %d]" % (
        lat_s, lat_n,
        lat_s_ind, lat_n_ind,
    ))
    wgt_cos_lat_V[0:lat_s_ind] = 0.0
    wgt_cos_lat_V[lat_n_ind:] = 0.0

    _wgt_cos_lat_T   = wgt_cos_lat_T[None, :]
    _wgt_cos_lat_V   = wgt_cos_lat_V[None, :]
    _wgt_cos_lat_Vm2 = wgt_cos_lat_V[None, 1:-1]

    verbose and print( "Cleared up cos_lat weighting. Result: coor['cos_lat'] = ", coor["cos_lat"])

    # Compute geometric factor
    Lw = coor["dx_T"][0, 0]
    Le = coor["dx_T"][1, 0]
    Lambda = Lw / Le

    H = np.sum(coor["dz_T"])
    verbose and print("Total depth H = ", H)
    
    wgt_mode1_dz_T = - (2 / H) * np.sin( coor["z_T"] * np.pi / H) * coor["dz_T"]
    wgt_mode1_dz_W = - (2 / H) * np.sin( coor["z_W"] * np.pi / H) * coor["dz_W"]
   
    # Reshape weighting
    _wgt_mode1_dz_T = wgt_mode1_dz_T[None, None, :]
    _wgt_mode1_dz_W = wgt_mode1_dz_W[None, None, :]
    _wgt_mode1_dz_Wm2 = wgt_mode1_dz_W[None, None, 1:-1]


 
    # Here, `s` stands for stability, defined as s=db/dz
    
    # These are old code to derive effective stability, considering the cancellation betwen
    # MOC and ZOC advection.

    # ZOC component
    w_ZOC = data["we"] 
    w_MOC = data["ww"] + data["we"]
   

 
    dbedz_W = ( data["be"][:, :, :-1] - data["be"][:, :, 1:] ) / coor["dz_W"][1:-1][None, None, :]
    dbwdz_W = ( data["bw"][:, :, :-1] - data["bw"][:, :, 1:] ) / coor["dz_W"][1:-1][None, None, :]
    dbwdy_V = ( data["bw"][:, 1:, :] - data["bw"][:, :-1, :] ) / dy_phy_V[1:-1][None, None, :]
    
    eff_dbdz_W = Lambda * dbedz_W + dbwdz_W
    new_data["mode1_dbdiffdt_dueto_ZOC_vert_adv"] = yavg(
        mode_integration(
            w_ZOC[:, :, 1:-1] * eff_dbdz_W,
            wgt = _wgt_mode1_dz_Wm2,
        ),
        wgt = _wgt_cos_lat_T,
    )

    dbdiffdt_dueto_ui = data["ui"] * (data["be"] - data["bw"]) / Lw
    new_data["mode1_dbdiffdt_dueto_ZOC_horz_adv"] = yavg(
        mode_integration(
            dbdiffdt_dueto_ui,
            wgt = _wgt_mode1_dz_T,
        ),
        wgt = _wgt_cos_lat_T,
    )
    
    new_data["mode1_dbdiffdt_dueto_ZOC"] = new_data["mode1_dbdiffdt_dueto_ZOC_vert_adv"] + new_data["mode1_dbdiffdt_dueto_ZOC_horz_adv"]

    # MOC component
    new_data["mode1_dbdiffdt_dueto_MOC_vert_adv"] = yavg(
        mode_integration(
            w_MOC[:, :, 1:-1] * dbwdz_W,
            wgt = _wgt_mode1_dz_Wm2,
        ),
        wgt = _wgt_cos_lat_T,
    )

    new_data["mode1_dbdiffdt_dueto_MOC_horz_adv"] = yavg(
        mode_integration(
            data["vw"][:, 1:-1, :] * dbwdy_V,
            wgt = _wgt_mode1_dz_Wm2,
        ),
        wgt = _wgt_cos_lat_Vm2,
    )
    
    new_data["mode1_dbdiffdt_dueto_MOC"] = new_data["mode1_dbdiffdt_dueto_MOC_vert_adv"] + new_data["mode1_dbdiffdt_dueto_MOC_horz_adv"]
    
    new_data["mode1_dbdiffdt_dueto_MOC"] *= 1e-6 #new_data["mode1_dbdiffdt_dueto_MOC_horz_adv"]
    
    # `SS` stands for "source and sink"
    SS_T = data["X_SS_"][:, 0, :, :, :]
    SS_S = data["X_SS_"][:, 1, :, :, :]
    
    SS_b_dueto_T = SS_T * T2b_factor
    SS_b_dueto_S = SS_S * S2b_factor
    
    dbdiffdt_dueto_SS_T = SS_b_dueto_T[:, 1, :, :] - SS_b_dueto_T[:, 0, :, :]
    dbdiffdt_dueto_SS_S = SS_b_dueto_S[:, 1, :, :] - SS_b_dueto_S[:, 0, :, :]
 
    new_data["mode1_dbdiffdt_dueto_SS_T"] = yavg(
        mode_integration(
            dbdiffdt_dueto_SS_T,
            wgt = _wgt_mode1_dz_T,
        ),
        wgt = _wgt_cos_lat_T,
    )

    new_data["mode1_dbdiffdt_dueto_SS_S"] = yavg(
        mode_integration(
            dbdiffdt_dueto_SS_S,
            wgt = _wgt_mode1_dz_T,
        ),
        wgt = _wgt_cos_lat_T,
    )
   
    # Diffusion
    VDIFU_T = data["X_VDIFU_"][:, 0, :, :, :]
    VDIFU_S = data["X_VDIFU_"][:, 1, :, :, :]
    HDIFU_T = data["X_HDIFU_"][:, 0, :, :, :]
    HDIFU_S = data["X_HDIFU_"][:, 1, :, :, :]

    VDIFU_b = TS2b(VDIFU_T, VDIFU_S)
    HDIFU_b = TS2b(HDIFU_T, HDIFU_S)

    dbdiffdt_dueto_VDIFU = Lambda * VDIFU_b[:, 1, :, :] - VDIFU_b[:, 0, :, :]
    dbdiffdt_dueto_HDIFU = Lambda * HDIFU_b[:, 1, :, :] - HDIFU_b[:, 0, :, :]

    new_data["mode1_dbdiffdt_dueto_VDIFU"] = yavg(
        mode_integration(
            dbdiffdt_dueto_VDIFU,
            wgt = _wgt_mode1_dz_T,
        ),
        wgt = _wgt_cos_lat_T,
    )
 
    new_data["mode1_dbdiffdt_dueto_HDIFU"] = yavg(
        mode_integration(
            dbdiffdt_dueto_HDIFU,
            wgt = _wgt_mode1_dz_T,
        ),
        wgt = _wgt_cos_lat_T,
    )

     
    # Psi
    f_co = 2 * omega * coor["sin_lat_T"]
    
    # psi_W = (Nx, Ny, Nz+1)
    f_W = np.ones_like(psi_W) * f_co[None, :, None]  # borrow the shape from psi_W
    mode1_f = yavg(
        mode_integration(
            f_W,
            wgt = _wgt_mode1_dz_W,
        ),
        wgt = wgt_cos_lat_T,
    )
    verbose and print("mode1_f = ", mode1_f) # (Ns,)
    
    new_data["mode1_psi"] = yavg(
        mode_integration(
            f_co[None, :, None] * psi_W,
            wgt = _wgt_mode1_dz_W,
        ),
        wgt = _wgt_cos_lat_T,
    ) / mode1_f

    #db_ew = data["be"] - data["bw_bnd"]
    #db_ew = np.sum(db_ew * mode_1_dz_weight_T[None, None, :], axis=2)
    #db_ew = np.average(db_ew, weights=coor["cos_lat"], axis=1) 
    #data["mode1_db_ew"] = db_ew
    
    
    # Convective mixing

    # notice that q is separated into T and S components. So z is now the fourth axis
    
    qw_b = TS2b(data["qw"][:, 0, :, :], data["qw"][:, 1, :, :])
    qe_b = TS2b(data["qe"][:, 0, :, :], data["qe"][:, 1, :, :])

    dq_b = qe_b - qw_b

    new_data["mode1_dbdiffdt_dueto_dq_b"] = yavg(
        mode_integration(
            dq_b,
            wgt = _wgt_mode1_dz_T,
        ),
        wgt = _wgt_cos_lat_T,
    )


    if merge is True:
        data.update(new_data)
        return data
    else:
        return new_data
