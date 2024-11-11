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

    if len(wgt.shape) != len(d.shape):
        _padded_idx = [  slice(None) if i == axis else None for i in range(len(d.shape)) ]
        _padded_wgt = wgt[*_padded_idx]
    else:
        _padded_wgt = wgt

    return np.sum(d * _padded_wgt, axis=axis)

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
# This function computes the exact terms for diagnostic relationship
def makeExtendedData(data, coor, lat_s = 40.0, lat_n = 68.0, verbose = False, merge = False):
   
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


    wgt_cos_lat_Vm2 = wgt_cos_lat_V[1:-1]
    #_wgt_cos_lat_T   = wgt_cos_lat_T[None, :]
    #_wgt_cos_lat_V   = wgt_cos_lat_V[None, :]
    #_wgt_cos_lat_Vm2 = wgt_cos_lat_V[None, 1:-1]

    verbose and print( "Cleared up cos_lat weighting. Result: wgt_cos_lat_V = ", wgt_cos_lat_V)

    # Compute geometric factor
    Lw = coor["dx_T"][0, 0]
    Le = coor["dx_T"][1, 0]
    Lambda = Lw / Le

    H = np.sum(coor["dz_T"])
    verbose and print("Total depth H = ", H)
    
    wgt_mode1_dz_T = - (2 / H) * np.sin( coor["z_T"] * np.pi / H) * coor["dz_T"]
    wgt_mode1_dz_W = - (2 / H) * np.sin( coor["z_W"] * np.pi / H) * coor["dz_W"]
    wgt_mode1_dz_Wm2 = wgt_mode1_dz_W[1:-1]
   
    # Reshape weighting
    #_wgt_mode1_dz_T = wgt_mode1_dz_T[None, None, :]
    #_wgt_mode1_dz_W = wgt_mode1_dz_W[None, None, :]
    #_wgt_mode1_dz_Wm2 = wgt_mode1_dz_W[None, None, 1:-1]


 
    # Here, `s` stands for stability, defined as s=db/dz
    
    # These are old code to derive effective stability, considering the cancellation betwen
    # MOC and ZOC advection.

    # ZOC component
    #w_ZOC = - data["we"] 
    #w_MOC = data["ww"] + data["we"]
    #dy_phy_V = R_earth * coor["dy_V"]

   
    #dTedz_W = ( Te[:, :, :-1] - Te[:, :,  1:] ) / coor["dz_W"][1:-1][None, None, :]
    #dTwdz_W = ( Tw[:, :, :-1] - Tw[:, :,  1:] ) / coor["dz_W"][1:-1][None, None, :]
    #dTwdy_V = ( Tw[:, 1:,  :] - Tw[:, :-1, :] ) / dy_phy_V[1:-1][None, :, None]
 
    #dSedz_W = ( Se[:, :, :-1] - Se[:, :,  1:] ) / coor["dz_W"][1:-1][None, None, :]
    #dSwdz_W = ( Sw[:, :, :-1] - Sw[:, :,  1:] ) / coor["dz_W"][1:-1][None, None, :]
    #dSwdy_V = ( Sw[:, 1:,  :] - Sw[:, :-1, :] ) / dy_phy_V[1:-1][None, :, None]
 
    #dbedz_W = ( data["be"][:, :, :-1] - data["be"][:, :, 1:] ) / coor["dz_W"][1:-1][None, None, :]
    #dbwdz_W = ( data["bw"][:, :, :-1] - data["bw"][:, :, 1:] ) / coor["dz_W"][1:-1][None, None, :]
    #dbwdy_V = ( data["bw"][:, 1:, :] - data["bw"][:, :-1, :] ) / dy_phy_V[1:-1][None, :, None]
    #eff_dbdz_W = Lambda * dbedz_W + dbwdz_W
    
    #eff_dTdz_W = Lambda * dTedz_W + dTwdz_W
    #eff_dSdz_W = Lambda * dSedz_W + dSwdz_W
    
    def getAdv(Xw, Xe, varname="X"):

        result = dict()

        dXedz_W = ( Xe[:, :, :-1] - Xe[:, :,  1:] ) / coor["dz_W"][1:-1][None, None, :]
        dXwdz_W = ( Xw[:, :, :-1] - Xw[:, :,  1:] ) / coor["dz_W"][1:-1][None, None, :]
        dXwdy_V = ( Xw[:, 1:,  :] - Xw[:, :-1, :] ) / dy_phy_V[1:-1][None, :, None]
        eff_dXdz_W = Lambda * dXedz_W + dXwdz_W
        
        result["d%sdiffdt_dueto_ZOC_vert_adv" % (varname,)] = yavg(
            mode_integration(
                w_ZOC[:, :, 1:-1] * eff_dXdz_W,
                wgt = wgt_mode1_dz_Wm2,
            ),
            wgt = wgt_cos_lat_T,
        )

        ui_dXdx = data["ui"] * (Xe - Xw) / Lw
        result["d%sdiffdt_dueto_ZOC_horz_adv" % (varname,)] = yavg(
            mode_integration(
                ui_dXdx,
                wgt = wgt_mode1_dz_T,
            ),
            wgt = wgt_cos_lat_T,
        )
        
        result["d%sdiffdt_dueto_ZOC" % (varname,)] = (
              result["d%sdiffdt_dueto_ZOC_vert_adv" % (varname,)] 
            + result["d%sdiffdt_dueto_ZOC_horz_adv" % (varname,)]
        )
 
        # MOC component
        result["d%sdiffdt_dueto_MOC_vert_adv" % (varname,)] = yavg(
            mode_integration(
                w_MOC[:, :, 1:-1] * dXwdz_W,
                wgt = wgt_mode1_dz_Wm2,
            ),
            wgt = wgt_cos_lat_T,
        )
        
        result["d%sdiffdt_dueto_MOC_horz_adv" % (varname,)] = yavg(
            mode_integration(
                data["vw"][:, 1:-1, :] * dXwdy_V,
                wgt = wgt_mode1_dz_T,
            ),
            wgt = wgt_cos_lat_Vm2,
        )
        
        result["d%sdiffdt_dueto_MOC" % (varname,)] = (
              result["d%sdiffdt_dueto_MOC_vert_adv" % (varname,)] 
            + result["d%sdiffdt_dueto_MOC_horz_adv" % (varname,)]
        )
        
        return result

    #Tw = data["X_"][:, 0, 0, :, :]
    #Te = data["X_"][:, 0, 1, :, :]
    #Sw = data["X_"][:, 1, 0, :, :]
    #Se = data["X_"][:, 1, 1, :, :]

    #advTEMP = getAdv(Xw=Tw, Xe=Te, varname="T")
    #advSALT = getAdv(Xw=Sw, Xe=Se, varname="S")
    #advBUOY = getAdv(Xw=data["bw"], Xe=data["be"], varname="b")
    
    new_data["mode1_dbdiffdt_dueto_ZOC_vert_adv"] = yavg(
        mode_integration(
            w_ZOC[:, :, 1:-1] * eff_dbdz_W,
            wgt = wgt_mode1_dz_Wm2,
        ),
        wgt = wgt_cos_lat_T,
    )

    dbdiffdt_dueto_ui = data["ui"] * (data["be"] - data["bw"]) / Lw
    new_data["mode1_dbdiffdt_dueto_ZOC_horz_adv"] = yavg(
        mode_integration(
            dbdiffdt_dueto_ui,
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_T,
    )
 
    new_data["mode1_dbdiffdt_dueto_ZOC"] = new_data["mode1_dbdiffdt_dueto_ZOC_vert_adv"] + new_data["mode1_dbdiffdt_dueto_ZOC_horz_adv"]


    new_data["mode1_dbdiffdt_dueto_ZOC_vert_adv"] = yavg(
        mode_integration(
            w_ZOC[:, :, 1:-1] * eff_dbdz_W,
            wgt = wgt_mode1_dz_Wm2,
        ),
        wgt = wgt_cos_lat_T,
    )

    dbdiffdt_dueto_ui = data["ui"] * (data["be"] - data["bw"]) / Lw
    new_data["mode1_dbdiffdt_dueto_ZOC_horz_adv"] = yavg(
        mode_integration(
            dbdiffdt_dueto_ui,
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_T,
    )
    
    new_data["mode1_dbdiffdt_dueto_ZOC"] = new_data["mode1_dbdiffdt_dueto_ZOC_vert_adv"] + new_data["mode1_dbdiffdt_dueto_ZOC_horz_adv"]

    # MOC component
    new_data["mode1_dbdiffdt_dueto_MOC_vert_adv"] = yavg(
        mode_integration(
            w_MOC[:, :, 1:-1] * dbwdz_W,
            wgt = wgt_mode1_dz_Wm2,
        ),
        wgt = wgt_cos_lat_T,
    )
    
    new_data["mode1_dbdiffdt_dueto_MOC_horz_adv"] = yavg(
        mode_integration(
            data["vw"][:, 1:-1, :] * dbwdy_V,
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_Vm2,
    )
    
    new_data["mode1_dbdiffdt_dueto_MOC"] = new_data["mode1_dbdiffdt_dueto_MOC_vert_adv"] + new_data["mode1_dbdiffdt_dueto_MOC_horz_adv"]
    """

    p = r"^d([a-zA-Z_0-9]+)(diffdt_dueto_[a-zA-Z_0-9]+)$"
    for d in [ advTEMP, advBUOY, advSALT  ]:
        
        tmp = dict()
        for k, v in d.items():

            match = re.match(p, k)

            if match:
                # Access the captured groups
                varname = match.group(1)
                suffix = match.group(2)
                print("Match: (%s, %s)" % (varname, suffix,))
            else:
                print("No match found for: ", k)

            new_varname = "mode1_db%s" % suffix
            
            if varname == "T":
                v = TS2b(v, 0.0)
                new_varname = "%s_T" % (new_varname,)
            elif varname == "S":
                v = TS2b(0.0, v)
                new_varname = "%s_S" % (new_varname,)

            tmp[new_varname] = v

        new_data.update(tmp)
    

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
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_T,
    )

    new_data["mode1_dbdiffdt_dueto_SS_S"] = yavg(
        mode_integration(
            dbdiffdt_dueto_SS_S,
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_T,
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
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_T,
    )
 
    new_data["mode1_dbdiffdt_dueto_HDIFU"] = yavg(
        mode_integration(
            dbdiffdt_dueto_HDIFU,
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_T,
    )

     
    # Psi
    f_co = 2 * omega * coor["sin_lat_T"]
    
    # psi_W = (Nx, Ny, Nz+1)
    f_W = np.ones_like(psi_W) * f_co[None, :, None]  # borrow the shape from psi_W
    mode1_f = yavg(
        mode_integration(
            f_W,
            wgt = wgt_mode1_dz_W,
        ),
        wgt = wgt_cos_lat_T,
    )

    if np.any( (mode1_f[1:] - mode1_f[:-1] ) != 0):
        raise Exception("mode1_f is not a constant in differencne Ns. Weird.")

    verbose and print("mode1_f = ", mode1_f[0]) # (Ns,)
    
    new_data["mode1_psi"] = yavg(
        mode_integration(
            f_co[None, :, None] * psi_W,
            wgt = wgt_mode1_dz_W,
        ),
        wgt = wgt_cos_lat_T,
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
            wgt = wgt_mode1_dz_T,
        ),
        wgt = wgt_cos_lat_T,
    )

    new_data["mode1_dbdiffdt_dueto_ZOCMOC"] = (
          new_data["mode1_dbdiffdt_dueto_ZOC"]
        + new_data["mode1_dbdiffdt_dueto_MOC"]
    )
    
    new_data["mode1_residue"] = (
          new_data["mode1_dbdiffdt_dueto_ZOC"]
        + new_data["mode1_dbdiffdt_dueto_MOC"]
        + new_data["mode1_dbdiffdt_dueto_dq_b"]
        + new_data["mode1_dbdiffdt_dueto_SS_T"]
        + new_data["mode1_dbdiffdt_dueto_SS_S"]
        + new_data["mode1_dbdiffdt_dueto_VDIFU"]
        + new_data["mode1_dbdiffdt_dueto_HDIFU"]
    )

    new_data["mode1_chi"] = yavg(
        mode_integration(
            data["chi"],
            wgt = wgt_mode1_dz_W,
        ),
        wgt = wgt_cos_lat_T,
    )




    if merge is True:
        data.update(new_data)
        return data
    else:
        return new_data
