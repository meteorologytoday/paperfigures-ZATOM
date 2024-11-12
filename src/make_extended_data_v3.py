from netCDF4 import Dataset
import os, re
import numpy as np

necessary_variables = [
    "X_ADV_MOC_",
    "X_ADV_ZOC_",
    "X_VDIFU_",
    "X_HDIFU_",
    "X_ZNLHDIFU_",
    "X_CVA_",
    "X_FRC_",
    "X_SS_",
] 
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
def makeExtendedData(data, coor, lat_s = 40.0, lat_n = 68.0, verbose = False, merge = False, mode=1):
  
    print("Computing extended data. Mode = %d" % (mode,)) 
    new_data = {}

    # Compute bw^*
    new_data["bw_bnd"] = data["be"] + 2 * ( data["bw"] - data["be"] )
    psi_W = (data["Psib"][:, :-1, :] + data["Psib"][:, 1:, :]) / 2



    # psi_W = (Nx, Ny, Nz+1)
    #f_W = np.ones_like(psi_W) * f_T[None, :, None]  # borrow the shape from psi_W
    
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
    
    wgt_mode_dz_T = - (2 / H) * np.sin( mode * coor["z_T"] * np.pi / H) * coor["dz_T"]
    wgt_mode_dz_W = - (2 / H) * np.sin( mode * coor["z_W"] * np.pi / H) * coor["dz_W"]
    #wgt_mode_dz_T = coor["dz_T"] / H
    #wgt_mode_dz_W = coor["dz_W"] / H
 

    wgt_mode_dz_Wm2 = wgt_mode_dz_W[1:-1]


    f_T = 2 * omega * coor["sin_lat_T"]
    f_T = f_T[None, :] # (s, y)
    invf_T = f_T**(-1)
    _const = 2 * (H / (mode * np.pi))**2.0
    def X2bdiff(varname):
        var_data = data[varname]
        var_TEMP = var_data[:, 0, :, :, :]
        var_SALT = var_data[:, 1, :, :, :]
 
        bdiff_dueto_TEMP = TS2b( var_TEMP[:, 1, :, :] - var_TEMP[:, 0, :, :], 0.0 )
        bdiff_dueto_SALT = TS2b( 0.0, var_SALT[:, 1, :, :] - var_SALT[:, 0, :, :] )
 
        bdiff_dueto_TEMP = _const * yavg(
            invf_T * mode_integration(
                bdiff_dueto_TEMP,
                wgt = wgt_mode_dz_T,
            ),
            wgt = wgt_cos_lat_T,
        )
 
        bdiff_dueto_SALT = _const * yavg(
            invf_T * mode_integration(
                bdiff_dueto_SALT,
                wgt = wgt_mode_dz_T,
            ),
            wgt = wgt_cos_lat_T,
        )
        
        return dict(
            TEMP = bdiff_dueto_TEMP,
            SALT = bdiff_dueto_SALT,
        ) 

    for component in [
        "ADV_MOC",
        "ADV_ZOC",
        "SS",
        "FRC",
        "VDIFU",
        "HDIFU",
        "ZNLHDIFU",
        "CVA",
    ]:
        varname = "X_%s_" % (component,)
        
        new_varname_ttl = "mode_dbdiffdt_dueto_%s" % (component,)
        new_varname_TEMP = "mode_dbdiffdt_dueto_%s_TEMP" % (component,)
        new_varname_SALT = "mode_dbdiffdt_dueto_%s_SALT" % (component,)

        tmp = X2bdiff(varname)
        new_data[new_varname_TEMP] = tmp["TEMP"]
        new_data[new_varname_SALT] = tmp["SALT"]
        new_data[new_varname_ttl]  = tmp["TEMP"] + tmp["SALT"]


    new_data["mode_psi"] = yavg(
        mode_integration(
            psi_W,
            wgt = wgt_mode_dz_W,
        ),
        wgt = wgt_cos_lat_T,
    )


    new_data["mode_dbdiffdt_dueto_ADV_ZOCMOC"] = (
          new_data["mode_dbdiffdt_dueto_ADV_ZOC"]
        + new_data["mode_dbdiffdt_dueto_ADV_MOC"]
    )
 
    new_data["mode_dbdiffdt_dueto_BGDIFU"] = (
          new_data["mode_dbdiffdt_dueto_VDIFU"]
        + new_data["mode_dbdiffdt_dueto_HDIFU"]
    )
    
    new_data["mode_dbdiffdt_sum"] = (
          new_data["mode_dbdiffdt_dueto_ADV_ZOC"]
        + new_data["mode_dbdiffdt_dueto_ADV_MOC"]
        + new_data["mode_dbdiffdt_dueto_CVA"]
        + new_data["mode_dbdiffdt_dueto_SS"]
        + new_data["mode_dbdiffdt_dueto_FRC"]
        + new_data["mode_dbdiffdt_dueto_VDIFU"]
        + new_data["mode_dbdiffdt_dueto_HDIFU"]
        + new_data["mode_dbdiffdt_dueto_ZNLHDIFU"]
    )



    # Diagnose chi tendency
    f_Tm2 = 2 * omega * coor["sin_lat_T"][1:-1]
    f_Tm2 = f_Tm2[None, :] # (s, y)
    invf_Tm2 = f_Tm2**(-1)
    _const = (H / (mode * np.pi))**2.0
    dy_Vm2 = coor["dy_V"][None, 1:-1, None] * R_earth # Transform into meters
    wgt_cos_lat_Tm2 = wgt_cos_lat_T[1:-1]
    def X2dbdy(varname, x_idx):
        var_data = data[varname]
        var_TEMP = var_data[:, 0, x_idx, :, :]
        var_SALT = var_data[:, 1, x_idx, :, :]
 
        b_dueto_TEMP = TS2b( var_TEMP, 0.0 )
        b_dueto_SALT = TS2b( 0.0, var_SALT )


        # Take y-derivative
        dbdy_dueto_TEMP = ( b_dueto_TEMP[:, 1:, :] - b_dueto_TEMP[:, :-1, :] ) / dy_Vm2
        dbdy_dueto_SALT = ( b_dueto_SALT[:, 1:, :] - b_dueto_SALT[:, :-1, :] ) / dy_Vm2
        
        # On Tm2 grid
        dbdy_dueto_TEMP = (dbdy_dueto_TEMP[:, 1:, :] + dbdy_dueto_TEMP[:, :-1, :] ) / 2.0
        dbdy_dueto_SALT = (dbdy_dueto_SALT[:, 1:, :] + dbdy_dueto_SALT[:, :-1, :] ) / 2.0
        

        dbdy_dueto_TEMP = _const * yavg(
            invf_Tm2 * mode_integration(
                dbdy_dueto_TEMP,
                wgt = wgt_mode_dz_T,
            ),
            wgt = wgt_cos_lat_Tm2,
        )
 
        dbdy_dueto_SALT = _const * yavg(
            invf_Tm2 * mode_integration(
                dbdy_dueto_SALT,
                wgt = wgt_mode_dz_T,
            ),
            wgt = wgt_cos_lat_Tm2,
        )
        
        return dict(
            TEMP = dbdy_dueto_TEMP,
            SALT = dbdy_dueto_SALT,
        ) 


    for component in [
        "ADV_MOC",
        "ADV_ZOC",
        "SS",
        "FRC",
        "VDIFU",
        "HDIFU",
        "ZNLHDIFU",
        "CVA",
    ]:
        varname = "X_%s_" % (component,)
        
        new_varname_ttl = "mode_dchidt_dueto_%s" % (component,)
        new_varname_TEMP = "mode_dchidt_dueto_%s_TEMP" % (component,)
        new_varname_SALT = "mode_dchidt_dueto_%s_SALT" % (component,)

        tmp = X2dbdy(varname, x_idx=1) # Eastern part
        new_data[new_varname_TEMP] = tmp["TEMP"]
        new_data[new_varname_SALT] = tmp["SALT"]
        new_data[new_varname_ttl]  = tmp["TEMP"] + tmp["SALT"]


    new_data["mode_dchidt_dueto_ADV_ZOCMOC"] = (
          new_data["mode_dchidt_dueto_ADV_ZOC"]
        + new_data["mode_dchidt_dueto_ADV_MOC"]
    )
 
    new_data["mode_dchidt_dueto_BGDIFU"] = (
          new_data["mode_dchidt_dueto_VDIFU"]
        + new_data["mode_dchidt_dueto_HDIFU"]
    )
    
    new_data["mode_dchidt_sum"] = (
          new_data["mode_dchidt_dueto_ADV_ZOC"]
        + new_data["mode_dchidt_dueto_ADV_MOC"]
        + new_data["mode_dchidt_dueto_CVA"]
        + new_data["mode_dchidt_dueto_SS"]
        + new_data["mode_dchidt_dueto_FRC"]
        + new_data["mode_dchidt_dueto_VDIFU"]
        + new_data["mode_dchidt_dueto_HDIFU"]
        + new_data["mode_dchidt_dueto_ZNLHDIFU"]
    )



    new_data["mode_chi"] = yavg(
        mode_integration(
            data["chi"],
            wgt = wgt_mode_dz_W,
        ),
        wgt = wgt_cos_lat_T,
    )




    if merge is True:
        data.update(new_data)
        return data
    else:
        return new_data
