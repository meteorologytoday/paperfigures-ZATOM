from netCDF4 import Dataset
import os, re
import numpy as np

"""
    This function load the scan data.
    It load every xxxx.nc file in the given folder and
    connect them.

    It also load the coordinate data. It selects the first
    xxxx.nc file.
"""
def makeExtendedData(data, coor):

    data["bw_bnd"] = data["be"] + 2 * ( data["bw"] - data["be"] )

    z1000_ind = np.argmin(np.abs(coor["z_T"] - 1000.0))
    y40_ind = np.argmin(np.abs(coor["y_T"] - 40.0))
    Lw = coor["dx_T"][0, 0]
    Le = coor["dx_T"][1, 0]

    b_mean = ( data["bw"] * Lw + data["be"] * Le ) / ( Lw + Le )
    data["s1000"] = ( b_mean[:, :, 0] - b_mean[:, :, z1000_ind] ) / ( coor["z_T"][z1000_ind] - coor["z_T"][0] )
    data["s1000"] = np.average( data["s1000"][:, :], weights=coor["cos_lat"], axis=1)
 
    data["s1000_hlat"] = ( b_mean[:, :, 0] - b_mean[:, :, z1000_ind] ) / ( coor["z_T"][z1000_ind] - coor["z_T"][0] )
    data["s1000_hlat"] = np.average( data["s1000_hlat"][:, y40_ind:], weights=coor["cos_lat"][y40_ind:], axis=1)
   
    b_mean_1000 = np.average( b_mean[:, :, 0:z1000_ind], weights=coor["dz_T"][0:z1000_ind], axis=2)
    data["db_ns"] = ( 
        np.average(b_mean_1000[:, :y40_ind], weights=coor["cos_lat"][:y40_ind], axis=1) - 
        np.average(b_mean_1000[:, y40_ind:], weights=coor["cos_lat"][y40_ind:], axis=1)
    )
    
    db_ew = data["be"] - data["bw_bnd"]
    db_ew = np.average(db_ew[:, :, 0:z1000_ind], weights=coor["dz_T"][0:z1000_ind], axis=2)
    db_ew = np.average(db_ew, weights=coor["cos_lat"], axis=1) 
    data["db_ew"] = db_ew

    # Convective mixing and convectivity
  
    be = data["be"]
    bw = data["bw"] 
    dwf_east = np.zeros(be.shape)
    dwf_west = np.zeros(be.shape)

    cvt_e = np.zeros( (len(data["Q"]),) )
    cvt_w = np.zeros( (len(data["Q"]),) )

    for j in range(len(data["Q"])):

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
                break

    return cnt / flags.size


