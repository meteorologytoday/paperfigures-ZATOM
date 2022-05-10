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
def loadScanData(folder, loaded_varnames, load_coor=True):

    filenames = os.listdir(folder)
    filenames.sort()
 
    valid_filenames = []
 
    data_by_file = []
    coor = {} if load_coor else None

    for (j, filename) in enumerate(filenames):
        m = re.search(r'^(?P<num>[0-9]{4}).nc$', filename)
        if m is not None:
            num = int(m.group('num'))
            
            if num == 0:
                print("Skip %s" % filename)
                continue

            valid_filenames.append(filename)

            d = {}
            d["filename"]  = filename
            with Dataset("%s/%s" % (folder, filename), "r") as f:
               
                if load_coor and ("y_V" not in coor):
                    for varname in ["y_V", "z_W", "y_T", "z_T"]:
                        coor[varname] = f.variables[varname][:]

                for varname in loaded_varnames:
                    d[varname] = f.variables[varname][:]
                    #print("Shape of variable %s : " % (varname,), d[varname].shape)

            data_by_file.append(d)

    # Processing data
    concat_dataset = {}
    for varname in loaded_varnames :
        tmp_data = None
        for i, d in enumerate(data_by_file):
            if i == 0:
                tmp_data = d[varname]
            else:
                tmp_data = np.concatenate((tmp_data, d[varname]), axis=0)
            
        concat_dataset[varname] = tmp_data


    if load_coor:
        coor["dz_T"] = coor["z_W"][:-1] - coor["z_W"][1:]
        coor["dy_T"] = coor["y_V"][1:]  - coor["y_V"][:-1]
        coor["dA_T"] = coor["dy_T"][:, None] * coor["dz_T"][None, :]
         
        coor["y_V"] *= 180.0/np.pi
        coor["y_T"] *= 180.0/np.pi
        coor["z_T"] *= -1 
        coor["z_W"] *= -1 
    
    return concat_dataset, coor






def detectRanges(arr):

    rngs = []
    vals = []
    
    N = len(arr)

    new_rng = True
    beg_i = 0

    detect_val = 0.0
    for i in range(N):

        if new_rng:
            new_rng = False
            detect_val = arr[beg_i]
        else:
            if arr[i] == detect_val:  # still the detected value
                pass # move on
            else: # chop
                new_rng = True
    
        # ending is a natural chop        
        if new_rng or i == N-1:
            vals.append(detect_val)
            rngs.append(slice(beg_i, i)) # remember this point is not included but to the next branch
            beg_i = i


    return vals, rngs
