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
def loadScanData(folder, loaded_varnames, load_coor=True, residue_threshold=np.nan):

    filenames = os.listdir(folder)
    filenames.sort()
 
    valid_filenames = []
 
    data_by_file = []
    coor = {} if load_coor else None

    load_anything = False

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
                    for varname in ["y_V", "z_W", "y_T", "z_T", "dx_T"]:
                        coor[varname] = f.variables[varname][:]
                        

                for varname in loaded_varnames:
                    d[varname] = f.variables[varname][:]
                    #print("Shape of variable %s : " % (varname,), d[varname].shape)
               
                attrs = f.ncattrs()
                if "xi" in attrs:
                    d["xi"] = np.zeros_like( d["res"] ) + f.getncattr("xi")
 
                
                if "Q" in attrs:
                    d["Q"] = np.zeros_like( d["res"] ) + f.getncattr("Q")
                
                

                
            data_by_file.append(d)
            
            load_anything = True

    if load_anything == False:
        raise Error("No data is detected.")


    loaded_varnames = loaded_varnames.copy()

    for varname in ["xi", "Q"]:
        if varname not in loaded_varnames:
            loaded_varnames.append(varname) 

    # Processing data
    concat_dataset = {}
    for varname in loaded_varnames:
        tmp_data = None
        for i, d in enumerate(data_by_file):
            if i == 0:
                tmp_data = d[varname]
            else:
                tmp_data = np.concatenate((tmp_data, d[varname]), axis=0)
            
        concat_dataset[varname] = tmp_data


    # Remove nan
    if np.isfinite(residue_threshold):
        valid_idx = concat_dataset["res"] < residue_threshold
        valid_idx[0] = True # A hard fix. We need to plot the first point. This should be better fixed in ZATOM code when producing output
        for varname in loaded_varnames:
            
            selector = [valid_idx,]
            for _ in range(len(concat_dataset[varname].shape)-1):
                selector.append(slice(None))
            
            concat_dataset[varname] = concat_dataset[varname][*selector]
    
    if load_coor:
        
        coor["dz_T"] = coor["z_W"][:-1] - coor["z_W"][1:]
        coor["dz_W"] = np.zeros(coor["z_W"].shape)
        coor["dz_W"][1:-1] = (coor["dz_T"][:-1] + coor["dz_T"][1:]) / 2.0
        coor["dz_W"][0] = coor["dz_T"][0]
        coor["dz_W"][-1] = coor["dz_T"][-1]

        coor["dy_T"] = coor["y_V"][1:]  - coor["y_V"][:-1]
        coor["dy_V"] = np.zeros(coor["y_V"].shape)
        coor["dy_V"][1:-1] = (coor["dy_T"][:-1] + coor["dy_T"][1:]) / 2.0
        coor["dy_V"][0] = coor["dy_T"][0]
        coor["dy_V"][-1] = coor["dy_T"][-1]

        coor["dA_T"] = coor["dy_T"][:, None] * coor["dz_T"][None, :]
        coor["cos_lat_T"] = np.cos(coor["y_T"])
        coor["cos_lat_V"] = np.cos(coor["y_V"])
        coor["sin_lat_T"] = np.sin(coor["y_T"])
        coor["sin_lat_V"] = np.sin(coor["y_V"])
        
        coor["y_V"] *= 180.0/np.pi
        coor["y_T"] *= 180.0/np.pi
        

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
