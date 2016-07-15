####################################################################################
# START OF LIBRARY                                                                 #
####################################################################################

import numpy as np
import math as m
from netCDF4 import Dataset

# WRITING A GENERIC 5D VARIABLE TO THE NETCDF OUTPUT FILE
def write5d( outfile, varname, data, dstr, d5, d4, d3, d2, d1 ):

    output = outfile.createVariable( varname, dstr, (d5, d4, d3, d2, d1,))

    output[:, :, :, :, :] = data

    return 

# WRITING A GENERIC 4D VARIABLE TO THE NETCDF OUTPUT FILE
def write4d( outfile, varname, data, dstr, d4, d3, d2, d1 ):

    output = outfile.createVariable( varname, dstr, (d4, d3, d2, d1,))

    output[:, :, :, :] = data

    return 


# WRITING A GENERIC 3D VARIABLE TO THE NETCDF OUTPUT FILE
def write3d( outfile, varname, data, dstr, d3, d2, d1 ):

    output = outfile.createVariable( varname, dstr, ("const t", d3, d2, d1,))

    output[:, :, :, :] = data

    return 


# CREATING THE NETCDF FILE
def create_ncfile( fout_name, tOut,nZ, nP, nR, ensemble, nSpecies ):

    # Generic dimensions
    outfile = Dataset( fout_name, "w", format="NETCDF4")
    outfile.createDimension("t_output",tOut+1)

    # Constant dimensions
    outfile.createDimension("const r",1)
    outfile.createDimension("const p",1)
    outfile.createDimension("const z",1)
    outfile.createDimension("const t",1)

    # Cylindrical dimensions
    outfile.createDimension("r index", nR)
    outfile.createDimension("phi index", nP)
    outfile.createDimension("z index", nZ)
    outfile.createDimension("r-ext index", nR+2)
    outfile.createDimension("phi-ext index", nP+2)
    outfile.createDimension("z-ext index", nZ+2)

    outfile.createDimension("bacteria ensemble",ensemble)
    outfile.createDimension("species", nSpecies)

    return outfile


