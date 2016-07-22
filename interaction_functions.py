import numpy as np
import time
import multiprocessing as mp
import math as m
from netCDF4 import Dataset

def f0( conc_00, n1 ):
    output = np.power( conc_00[n1,0,:,:,:] , 1) * np.power( conc_00[n1,1,:,:,:] , 1) * (-1.000000) + np.power( conc_00[n1,2,:,:,:] , 1) * (2.000000)
    return output

def f1( conc_00, n1 ):
    output = np.power( conc_00[n1,0,:,:,:] , 1) * np.power( conc_00[n1,1,:,:,:] , 1) * (1.000000) + np.power( conc_00[n1,2,:,:,:] , 1) * (-2.000000)
    return output

def f2( conc_00, n1 ):
    output = np.power( conc_00[n1,0,:,:,:] , 1) * np.power( conc_00[n1,1,:,:,:] , 1) * (-1.000000) + np.power( conc_00[n1,2,:,:,:] , 1) * (2.000000)
    return output

