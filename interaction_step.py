import numpy as np
import time
import multiprocessing as mp
import math as m
from netCDF4 import Dataset

####################################################################################
# INTERACTION STEP                                                                 #
####################################################################################
def interaction_step( conc_ensemble_00, interact_coeff, interact_arr, nSpecies, interact_dim1, n1, n21 ):

    # Calling and reshaping the array from the buffer
    interact0 = np.frombuffer( interact_arr.get_obj() )
    interact  = interact0.reshape(interact_dim1)

    holder = np.array([ conc_ensemble_00[n1, n21] * interact_coeff[n21, n22]    \
               for n22 in np.arange(nSpecies) ])

    interact[n1,n21,:,:,:] = np.sum(holder, axis=0)

    return


