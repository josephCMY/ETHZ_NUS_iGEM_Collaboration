import numpy as np
import interaction_functions as ifun

####################################################################################
# INTERACTION STEP                                                                 #
####################################################################################
def interaction_step( conc_00, interact_arr, n1, n2, fInd, interact_dim1 ):

    # Calling and reshaping the array from the buffer
    interact0 = np.frombuffer( interact_arr.get_obj() )
    interact  = interact0.reshape(interact_dim1)

    exec "output = ifun.f%d(conc_00, n1)" % fInd
    
    interact[n1,n2,:,:,:] = output[:,:,:]

    return


