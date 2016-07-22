import numpy as np
import time
import math as m
import multiprocessing as mp
import diffusion_step as diff
from interaction_step import interaction_step as iStep
from netCDF4 import Dataset


####################################################################################
# EULER FORWARD STEPPING INTEGRATOR [ ERROR TERM: O(dt) ]                          #
####################################################################################

def euler_forward_integrator( conc_all, nR, nP, nZ, t00, cyl_rL, cyl_rR, cyl_r0, cyl_pL, cyl_pR, cyl_p0, cyl_zL, cyl_zR, cyl_z0, cyl_rA1, cyl_rA2, cyl_rA3, cyl_rB1, cyl_rB2, cyl_rB3, dp, dr, dz, store_tt_arr, tt_marker, store_tt, nT, cyl_invR2, D, dt, tOut, ensemble, nSpecies, pore_z, pore_p, outer_concentration, pore, iFunc_ind):

    conc_00 = np.zeros([ensemble, nSpecies, nZ+2, nP+2, nR+2],dtype="f4")

    conc_00[ :, :, cyl_z0, cyl_p0, cyl_r0 ] = conc_all[ 0, : , : , : , :, : ]

    print "%.3f seconds --- EULER FORWARD STEP INTEGRATOR ACTIVATING" \
    % (time.time() - t00)

    # EULER FORWARD STEPPING TIME INTEGRATOR. THIS IS THE SIMPLEST AND MOST CRUDE
    # OF SCHEMES
    tstep = 0 
    t11 = time.time()

    # Generating the shared memory for the diffusion step
    diff_dim0 = ensemble*nSpecies*(nZ+2)*(nP+2)*(nR+2)
    diff_dim1 = (ensemble, nSpecies, nZ+2, nP+2, nR+2)
    diff_arr = mp.Array('d', diff_dim0) 

    # Generating the shared memory for the interaction step
    interact_dim0 = ensemble*nSpecies*(nZ+2)*(nP+2)*(nR+2)
    interact_dim1 = (ensemble, nSpecies, nZ+2, nP+2, nR+2)
    interact_arr  = mp.Array('d', interact_dim0)

    shp = np.shape(conc_00[0,:,:,:,:])
    blank = np.zeros(shp, dtype="f8")
 
    while tstep <= nT:

        # Handling the diffusion processes
        processes = [mp.Process( target = diff.diffusion_step\
                            , args = ( conc_00, nR, nP, nZ, t00, cyl_rL, cyl_rR   \
                                      , cyl_r0, cyl_pL, cyl_pR, cyl_p0, cyl_zL    \
                                      , cyl_zR, cyl_z0, cyl_rA1, cyl_rA2, cyl_rA3 \
                                      , cyl_rB1, cyl_rB2, cyl_rB3, dp, dr, dz     \
                                      , store_tt_arr, tt_marker, store_tt, nT, n1 \
                                      , cyl_invR2, D, dt, tOut, diff_arr, n2      \
                                      , diff_dim1, pore_z[n1][n2], pore_p[n1][n2]  \
                                      , outer_concentration[n2], pore[n2] )\
                          )\
                      for n1 in np.arange(ensemble) for n2 in np.arange(nSpecies)]

        for p in processes:
            p.start()

        for p in processes:
            p.join(timeout=None)

        diffuse0 = np.frombuffer(diff_arr.get_obj())
        diffuse = diffuse0.reshape(diff_dim1)
        
        # Handling the species interaction step
        processes = [mp.Process( target =  iStep
                                , args = ( conc_00, interact_arr, n1		\
                                          , iFunc_ind[fInd], fInd, interact_dim1))\
                     for n1 in np.arange(ensemble) \
                     for fInd in np.arange(len(iFunc_ind))]

        for p in processes:
            p.start()

        for p in processes:
            p.join(timeout=None)


        interact0 = np.frombuffer( interact_arr.get_obj() )
        interact  = interact0.reshape(interact_dim1)

#        print np.sum(interact)


        conc_11 = conc_00 + diffuse*dt + interact*dt

        if tstep == store_tt:
            print ""
            print "    %.3f seconds since start of integration" \
                   % (time.time() - t11)
            print "    Storing result of tstep %d" % tstep


            conc_all[tt_marker+1,:,:,:,:,:] = conc_11[:,:,cyl_z0, cyl_p0, cyl_r0]

            if tt_marker < tOut - 1 :
                tt_marker = tt_marker + 1
                store_tt = store_tt_arr[ tt_marker ]

        tstep += 1
        conc_00 = conc_11

    return conc_all
