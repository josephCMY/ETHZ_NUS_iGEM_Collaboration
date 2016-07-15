####################################################################################
# START OF LIBRARY                                                                 #
####################################################################################

import numpy as np
import math as m
import time
import writeout as wo

####################################################################################
# GENERATING GLOBAL FIELDS                                                         #
####################################################################################

def initialize(outfile, Lr, Lz, nZ, nP, nR, initial_concentration, ensemble, tOut, nSpecies):

    print "    Preparing coordinate grid."
    cylR_int, cylR_ext = create_cylinder_coordR( Lr, nZ, nP, nR )
    cylP_int, cylP_ext = create_cylinder_coordP(     nZ, nP, nR )
    cylZ_int, cylZ_ext = create_cylinder_coordZ( Lz, nZ, nP, nR )

    print "    Writing out the coordinate grid."
    wo.write3d( outfile, "cylR_int", cylR_int, "f8", "z index", "phi index"\
               , "r index")
    wo.write3d( outfile, "cylP_int", cylP_int, "f8", "z index", "phi index"\
               , "r index")
    wo.write3d( outfile, "cylZ_int", cylR_int, "f8", "z index", "phi index"\
               , "r index")

    print "    Initializing staggering grids."
    cyl_rL, cyl_rR, cyl_r0, cyl_pL, cyl_pR, cyl_p0, cyl_zL, cyl_zR, cyl_z0 \
    = create_stagger_indices( nR, nP, nZ )

    print "    Preparing concentration storage grid."
    conc_all_container = np.zeros([tOut+1, ensemble, nSpecies, nZ, nP, nR], dtype = "f8")


    print "    Inputting initial concentration field."  
    initial_concentration1 = np.array(initial_concentration)

    for i in np.arange(ensemble):
        for j in np.arange(nSpecies):
            conc_all_container[0,i,j,:,:,:] = initial_concentration1[i,j]

    print "    Preparing coefficients for calculations."
    cyl_invR1 = np.ones([nZ, nP, nR], dtype = "f8")/cylR_int
    cyl_invR2 = cyl_invR1*cyl_invR1

    cyl_rLL \
    = cylR_ext[ cyl_z0, cyl_p0, cyl_r0 ] - cylR_ext[ cyl_z0, cyl_p0, cyl_rL ]
    cyl_rRL \
    = cylR_ext[ cyl_z0, cyl_p0, cyl_rR ] - cylR_ext[ cyl_z0, cyl_p0, cyl_r0 ]
    cyl_rCo \
    = cyl_rLL + cyl_rRL

    cyl_rA1 = 2 *np.reciprocal( cyl_rLL * cyl_rCo )
    cyl_rA2 = 2 *np.reciprocal( cyl_rRL * cyl_rCo )
    cyl_rA3 = 2 *np.reciprocal( cyl_rRL * cyl_rLL )

    cyl_rB1 = (cyl_rRL / ( cyl_rLL * ( cyl_rLL + cyl_rRL ) )) * cyl_invR1
    cyl_rB2 = (cyl_rLL / ( cyl_rRL * ( cyl_rLL + cyl_rRL ) )) * cyl_invR1
    cyl_rB3 = (( cyl_rLL - cyl_rRL ) / ( cyl_rLL * cyl_rRL )) * cyl_invR1

    return cyl_rA1, cyl_rA2, cyl_rA3, cyl_rB1, cyl_rB2, cyl_rB3, cyl_invR2, conc_all_container, cyl_rL, cyl_rR, cyl_r0, cyl_pL, cyl_pR, cyl_p0, cyl_zL, cyl_zR, cyl_z0


####################################################################################
# SETTING UP INITIAL CONCENTRATION FIELDS                                          #
####################################################################################

def initial_bacteria_concentration( ensemble, nZ, nP, nR, conc_ini ):
    conc_00 = np.zeros( [ensemble, nZ, nP, nR], dtype = "f4")
    for i in np.arange(ensemble):
        conc_00[i,:,:,:] = conc_ini[i]

    return conc_00


####################################################################################
# CREATION OF THE CYLINDRICAL COMPONENT COORDINATE GRIDS                           #
####################################################################################

def create_cylinder_coordR( Lr, nZ, nP, nR):
    coordR_int = np.ones( [  nZ  ,   nP  ,   nR  ], dtype = "f8" )
    coordR_ext = np.ones( [nZ + 2, nP + 2, nR + 2], dtype = "f8" )

    dr = float(Lr)/float(nR)

    for i in np.arange(nR):
        coordR_int[:,:,i  ] = dr * ( i + 1 )
        coordR_ext[:,:,i+1] = dr * ( i + 1 )

    coordR_ext[:,:,0   ] = dr * -1.0
    coordR_ext[:,:,nR+1] = dr * (nR + 1)

    return coordR_int, coordR_ext



def create_cylinder_coordP( nZ, nP, nR):
    coordP_int = np.ones( [  nZ  ,   nP  ,   nR  ], dtype = "f8" )
    coordP_ext = np.ones( [nZ + 2, nP + 2, nR + 2], dtype = "f8" )

    dp = 2.0*m.pi / float(nP)

    for j in np.arange(nP):
        coordP_int[:,j  ,:] = dp * j
        coordP_ext[:,j+1,:] = dp * j

    coordP_int[:,0   ,:] = -dp
    coordP_ext[:,nP+1,:] = 2*m.pi

    return coordP_int, coordP_ext



def create_cylinder_coordZ( Lz, nZ, nP, nR):
    coordZ_int = np.ones( [  nZ  ,   nP  ,   nR  ], dtype = "f8" )
    coordZ_ext = np.ones( [nZ + 2, nP + 2, nR + 2], dtype = "f8" )

    dz = float(Lz)/float(nZ)

    for k in np.arange(nZ):
        coordZ_int[k  ,:,:] = dz * k
        coordZ_ext[k+1,:,:] = dz * k

    coordZ_int[0   ,:,:] = -dz
    coordZ_ext[nZ+1,:,:] =  Lz

    return coordZ_int, coordZ_ext




####################################################################################
# SETTING UP INDICE FUNCTION FOR BOTH TYPES OF GRIDS                               #
####################################################################################

def create_stagger_indices( nR, nP, nZTH ):

    shp = [nZTH, nP, nR]

    rL = np.ones(shp, dtype="i4")
    rR = np.ones(shp, dtype="i4")
    r0 = np.ones(shp, dtype="i4")

    for i in np.arange(shp[2], dtype="i4"):
        rL[:,:,i] = i
        rR[:,:,i] = i + 2
        r0[:,:,i] = i + 1


    pL = np.ones(shp, dtype="i4")
    pR = np.ones(shp, dtype="i4")
    p0 = np.ones(shp, dtype="i4")

    for j in np.arange(shp[1], dtype="i4"):
        pL[:,j,:] = j  
        pR[:,j,:] = j + 2
        p0[:,j,:] = j + 1


    zL = np.ones(shp, dtype="i4")
    zR = np.ones(shp, dtype="i4")
    z0 = np.ones(shp, dtype="i4")

    for k in np.arange(shp[0], dtype="i4"):
        zL[k,:,:] = k 
        zR[k,:,:] = k + 2
        z0[k,:,:] = k + 1

    return rL, rR, r0, pL, pR, p0, zL, zR, z0


