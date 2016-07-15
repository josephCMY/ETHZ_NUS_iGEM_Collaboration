import numpy as np
import time
import multiprocessing as mp
import math as m
from netCDF4 import Dataset

####################################################################################
# DIFFUSION STEP                                                                   #
####################################################################################
def diffusion_step( conc_00, nR, nP, nZ, t00, cyl_rL, cyl_rR, cyl_r0, cyl_pL, cyl_pR, cyl_p0, cyl_zL, cyl_zR, cyl_z0, cyl_rA1, cyl_rA2, cyl_rA3, cyl_rB1, cyl_rB2, cyl_rB3, dp, dr, dz, store_tt_arr, tt_marker, store_tt, nT, nn, cyl_invR2, D, dt, tOut, diff_arr, n2, diff_dim1, pore_z, pore_p, outer_concentration, nPores):


    # BRINGING IN PERIODIC AND NO-PENETRATION BOUNDARIES
    conc_01 = cyl_ext_points( conc_00[nn,n2,:,:,:], nR, nP, nZ, pore_z, pore_p   \
                             , outer_concentration, nPores)
    
    
    diffuse0 = np.frombuffer(diff_arr.get_obj())
    diffuse = diffuse0.reshape(diff_dim1)
    # DETERMINING THE RATE OF CHANGE DUE TO DIFFUSION	
    diffuse[nn,n2,cyl_z0,cyl_p0,cyl_r0] = cyl_laplacian( conc_01, cyl_rL, cyl_r0 \
                              , cyl_rR   \
                              , cyl_pL, cyl_p0, cyl_pR, cyl_zL, cyl_z0    \
                              , cyl_zR, dr, dp, dz, cyl_rA1, cyl_rA2      \
                              , cyl_rA3, cyl_rB1, cyl_rB2, cyl_rB3        \
                              , cyl_invR2 )*D[n2]

    return


####################################################################################
# Cylinder laplacian                                                               #
####################################################################################

# LAPLACIAN FOR CYLINDER
def cyl_laplacian( c, rL, r0, rR, pL, p0, pR, zL, z0, zR, dr\
                  , dp, dz, rA1, rA2, rA3, rB1, rB2, rB3, invR2):

# DO NOT ERASE THIS COMMENTARY (FOR MODELLER'S REFERENCE)
#    diffR1 = c[z0,p0,rR]*rB2 - rB1*c[z0,p0,rL] - rB3*c[z0,p0,r0]
#    diffR2 = rA1*c[z0,p0,rL] + rA2*c[z0,p0,rR] - rA3*c[z0,p0,r0]
#    diffP2 = (invR2/(dp*dp)) * ( c[z0,pL,r0] + c[z0,pR,r0] - 2 * c[z0,p0,r0] )
#    diffZ2 = ( c[zL,p0,r0] + c[zR,p0,r0] - 2 * c[z0,p0,r0] ) / (dz*dz)

    lpc = c[z0,p0,rL]*(rA1-rB1) + c[z0,p0,rR]*(rB2+rA2)     \
          + (invR2/(dp*dp)) * ( c[z0,pL,r0] + c[z0,pR,r0] ) \
          + ( c[zL,p0,r0] + c[zR,p0,r0]) / (dz*dz)          \
          + c[z0,p0,r0]*( -rB3-rA3 - 2*invR2/(dp*dp) - 2/(dz*dz))



    return lpc


# GENERATING THE CONCENTRATION OF THE EXTERIOR POINT FOR CYLINDER
def cyl_ext_points( arr , nR , nP , nZ, pore_z0, pore_p0, outer_conc, nPores ):
    out = arr[:,:,:]

    # Generating the phi-direction periodic boundaries
    out[:   , 0   , :   ] = out[:   , nP  , :   ]
    out[:   , nP+1, :   ] = out[:   , 1   , :   ]

    # Preparing the r-direction no-penetration boundary
    out[:   , :   , nR+1] = out[:   , :   , nR  ]

    # Preparing the z-direction no-penetration boundary
    out[0   , :   , :   ] = out[1   , :   , :   ]
    out[nZ+1, :   , :   ] = out[nZ  , :   , :   ]

    # including the pores
    if nPores!=0:
        pore_r0 = np.ones(nPores, dtype="i4") * (nR+1)
        out[pore_z0, pore_p0, pore_r0] = outer_conc

    # Setting up the r-direction cross-linking boundary
    out = cyl_cross_linker(out, nP)

    return out


def cyl_cross_linker(arr, nP):
    out = arr[:,:,:]
    for j in np.arange(nP/2):
        out[:,j+1,0] = arr[:,j+1+nP/2,1]
        out[:,j+1+nP/2,0] = arr[:,j+1,1]
    return out


