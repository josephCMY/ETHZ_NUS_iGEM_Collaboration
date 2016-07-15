####################################################################################
# FUNCTION TO GENERATE A RANDOM ARRAY OF LOCATIONS FOR THE HOLES                   #
####################################################################################
import numpy as np


def pore_generator(ensemble, nSpecies, pore, nZ, nP):
    pore_z = (np.arange(ensemble)).tolist()
    pore_p = (np.arange(ensemble)).tolist()

    for n1 in np.arange(ensemble):
        pore_z[n1] = (np.arange(nSpecies)*0).tolist()
        pore_p[n1] = (np.arange(nSpecies)*0).tolist()

        for n2 in np.arange(nSpecies):
            pore_n = pore[n2]
            if pore_n!=0:
                pore_z[n1][n2] = np.random.random_integers(1,high=nZ, size= pore_n)
                pore_p[n1][n2] = np.random.random_integers(1,high=nP, size= pore_n)

            else:
                pore_z[n1][n2] = 0
                pore_p[n1][n2] = 0

    return pore_z, pore_p
