####################################################################################
# MAIN RUNNING SCRIPT                                                              #
####################################################################################

print ""
print "########################################################"
print "# BACTERIAL CYLINDER REACTION-DIFFUSION MODEL          #"
print "########################################################"
print ""

import time

print "%.3f seconds --- Starting up simulation" % (0.0)
print ""
print "%.3f seconds --- Importing requisite libraries" % (0.0)
t00 = time.time()

import numpy as np
import math as m
from netCDF4 import Dataset

import initialize as ini
import solvers    as sol
import writeout   as wo
import pore_generator as po
print ""
print "%.3f seconds --- Libraries imported" % (time.time() - t00)


####################################################################################
# SECTION 0: NAMELIST VARIABLES FOR USE                                            #
####################################################################################
print ""
print "%.3f seconds --- Reading namelist" % (time.time() - t00)

execfile("namelist.input")

print ""
print "%.3f seconds --- Generating pore locations." % (time.time() - t00)

pore_z, pore_p = po.pore_generator(ensemble, nSpecies, pore, nZ, nP)

print "%.3f seconds --- Saving pore locations." % (time.time() - t00)

pore_string = ""

for n1 in np.arange(ensemble):
    for n2 in np.arange(nSpecies):
        pore_string += "ensemble: " + str(n1) + " , species: " + names[n2] + "\n"
        pore_string += "z-indices: " + str(pore_z[n1][n2]) + "\n"
        pore_string += "p-indices: " + str(pore_p[n1][n2]) + "\n"
        pore_string += "\n"        

textfile = open("pore_information.txt","w")
textfile.write(pore_string)
textfile.close()

print ""
print "%.3f seconds --- Creating output netCDF file." % (time.time() - t00)

outfile = wo.create_ncfile( fout_name , tOut, nZ, nP, nR, ensemble, nSpecies )

####################################################################################
# SECTION 1: INITIALIZING THE SIMULATION                                           #
####################################################################################

print ""
print "%.3f seconds --- Initializing requisite grids" % (time.time() - t00)

cyl_rA1, cyl_rA2, cyl_rA3, cyl_rB1, cyl_rB2, cyl_rB3, cyl_invR2, conc_all, cyl_rL\
, cyl_rR, cyl_r0, cyl_pL, cyl_pR, cyl_p0, cyl_zL, cyl_zR, cyl_z0 \
= ini.initialize(outfile, Lr, Lz, nZ, nP, nR, initial_concentration, ensemble, tOut, nSpecies)


dp = 2.0 * m.pi/nP
dr = Lr/nR
dz = Lz/nZ

dnT = nT/tOut
store_tt_arr = np.arange( 1 , tOut + 1 , dtype="i16" ) * dnT
tt_marker = 0
store_tt = store_tt_arr[tt_marker]

####################################################################################
# SECTION 2: RUNNING THE SIMULATION INTEGRATOR                                     #
####################################################################################

print ""
print "%.3f seconds --- Proceeding with time integration" % (time.time() - t00)

conc_all = sol.euler_forward_integrator(conc_all, nR, nP, nZ, t00, cyl_rL, cyl_rR, cyl_r0, cyl_pL, cyl_pR, cyl_p0, cyl_zL, cyl_zR, cyl_z0, cyl_rA1, cyl_rA2, cyl_rA3, cyl_rB1, cyl_rB2, cyl_rB3, dp, dr, dz, store_tt_arr, tt_marker, store_tt, nT, cyl_invR2, D, dt, tOut, ensemble, nSpecies, pore_z, pore_p, outer_concentration, pore, interact_coeff)

print ""
print "%.3f seconds --- Time integration completed" % (time.time() - t00)

print ""
print "%.3f seconds --- Storing results" % (time.time() - t00)

print "    Storing results from all concentration fields"

for i in np.arange(ensemble):
    for j in np.arange(nSpecies):
        varname = "ensemble member: " + str(i) + ", species: " + names[j]
        wo.write4d( outfile, varname, conc_all[:,i,j,:,:,:], "f4", "t_output", "z index", "phi index", "r index" )


print "    Storing results from averaged fields"

for i in np.arange(nSpecies):
    varname = "averaged " + names[i] + " fields"
    variable = np.mean( conc_all[:,:,i,:,:,:], axis = 1)
    wo.write4d( outfile, varname, variable, "f4", "t_output", "z index", "phi index", "r index" )


print ""
print "%.3f seconds --- Simulation completed" % (time.time() - t00)


