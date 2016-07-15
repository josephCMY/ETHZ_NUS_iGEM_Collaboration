import numpy as np
from math import pi
import scipy as sc
from matplotlib import pyplot as plt
from netCDF4 import Dataset

execfile('namelist.input')
print fout_name
f_nc = Dataset(fout_name,"r")

vlist =  [var for var in f_nc.variables]

dp = 2*pi/nP
dz = Lz/nZ
dr = Lr/nR

vol = Lz*pi*(Lr*Lr)

storage = np.zeros([tOut+1, nSpecies, ensemble])

rCoord0 = f_nc["cylR_int"]
rCoord1 = np.zeros([tOut+1, nZ, nP, nR+1],dtype="f8")

for k in np.arange(tOut +1):
    rCoord1[k,:,:,1:(nR+1)] = rCoord0[0,:,:,:]

rCoord = np.zeros([tOut+1, nZ, nP+1, nR+1],dtype="f8")
rCoord[:,:,1:(nP+1),:] = rCoord1[:,:,:,:]
rCoord[:,:,0,:] = rCoord1[:,:,0,:]

vol = np.trapz(rCoord,dx=dr,axis=-1) #x=rCoord,  axis=-1)
vol = np.trapz(vol, dx=dp, axis=-1)
vol = np.trapz(vol, dx=dz, axis=-1)
vol = np.mean(vol)

for i in np.arange(nSpecies):
    for j in np.arange(ensemble):
        varname = "ensemble member: " + str(j) + ", species: " + names[i]
        conc0 = f_nc[varname]
        conc01 = np.zeros([tOut+1, nZ, nP, nR+1], dtype="f8")
        conc1 = np.zeros([tOut+1, nZ, nP+1, nR+1], dtype="f8")

        conc01[:,:,:,1:(nR+1)] = conc0[:,:,:,:]
        conc1[:,:,1:(nP+1),:] = conc01[:,:,:,:]
        conc1[:,:,0,:] = conc01[:,:,nP-1,:]
        conc1 = conc1*rCoord
        
        tot = np.trapz(conc1,dx=dr,axis=-1) #x=rCoord,  axis=-1)
        tot = np.trapz(tot, dx=dp, axis=-1)
        tot = np.trapz(tot, dx=dz, axis=-1)
        aveC = tot/vol
        storage[:,i,j] = aveC


xAxis = np.arange(91)
for i in np.arange(nSpecies):
    plt.plot(xAxis, np.mean(storage[:,i,:], axis=-1), 'k')


colors = ['g','y','c','r','b']

for i in np.arange(nSpecies):
    for j in np.arange(ensemble):
        plt.plot(xAxis, storage[:,i,j], colors[i])

plt.show()
