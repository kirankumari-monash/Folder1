#!/usr/bin/puthon 
import os
#import pandas as pd 
import numpy as np
#import netCDF4 as nc4
from scipy.io import netcdf
r_ens = 0
for traj in range(5):
  f = netcdf.NetCDFFile('net_%s.nc' %(traj+101),'r')
#########################################################
#####   READING FROM NetCDF
#########################################################

  Nsample = f.dimensions['No_of_samples']
#  print Nsample
  Nsample = Nsample -1
  Ndim = f.dimensions['Ndim']
#  print Ndim
  Ndim = Ndim -1
#  print Ndim
  NBeads = f.dimensions['NBeads']
#  print NBeads
  NBeads = NBeads - 1
  timenet = f.variables['Time']
  time = timenet[:]*1
#  print time
#  print time[1]
  confinet = f.variables['configuration']
  confi = confinet[:]*1
#  print confi
#  print confi[0, 0, :]
#  print "changing"
#  print confi[NBeads, Ndim, Nsample]

  gradnet = f.variables['Gradient']
  grad = gradnet[:]*1
 
##########################################################
#####    PROPERTY CALCULATION
##########################################################
################ Re2  #####################################
  dist= np.zeros((Ndim, Nsample))
  r_sqs= np.zeros((Nsample))
  rg_sqs = np.zeros((NBeads, Ndim, Nsample))

  dist = confi[NBeads,:,:]-confi[0,:,:]
  dist = np.square(dist)
  r_sqs = np.sum(dist[:],axis=0)
  print 're2_%s'%(traj+1), r_sqs

###############   Rg #####################################


  r_ens = r_ens+ r_sqs
##########################################################
######### Ensemble averaging 
##########################################################
  

r_ens = r_ens/(traj+1)
print 're2_ensemble', r_ens

