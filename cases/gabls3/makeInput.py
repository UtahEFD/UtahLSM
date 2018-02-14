#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

##########################################
# Manual entry of soil data for 20121024 #
##########################################

# soil sensor depths and interpolation levels
measZ = [0.0,0.005,0.01,0.02,0.04,0.06,0.08,0.12,0.20,0.30,0.50,1.00,2.00]
nSoil = len(measZ)

# soil temperature
measSoilTemp = [287.450,287.250,288.677,289.777,290.335,291.340,291.080,291.181,290.520,289.958,288.689,285.350,283.150]

# soil moisture
measSoilMois = [0.24,0.247,0.269,0.274,0.310,0.330,0.330,0.360,0.450,0.470,0.470,0.570,0.570]

# soil properties
# Lamont - 0-19: silty loam, 20-34: clay, 35-60: clay loam
porosity     =  np.array([0.47,0.47,0.47,0.47,0.47,0.47,0.673,0.863,0.863,0.863,0.863,0.863,0.863])
satPotential =  np.array([0.405,0.405,0.405,0.405,0.405,0.405,0.381,0.356,0.356,0.356,0.356,0.356,0.356])
satHydrCond  =  np.array([1.3,1.3,1.3,1.3,1.3,1.3,4.65,8.0,8.0,8.0,8.0,8.0,8.0])
soilExponent =  np.array([11.4,11.4,11.4,11.4,11.4,11.4,9.575,7.75,7.75,7.75,7.75,7.75,7.75])
heatCapSoil  =  np.array([1.090,1.090,1.090,1.090,1.090,1.090,0.965,0.840,0.840,0.840,0.840,0.840,0.840])

# proper units and sign
satPotential = -satPotential
satHydrCond  = satHydrCond * 1e-6
heatCapSoil  = heatCapSoil * 1e6

# write output
of = open('inputSoil.dat','w')
os = '{0:^15} {1:^15s} {2:^15s} {3:^15s} {4:^19s} {5:^11s} {6:^16s} {7:^15s}\n'
os = os.format('soil_z','soil_T','soil_q','porosity','psi_nsat','K_nsat','b','Ci')
of.write(os)
for z in range(nSoil):
	os = "{0:15.8E}  {1:2.8E}  {2:2.8E}  {3:2.8E}  {4:2.8E}  {5:2.8E}  {6:2.8E}  {7:2.8E}\n"
	os = os.format(measZ[z],measSoilTemp[z],measSoilMois[z],porosity[z],satPotential[z],satHydrCond[z],soilExponent[z],heatCapSoil[z])
	of.write(os)
of.close()

################################################
# Read met tower data for u,v,T,q for 20121024 #
################################################

# open MET data
met = np.array([float(d) for d in open('utahLES.dat').read().split()])

# grab variables
tt = met[0::7]
uc = met[1::7]
vc = met[2::7]
pt = met[4::7]
qv = met[5::7]

##############################
# Write all time series data #
##############################
of = open('inputMetr.dat','w')
os = '{0:^15s} {1:^15s} {2:^15s} {3:^15s}\n'.format('atm_u','atm_v','atm_T','atm_q')
of.write(os)
for t in range(len(tt)):
	os = "{0:15.8E}  {1:15.8E}  {2:15.8E}  {3:15.8E}\n".format(uc[t], vc[t], pt[t], qv[t])
	of.write(os)
of.close()