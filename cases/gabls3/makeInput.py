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

# soil type from USDA 11-category + peat
#  1 = sand
#  2 = loamy sand
#  3 = sandy loam
#  4 = silty loam
#  5 = loam
#  6 = sandy clay loam
#  7 = silty clay loam
#  8 = clay loam
#  9 = sandy clay
# 10 = silty clay
# 11 = clay
# 12 = peat
soilType = np.array([11,11,11,11,11,11,12,12,12,12,12,12,12])

# write output
of = open('inputSoil.dat','w')
os = '{0:^15} {1:^15s} {2:^15s} {3:^15s}\n'
os = os.format('soil_z','soil_type','soil_T','soil_q')
of.write(os)
for z in range(nSoil):
	os = "{0:15.8E}  {1:2.8E}  {2:2.8E}  {3:2.8E}\n"
	os = os.format(measZ[z],soilType[z],measSoilTemp[z],measSoilMois[z])
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