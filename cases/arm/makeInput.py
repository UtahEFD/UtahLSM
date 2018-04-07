#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
import os,sys
import pylab as pl

####################################
# UtahLSM test case: June 2016     #
# ARM SGP site in Lamont, Oklahoma #
####################################

# soil data taken from STAMP manual
# https://www.arm.gov/publications/tech_reports/handbooks/stamp_handbook.pdf
# soil levels are 5cm, 10cm, 20cm, 50cm, 75/100cm
# soil types are SiL, SiL, C, CL, CL

# soil sensor depths for STAMP
measZ = [0.0, 0.05, 0.10, 0.20, 0.50, 0.75, 1.0]
nSoil = len(measZ)

# soil temperature
measSoilTemp = [301.95,301.95,300.45,297.52,293.52,291.90,291.95]

# soil moisture
measSoilMois = [0.1115,0.1115,0.1707,0.3603,0.3276,0.1549,20.8900]

# soil type from USDA 11-category + peat
# Lamont - 0-19: silty loam, 20-34: clay, 35-60: clay loam
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

soilType     = [4,4,4,4,11,8,8,8]

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
met = nc.MFDataset('observations/sgpmetE13.b1.*.cdf')

# grab variables
tt = met.variables['time'][:]
ws = met.variables['wspd_vec_mean'][:]
wd = met.variables['wdir_vec_mean'][:]
pt = met.variables['temp_mean'][:] + 273.15
pa = met.variables['atmos_pressure'][:]*10
pv = met.variables['vapor_pressure_mean'][:]*10

# fix times to be continuous across days
#for ttt in range(1,31):
#    ti = ttt*1440
 #   tt[ti::] = tt[ti::]+tt[ti-1]+60

# compute wind components
uc = -ws*np.sin(wd*np.pi/180)
vc = -ws*np.cos(wd*np.pi/180)
print(len(uc))
# compute water vapor mixing ratio
qv = 0.622 * pv / (pa-pv)

########################################
# Read net radiation data for 20121024 #
########################################

# open CO2FLX data
data = nc.MFDataset('observations/sgpseb*')

# grab variables
tt   = data.variables['time_offset'][:]
rNet = data.variables['net_radiation'][:]

nt   = len(tt)
tt   = np.arange(0,nt*1800,1800)
tt1m = np.arange(0,30*86400,60)
# interpolate from 30-minute to 1-minute frequency to match MET
#tt1m   = np.arange(tt[0],tt[-1]+1800,60)
rNet1m = np.interp(tt1m,tt,rNet)

tc = 0
badT = []
for r in pt:
	if (str(r)=='--'):
		badT.append(tc)
	tc+=1
for b in badT:
	pt[b] = pt[badT[0]-1]

tc = 0
badQ = []
for r in qv:
	if (str(r)=='--'):
		badQ.append(tc)
	tc+=1

for b in badQ:
	qv[b] = qv[badQ[0]-1]

##############################
# Write all time series data #
##############################
of = open('inputMetr.dat','w')
os1 = '{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s}\n'.format('atm_u','atm_v','atm_T','atm_q', 'R_net')
of.write(os1)
for t in range(len(tt1m)):
	os1 = "{0:15.8E}  {1:2.8E}  {2:2.8E}  {3:2.8E}  {4:2.8E}\n".format(uc[t], vc[t], pt[t], qv[t], rNet1m[t])
	of.write(os1)
of.close()