#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

##########################################
# Manual entry of soil data for 20121024 #
##########################################

# soil sensor depths and interpolation levels
measZ = [0.0, 0.025, 0.05, 0.15, 0.25, 0.35, 0.60]
nSoil = len(measZ)

# soil temperature
measSoilTemp = [295.68, 295.68, 295.52, 294.75, 293.69, 293.35, 292.66]

# soil moisture
measSoilMois = [0.135, 0.135, 0.2535, 0.2528, 0.3666, 0.3272, 0.3230]

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
met = nc.MFDataset('sgpmetE13.b1.*.cdf')

# grab variables
tt = met.variables['time'][:]
ws = met.variables['wspd_vec_mean'][:]
wd = met.variables['wdir_vec_mean'][:]
pt = met.variables['temp_mean'][:] + 273.15
pa = met.variables['atmos_pressure'][:]
pv = met.variables['vapor_pressure_mean'][:]

# fix times to be continuous across days
tt[1440::] = tt[1440::]+tt[1439]+60

# compute wind components
uc = -ws*np.sin(wd*np.pi/180)
vc = -ws*np.cos(wd*np.pi/180)

# compute water vapor mixing ratio
qv = 0.622 * pv / (pa-pv)

########################################
# Read net radiation data for 20121024 #
########################################

# open CO2FLX data
data = nc.MFDataset('sgp30co2flx4mmetC1.b1.*')

# grab variables
tt   = data.variables['time_offset'][:]
lwD  = data.variables['r_down_long_hemisp'][:]
lwU  = data.variables['r_up_long_hemisp'][:]
swD  = data.variables['r_down_short_hemisp'][:]
swU  = data.variables['r_up_short_hemisp'][:]
rNet = data.variables['r_net'][:]
G    = data.variables['mean_g_soil'][:]
H    = data.variables['h'][:]
LE   = data.variables['le'][:]

# fix times to be continuous across days
tt[48::] = tt[48::]+tt[47]+1800
tt=tt-tt[0]

# interpolate from 30-minute to 1-minute frequency to match MET
tt1m   = np.arange(tt[0],tt[-1]+1800,60)
lwD1m  = np.interp(tt1m,tt,lwD)
lwU1m  = np.interp(tt1m,tt,lwU)
swD1m  = np.interp(tt1m,tt,swD)
swU1m  = np.interp(tt1m,tt,swU)
rNet1m = np.interp(tt1m,tt,rNet)
G1m    = np.interp(tt1m,tt,G)
H1m    = np.interp(tt1m,tt,H)
LE1m   = np.interp(tt1m,tt,LE)

##############################
# Write all time series data #
##############################
of = open('inputMetr.dat','w')
os = '{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s}\n'.format('atm_u','atm_v','atm_T','atm_q', 'R_net')
of.write(os)
for t in range(len(tt1m)):
	os = "{0:15.8E}  {1:2.8E}  {2:2.8E}  {3:2.8E}  {4:2.8E}\n".format(uc[t], vc[t], pt[t], qv[t], rNet1m[t])
	of.write(os)
of.close()