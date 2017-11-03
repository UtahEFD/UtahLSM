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
measSoilTemp = [297.18, 295.68, 295.52, 294.75, 293.69, 293.35, 292.66]

# soil moisture
measSoilMois = [0.078, 0.135, 0.2535, 0.2528, 0.3666, 0.3272, 0.3230]

# soil properties
# Lamont - 0-19: silty loam, 20-34: clay, 35-60: clay loam
porosity     =  np.array([0.485,0.485,0.485,0.485,0.482,0.476,0.476])
satPotential =  np.array([0.786,0.786,0.786,0.786,0.405,0.630,0.630])
satHydrCond  =  np.array([7.2,7.2,7.2,7.2,1.3,2.5,2.5])
soilExponent =  np.array([5.3,5.3,5.3,5.3,11.4,8.52,8.52])
heatCapSoil  =  np.array([1.27,1.27,1.27,1.27,1.09,1.23,1.23])

# proper units and sign
satPotential = -satPotential
satHydrCond  = satHydrCond * 1e-6
heatCapSoil  = heatCapSoil * 1e6

# write output
of = open('inputSoil.dat','w')
os = '{0:^15} {1:^15s} {2:^15s} {3:^15s} {4:^19s} {5:^11s} {6:^16s} {7:^15s}\n'
os = os.format('zSoil','tempSoil','moisSoil','porosity','moisPotSat','hydCondSat','b','Ci')
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
os = '{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s}\n'.format('uComp','vComp','pTemp','wvMix', 'rNet')
of.write(os)
for t in range(len(tt1m)):
	os = "{0:15.8E}  {1:2.8E}  {2:2.8E}  {3:2.8E}  {4:2.8E}\n".format(uc[t], vc[t], pt[t], qv[t], rNet1m[t])
	of.write(os)
of.close()