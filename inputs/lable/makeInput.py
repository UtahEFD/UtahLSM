#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

##########################################
# Manual entry of soil data for 20121024 #
##########################################

# soil sensor depths and interpolation levels
measZ = [0.0, 0.05, 0.15, 0.25, 0.35, 0.60]
#initZ = [0.0, 0.05, 0.025, 0.05, 0.075, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 1.0, 2.0]
nSoil = len(measZ)

# soil temperature
measSoilTemp = [295.05, 293.09, 293.28, 293.22, 293.15, 292.63]
#initSoilTemp = np.interp(initZ, measZ, measSoilTemp)

# soil moisture
measSoilMois = [0.23, 0.25, 0.26, 0.27, 0.33, 0.34]
#initSoilMois = np.interp(initZ, measZ, measSoilMois)

# soil properties
# Lamont - 0-19: silty loam, 20-34: clay, 35-60: clay loam
porosity     = np.zeros((nSoil,1))
satPotential = np.zeros((nSoil,1))
satHydrCond  = np.zeros((nSoil,1))
soilExponent = np.zeros((nSoil,1))
heatCapSoil  = np.zeros((nSoil,1))
profileTemp  = np.zeros((nSoil,1))
profileMois  = np.zeros((nSoil,1))

# if interpolating to more depths
#porosity[:,0]     =  [0.485,0.485, 0.485,0.485,0.485,0.485,0.482,0.482,0.482,0.476,0.476,0.476,0.476,0.476,0.476,0.476,0.476]
#satPotential[:,0] =  [0.786,0.786,0.786,0.786,0.786,0.786,0.405,0.405,0.405,0.630,0.630,0.630,0.630,0.630,0.630,0.630,0.630]
#satHydrCond[:,0]  =  [7.2,7.2,7.2,7.2,7.2,7.2,1.3,1.3,1.3,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5]
#soilExponent[:,0] =  [5.3,5.3,5.3,5.3,5.3,5.3,11.4,11.4,11.4,8.52,8.52,8.52,8.52,8.52,8.52,8.52,8.52]
#heatCapSoil[:,0]  =  [1.27,1.27,1.27,1.27,1.27,1.27,1.09,1.09,1.09,1.23,1.23,1.23,1.23,1.23,1.23,1.23,1.23]
#porosity[:,0]     =  [0.485,0.485, 0.485,0.485,0.485,0.485,0.482,0.482,0.482,0.476,0.476,0.476,0.476,0.476,0.476,0.476,0.476]
#satPotential[:,0] =  [0.786,0.786,0.786,0.786,0.786,0.786,0.405,0.405,0.405,0.630,0.630,0.630,0.630,0.630,0.630,0.630,0.630]
#satHydrCond[:,0]  =  [7.2,7.2,7.2,7.2,7.2,7.2,1.3,1.3,1.3,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5]
#soilExponent[:,0] =  [5.3,5.3,5.3,5.3,5.3,5.3,11.4,11.4,11.4,8.52,8.52,8.52,8.52,8.52,8.52,8.52,8.52]
#heatCapSoil[:,0]  =  [1.27,1.27,1.27,1.27,1.27,1.27,1.09,1.09,1.09,1.23,1.23,1.23,1.23,1.23,1.23,1.23,1.23]

#profileTemp[:,0]  = initSoilTemp
#profileMois[:,0]  = initSoilMois

porosity[:,0]     =  [0.485,0.485,0.485,0.482,0.476,0.476]
satPotential[:,0] =  [0.786,0.786,0.786,0.405,0.630,0.630]
satHydrCond[:,0]  =  [7.2,7.2,7.2,1.3,2.5,2.5]
soilExponent[:,0] =  [5.3,5.3,5.3,11.4,8.52,8.52]
heatCapSoil[:,0]  =  [1.27,1.27,1.27,1.09,1.23,1.23]

profileTemp[:,0]  = measSoilTemp
profileMois[:,0]  = measSoilMois


# proper units and sign
satPotential = -satPotential
satHydrCond  = satHydrCond * 1e-6
heatCapSoil  = heatCapSoil * 1e6

# soil properties array
soilOut = np.r_[porosity.T, satPotential.T, satHydrCond.T, soilExponent.T, heatCapSoil.T]

# write input files
#np.savetxt('soilLevels.ini',       initZ,       fmt='%.7e',  delimiter='\t', newline='\n')
np.savetxt('soilLevels.ini',       measZ,       fmt='%.7e',  delimiter='\t', newline='\n')
np.savetxt('soilTypeParams.ini',   soilOut,     fmt='%14.7e',delimiter='\t', newline='\n')
np.savetxt('soilTemperature.ini',  profileTemp, fmt='%.7e',  delimiter='\t', newline='\n')
np.savetxt('soilMoisture.ini',     profileMois, fmt='%.7e',  delimiter='\t', newline='\n')

################################################
# Read met tower data for u,v,T,q for 20121024 #
################################################

# open MET data
met = nc.MFDataset('sgpmetE13.b1.*.cdf')

# time indices (12UTC-12UTC)
tstart = 0
tfinal = 2880

# grab variables
tt = met.variables['time'][tstart:tfinal+1]
ws = met.variables['wspd_vec_mean'][tstart:tfinal+1]
wd = met.variables['wdir_vec_mean'][tstart:tfinal+1]
pt = met.variables['temp_mean'][tstart:tfinal+1] + 273.15
pa = met.variables['atmos_pressure'][tstart:tfinal+1]
pv = met.variables['vapor_pressure_mean'][tstart:tfinal+1]

# fix times to be continuous across days
tt[1440::] = tt[1440::]+tt[1439]+60

# compute wind components
uc = -ws*np.sin(wd*np.pi/180)
vc = -ws*np.cos(wd*np.pi/180)

# compute water vapor mixing ratio
qv = 0.622 * pv / (pa-pv)

# interpolate between times to increase input times
#dt     = 10
#ttNew  = np.arange(tt[0],tt[-1]+1,dt) 
#nTimes = len(ttNew)
#
#ucNew = np.interp(ttNew, tt, uc)
#vcNew = np.interp(ttNew, tt, vc)
#qvNew = np.interp(ttNew, tt, qv)
#ptNew = np.interp(ttNew, tt, pt)
#
ttNew = tt
nTimes = len(tt)
ucNew = uc
vcNew = vc
qvNew = qv
ptNew = pt

of = open('timeSeriesMET.dat','w')
for t in range(len(ttNew)):
	os = "{0:.6E}\t{1:13.6E}\t{2:2.6E}\t{3:2.6E}\t{4:2.6E}\t{5:2.6E}\t{6:2.6E}\n".format(ttNew[t], ucNew[t], vcNew[t], 0, ptNew[t], qvNew[t],0)
	of.write(os)
of.close()

########################################
# Read net radiation data for 20121024 #
########################################

# open CO2FLX data
data = nc.MFDataset('sgp30co2flx4mmetC1.b1.*')

# time indices
tstart = 0
tfinal = 96

# grab variables
tt   = data.variables['time_offset'][tstart:tfinal+1]
lwD  = data.variables['r_down_long_hemisp'][tstart:tfinal+1]
lwU  = data.variables['r_up_long_hemisp'][tstart:tfinal+1]
swD  = data.variables['r_down_short_hemisp'][tstart:tfinal+1]
swU  = data.variables['r_up_short_hemisp'][tstart:tfinal+1]
rNet = data.variables['r_net'][tstart:tfinal+1]
G    = data.variables['mean_g_soil'][tstart:tfinal+1]
H    = data.variables['h'][tstart:tfinal+1]
LE   = data.variables['le'][tstart:tfinal+1]

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

of = open('timeSeriesRAD.dat','w')
for t in range(len(tt1m)):
	os = "{0:.6E}\t{1:13.6E}\n".format(tt1m[t], rNet1m[t])
	of.write(os)
of.close()