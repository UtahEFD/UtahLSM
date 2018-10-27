#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
import pylab as pl

##########################################
# Manual entry of soil data for 20121024 #
##########################################

# soil sensor depths and interpolation levels
measZ = [0.0, 0.025, 0.05, 0.15, 0.25, 0.35, 0.60]
nSoil = len(measZ)

# soil temperature
measSoilTemp = [299.52, 295.68, 295.52, 294.75, 293.69, 293.35, 292.66]

# soil moisture
measSoilMois = [0.2448, 0.2448, 0.2535, 0.2528, 0.3666, 0.3272, 0.3230]

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
rh = met.variables['rh_mean'][:]

# compute wind components
uc = -ws*np.sin(wd*np.pi/180)
vc = -ws*np.cos(wd*np.pi/180)
print(len(uc))
# compute water vapor mixing ratio
qv = 0.622 * pv / (pa-pv) * 1.35
es = 6.1078 * np.exp(17.269*(pt-273.15)/(pt-35.86))
qs = 0.622 * (es / (pa-0.378*es))

qv2 = rh/100 * qs

sl = np.arange(0,0.031,0.001)
tt1 = tt
########################################
# Read net radiation data for 20121024 #
########################################

# open CO2FLX data
data = nc.MFDataset('observations/sgp30co2flx4mmetC1.b1.*')

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
qv3  = data.variables['mean_q'][:]*18.0153/1000/1.2/1000


# construct time
nt   = len(tt)
tt   = np.arange(0,nt*1800,1800)
tt1m = np.arange(0,31*86400,60)

# interpolate from 30-minute to 1-minute frequency to match MET
lwD1m  = np.interp(tt1m,tt,lwD)
lwU1m  = np.interp(tt1m,tt,lwU)
swD1m  = np.interp(tt1m,tt,swD)
swU1m  = np.interp(tt1m,tt,swU)
rNet1m = np.interp(tt1m,tt,rNet)
G1m    = np.interp(tt1m,tt,G)
H1m    = np.interp(tt1m,tt,H)
LE1m   = np.interp(tt1m,tt,LE)

tc = 0
badU = []
for r in ws:
	if (r==0):
		badU.append(tc)
	tc+=1

for b in badU:
	uc[b] = uc[badU[0]-1]
	vc[b] = vc[badU[0]-1]

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
	
tc = 0
badQ = []
for r in qv:
	if (str(r)=='--'):
		badQ.append(tc)
	tc+=1

for b in badQ:
	qv3[b] = qv3[badQ[0]-1]


nt   = len(tt1)
tt1  = np.arange(0,nt*60,60)
nt   = len(tt)
tt2  = np.arange(0,nt*1800,1800)
#pl.plot(tt1/3600/24,ws)
#pl.plot(tt2/3600/24,qv3)
#pl.plot(sl,sl,color='k')
#pl.show()

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
