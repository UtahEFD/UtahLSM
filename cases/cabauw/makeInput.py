#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
import sys
import pylab as pl

##########################################
# Manual entry of soil data for 20121024 #
##########################################

# soil temperature data
obs   = nc.Dataset('observations/cesar_soil_heat_lb1_t10_v1.0_201709.nc')
t_obs = obs.variables['time'][:]*3600.
st_00 = obs.variables['TS00'][:]+273.15
st_02 = obs.variables['TS02'][:]+273.15
st_04 = obs.variables['TS04'][:]+273.15
st_06 = obs.variables['TS06'][:]+273.15
st_08 = obs.variables['TS08'][:]+273.15
st_12 = obs.variables['TS12'][:]+273.15
st_20 = obs.variables['TS20'][:]+273.15
st_30 = obs.variables['TS30'][:]+273.15
st_50 = obs.variables['TS50'][:]+273.15
st_ob = np.array([st_00,st_02,st_04,st_04,st_08,st_12,st_20,st_30,st_50])
z_obs = np.array([0.00,0.02,0.04,0.06,0.08,0.12,0.20,0.30,0.50])
nsoil = len(z_obs)
ntime = len(t_obs)

# estimate surface moisture
obs   = nc.Dataset('observations/cesar_surface_meteo_lc1_t10_v1.0_201709.nc')
rh_00 = obs.variables['RH002'][:]/100
psi_s = -.405
eta_s = .482
b     = 11.4
g     = 9.81
Rv    = 416.5
sm_00 = eta_s * ( (g*psi_s)/(Rv*st_00*np.log(rh_00)) )**(1/b)

# soil moisture data
obs   = nc.Dataset('observations/cesar_soil_water_lb1_t10_v1.1_201709.nc')
sm_05 = obs.variables['TH05'][:]
sm_19 = obs.variables['TH19'][:]
sm_33 = obs.variables['TH33'][:]
sm_40 = obs.variables['TH40'][:]
sm_56 = obs.variables['TH56'][:]
sm_ob = np.array([sm_00,sm_05,sm_19,sm_33,sm_40,sm_56])
z_obm = np.array([0.0,0.05,0.19,0.33,0.40,0.56])


# interpolate soil moisture to temperature grid
sm_oi = np.zeros((nsoil,ntime))
for t in range(ntime):
	sm_oi[:,t] = np.interp(z_obs,z_obm,sm_ob[:,t])

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
stype = np.full((nsoil),11)

# write output
of = open('inputSoil.dat','w')
os = '{0:^15} {1:^15s} {2:^15s} {3:^15s}\n'
os = os.format('soil_z','soil_type','soil_T','soil_q')
of.write(os)
for z in range(nsoil):
	os = "{0:15.8E}  {1:2.8E}  {2:2.8E}  {3:2.8E}\n"
	os = os.format(z_obs[z],stype[z],st_ob[z,0],sm_oi[z,0])
	of.write(os)
of.close()

print(ntime)

###################################
# Read met tower data for u,v,T,q #
###################################

# open MET data
met = nc.MFDataset('observations/cesar_surface_meteo_lc1_t10_v1.0_201709.nc')

# grab variables
ws = met.variables['F010'][:]
wd = met.variables['D010'][:]
pt = met.variables['TA002'][:]
pa = met.variables['P0'][:]
qv = met.variables['Q002'][:]/1000

# compute wind components
uc = -ws*np.sin(wd*np.pi/180)
vc = -ws*np.cos(wd*np.pi/180)

#################################
# Read radiation data for R_net #
#################################
rad = nc.MFDataset('observations/cesar_surface_radiation_lc1_t10_v1.0_201709.nc')
swu = rad.variables['SWU'][:]
swd = rad.variables['SWD'][:]
lwu = rad.variables['LWU'][:]
lwd = rad.variables['LWD'][:]
net = swd-swu+lwd-lwu

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

##############################
# Write all time series data #
##############################
of = open('inputMetr.dat','w')
os = '{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s}\n'.format('atm_u','atm_v','atm_T','atm_q', 'R_net')
of.write(os)
for t in range(ntime):
	os = "{0:15.8E}  {1:2.8E}  {2:2.8E}  {3:2.8E}  {4:2.8E}\n".format(uc[t], vc[t], pt[t], qv[t], net[t])
	of.write(os)
of.close()
