#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
import os,sys,json
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
z_obs = [0.0, 0.05, 0.10, 0.20, 0.50, 0.75, 1.0]
nsoil = len(z_obs)

# soil temperature
st_ob = [301.95,301.95,300.45,297.52,293.52,291.90,291.95]

# soil moisture
sm_ob = [0.1115,0.1115,0.1707,0.3603,0.3276,0.1549,.2089]

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
stype     = [4,4,4,4,11,8,8,8]

################################################
# Read met tower data for u,v,T,q for 20121024 #
################################################

# open MET data
met = nc.MFDataset('observations/sgpmetE13.b1.*.cdf')

# grab variables
tm = met.variables['time'][:]
ws = met.variables['wspd_vec_mean'][:]
wd = met.variables['wdir_vec_mean'][:]
pt = met.variables['temp_mean'][:] + 273.15
pa = met.variables['atmos_pressure'][:]*10
pv = met.variables['vapor_pressure_mean'][:]*10

# compute wind components
uc = -ws*np.sin(wd*np.pi/180)
vc = -ws*np.cos(wd*np.pi/180)

# compute water vapor mixing ratio
qs = 0.622 * pv / (pa-pv)

# time dimension
dt    = tm[1] - tm[0]
ntime = len(tm)
t_utc = (tm-86400)%86400

########################################
# Read net radiation data for 20121024 #
########################################

# open CO2FLX data
data = nc.MFDataset('observations/sgpseb*')

# grab variables
tm2  = data.variables['time_offset'][:]
rNet = data.variables['net_radiation'][:]

# construct time
nt    = len(tm2)
dt2   = tm2[1]-tm2[0]
tm2   = np.arange(0,nt*dt2,dt2)

# interpolate from 30-minute to 1-minute frequency to match MET
rNet1m = np.interp(tm,tm2,rNet)

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
for r in qs:
	if (str(r)=='--'):
		badQ.append(tc)
	tc+=1

for b in badQ:
	qs[b] = qs[badQ[0]-1]

##############################
# Write all time series data #
##############################
metr = {}
metr['time'] = {}
metr['data'] = {}

metr['time']['ntime'] = ntime
metr['time']['tstep'] = float(dt)
metr['data']['atm_U'] = ws.tolist()
metr['data']['atm_T'] = pt.tolist()
metr['data']['atm_q'] = qs.tolist()
metr['data']['atm_p'] = pa.tolist()
metr['data']['R_net'] = rNet1m.tolist()
with open('inputOffline.json', 'w') as outfile:  
    json.dump(metr,outfile,indent=4)

########################
# Settings for UtahLSM #
########################
namelist = {}
namelist['grid'] = {}
namelist['length'] = {}
namelist['soil'] = {}
namelist['radiation'] = {}
namelist['output'] = {}

# grid section
namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1

# length scale section
namelist['length']['z_o'] = 0.0500
namelist['length']['z_t'] = 0.0005
namelist['length']['z_m'] = 10.0 
namelist['length']['z_s'] = 2.0

# soil section
namelist['soil']['nsoil']     = nsoil
namelist['soil']['param']     = 3
namelist['soil']['model']     = 2
namelist['soil']['soil_z']    = z_obs
namelist['soil']['soil_type'] = stype
namelist['soil']['soil_T']    = st_ob
namelist['soil']['soil_q']    = sm_ob

# radiation section
namelist['radiation']['utc_start']  = t_utc[0] 
namelist['radiation']['comp_rad']   = 0
namelist['radiation']['albedo']     = 0.25
namelist['radiation']['emissivity'] = 0.96
namelist['radiation']['latitude']   = 36.6906
namelist['radiation']['longitude']  = 97.5564
namelist['radiation']['julian_day'] = 153

# output section
namelist['output']['save'] = 1
namelist['output']['fields'] = ['all']

with open('inputLSM.json', 'w') as outfile:  
    json.dump(namelist,outfile,indent=4)