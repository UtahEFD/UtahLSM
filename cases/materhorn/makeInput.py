#!/usr/bin/env python
import sys, json
import netCDF4 as nc
import numpy as np
import pylab as pl
import datetime as dt

#############
# Constants #
#############
Lv = 2.5E6
Rd = 287
Rv = 461
T0 = 273.15
E0 = 6.11

######################################
# Data from Materhorn Field Campaign #
######################################

# materhorn data files
obs = open('observations/efs_sage_iop5_5min.dat','r').read().split()[20::]
lem = open('observations/lems_iop5_10sec.dat','r').read().split()[17::]

# atmospheric data
pa    = np.array(obs[5::20]).astype(float)*100.0
ws    = np.array(obs[6::20]).astype(float)
wd    = np.array(obs[7::20]).astype(float)
pt    = np.array(obs[8::20]).astype(float)+T0
rh    = np.array(obs[9::20]).astype(float)
net   = np.array(obs[10::20]).astype(float)
ntime = len(ws)

# convert relative humidity to mixing ratio
es = E0*np.exp(Lv/Rv * (1/T0 - 1/pt))
e  = rh*es/100
sh = e * Rd / (Rv * (pa/100 - e))
qs = sh / (sh + 1)

# soil moisture
sm_02 = np.float(lem[15])
sm_07 = np.float(obs[11])
sm_25 = np.float(lem[13])
sm_00 = sm_02
sm_ob = np.array([sm_00,sm_02,sm_07,sm_25])
z_obm = np.array([0.000,0.025,0.075,0.25])

# soil temperature
st_01 = np.float(obs[12])+T0
st_02 = np.float(obs[13])+T0
st_05 = np.float(obs[14])+T0
st_07 = np.float(obs[15])+T0
st_10 = np.float(obs[16])+T0
st_15 = np.float(obs[17])+T0
st_25 = np.float(obs[18])+T0
st_70 = np.float(obs[19])+T0
st_00 = st_01
st_ob = np.array([st_00,st_01,st_02,st_05,st_07,st_10,st_15,st_25,st_70])
z_obs = np.array([0.000,0.01,0.025,0.050,0.075,0.100,0.150,0.250,0.700])
nsoil = len(z_obs)

# interpolate soil moisture to temperature grid
sm_oi = np.interp(z_obs,z_obm,sm_ob)

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
stype = np.full((nsoil),4)

##############################
# Write all time series data #
##############################
metr = {}
metr['time'] = {}
metr['data'] = {}

metr['time']['ntime'] = ntime 
metr['time']['tstep'] = 300.0
metr['data']['atm_U'] = ws.tolist()
metr['data']['atm_T'] = pt.tolist()
metr['data']['atm_q'] = qs.tolist()
metr['data']['atm_p'] = pa.tolist()
metr['data']['R_net'] = net.tolist()
with open('inputOffline.json', 'w') as outfile:  
    json.dump(metr,outfile,indent=4)

########################
# Settings for UtahLSM #
########################
namelist = {}
namelist['time'] = {}
namelist['grid'] = {}
namelist['length'] = {}
namelist['soil'] = {}
namelist['radiation'] = {}
namelist['output'] = {}

# grid section
namelist['time']['step_seb'] = 10
namelist['time']['step_dif'] = 10

# grid section
namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1

# length scale section
namelist['length']['z_o'] = 0.15
namelist['length']['z_t'] = 0.0015
namelist['length']['z_m'] = 10.0 
namelist['length']['z_s'] = 2.0

# soil section
namelist['soil']['nsoil']     = nsoil
namelist['soil']['param']     = 3
namelist['soil']['model']     = 2
namelist['soil']['soil_z']    = z_obs.tolist()
namelist['soil']['soil_type'] = stype.tolist()
namelist['soil']['soil_T']    = st_ob.tolist()
namelist['soil']['soil_q']    = sm_oi.tolist()

# radiation section
namelist['radiation']['utc_start']  = 12.0 
namelist['radiation']['comp_rad']   = 0
namelist['radiation']['albedo']     = 0.33
namelist['radiation']['emissivity'] = 0.99
namelist['radiation']['latitude']   = 40.121360
namelist['radiation']['longitude']  = -113.129070
namelist['radiation']['julian_day'] = 133

# output section
namelist['output']['save'] = 1
namelist['output']['fields'] = ['all']

with open('inputLSM.json', 'w') as outfile:  
    json.dump(namelist,outfile,indent=4)