#!/usr/bin/env python
# 
# UtahLSM
# 
# Copyright (c) 2017–2023 Jeremy A. Gibbs
# Copyright (c) 2017–2023 Rob Stoll
# Copyright (c) 2017–2023 Eric Pardyjak
# Copyright (c) 2017–2023 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

import json
import numpy as np
import netCDF4 as nc
import time 

#############
# Constants #
#############
Lv = 2.5E6
Rd = 287
Rv = 461
T0 = 273.15
E0 = 6.11
dt = 300.0

######################################
# Data from Materhorn Field Campaign #
######################################

# materhorn data files
obs = open('observations/efs_sage_iop5_5min.dat','r').read().split()[20::]
lem = open('observations/lems_iop5_10sec.dat','r').read().split()[17::]

# atmospheric data
pa    = np.array(obs[5::20]).astype(float) * 100.0
ws    = np.array(obs[6::20]).astype(float)
wd    = np.array(obs[7::20]).astype(float)
pt    = np.array(obs[8::20]).astype(float) + T0
rh    = np.array(obs[9::20]).astype(float)
net   = np.array(obs[10::20]).astype(float)
ntime = len(ws)

# convert relative humidity to mixing ratio
es = E0 * np.exp(Lv / Rv * (1 / T0 - 1 / pt))
e  = rh * es / 100
sh = e * Rd / (Rv * (pa / 100 - e))
qs = sh / (sh + 1)

# soil moisture
sm_02 = float(lem[15])
sm_07 = float(obs[11])
sm_25 = float(lem[13])
sm_00 = sm_02
sm_ob = np.array([sm_00,sm_02,sm_07,sm_25])
z_obm = np.array([0.000,0.025,0.075,0.25])

# soil temperature
st_01 = float(obs[12]) + T0
st_02 = float(obs[13]) + T0
st_05 = float(obs[14]) + T0
st_07 = float(obs[15]) + T0
st_10 = float(obs[16]) + T0
st_15 = float(obs[17]) + T0
st_25 = float(obs[18]) + T0
st_70 = float(obs[19]) + T0
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

# initialization file
init             = nc.Dataset('lsm_init.nc','w')
init.description = "UtahLSM input file"
init.source      = "Jeremy A. Gibbs"
init.history     = "Created " + time.ctime(time.time())

# add dimensions
init.createDimension('z', nsoil)

# add variables
init_z = init.createVariable("soil_z", "f4", ("z",))
init_z.long_name = "z-distance"
init_z.units = "m"
init_T = init.createVariable("soil_T", "f4", ("z"))
init_T.long_name = "soil temperature"
init_T.units = "T"
init_q = init.createVariable("soil_q", "f4", ("z"))
init_q.long_name = "soil moisture"
init_q.units = "g g-1"
init_i = init.createVariable("soil_type", "f4", ("z"))
init_i.long_name = "soil type"
init_i.units = ""

# write initial data
init_z[:] = z_obs
init_T[:] = st_ob
init_q[:] = sm_oi
init_i[:] = stype

# close file
init.close()

##############################
# Write all time-series data #
##############################

# time-series file
metr             = nc.Dataset('lsm_offline.nc','w')
metr.description = "UtahLSM input file for offline run"
metr.source      = "Jeremy A. Gibbs"
metr.history     = "Created " + time.ctime(time.time())

# add dimensions
metr.createDimension('t', ntime)

# add variables
metr_s = metr.createVariable("tstep", "f4", ())
metr_s.long_name = "time step for input offline data"
metr_s.units = "s"
metr_u = metr.createVariable("atm_U", "f4", ("t"))
metr_u.long_name = "wind speed"
metr_u.units = "m s-1"
metr_t = metr.createVariable("atm_T", "f4", ("t"))
metr_t.long_name = "temperature"
metr_t.units = "K"
metr_q = metr.createVariable("atm_q", "f4", ("t"))
metr_q.long_name = "mixing ratio"
metr_q.units = "g g-1"
metr_p = metr.createVariable("atm_p", "f4", ("t"))
metr_p.long_name = "pressure"
metr_p.units = "hPa"
metr_r = metr.createVariable("R_net", "f4", ("t"))
metr_r.long_name = "net radiation"
metr_r.units = "W m-2"

# write time-series data
metr_s[:] = float(dt)
metr_u[:] = ws
metr_t[:] = pt
metr_q[:] = qs
metr_p[:] = pa
metr_r[:] = net

# close file
metr.close()

########################
# Settings for UtahLSM #
########################
namelist = {}
namelist['time'] = {}
namelist['grid'] = {}
namelist['surface'] = {}
namelist['soil'] = {}
namelist['radiation'] = {}
namelist['output'] = {}

# time section
namelist['time']['step_seb']   = 1
namelist['time']['step_dif']   = 1
namelist['time']['utc_start']  = 12.0
namelist['time']['julian_day'] = 133

# grid section
namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1

# surface scale section
namelist['surface']['z_o']        = 0.15
namelist['surface']['z_t']        = 0.0015
namelist['surface']['z_m']        = 10.0
namelist['surface']['z_s']        = 2.0
namelist['surface']['albedo']     = 0.33
namelist['surface']['emissivity'] = 0.99

# soil section
namelist['soil']['nsoil'] = nsoil
namelist['soil']['param'] = 3
namelist['soil']['model'] = 2

# radiation section
namelist['radiation']['comp_rad']  = 0
namelist['radiation']['latitude']  = 40.121360
namelist['radiation']['longitude'] = -113.129070

# output section
namelist['output']['save']   = 1
namelist['output']['fields'] = ['all']

with open('lsm_namelist.json', 'w') as outfile:
    json.dump(namelist,outfile,indent=4)
