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

##########################################
# Manual entry of soil data for 20121024 #
##########################################

# soil sensor depths and interpolation levels
z_obs = [0.0, 0.025, 0.05, 0.15, 0.25, 0.35, 0.60]
nsoil = len(z_obs)

# soil temperature
st_ob = [299.52, 295.68, 295.52, 294.75, 293.69, 293.35, 292.66]

# soil moisture
sm_ob = [0.2448, 0.2448, 0.2535, 0.2528, 0.3666, 0.3272, 0.3230]

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
stype     = [4,4,4,4,11,8,8]

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
init_q[:] = sm_ob
init_i[:] = stype

# close file
init.close()

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
pa = met.variables['atmos_pressure'][:] * 10
pv = met.variables['vapor_pressure_mean'][:] * 10
rh = met.variables['rh_mean'][:]

# compute wind components
uc = -ws * np.sin(wd * np.pi / 180)
vc = -ws * np.cos(wd * np.pi / 180)

# compute water vapor mixing ratio
qs = 0.622 * pv / (pa - pv) * 1.35

# time dimension
dt    = tm[1] - tm[0]
ntime = len(tm)
t_utc = (tm - 86400) % 86400

########################################
# Read net radiation data for 20121024 #
########################################

# open CO2FLX data
data = nc.MFDataset('observations/sgp30co2flx4mmetC1.b1.*')

# grab variables
tm2  = data.variables['time_offset'][:]
lwD  = data.variables['r_down_long_hemisp'][:]
lwU  = data.variables['r_up_long_hemisp'][:]
swD  = data.variables['r_down_short_hemisp'][:]
swU  = data.variables['r_up_short_hemisp'][:]
rNet = data.variables['r_net'][:]
G    = data.variables['mean_g_soil'][:]
H    = data.variables['h'][:]
LE   = data.variables['le'][:]

# construct time
nt    = len(tm2)
dt2   = tm2[1] - tm2[0]
tm2    = np.arange(0,nt * dt2,dt2)

# interpolate from 30-minute to 1-minute frequency to match MET
lwD1m  = np.interp(tm,tm2,lwD)
lwU1m  = np.interp(tm,tm2,lwU)
swD1m  = np.interp(tm,tm2,swD)
swU1m  = np.interp(tm,tm2,swU)
rNet1m = np.interp(tm,tm2,rNet)
G1m    = np.interp(tm,tm2,G)
H1m    = np.interp(tm,tm2,H)
LE1m   = np.interp(tm,tm2,LE)

tc = 0
badU = []
for r in ws:
    if (r == 0):
        badU.append(tc)
    tc += 1

for b in badU:
    uc[b] = uc[badU[0] - 1]
    vc[b] = vc[badU[0] - 1]

tc = 0
badT = []
for r in pt:
    if (str(r) == '--'):
        badT.append(tc)
    tc += 1
for b in badT:
    pt[b] = pt[badT[0] - 1]

tc = 0
badQ = []
for r in qs:
    if (str(r) == '--'):
        badQ.append(tc)
    tc += 1

for b in badQ:
    qs[b] = qs[badQ[0] - 1]

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
metr_r[:] = rNet1m

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
namelist['time']['utc_start']  = t_utc[0]
namelist['time']['julian_day'] = 298

# grid section
namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1

# surface section
namelist['surface']['z_o']        = 0.0500
namelist['surface']['z_t']        = 0.0005
namelist['surface']['z_m']        = 10.0
namelist['surface']['z_s']        = 2.0
namelist['surface']['albedo']     = 0.25
namelist['surface']['emissivity'] = 0.96

# soil section
namelist['soil']['nsoil'] = nsoil
namelist['soil']['param'] = 3
namelist['soil']['model'] = 2

# radiation section
namelist['radiation']['comp_rad']  = 0
namelist['radiation']['latitude']  = 36.6906
namelist['radiation']['longitude'] = 97.5564

# output section
namelist['output']['save']   = 1
namelist['output']['fields'] = ['all']

with open('lsm_namelist.json', 'w') as outfile:
    json.dump(namelist,outfile,indent=4)