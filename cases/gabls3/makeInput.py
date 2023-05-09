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
# Soil data taken from Cesar Observatory #
##########################################

# measurements are saved 10 minutes
# so there are 144 per day
date = 2
tidx = (date - 1) * 144

# soil temperature data
obs   = nc.Dataset('observations/cesar_soil_heat_lb1_t10_v1.0_200607.nc')
t_obs = obs.variables['time'][tidx] * 3600.
st_00 = obs.variables['TS00'][tidx] + 273.15
st_02 = obs.variables['TS02'][tidx] + 273.15
st_04 = obs.variables['TS04'][tidx] + 273.15
st_06 = obs.variables['TS06'][tidx] + 273.15
st_08 = obs.variables['TS08'][tidx] + 273.15
st_12 = obs.variables['TS12'][tidx] + 273.15
st_20 = obs.variables['TS20'][tidx] + 273.15
st_30 = obs.variables['TS30'][tidx] + 273.15
st_50 = obs.variables['TS50'][tidx] + 273.15
st_ob = np.array([st_00,st_02,st_04,st_04,st_08,st_12,st_20,st_30,st_50])
z_obs = np.array([0.00,0.02,0.04,0.06,0.08,0.12,0.20,0.30,0.50])
nsoil = len(z_obs)

# soil moisture data (10-minute)
obs   = nc.Dataset('observations/cesar_soil_water_lb1_t10_v1.1_200607.nc')
sm_03 = obs.variables['TH03'][tidx]
sm_08 = obs.variables['TH08'][tidx]
sm_20 = obs.variables['TH20'][tidx]

# soil moisture data (daily average from prior day)
# Julian day for July 2, 2006 is 183, so day prior is 182
jdi   = 181
obs   = nc.Dataset('observations/cesar_tdr_soilmoisture_la1_t1d_v1.0_2006.nc')
sm_05 = np.mean([obs.variables['SM1'][jdi],obs.variables['SM7'][jdi], obs.variables['SM13'][jdi],obs.variables['SM19'][jdi]])
sm_15 = np.mean([obs.variables['SM2'][jdi],obs.variables['SM8'][jdi], obs.variables['SM14'][jdi],obs.variables['SM20'][jdi]])
sm_30 = np.mean([obs.variables['SM3'][jdi],obs.variables['SM9'][jdi], obs.variables['SM15'][jdi],obs.variables['SM21'][jdi]])
sm_45 = np.mean([obs.variables['SM4'][jdi],obs.variables['SM10'][jdi],obs.variables['SM16'][jdi],obs.variables['SM22'][jdi]])
sm_60 = np.mean([obs.variables['SM5'][jdi],obs.variables['SM11'][jdi],obs.variables['SM17'][jdi],obs.variables['SM23'][jdi]])
sm_73 = np.mean([obs.variables['SM6'][jdi],obs.variables['SM12'][jdi],obs.variables['SM18'][jdi],obs.variables['SM24'][jdi]])

sm_ob = np.array([sm_05,sm_05,sm_15,sm_30,sm_45,sm_60,sm_73])
z_obm = np.array([0.0,0.05,0.15,0.30,0.45,0.60,0.725])

# interpolate soil moisture and temperature to regular grid
z_int = np.arange(0,0.51,0.05)
nsoil = len(z_int)
st_oi = np.interp(z_int,z_obs,st_ob)
sm_oi = np.interp(z_int,z_obm,sm_ob)

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
# 13 = B11
# 14 = O12
# 15 = O16
stype = np.full((nsoil),11)
# stype = np.array([11,11,11,11,11,11,11,11,11,11,12])
# stype = np.array([13,13,13,13,13,13,14,14,15])

# initialization file
init             = nc.Dataset('lsm_init.nc','w')
init.description = "UtahLSM input file"
init.source      = "Jeremy A. Gibbs"
init.history     = "Created " + time.ctime(time.time())

# add dimensions
init.createDimension('z', nsoil)

# add variables
init_z = init.createVariable("soil_z", "f8", ("z",))
init_z.long_name = "z-distance"
init_z.units = "m"
init_T = init.createVariable("soil_T", "f8", ("z"))
init_T.long_name = "soil temperature"
init_T.units = "T"
init_q = init.createVariable("soil_q", "f8", ("z"))
init_q.long_name = "soil moisture"
init_q.units = "g g-1"
init_i = init.createVariable("soil_type", "f8", ("z"))
init_i.long_name = "soil type"
init_i.units = ""

# write initial data
init_z[:] = z_int
init_T[:] = st_oi
init_q[:] = sm_oi
init_i[:] = stype

# close file
init.close()

###################################
# Read met tower data for u,v,T,q #
###################################

# end time is 9 hours, so 6*9 = 54
tend  = tidx + 55

# open MET data
met = nc.MFDataset('observations/cesar_surface_meteo_lc1_t10_v1.0_200607.nc')

# grab variables
tm = met.variables['time'][tidx:tend] * 3600.
ws = met.variables['F010'][tidx:tend]
wd = met.variables['D010'][tidx:tend]
pt = met.variables['TA002'][tidx:tend]
pa = met.variables['P0'][tidx:tend]
qs = met.variables['Q002'][tidx:tend] / 1000

# compute wind components
uc = -ws * np.sin(wd * np.pi / 180)
vc = -ws * np.cos(wd * np.pi / 180)

# time dimension
dt    = tm[1] - tm[0]
ntime = len(tm)
t_utc = (tm - 86400) % 86400

#################################
# Read radiation data for R_net #
#################################
rad = nc.MFDataset('observations/cesar_surface_radiation_lc1_t10_v1.0_200607.nc')
swu = rad.variables['SWU'][tidx:tend]
swd = rad.variables['SWD'][tidx:tend]
lwu = rad.variables['LWU'][tidx:tend]
lwd = rad.variables['LWD'][tidx:tend]
net = swd - swu + lwd - lwu

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
metr_s = metr.createVariable("tstep", "f8", ())
metr_s.long_name = "time step for input offline data"
metr_s.units = "s"
metr_u = metr.createVariable("atm_U", "f8", ("t"))
metr_u.long_name = "wind speed"
metr_u.units = "m s-1"
metr_t = metr.createVariable("atm_T", "f8", ("t"))
metr_t.long_name = "temperature"
metr_t.units = "K"
metr_q = metr.createVariable("atm_q", "f8", ("t"))
metr_q.long_name = "mixing ratio"
metr_q.units = "g g-1"
metr_p = metr.createVariable("atm_p", "f8", ("t"))
metr_p.long_name = "pressure"
metr_p.units = "hPa"
metr_r = metr.createVariable("R_net", "f8", ("t"))
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
namelist['time']['utc_start']  = t_utc[0]
namelist['time']['julian_day'] = 183

# grid section
namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1

# surface section
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
namelist['radiation']['latitude']  = 51.9711
namelist['radiation']['longitude'] = -4.9267

# output section
namelist['output']['save']   = 1
namelist['output']['fields'] = ['all']

with open('lsm_namelist.json', 'w') as outfile:
    json.dump(namelist, outfile, indent=4)
