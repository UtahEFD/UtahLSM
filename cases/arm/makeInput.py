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

# compute wind components
uc = -ws * np.sin(wd * np.pi / 180)
vc = -ws * np.cos(wd * np.pi / 180)

# compute water vapor mixing ratio
qs = 0.622 * pv / (pa - pv)

# time dimension
dt    = tm[1] - tm[0]
ntime = len(tm)
t_utc = (tm - 86400) % 86400

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
dt2   = tm2[1] - tm2[0]
tm2   = np.arange(0,nt * dt2,dt2)

# interpolate from 30-minute to 1-minute frequency to match MET
rNet1m = np.interp(tm,tm2,rNet)

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
namelist['length'] = {}
namelist['soil'] = {}
namelist['radiation'] = {}
namelist['output'] = {}

# time section
namelist['time']['step_seb'] = 1
namelist['time']['step_dif'] = 1

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

with open('lsm_namelist.json', 'w') as outfile:
    json.dump(namelist,outfile,indent=4)
