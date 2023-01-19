#!/usr/bin/env python
import json
import numpy as np
import netCDF4 as nc
import time 

##########################################
# Soil data taken from Cesar Observatory #
##########################################

# soil temperature data
obs   = nc.Dataset('observations/cesar_soil_heat_lb1_t10_v1.0_201709.nc')
t_obs = obs.variables['time'][:] * 3600.
st_00 = obs.variables['TS00'][0] + 273.15
st_02 = obs.variables['TS02'][0] + 273.15
st_04 = obs.variables['TS04'][0] + 273.15
st_06 = obs.variables['TS06'][0] + 273.15
st_08 = obs.variables['TS08'][0] + 273.15
st_12 = obs.variables['TS12'][0] + 273.15
st_20 = obs.variables['TS20'][0] + 273.15
st_30 = obs.variables['TS30'][0] + 273.15
st_50 = obs.variables['TS50'][0] + 273.15
st_ob = np.array([st_00,st_02,st_04,st_04,st_08,st_12,st_20,st_30,st_50])
z_obs = np.array([0.00,0.02,0.04,0.06,0.08,0.12,0.20,0.30,0.50])
nsoil = len(z_obs)
ntime = len(t_obs)

# estimate surface moisture
obs   = nc.Dataset('observations/cesar_surface_meteo_lc1_t10_v1.0_201709.nc')
rh_00 = obs.variables['RH002'][0] / 100
psi_s = -.405
eta_s = .482
b     = 11.4
g     = 9.81
Rv    = 416.5
sm_00 = eta_s * ((g * psi_s) / (Rv * st_00 * np.log(rh_00)))**(1 / b)

# soil moisture data
obs   = nc.Dataset('observations/cesar_soil_water_lb1_t10_v1.1_201709.nc')
sm_05 = obs.variables['TH05'][0]
sm_19 = obs.variables['TH19'][0]
sm_33 = obs.variables['TH33'][0]
sm_40 = obs.variables['TH40'][0]
sm_56 = obs.variables['TH56'][0]
sm_ob = np.array([sm_00,sm_05,sm_19,sm_33,sm_40,sm_56])
z_obm = np.array([0.0,0.05,0.19,0.33,0.40,0.56])

# interpolate soil moisture to temperature grid
sm_oi    = np.zeros(nsoil)
sm_oi[:] = np.interp(z_obs,z_obm,sm_ob[:])

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

###################################
# Read met tower data for u,v,T,q #
###################################

# open MET data
met = nc.MFDataset('observations/cesar_surface_meteo_lc1_t10_v1.0_201709.nc')

# grab variables
tm = met.variables['time'][:] * 3600.
ws = met.variables['F010'][:]
wd = met.variables['D010'][:]
pt = met.variables['TA002'][:]
pa = met.variables['P0'][:]
qs = met.variables['Q002'][:] / 1000

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
rad = nc.MFDataset('observations/cesar_surface_radiation_lc1_t10_v1.0_201709.nc')
swu = rad.variables['SWU'][:]
swd = rad.variables['SWD'][:]
lwu = rad.variables['LWU'][:]
lwd = rad.variables['LWD'][:]
net = swd - swu + lwd - lwu

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

# time section
namelist['time']['step_seb'] = 1
namelist['time']['step_dif'] = 1

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
namelist['soil']['param']     = 1
namelist['soil']['model']     = 2

# radiation section
namelist['radiation']['utc_start']  = t_utc[0]
namelist['radiation']['comp_rad']   = 0
namelist['radiation']['albedo']     = 0.33
namelist['radiation']['emissivity'] = 0.99
namelist['radiation']['latitude']   = 51.9711
namelist['radiation']['longitude']  = -4.9267
namelist['radiation']['julian_day'] = 244

# output section
namelist['output']['save'] = 1
namelist['output']['fields'] = ['all']

with open('lsm_namelist.json', 'w') as outfile:
    json.dump(namelist,outfile,indent=4)
