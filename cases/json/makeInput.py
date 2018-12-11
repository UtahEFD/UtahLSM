#!/usr/bin/env python
import sys, json
import netCDF4 as nc
import numpy as np

##########################################
# Soil data taken from Cesar Observatory #
##########################################

# measurements are saved 10 minutes
# so there are 144 per day
date = 2
tidx = (date-1)*144

# soil temperature data
obs   = nc.Dataset('observations/cesar_soil_heat_lb1_t10_v1.0_200607.nc')
t_obs = obs.variables['time'][tidx]*3600.
st_00 = obs.variables['TS00'][tidx]+273.15
st_02 = obs.variables['TS02'][tidx]+273.15
st_04 = obs.variables['TS04'][tidx]+273.15
st_06 = obs.variables['TS06'][tidx]+273.15
st_08 = obs.variables['TS08'][tidx]+273.15
st_12 = obs.variables['TS12'][tidx]+273.15
st_20 = obs.variables['TS20'][tidx]+273.15
st_30 = obs.variables['TS30'][tidx]+273.15
st_50 = obs.variables['TS50'][tidx]+273.15
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
stype = np.full((nsoil),11)

###################################
# Read met tower data for u,v,T,q #
###################################

# end time is 9 hours, so 6*9 = 54
tend  = tidx+55

# open MET data
met = nc.MFDataset('observations/cesar_surface_meteo_lc1_t10_v1.0_200607.nc')

# grab variables
tm = met.variables['time'][tidx:tend]*3600.
ws = met.variables['F010'][tidx:tend]
wd = met.variables['D010'][tidx:tend]
pt = met.variables['TA002'][tidx:tend]
pa = met.variables['P0'][tidx:tend]
qs = met.variables['Q002'][tidx:tend]/1000

# compute wind components
uc = -ws*np.sin(wd*np.pi/180)
vc = -ws*np.cos(wd*np.pi/180)

# time dimension
dt    = tm[1] - tm[0]
ntime = len(tm)
t_utc = (tm-86400)%86400

#################################
# Read radiation data for R_net #
#################################
rad = nc.MFDataset('observations/cesar_surface_radiation_lc1_t10_v1.0_200607.nc')
swu = rad.variables['SWU'][tidx:tend]
swd = rad.variables['SWD'][tidx:tend]
lwu = rad.variables['LWU'][tidx:tend]
lwd = rad.variables['LWD'][tidx:tend]
net = swd-swu+lwd-lwu

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
metr['data']['R_net'] = pt.tolist()
with open('inputOffline.json', 'w') as outfile:  
    json.dump(metr,outfile,indent=4)

########################
# Settings for UtahLSM #
########################
namelist = {}
namelist['length'] = {}
namelist['soil'] = {}
namelist['radiation'] = {}

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
namelist['radiation']['utc_start']  = t_utc[0] 
namelist['radiation']['comp_rad']   = 0
namelist['radiation']['albedo']     = 0.33
namelist['radiation']['emissivity'] = 0.99
namelist['radiation']['latitude']   = 51.9711
namelist['radiation']['longitude']  = -4.9267
namelist['radiation']['julian_day'] = 183
with open('inputLSM.json', 'w') as outfile:  
    json.dump(namelist,outfile,indent=4)