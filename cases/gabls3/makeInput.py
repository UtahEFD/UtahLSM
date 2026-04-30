#!/usr/bin/env python
#
# UtahLSM
#
# Copyright (c) 2017-2025 Jeremy A. Gibbs
# Copyright (c) 2017-2025 Rob Stoll
# Copyright (c) 2017-2025 Eric Pardyjak
# Copyright (c) 2017-2025 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#

import json
import time

import netCDF4 as nc
import numpy as np


##########################################
# GABLS3 stable subset (2006-07-02 00 UTC -> 09 UTC)
##########################################

# Cabauw observations are saved every 10 minutes, so 144 per day.
# We run the 9 h stable subset that aligns with the LES studies of GABLS3
# stable conditions, not the full 24 h diurnal case.
steps_per_day = 144
steps_per_hour = 6
start_day = 2
start_hour_utc = 0
case_duration_hours = 9
tidx = ((start_day - 1) * steps_per_day) + (start_hour_utc * steps_per_hour)
tend = tidx + (case_duration_hours * steps_per_hour)

# Soil column: 61 uniform layers from 0 to 0.6 m. Matches lsm_namelist.json
# nz=61 and the uniform-spacing requirement enforced in Input._validate_uniform_soil_z.
dz = 0.01
soil_depth = 0.6
z_int = np.linspace(0.0, soil_depth, int(round(soil_depth / dz)) + 1)
nsoil = len(z_int)


##########################################
# Soil temperature from CESAR observations
##########################################

# Cabauw soil-heat observations at 0, 2, 4, 6, 8, 12, 20, 30, 50 cm. The
# 0.5 m sensor is the deepest available; values below 0.5 m are held at the
# 0.5 m reading rather than extrapolated.
soil_t_obs = nc.Dataset('observations/cesar_soil_heat_lb1_t10_v1.0_200607.nc')
z_obs_T_full = np.array([0.00, 0.02, 0.04, 0.06, 0.08, 0.12, 0.20, 0.30, 0.50])
ts_names = ['TS00', 'TS02', 'TS04', 'TS06', 'TS08', 'TS12', 'TS20', 'TS30', 'TS50']
# CESAR sensors occasionally drop out at individual time steps; skip any
# masked sensor before interpolating so np.interp doesn't propagate NaN.
st_vals = []
z_vals = []
for name, depth in zip(ts_names, z_obs_T_full):
    raw = soil_t_obs.variables[name][tidx]
    if np.ma.is_masked(raw):
        continue
    st_vals.append(float(raw) + 273.15)
    z_vals.append(depth)
soil_t_obs.close()
if not st_vals:
    raise RuntimeError(
        f"All soil-temperature sensors are masked at tidx={tidx}; cannot "
        f"build an initial profile."
    )
z_obs_T = np.array(z_vals)
st_ob = np.array(st_vals)
st_oi = np.interp(z_int, z_obs_T, st_ob)


##########################################
# Soil moisture from CESAR observations
##########################################

# Anchor the upper column with the 10-minute Campbell-calibrated TH probes
# (TH03/TH08/TH20) since they capture the actual start-time state, not the
# slowly-varying daily mean. The EB-field probes (TH05/TH19/TH33/...) are
# all flagged -9999 at this time so we cannot use them. Below the deepest
# valid TH probe (0.20 m) we splice in the previous-day TDR daily-mean
# profile to fill the deeper column. The TDR network shows large day-to-day
# swings, so it is treated as a fallback below the start-state TH probes
# rather than the primary anchor.
sm10 = nc.Dataset('observations/cesar_soil_water_lb1_t10_v1.1_200607.nc')
sm_03 = float(sm10.variables['TH03'][tidx])
sm_08 = float(sm10.variables['TH08'][tidx])
sm_20 = float(sm10.variables['TH20'][tidx])
sm10.close()

tdr = nc.Dataset('observations/cesar_tdr_soilmoisture_la1_t1d_v1.0_2006.nc')
jdi = 181  # 0-indexed prior day for July 2
sm_30 = float(np.mean([tdr.variables[f'SM{i}'][jdi] for i in (3, 9, 15, 21)]))
sm_45 = float(np.mean([tdr.variables[f'SM{i}'][jdi] for i in (4, 10, 16, 22)]))
sm_60 = float(np.mean([tdr.variables[f'SM{i}'][jdi] for i in (5, 11, 17, 23)]))
sm_73 = float(np.mean([tdr.variables[f'SM{i}'][jdi] for i in (6, 12, 18, 24)]))
tdr.close()

# Duplicate TH03 at z=0 so the surface knot reflects the start-time state.
z_obs_q = np.array([0.00, 0.03, 0.08, 0.20, 0.30, 0.45, 0.60, 0.725])
sm_ob = np.array([sm_03, sm_03, sm_08, sm_20, sm_30, sm_45, sm_60, sm_73])
sm_oi = np.interp(z_int, z_obs_q, sm_ob)

# Soil type: clay over peat. TH20=0.55 already exceeds cosby clay porosity
# (0.468), so the clay/peat boundary must sit above 0.20 m. We put it at
# 0.15 m: clay holds the upper column where θ <= 0.43, peat absorbs the
# wetter 0.20 m+ values. The Jarvis canopy roots span both layers (the
# upper clay has cosby wilt=0.220, the deeper peat has wilt=0.396); the
# wet TH20 anchor keeps the root-weighted moisture comfortably above the
# root-weighted wilting so f4 stays positive.
clay_peat_boundary = 0.15
stype = np.where(z_int <= clay_peat_boundary, 'clay', 'peat').astype('U8')


#######################
# Initialization file #
#######################

init = nc.Dataset('lsm_init.nc', 'w')
init.description = (
    "UtahLSM input file for the GABLS3 stable subset (2006-07-02 00 UTC -> "
    "09 UTC). Soil temperature and moisture are taken from the Cabauw CESAR "
    "observations (cesar_soil_heat_lb1, cesar_soil_water_lb1, "
    "cesar_tdr_soilmoisture_la1)."
)
init.source = "Jeremy A. Gibbs"
init.history = "Created " + time.ctime(time.time())

init.createDimension('z', nsoil)

init_z = init.createVariable("soil_z", "f8", ("z",))
init_z.long_name = "z-distance"
init_z.units = "m"
init_T = init.createVariable("soil_T", "f8", ("z",))
init_T.long_name = "soil temperature"
init_T.units = "K"
init_q = init.createVariable("soil_q", "f8", ("z",))
init_q.long_name = "soil moisture"
init_q.units = "m3 m-3"
init_i = init.createVariable("soil_type", "str", ("z",))
init_i.long_name = "soil type"
init_i.units = ""

init_z[:] = z_int
init_T[:] = st_oi
init_q[:] = sm_oi
init_i[:] = stype
init.close()


###################################
# Read met tower data for offline #
###################################

# UtahLSM offline mode expects near-surface forcing and net radiation, so we
# retain the Cabauw observation time series over the official 24 h GABLS3
# window even though the parent SCM case is forced geostrophically.
met = nc.MFDataset('observations/cesar_surface_meteo_lc1_t10_v1.0_200607.nc')
tm = np.asarray(met.variables['time'][tidx:tend], dtype=float) * 3600.0
ws = np.asarray(met.variables['F010'][tidx:tend], dtype=float)
pt = np.asarray(met.variables['TA002'][tidx:tend], dtype=float)
pa = np.asarray(met.variables['P0'][tidx:tend], dtype=float) * 100.0
qs = np.asarray(met.variables['Q002'][tidx:tend], dtype=float) / 1000.0
met.close()

# CESAR's time variable is stored as float32, so tm[1]-tm[0] comes out
# as 599.9977 s instead of an exact 600 s. Round to the nearest integer
# second so the model's runtime accumulator stays on round-number times
# (32400 s, not 32399.9 s).
dt = float(round(tm[1] - tm[0]))
ntime = len(tm)
t_utc = (np.round(tm) % 86400).astype(float)

rad = nc.MFDataset('observations/cesar_surface_radiation_lc1_t10_v1.0_200607.nc')
swu = np.asarray(rad.variables['SWU'][tidx:tend], dtype=float)
swd = np.asarray(rad.variables['SWD'][tidx:tend], dtype=float)
lwu = np.asarray(rad.variables['LWU'][tidx:tend], dtype=float)
lwd = np.asarray(rad.variables['LWD'][tidx:tend], dtype=float)
net = swd - swu + lwd - lwu
rad.close()


##############################
# Write all time-series data #
##############################

metr = nc.Dataset('lsm_offline.nc', 'w')
metr.description = "UtahLSM input file for offline run"
metr.source = "Jeremy A. Gibbs"
metr.history = "Created " + time.ctime(time.time())

metr.createDimension('t', ntime)

metr_s = metr.createVariable("tstep", "f8", ())
metr_s.long_name = "time step for input offline data"
metr_s.units = "s"
metr_u = metr.createVariable("atm_U", "f8", ("t",))
metr_u.long_name = "wind speed"
metr_u.units = "m s-1"
metr_t = metr.createVariable("atm_T", "f8", ("t",))
metr_t.long_name = "temperature"
metr_t.units = "K"
metr_q = metr.createVariable("atm_q", "f8", ("t",))
metr_q.long_name = "mixing ratio"
metr_q.units = "g g-1"
metr_p = metr.createVariable("atm_p", "f8", ("t",))
metr_p.long_name = "pressure"
metr_p.units = "Pa"
metr_r = metr.createVariable("R_net", "f8", ("t",))
metr_r.long_name = "net radiation"
metr_r.units = "W m-2"

metr_s[:] = dt
metr_u[:] = ws
metr_t[:] = pt
metr_q[:] = qs
metr_p[:] = pa
metr_r[:] = net
metr.close()


########################
# Settings for UtahLSM #
########################

namelist = {}
namelist['general'] = {}
namelist['numerics'] = {}
namelist['numerics']['iterations'] = {}
namelist['numerics']['tolerances'] = {}
namelist['time'] = {}
namelist['grid'] = {}
namelist['surface'] = {}
namelist['soil'] = {}
namelist['canopy'] = {}
namelist['radiation'] = {}
namelist['output'] = {}

namelist['general']['log_level'] = "info"

namelist['numerics']['heat_diffusion_back_weight'] = float(0.5)
namelist['numerics']['warm_start_turbulence'] = True
namelist['numerics']['initialize_surface_temperature_from_seb'] = True
namelist['numerics']['iterations']['sfc_flux'] = 100
namelist['numerics']['iterations']['seb_bracket'] = 100
namelist['numerics']['iterations']['seb_root'] = 100
namelist['numerics']['iterations']['smb_flux'] = 100
namelist['numerics']['iterations']['moisture_picard'] = 100
namelist['numerics']['iterations']['coupling'] = 20
namelist['numerics']['tolerances']['sfc_flux'] = 1e-3
namelist['numerics']['tolerances']['seb_root'] = 1e-6
namelist['numerics']['tolerances']['smb_flux'] = 1e-3
namelist['numerics']['tolerances']['moisture_picard'] = 1e-5
namelist['numerics']['tolerances']['moisture_bounds'] = 1e-3
namelist['numerics']['tolerances']['coupling_temp'] = 1e-2
namelist['numerics']['tolerances']['coupling_mois'] = 1e-5

namelist['time']['utc_start'] = int(round(t_utc[0]))
namelist['time']['utc_year'] = 2006
namelist['time']['julian_day'] = 183

namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1
namelist['grid']['nz'] = nsoil

# z_o = 0.03 m is the *local* roughness length for Cabauw grass. The 0.15 m
# value sometimes used in mesoscale / SCM setups represents an effective
# fetch-aware roughness for the surrounding terrain, not the local patch.
# We are running locally, so use the local value.
namelist['surface']['z_o'] = float(0.03)
namelist['surface']['z_t'] = float(0.0015)
namelist['surface']['z_m'] = float(10.0)
namelist['surface']['z_s'] = float(2.0)
namelist['surface']['albedo'] = float(0.33)
namelist['surface']['emissivity'] = float(0.99)
namelist['surface']['model'] = 1
namelist['surface']['psi_stable'] = "beljaars-holtslag"
namelist['surface']['zeta_max'] = float(1.0)
# Beljaars (1995) recommends 0.5-1.0 m/s gustiness under stable conditions;
# 2.0 m/s overdrives u* and prevents nocturnal decoupling. We apply it only
# under stable stratification so it doesn't double-count daytime convective
# wind variance, which is already handled by mean wind.
namelist['surface']['gustiness'] = float(2.0)
namelist['surface']['gustiness_stable_only'] = True

# Cosby et al. (1984) gives lower θ_wilt for clay (0.220 vs Clapp-Hornberger
# 0.287) which is closer to the observed Cabauw root-zone moisture and lets
# the Jarvis canopy actually transpire. Peat parameters are common across
# the three packaged datasets (the peat node here is for the deeper Cabauw
# layer, not a true bog).
namelist['soil']['properties'] = "cosby"
namelist['soil']['model'] = 1

# These Jarvis parameters remain a UtahLSM-specific canopy choice. Only LAI and
# vegetation fraction come directly from the published GABLS3 specification.
namelist['canopy']['model'] = "jarvis"
namelist['canopy']['lai'] = float(2.0)
namelist['canopy']['veg_fraction'] = float(1.0)
# Cabauw grass roots are concentrated in the upper 0.4 m of clay. The clay
# layer is now thick enough (0-0.4 m) to fully contain the root profile, so
# we use the original Jackson β=0.943 (typical for grasslands) without
# extending into peat.
namelist['canopy']['rooting_depth'] = float(0.4)
namelist['canopy']['beta'] = float(0.943)
namelist['canopy']['rs_min'] = float(40.0)
namelist['canopy']['rs_max'] = float(5000.0)
namelist['canopy']['rg_half'] = float(100.0)
namelist['canopy']['vpd_coef'] = float(1.0e-4)
namelist['canopy']['t_opt'] = float(298.0)
namelist['canopy']['t_coef'] = float(1.6e-3)
# In-canopy aerodynamic resistance for ground heat transport [s/m].
# Used in series with the top-cell soil conductive resistance and
# scaled by veg_fraction. Tuned against the Cabauw eddy-covariance
# day-2 fluxes: r_ground=300 s/m matches the observed H/LE/G partition
# and the TS00 surface temperature evolution to within ~0.1 K. Lies
# in the upper end of the Choudhury & Monteith (1988) range for a
# closed-canopy grassland.
namelist['canopy']['r_ground'] = float(300.0)

namelist['radiation']['model'] = 0
namelist['radiation']['latitude'] = float(51.9711)
namelist['radiation']['longitude'] = float(-4.9267)

namelist['output']['save'] = True
namelist['output']['fields'] = ['all']

with open('lsm_namelist.json', 'w', encoding='utf-8') as outfile:
    json.dump(namelist, outfile, indent=4)
