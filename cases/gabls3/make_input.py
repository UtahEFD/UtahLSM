#!/usr/bin/env python
#
# UtahLSM
#
# Copyright (c) 2017-2026 Jeremy A. Gibbs
# Copyright (c) 2017-2026 Rob Stoll
# Copyright (c) 2017-2026 Eric Pardyjak
# Copyright (c) 2017-2026 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#

from __future__ import annotations

import json
import time
from typing import Any

import netCDF4 as nc
import numpy as np
from numpy.typing import NDArray


####################################################
# GABLS3 stable subset (2006-07-02 00 UTC -> 09 UTC)
####################################################

# Soil column: 61 uniform layers from 0 to 0.6 m.
dz: float = 0.01
soil_depth: float = 0.6
soil_zlev_int: NDArray[np.float64] = np.linspace(0.0, soil_depth, int(round(soil_depth / dz)) + 1)
nsoil: int = len(soil_zlev_int)

##########################################
# Soil temperature from CESAR observations
##########################################

# Cabauw soil-heat observations at 0, 2, 4, 6, 8, 12, 20, 30, 50 cm.
soil_temp_dat: nc.Dataset = nc.Dataset('observations/gabls3_soil_heat.nc')
soil_temp_lev: NDArray[np.float64] = np.array([0.00, 0.02, 0.04, 0.06, 0.08, 0.12, 0.20, 0.30, 0.50])
soil_temp_var: list[str] = ['TS00', 'TS02', 'TS04', 'TS06', 'TS08', 'TS12', 'TS20', 'TS30', 'TS50']

# Skip any masked sensors before interpolating.
soil_temp_obs = np.ma.stack(
    [soil_temp_dat.variables[name][0] for name in soil_temp_var]
).astype(float)
soil_temp_dat.close()
valid_T: NDArray[np.bool_] = ~np.ma.getmaskarray(soil_temp_obs)
if not np.any(valid_T):
    raise RuntimeError(
        "All soil-temperature sensors are masked; cannot build an initial profile."
    )
soil_temp_lev: NDArray[np.float64] = soil_temp_lev[valid_T]
soil_temp_obs: NDArray[np.float64] = (
    np.asarray(soil_temp_obs[valid_T], dtype=float) + 273.15
)
soil_temp_ini: NDArray[np.float64] = np.interp(soil_zlev_int, soil_temp_lev, soil_temp_obs)

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
soil_mois_thc: nc.Dataset = nc.Dataset('observations/gabls3_soil_mois_thc.nc')
soil_mois_z03: float = float(soil_mois_thc.variables['TH03'][0])
soil_mois_z08: float = float(soil_mois_thc.variables['TH08'][0])
soil_mois_z20: float = float(soil_mois_thc.variables['TH20'][0])
soil_mois_thc.close()

soil_mois_tdr: nc.Dataset = nc.Dataset('observations/gabls3_soil_mois_tdr.nc')
soil_mois_z30: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (3,9,15,21)]))
soil_mois_z45: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (4,10,16,22)]))
soil_mois_z60: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (5,11,17,23)]))
soil_mois_z73: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (6,12,18,24)]))
soil_mois_tdr.close()

# Duplicate TH03 at z=0 so the surface knot reflects the start-time state.
soil_mois_lev: NDArray[np.float64] = np.array([0.00, 0.03, 0.08, 0.20, 0.30, 0.45, 0.60, 0.725])
soil_mois_obs: NDArray[np.float64] = np.array([soil_mois_z03, soil_mois_z03, soil_mois_z08, 
    soil_mois_z20, soil_mois_z30, soil_mois_z45, soil_mois_z60, soil_mois_z73])
soil_mois_ini: NDArray[np.float64] = np.interp(soil_zlev_int, soil_mois_lev, soil_mois_obs)

# Soil type: clay over peat. TH20=0.55 already exceeds cosby clay porosity
# (0.468), so the clay/peat boundary must sit above 0.20 m. We put it at
# 0.15 m: clay holds the upper column where θ <= 0.43, peat absorbs the
# wetter 0.20 m+ values. The Jarvis canopy roots span both layers (the
# upper clay has cosby wilt=0.220, the deeper peat has wilt=0.396); the
# wet TH20 anchor keeps the root-weighted moisture comfortably above the
# root-weighted wilting so f4 stays positive.
clay_peat_boundary: float = 0.15
stype: NDArray[np.str_] = np.where(soil_zlev_int <= clay_peat_boundary, 'clay', 'peat').astype('U8')

#######################
# Initialization file #
#######################

init = nc.Dataset('lsm_init.nc', 'w')
init.description = (
    "UtahLSM input file for the GABLS3 stable subset (2006-07-02 00 UTC -> "
    "09 UTC). Soil temperature and moisture are taken from the Cabauw CESAR "
    "observations."
)
init.source = "Jeremy A. Gibbs"
init.history = "Created " + time.ctime(time.time())

init.createDimension('z', nsoil)

init_z = init.createVariable("soil_z", "f8", ("z",))  # type: ignore[assignment]
init_z.long_name = "z-distance"
init_z.units = "m"
init_T = init.createVariable("soil_T", "f8", ("z",))  # type: ignore[assignment]
init_T.long_name = "soil temperature"
init_T.units = "K"
init_q = init.createVariable("soil_q", "f8", ("z",))  # type: ignore[assignment]
init_q.long_name = "soil moisture"
init_q.units = "m3 m-3"
init_i = init.createVariable("soil_type", "str", ("z",))  # type: ignore[assignment]
init_i.long_name = "soil type"
init_i.units = ""

init_z[:] = soil_zlev_int
init_T[:] = soil_temp_ini
init_q[:] = soil_mois_ini
init_i[:] = stype
init.close()

###################################
# Read met tower data for offline #
###################################

# UtahLSM offline mode expects near-surface forcing and net radiation, so we
# retain the Cabauw observation time series over the official 24 h GABLS3
# window even though the parent SCM case is forced geostrophically.
met: nc.MFDataset = nc.MFDataset('observations/gabls3_surf_metr.nc')
tm: NDArray[np.float64] = np.asarray(met.variables['time'][:], dtype=float) * 3600.0
ws: NDArray[np.float64] = np.asarray(met.variables['F010'][:], dtype=float)
pt: NDArray[np.float64] = np.asarray(met.variables['TA002'][:], dtype=float)
pa: NDArray[np.float64] = np.asarray(met.variables['P0'][:], dtype=float) * 100.0
qs: NDArray[np.float64] = np.asarray(met.variables['Q002'][:], dtype=float) / 1000.0
met.close()

# CESAR's time variable is stored as float32, so tm[1]-tm[0] comes out
# as 599.9977 s instead of an exact 600 s. Round to the nearest integer
# second so the model's runtime accumulator stays on round-number times
# (32400 s, not 32399.9 s).
dt: float = float(round(tm[1] - tm[0]))
ntime: int = len(tm)
t_utc: NDArray[np.float64] = (np.round(tm) % 86400).astype(float)

rad: nc.MFDataset = nc.MFDataset('observations/gabls3_surf_radn.nc')
swu: NDArray[np.float64] = np.asarray(rad.variables['SWU'][:], dtype=float)
swd: NDArray[np.float64] = np.asarray(rad.variables['SWD'][:], dtype=float)
lwu: NDArray[np.float64] = np.asarray(rad.variables['LWU'][:], dtype=float)
lwd: NDArray[np.float64] = np.asarray(rad.variables['LWD'][:], dtype=float)
rad.close()

flux: nc.MFDataset = nc.MFDataset('observations/gabls3_surf_flux.nc')
shf_obs: NDArray[np.float64] = np.asarray(flux.variables['H'][:], dtype=float)
lhf_obs: NDArray[np.float64] = np.asarray(flux.variables['LE'][:], dtype=float)
flux.close()

soil_heat: nc.MFDataset = nc.MFDataset('observations/gabls3_soil_heat.nc')
ghf_obs: NDArray[np.float64] = np.asarray(soil_heat.variables['FG0'][:], dtype=float)
soil_heat.close()

# Prescribed unresolved SEB storage/closure term. Positive values remove
# energy from the modeled H/LE/G partition:
#   Rn - H - LE - G - seb_storage = 0
# Here Rn comes from the independent radiation components and G is the
# Fourier-extrapolated 0 cm soil heat flux.
seb_storage: NDArray[np.float64] = lwd - lwu + swd - swu - shf_obs - lhf_obs - ghf_obs

##############################
# Write all time-series data #
##############################

metr = nc.Dataset('lsm_offline.nc', 'w')
metr.description = "UtahLSM input file for offline run"
metr.source = "Jeremy A. Gibbs"
metr.history = "Created " + time.ctime(time.time())

metr.createDimension('t', ntime)

metr_s = metr.createVariable("tstep", "f8", ())  # type: ignore[assignment]
metr_s.long_name = "time step for input offline data"
metr_s.units = "s"
metr_u = metr.createVariable("atm_U", "f8", ("t",))  # type: ignore[assignment]
metr_u.long_name = "wind speed"
metr_u.units = "m s-1"
metr_t = metr.createVariable("atm_T", "f8", ("t",))  # type: ignore[assignment]
metr_t.long_name = "temperature"
metr_t.units = "K"
metr_q = metr.createVariable("atm_q", "f8", ("t",))  # type: ignore[assignment]
metr_q.long_name = "mixing ratio"
metr_q.units = "g g-1"
metr_p = metr.createVariable("atm_p", "f8", ("t",))  # type: ignore[assignment]
metr_p.long_name = "pressure"
metr_p.units = "Pa"
metr_swd = metr.createVariable("sw_in", "f8", ("t",))  # type: ignore[assignment]
metr_swd.long_name = "downwelling shortwave radiation"
metr_swd.units = "W m-2"
metr_swu = metr.createVariable("sw_out", "f8", ("t",))  # type: ignore[assignment]
metr_swu.long_name = "upwelling (reflected) shortwave radiation"
metr_swu.units = "W m-2"
metr_lwd = metr.createVariable("lw_in", "f8", ("t",))  # type: ignore[assignment]
metr_lwd.long_name = "downwelling longwave radiation"
metr_lwd.units = "W m-2"
metr_lwu = metr.createVariable("lw_out", "f8", ("t",))  # type: ignore[assignment])
metr_lwu.long_name = "upwelling (emitted) longwave radiation"
metr_lwu.units = "W m-2"
metr_sto = metr.createVariable("seb_storage", "f8", ("t",))  # type: ignore[assignment]
metr_sto.long_name = "prescribed surface energy storage or closure term"
metr_sto.units = "W m-2"
metr_sto.comment = "Computed as LWD-LWU+SWD-SWU-H-LE-FG0 from CESAR observations"

metr_s[:] = dt
metr_u[:] = ws
metr_t[:] = pt
metr_q[:] = qs
metr_p[:] = pa
metr_swd[:] = swd
metr_swu[:] = swu
metr_lwd[:] = lwd
metr_lwu[:] = lwu
metr_sto[:] = seb_storage
metr.close()


########################
# Settings for UtahLSM #
########################

namelist: dict[str, Any] = {}
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
namelist['surface']['z_t'] = float(0.003)
namelist['surface']['z_m'] = float(10.0)
namelist['surface']['z_s'] = float(2.0)
namelist['surface']['albedo'] = float(0.22)
namelist['surface']['emissivity'] = float(0.99)
namelist['surface']['model'] = 1
namelist['surface']['psi_stable'] = "beljaars-holtslag"
namelist['surface']['zeta_max'] = float(1.0)
# Beljaars (1995) recommends 0.5-1.0 m/s gustiness under stable conditions;
# 2.0 m/s overdrives u* and prevents nocturnal decoupling. We apply it only
# under stable stratification so it doesn't double-count daytime convective
# wind variance, which is already handled by mean wind.
namelist['surface']['gustiness'] = float(1.0)
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
namelist['canopy']['rs_min'] = float(75.0)
namelist['canopy']['rs_max'] = float(5000.0)
namelist['canopy']['rg_half'] = float(30.0)
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
namelist['canopy']['r_ground'] = float(200.0)

namelist['radiation']['model'] = 0
namelist['radiation']['latitude'] = float(51.9711)
namelist['radiation']['longitude'] = float(4.9267)

namelist['output']['save'] = True
namelist['output']['fields'] = ['all']

with open('lsm_namelist.json', 'w', encoding='utf-8') as outfile:
    json.dump(namelist, outfile, indent=4)
