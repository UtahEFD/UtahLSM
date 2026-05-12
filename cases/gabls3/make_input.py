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

"""Generates NetCDF initial condition and namelist input for the GABLS3 case."""

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
soil_zlev_int: NDArray[np.float64] = np.linspace(0.0, soil_depth, round(soil_depth / dz) + 1)
nsoil: int = len(soil_zlev_int)

##########################################
# Soil temperature from CESAR observations
##########################################

# Cabauw soil-heat observations at 0, 2, 4, 6, 8, 12, 20, 30, 50 cm.
soil_temp_dat: nc.Dataset = nc.Dataset('observations/gabls3_soil_heat.nc')
soil_temp_lev: NDArray[np.float64] = np.array([0.00,0.02,0.04,0.06,0.08,0.12,0.20,0.30,0.50])
soil_temp_var: list[str] = ['TS00','TS02','TS04','TS06','TS08','TS12','TS20','TS30','TS50']

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
soil_temp_lev = soil_temp_lev[valid_T]
soil_temp_obs = (np.asarray(soil_temp_obs[valid_T], dtype=float) + 273.15)
soil_temp_ini: NDArray[np.float64] = np.interp(soil_zlev_int, soil_temp_lev, soil_temp_obs)

#######################################
# Soil moisture from CESAR observations
#######################################

# TH Campbell probes: TR-384 Section 20 validates 3 cm and 8 cm as reliable, but 20 cm is flagged
#   unreliable (can read above clay porosity).
#   see: https://cdn.knmi.nl/knmi/pdf/bibliotheek/knmipubTR/TR384.pdf
soil_mois_thc: nc.Dataset = nc.Dataset('observations/gabls3_soil_mois_thc.nc')
soil_mois_z03: float = soil_mois_thc.variables['TH03'][0]
soil_mois_z08: float = soil_mois_thc.variables['TH08'][0]
soil_mois_thc.close()

# TDR network (previous-day daily means) anchors the profile from 30 cm down. The 15 cm TDR knot
#   is omitted: its value would create an artificial decrease below TH08 due to inter-instrument
#   calibration offsets. Linear interpolation fills the 8–30 cm gap.
soil_mois_tdr: nc.Dataset = nc.Dataset('observations/gabls3_soil_mois_tdr.nc')
soil_mois_z30: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (3,9,15,21)]))
soil_mois_z45: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (4,10,16,22)]))
soil_mois_z60: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (5,11,17,23)]))
soil_mois_z73: float = float(np.mean([soil_mois_tdr.variables[f'SM{i}'][0] for i in (6,12,18,24)]))
soil_mois_tdr.close()

# Surface pinned to TH03; TH08 anchors shallow clay; TDR from 30 cm down.
soil_mois_lev: NDArray[np.float64] = np.array([0.00, 0.03, 0.08, 0.30, 0.45, 0.60, 0.725])
soil_mois_obs: NDArray[np.float64] = np.array([soil_mois_z03, soil_mois_z03, soil_mois_z08,
    soil_mois_z30, soil_mois_z45, soil_mois_z60, soil_mois_z73])
soil_mois_ini: NDArray[np.float64] = np.interp(soil_zlev_int, soil_mois_lev, soil_mois_obs)

# Soil type: Cabauw places clay in 0–18 cm, a clay-peat mix from 18–60 cm, and heavier peat below.
#   Boundary at 0.18 m. Cabauw subsoil is moderately decomposed Holocene peat -> Letts hemic.
clay_peat_boundary: float = 0.18
stype: NDArray[np.str_] = np.where(soil_zlev_int <= clay_peat_boundary,
                                   'clay',
                                   'peat_hemic').astype('U16')

#####################
# Initialization file
#####################

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

#################################
# Read met tower data for offline
#################################

# UtahLSM offline mode expects near-surface forcing and net radiation.
met: nc.MFDataset = nc.MFDataset('observations/gabls3_surf_metr.nc')
tm: NDArray[np.float64] = np.asarray(met.variables['time'][:], dtype=float) * 3600.0
ws: NDArray[np.float64] = np.asarray(met.variables['F010'][:], dtype=float)
pt: NDArray[np.float64] = np.asarray(met.variables['TA002'][:], dtype=float)
pa: NDArray[np.float64] = np.asarray(met.variables['P0'][:], dtype=float) * 100.0
qs: NDArray[np.float64] = np.asarray(met.variables['Q002'][:], dtype=float) / 1000.0
met.close()

dt: float = round(tm[1] - tm[0])
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

# Prescribed unresolved SEB storage/closure term: Rn - H - LE - G - seb_storage = 0
#   Rn comes from the independent radiation components and G is the Fourier-extrapolated
#   0-cm soil heat flux.
seb_storage: NDArray[np.float64] = lwd - lwu + swd - swu - shf_obs - lhf_obs - ghf_obs

############################
# Write all time-series data
############################

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
metr_lwd = metr.createVariable("lw_in", "f8", ("t",))  # type: ignore[assignment]
metr_lwd.long_name = "downwelling longwave radiation"
metr_lwd.units = "W m-2"
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
metr_lwd[:] = lwd
metr_sto[:] = seb_storage
metr.close()

######################
# Settings for UtahLSM
######################

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

# general settings
namelist['general']['log_level'] = "info"

# numerics settings
namelist['numerics']['heat_diffusion_back_weight'] = 0.5
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

# time settings
namelist['time']['utc_start'] = round(t_utc[0])
namelist['time']['utc_year'] = 2006
namelist['time']['julian_day'] = 183

# grid settings
namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1
namelist['grid']['nz'] = nsoil

# surface settings
namelist['surface']['z_o'] = 0.03
namelist['surface']['z_t'] = 0.003
namelist['surface']['z_m'] = 10.0
namelist['surface']['z_s'] = 2.0
namelist['surface']['albedo'] = 0.23
namelist['surface']['emissivity'] = 0.99
namelist['surface']['model'] = "most"
namelist['surface']['psi_stable'] = "beljaars-holtslag"
namelist['surface']['zeta_max'] = 1.0
namelist['surface']['gustiness'] = 1.0
namelist['surface']['gustiness_stable_only'] = True

# soil settings
namelist['soil']['properties'] = "cosby_letts"
namelist['soil']['model'] = "van-genuchten"
namelist['soil']['thermal_conductivity_model'] = "johansen"

# canopy settings
namelist['canopy']['model'] = "jarvis"
namelist['canopy']['lai'] = 2.0
namelist['canopy']['veg_fraction'] = 1.0
namelist['canopy']['rooting_depth'] = 0.4
namelist['canopy']['beta'] = 0.943
namelist['canopy']['rs_min'] = 100.0
namelist['canopy']['rs_max'] = 5000.0
namelist['canopy']['rg_half'] = 50.0
namelist['canopy']['vpd_coef'] = 1.0e-4
namelist['canopy']['t_opt'] = 298.0
namelist['canopy']['t_coef'] = 1.6e-3
namelist['canopy']['r_ground'] = 200.0

# radiation settings
namelist['radiation']['model'] = "forcing"
namelist['radiation']['latitude'] = 51.9711
namelist['radiation']['longitude'] = 4.9267

# output settings
namelist['output']['save'] = True
namelist['output']['fields'] = ['all']

# write the namelist
with open('lsm_namelist.json', 'w', encoding='utf-8') as outfile:
    json.dump(namelist, outfile, indent=4)
