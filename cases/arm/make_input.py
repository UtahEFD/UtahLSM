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

"""Generates NetCDF initial condition and namelist input for the ARM case."""

from __future__ import annotations

import json
import time
from typing import Any

import netCDF4 as nc
import numpy as np
from numpy.typing import NDArray


def _clean_variable(
    dataset: nc.Dataset,
    name: str,
    *,
    min_value: float | None = None,
    max_value: float | None = None,
) -> NDArray[np.float64]:
    """Return a variable with missing, failed-QC, and out-of-range values masked."""
    variable = dataset.variables[name]
    values = np.asarray(variable[:], dtype=float)
    valid = np.isfinite(values)

    for attr in ("missing_value", "_FillValue"):
        if hasattr(variable, attr):
            valid &= values != float(getattr(variable, attr))

    if min_value is not None:
        valid &= values >= min_value
    if max_value is not None:
        valid &= values <= max_value

    qc_name = f"qc_{name}"
    if qc_name in dataset.variables:
        qc = np.asarray(dataset.variables[qc_name][:])
        valid &= qc == 0

    cleaned = values.astype(float)
    cleaned[~valid] = np.nan
    return cleaned

def _interp_clean(
    source_time: NDArray[np.float64],
    source_values: NDArray[np.float64],
    target_time: NDArray[np.float64],
    name: str,
) -> NDArray[np.float64]:
    """Linearly interpolate finite source values onto target_time."""
    valid = np.isfinite(source_values)
    if np.count_nonzero(valid) < 2:
        raise RuntimeError(f"{name} has fewer than two valid samples after QC.")
    if target_time[0] < source_time[valid][0] or target_time[-1] > source_time[valid][-1]:
        raise RuntimeError(
            f"{name} valid samples do not cover the requested forcing period."
        )
    return np.interp(target_time, source_time[valid], source_values[valid])

def _read_interp(
    dataset: nc.Dataset,
    name: str,
    target_time: NDArray[np.float64],
    *,
    min_value: float | None = None,
    max_value: float | None = None,
) -> NDArray[np.float64]:
    """Read a QC-cleaned variable and interpolate it to the forcing timeline."""
    source_time = np.asarray(dataset.variables['time'][:], dtype=float)
    values = _clean_variable(dataset, name, min_value=min_value, max_value=max_value)
    return _interp_clean(source_time, values, target_time, name)


def _initial_depth_mean(
    dataset: nc.Dataset,
    variable_prefix: str,
    sites: tuple[str, ...],
    *,
    min_value: float | None = None,
    max_value: float | None = None,
) -> NDArray[np.float64]:
    """Return the first-time mean profile across sites after QC filtering."""
    stack = np.stack([
        _clean_variable(
            dataset,
            f"{variable_prefix}_{site}",
            min_value=min_value,
            max_value=max_value,
        )[0]
        for site in sites
    ])
    valid = np.isfinite(stack)
    counts = np.sum(valid, axis=0)
    totals = np.nansum(stack, axis=0)
    return np.divide(
        totals,
        counts,
        out=np.full(stack.shape[1], np.nan),
        where=counts > 0,
    )


##########################################################################
# ARM SGP site in Lamont, Oklahoma (2017-06-17 18 UTC -> 2017-06-18 18 UTC
##########################################################################

# Soil column: uniform layers through 1 m. Observations below 50 cm are not
# trusted for initialization, but stopping the model at 50 cm makes the 50-cm
# validation level a lower-boundary node. Extending below the trusted profile
# keeps 50 cm prognostic while holding the deepest observed T/q values below it.
dz: float = 0.01
soil_depth: float = 1.0
rooting_depth: float = 0.5
soil_zlev_int: NDArray[np.float64] = np.linspace(0.0, soil_depth, round(soil_depth / dz) + 1)
nsoil: int = len(soil_zlev_int)

############################################################################
# Soil temperature and moisture data taken from STAMP instrument
# https://www.arm.gov/publications/tech_reports/handbooks/stamp_handbook.pdf
############################################################################

# soil levels are 5cm, 10cm, 20cm, 50cm, 75cm, 100cm
# soil types are SiL, SiL, C, CL, CL, CL. The 75/100-cm data are excluded
# because they are frequently flagged and are not representative enough for
# this one-day initialization. The texture break is represented halfway
# between the 10-cm silty-loam and 20-cm clay sensors; the moisture profile
# keeps the observed 10-cm and 20-cm values as separate layer plateaus instead
# of linearly smearing the sharp textural storage jump.
soil_texture_break: float = 0.15

soil_data: nc.Dataset = nc.Dataset('observations/arm_soil_data.nc')
soil_data_lev_raw: NDArray[np.float64] = (
    np.asarray(soil_data.variables['depth'][:], dtype=float) / 100.0
)
sites = ("west", "east", "south")

# Average temperature across the three measurement sites
soil_temp_obs_raw: NDArray[np.float64] = _initial_depth_mean(
    soil_data, "soil_temperature", sites, min_value=-50.0, max_value=80.0
)[:-2] + 273.15
valid_temp = np.isfinite(soil_temp_obs_raw)
if np.count_nonzero(valid_temp) < 2:
    raise RuntimeError("ARM soil temperature profile has fewer than two valid depths.")
soil_temp_lev = np.insert(soil_data_lev_raw[:-2][valid_temp], 0, 0.0)
soil_temp_obs = np.insert(
    soil_temp_obs_raw[valid_temp], 0, soil_temp_obs_raw[valid_temp][0]
)
soil_temp_ini: NDArray[np.float64] = np.interp(
    soil_zlev_int, soil_temp_lev, soil_temp_obs
)

soil_mois_obs_raw: NDArray[np.float64] = _initial_depth_mean(
    soil_data, "soil_specific_water_content", sites, min_value=0.0, max_value=100.0
)[:-2] / 100.0
valid_mois = np.isfinite(soil_mois_obs_raw)
if np.count_nonzero(valid_mois) < 2:
    raise RuntimeError("ARM soil moisture profile has fewer than two valid depths.")
soil_mois_lev = np.insert(soil_data_lev_raw[:-2][valid_mois], 0, 0.0)
soil_mois_obs = np.insert(
    soil_mois_obs_raw[valid_mois], 0, soil_mois_obs_raw[valid_mois][0]
)
soil_mois_ini: NDArray[np.float64] = np.interp(
    soil_zlev_int, soil_mois_lev, soil_mois_obs
)
if np.any(np.isclose(soil_mois_lev, 0.10)) and np.any(
    np.isclose(soil_mois_lev, 0.20)
):
    soil_mois_05 = float(soil_mois_obs[np.argmin(np.abs(soil_mois_lev - 0.05))])
    soil_mois_10 = float(soil_mois_obs[np.argmin(np.abs(soil_mois_lev - 0.10))])
    soil_mois_20 = float(soil_mois_obs[np.argmin(np.abs(soil_mois_lev - 0.20))])
    soil_mois_50 = float(soil_mois_obs[np.argmin(np.abs(soil_mois_lev - 0.50))])
    shallow = soil_zlev_int <= 0.10
    transition = np.logical_and(
        soil_zlev_int > 0.10, soil_zlev_int < soil_texture_break
    )
    clay_plateau = np.logical_and(
        soil_zlev_int >= soil_texture_break, soil_zlev_int <= 0.20
    )
    deeper = soil_zlev_int > 0.20
    soil_mois_ini[shallow] = np.interp(
        soil_zlev_int[shallow],
        np.array([0.0, 0.05, 0.10]),
        np.array([soil_mois_05, soil_mois_05, soil_mois_10]),
    )
    soil_mois_ini[transition] = soil_mois_10
    soil_mois_ini[clay_plateau] = soil_mois_20
    soil_mois_ini[deeper] = np.interp(
        soil_zlev_int[deeper],
        np.array([0.20, 0.50]),
        np.array([soil_mois_20, soil_mois_50]),
    )

stype: NDArray[np.str_] = np.empty(nsoil,dtype='S15')
stype[soil_zlev_int<soil_texture_break] = 'silty_loam'
stype[
    np.logical_and(soil_zlev_int >= soil_texture_break, soil_zlev_int < 0.50)
] = 'clay'
stype[soil_zlev_int>=0.5] = 'clay_loam'

# initialization file
init             = nc.Dataset('lsm_init.nc','w')
init.description = "UtahLSM input file"
init.source      = "Jeremy A. Gibbs"
init.history     = "Created " + time.ctime(time.time())

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

# write initial data
init_z[:] = soil_zlev_int
init_T[:] = soil_temp_ini
init_q[:] = soil_mois_ini
init_i[:] = stype
init.close()

#################################
# Read met tower data for offline
#################################

# UtahLSM offline mode expects 1-minute near-surface forcing. The ARM radiation
# and turbulent-flux observations are 30-minute averages, so those fields are
# QC-filtered and linearly interpolated onto the MET timeline.
met: nc.Dataset = nc.Dataset('observations/arm_surf_metr.nc')
tm: NDArray[np.float64] = np.asarray(met.variables['time'][:], dtype=float)
ws: NDArray[np.float64] = _read_interp(met, 'wspd_vec_mean', tm, min_value=0.0)
pt: NDArray[np.float64] = _read_interp(met, 'temp_mean', tm) + 273.15
pa_kpa: NDArray[np.float64] = _read_interp(met, 'atmos_pressure', tm, min_value=0.0)
pv_kpa: NDArray[np.float64] = _read_interp(met, 'vapor_pressure_mean', tm, min_value=0.0)

# UtahLSM casts precipitation as water mass flux (kg m^-2 s^-1). ARM PWD is
# mm/hr, and 1 mm water over 1 m^2 is 1 kg m^-2.
pr_mm_hr: NDArray[np.float64] = _read_interp(
    met, 'pwd_precip_rate_mean_1min', tm, min_value=0.0
)
pr: NDArray[np.float64] = pr_mm_hr / 3600.0
met.close()

# Pressure is stored in Pa in lsm_offline.nc. Specific humidity is kg/kg.
pa: NDArray[np.float64] = pa_kpa * 1000.0
qs: NDArray[np.float64] = 0.622 * pv_kpa / (pa_kpa - pv_kpa)

dt: float = round(tm[1] - tm[0])
ntime: int = len(tm)
t_utc: NDArray[np.float64] = (np.round(tm) % 86400).astype(float)

###########################################
# Read radiation and flux data for offline
###########################################

radn: nc.Dataset = nc.Dataset('observations/arm_surf_radn.nc')
swd: NDArray[np.float64] = _read_interp(radn, 'down_short_hemisp', tm, min_value=0.0)
swu: NDArray[np.float64] = _read_interp(radn, 'up_short_hemisp', tm, min_value=0.0)
lwd: NDArray[np.float64] = _read_interp(radn, 'down_long', tm, min_value=0.0)
lwu: NDArray[np.float64] = _read_interp(radn, 'up_long', tm, min_value=0.0)

# ARM surface_soil_heat_flux_avg is positive upward. UtahLSM's SEB closure term
# below uses G positive into the ground.
ghf_up: NDArray[np.float64] = _read_interp(radn, 'surface_soil_heat_flux_avg', tm)
ghf_obs: NDArray[np.float64] = -ghf_up
radn.close()

flux: nc.Dataset = nc.Dataset('observations/arm_surf_flux.nc')
shf_obs: NDArray[np.float64] = _read_interp(flux, 'corrected_sensible_heat_flux', tm)
lhf_obs: NDArray[np.float64] = _read_interp(flux, 'corrected_latent_heat_flux', tm)
flux.close()

# Prescribed unresolved SEB storage/closure term:
#   Rn - H - LE - G - seb_storage = 0
# where Rn = SWD - SWU + LWD - LWU, H and LE are positive upward, and G is
# positive downward into the ground.
rnet_obs: NDArray[np.float64] = swd - swu + lwd - lwu
seb_storage: NDArray[np.float64] = rnet_obs - shf_obs - lhf_obs - ghf_obs

##############################
# Write all time-series data #
##############################

# time-series file
metr: nc.Dataset = nc.Dataset('lsm_offline.nc','w')
metr.description = "UtahLSM input file for offline run"
metr.source      = "Jeremy A. Gibbs"
metr.history     = "Created " + time.ctime(time.time())

# add dimensions
metr.createDimension('t', ntime)

# add variables
metr_step = metr.createVariable("tstep", "f8", ())
metr_step.long_name = "time step for input offline data"
metr_step.units = "s"
metr_wind = metr.createVariable("atm_U", "f8", ("t"))
metr_wind.long_name = "wind speed"
metr_wind.units = "m s-1"
metr_temp = metr.createVariable("atm_T", "f8", ("t"))
metr_temp.long_name = "temperature"
metr_temp.units = "K"
metr_mixr = metr.createVariable("atm_q", "f8", ("t"))
metr_mixr.long_name = "specific humidity"
metr_mixr.units = "kg kg-1"
metr_pres = metr.createVariable("atm_p", "f8", ("t"))
metr_pres.long_name = "pressure"
metr_pres.units = "Pa"
metr_rain = metr.createVariable("precip", "f8", ("t"))
metr_rain.long_name = "precipitation (water mass flux)"
metr_rain.units = "kg m-2 s-1"
metr_swd = metr.createVariable("sw_in", "f8", ("t"))
metr_swd.long_name = "downwelling shortwave radiation"
metr_swd.units = "W m-2"
metr_lwd = metr.createVariable("lw_in", "f8", ("t"))
metr_lwd.long_name = "downwelling longwave radiation"
metr_lwd.units = "W m-2"
metr_sto = metr.createVariable("seb_storage", "f8", ("t"))
metr_sto.long_name = "prescribed surface energy storage or closure term"
metr_sto.units = "W m-2"
metr_sto.comment = (
    "Computed as SWD-SWU+LWD-LWU-H-LE-G from SEBS and ECOR observations; "
    "G is positive into the ground."
)
metr_rnet = metr.createVariable("R_net_obs", "f8", ("t"))
metr_rnet.long_name = "observed net radiation"
metr_rnet.units = "W m-2"
metr_shf = metr.createVariable("H_obs", "f8", ("t"))
metr_shf.long_name = "observed sensible heat flux"
metr_shf.units = "W m-2"
metr_shf.positive = "up"
metr_lhf = metr.createVariable("LE_obs", "f8", ("t"))
metr_lhf.long_name = "observed latent heat flux"
metr_lhf.units = "W m-2"
metr_lhf.positive = "up"
metr_ghf = metr.createVariable("G_obs", "f8", ("t"))
metr_ghf.long_name = "observed ground heat flux"
metr_ghf.units = "W m-2"
metr_ghf.positive = "down"

# write time-series data
metr_step[:] = float(dt)
metr_wind[:] = ws
metr_temp[:] = pt
metr_mixr[:] = qs
metr_pres[:] = pa
metr_rain[:] = pr
metr_swd[:] = swd
metr_lwd[:] = lwd
metr_sto[:] = seb_storage
metr_rnet[:] = rnet_obs
metr_shf[:] = shf_obs
metr_lhf[:] = lhf_obs
metr_ghf[:] = ghf_obs

# close file
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

# general section
namelist['general']['log_level'] = "info"

# numerics section
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

# time section
namelist['time']['utc_start'] = round(float(t_utc[0]))
namelist['time']['utc_year'] = 2017
namelist['time']['julian_day'] = 168

# grid section
namelist['grid']['nx'] = 1
namelist['grid']['ny'] = 1
namelist['grid']['nz'] = nsoil

# surface section
namelist['surface']['z_o']        = float(0.0500)
namelist['surface']['z_t']        = float(0.0005)
namelist['surface']['z_m']        = float(10.0)
namelist['surface']['z_s']        = float(2.0)
namelist['surface']['albedo']     = float(0.19)
namelist['surface']['emissivity'] = float(0.96)
namelist['surface']['model']      = "most"
namelist['surface']['psi_stable'] = "beljaars-holtslag"
namelist['surface']['zeta_max'] = 1.0
namelist['surface']['gustiness'] = 1.0
namelist['surface']['gustiness_stable_only'] = True

# soil section
namelist['soil']['properties'] = "clapp-hornberger"
namelist['soil']['model'] = "van-genuchten"
namelist['soil']['thermal_conductivity_model'] = "johansen"

# canopy settings: ARM SGP C1 in late June. Observed albedo of ~0.19
# indicates active green pasture rather than wheat stubble, so values
# lean to the warm-season pasture end (LAI 2.5, near-full cover, C4
# t_opt). rooting_depth is held at the previous effective 0.5 m value so
# extending the soil column does not also change transpiration access.
namelist['canopy']['model'] = "jarvis"
namelist['canopy']['lai'] = 2.5
namelist['canopy']['veg_fraction'] = 0.95
namelist['canopy']['rooting_depth'] = min(rooting_depth, soil_depth)
namelist['canopy']['beta'] = 0.943
namelist['canopy']['rs_min'] = 130.0
namelist['canopy']['rs_max'] = 5000.0
namelist['canopy']['rg_half'] = 50.0
namelist['canopy']['vpd_coef'] = 1.0e-4
namelist['canopy']['t_opt'] = 300.0
namelist['canopy']['t_coef'] = 1.6e-3
# ARM flux partitioning with this canopy is over-coupled thermally at the
# soil surface: G is too large and H is too small. A stronger in-canopy
# ground resistance keeps the soil heat exchange closer to observed G0.
namelist['canopy']['r_ground'] = 400.0

# radiation section
namelist['radiation']['model']     = "forcing"
namelist['radiation']['latitude']  = float(36.6049)
namelist['radiation']['longitude'] = float(-97.4856)

# output section
namelist['output']['save']   = True
namelist['output']['fields'] = ['all']

with open('lsm_namelist.json', 'w', encoding='utf-8') as outfile:
    json.dump(namelist,outfile,indent=4)
