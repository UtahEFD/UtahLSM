#
# UtahLSM
#
# Copyright (c) 2017–2026 Jeremy A. Gibbs
# Copyright (c) 2017–2026 Rob Stoll
# Copyright (c) 2017–2026 Eric Pardyjak
# Copyright (c) 2017–2026 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#

"""Core Data Models for UtahLSM.

This module defines the core data structures used throughout the UtahLSM model.
These structures are implemented as Python `dataclasses` to provide a clear
and robust way to manage the model's state and configuration. The module is
divided into two main sections:

1.  **State Data Models**: Represent the physical state of different components
    of the model at a given time (e.g., `AtmosphericState`, `SoilState`).
2.  **Configuration Data Models**: Hold the static parameters and settings
    loaded from the `lsm_namelist.json` file (e.g., `TimeConfig`, `GridConfig`).
"""

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from ._types import FloatOrArray, FloatOrArrayLike

# --- State Data Models ---

@dataclass
class AtmosphericState:
    """Holds the state of the atmosphere at a given time step.

    This data is provided by an external forcing file or a coupled
    atmospheric model, and represents the near-surface atmospheric
    conditions driving the land-surface model.

    Attributes:
        wind_speed: Wind speed [m/s] (scalar or per-column array).
        temperature: Air temperature [K] (scalar or per-column array).
        specific_humidity: Specific humidity [kg/kg] (scalar or per-column array).
        pressure: Atmospheric pressure [Pa] (scalar or per-column array).
        sw_in: Downwelling shortwave radiation [W/m^2] (forcing input).
        lw_in: Downwelling longwave radiation [W/m^2] (forcing input).
        sw_out: Upwelling (reflected) shortwave radiation [W/m^2]
            (diagnostic; written by the radiation/SEB solvers from T_s
            and the surface optical properties).
        lw_out: Upwelling (emitted) longwave radiation [W/m^2]
            (diagnostic; written by the radiation/SEB solvers from T_s
            and the surface optical properties).
        radiation_net: Net radiation [W/m^2] (diagnostic; equals
            ``sw_in - sw_out + lw_in - lw_out`` from the converged T_s).
        seb_storage: Prescribed surface-energy storage/closure term [W/m^2].
            The SEB is solved as ``radiation_net - H - LE - G
            - seb_storage = 0``. Defaults to zero for energy-conserving
            model runs.
    """
    wind_speed: FloatOrArray = 0.0
    temperature: FloatOrArray = 0.0
    specific_humidity: FloatOrArray = 0.0
    pressure: FloatOrArray = 0.0
    sw_in: FloatOrArray = 0.0
    lw_in: FloatOrArray = 0.0
    sw_out: FloatOrArray = 0.0
    lw_out: FloatOrArray = 0.0
    radiation_net: FloatOrArray = 0.0
    seb_storage: FloatOrArray = 0.0

@dataclass
class SoilState:
    """Holds the prognostic state of the soil column.

    This class represents the vertical profile of temperature and moisture
    within the soil, which evolves over time by the model's diffusion solvers.

    Attributes:
        temperature: Soil temperature profile [K] (nz or nz-by-ncol).
        moisture: Soil moisture profile [m^3/m^3] (nz or nz-by-ncol).
        type: Soil type name for each layer (string, e.g., 'clay', 'sand',
            'b11'). Names are lowercase and must match keys in the loaded
            soil properties dataset.
    """
    temperature: NDArray[np.float64] = field(
        default_factory=lambda: np.array([]))
    moisture: NDArray[np.float64] = field(
        default_factory=lambda: np.array([]))
    type: NDArray[np.object_] = field(
        default_factory=lambda: np.array([], dtype=object))

@dataclass
class SurfaceFluxes:
    """Holds all surface flux quantities.

    Attributes:
        kinematic_heat: Kinematic heat flux (w'T') [K m/s].
        kinematic_moisture: Kinematic moisture flux (w'q') [kg/kg m/s].
        sensible_heat: Sensible heat flux [W/m^2].
        latent_heat: Latent heat flux [W/m^2].
        ground_heat: Ground heat flux [W/m^2].
    """
    kinematic_heat: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    kinematic_moisture: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    sensible_heat: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    latent_heat: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    ground_heat: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))

@dataclass
class TurbulenceScales:
    """Holds Monin-Obukhov turbulence scales.

    Attributes:
        friction_velocity: Friction velocity (u*) [m/s].
        obukhov_length: Obukhov length (L) [m].
    """
    friction_velocity: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    obukhov_length: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))

@dataclass
class SurfaceState:
    """Holds the diagnostic state of the land surface at a given time step.

    These variables are calculated by the model and represent the interaction
    between the soil, the surface, and the atmosphere.

    Attributes:
        temperature: Radiative skin temperature [K] (scalar or per-column
            array). This is the temperature that closes the SEB and
            drives radiative emission and turbulent fluxes.
        soil_top_temperature: Temperature at the top of the mineral soil
            column [K]. Equal to ``temperature`` when there is no
            in-canopy thermal resistance; otherwise lower (cooler) by
            ``G * r_canopy_thermal``. Used as the Dirichlet BC for the
            heat-diffusion solver.
        moisture: Surface moisture content [kg/kg] (scalar or per-column array).
        specific_humidity: Surface-air specific humidity [kg/kg]
            (scalar or per-column array).
        fluxes: A dataclass containing all surface fluxes.
        turbulence: A dataclass containing turbulence scales.
        seb_residual: Final post-solve surface energy budget residual
            ``Rn - H - LE - G - seb_storage`` [W/m^2], one per column.
            Diagnostic for monitoring thermal drift over long integrations;
            the magnitude reflects how cleanly the SEB closed at the
            converged Obukhov length.
        air_density: Moist-air density [kg/m^3] cached once per timestep in
            ``_load_atm_state`` from the atmospheric state.
    """
    temperature: FloatOrArray = 0.0
    soil_top_temperature: FloatOrArray = 0.0
    moisture: FloatOrArray = 0.0
    specific_humidity: FloatOrArray = 0.0
    fluxes: SurfaceFluxes = field(default_factory=SurfaceFluxes)
    turbulence: TurbulenceScales = field(default_factory=TurbulenceScales)
    seb_residual: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    air_density: FloatOrArray = 0.0

@dataclass
class CanopyState:
    """Holds canopy / vegetation diagnostics for a single time step.

    All per-column quantities carry shape ``(ncol,)``; per-layer-per-
    column quantities carry shape ``(nz, ncol)``.

    Attributes:
        resistance: Bulk stomatal + cuticular resistance r_s [s/m].
        theta_root: Root-zone-weighted soil moisture [m^3/m^3].
        transpiration: Canopy transpiration flux [kg/m^2/s].
        wet_evaporation: Wet-canopy vapor exchange [kg/m^2/s]. Positive
            values evaporate canopy water storage; negative values are
            dewfall/condensation into storage.
        evap_soil: Bare-soil evaporation flux [kg/m^2/s]. Sum with
            ``transpiration`` and ``wet_evaporation`` matches total E.
        water_storage: Canopy intercepted/dew water storage [kg/m^2].
        water_capacity: Maximum canopy water storage [kg/m^2].
        latent_veg: Latent heat flux contribution from dry canopy
            transpiration [W/m^2].
        latent_wet: Latent heat flux contribution from wet canopy
            evaporation/dewfall [W/m^2].
        latent_soil: Latent heat flux contribution from bare soil [W/m^2].
        root_uptake: Per-layer root extraction rate [m^3 water /
            m^3 soil / s], shape (nz, ncol). Negative sign in the
            moisture budget (sink from soil to canopy).
    """
    resistance: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    theta_root: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    transpiration: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    wet_evaporation: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    evap_soil: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    water_storage: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    water_capacity: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    latent_veg: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    latent_wet: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    latent_soil: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1))
    root_uptake: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros((1, 1)))

@dataclass
class SolverState:
    """Holds intermediate variables used by the numerical solvers.

    This is a convenience class to store values that are calculated in one
    part of a solver and needed in another, avoiding recalculation and
    clarifying the data flow within complex numerical schemes.

    Attributes:
        conductivity_thermal_mid: Thermal conductivity at the midpoint between
            the top two soil layers [W/m/K].
        diffusion_e: Sub-diagonal buffer for the tridiagonal diffusion solver,
            shape (nz - 1, ncol). Pre-allocated once to avoid per-timestep
            heap churn; reused by both heat and moisture linear solves.
        diffusion_f: Main-diagonal buffer, shape (nz - 1, ncol).
        diffusion_g: Super-diagonal buffer, shape (nz - 1, ncol).
        diffusion_r: Right-hand-side buffer, shape (nz - 1, ncol).
    """
    conductivity_thermal_mid: FloatOrArray = 0.0
    diffusion_e: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(0))
    diffusion_f: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(0))
    diffusion_g: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(0))
    diffusion_r: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(0))

@dataclass(frozen=True)
class ForcingData:
    """Represents the entire time-series of meteorological forcing data.

    This class is the forcing data provided by an offline file and is not
    used if the land-surface model is driven by a coupled atmospheric model.

    Attributes:
        ntime: The number of time steps in the forcing data.
        tstep: The time step interval [s].
        atmos: A list of `AtmosphericState` objects, one for each time step.
    """
    ntime: int
    tstep: float
    # A list of atmospheric states, one for each timestep
    atmos: list[AtmosphericState]

# --- Configuration Data Models ---

@dataclass(frozen=True)
class GeneralConfig:
    """General simulation settings.

    Attributes:
        log_level: Logging level for the simulation (e.g., 'info', 'debug').
    """
    log_level: str

@dataclass(frozen=True)
class IterationsConfig:
    """Maximum iterations for looping procedures.

    Several fields require an iterative approach to solve. This
    dataclass sets a maximum number of iterations to reach convergence.

    Attributes:
        sfc_flux: iterations to solve Obukhov length.
        seb_bracket: iterations to find root brackets.
        seb_root: iterations to find seb root.
        smb_flux: iterations to solve the surface moisture budget root.
        moisture_picard: iterations for the mixed-form soil moisture Picard
            solve.
        coupling: iterations to solve coupled system.
    """
    sfc_flux: int
    seb_bracket: int
    seb_root: int
    smb_flux: int
    moisture_picard: int
    coupling: int

@dataclass(frozen=True)
class TolerancesConfig:
    """Convergence criteria for fields requiring an iterative solution.

    Several fields require an iterative approach to solve. This
    dataclass sets a tolerance needed to achieve convergence.

    Attributes:
        sfc_flux: tolerance for Obukhov length.
        seb_root: tolerance for seb root.
        smb_flux: tolerance for SMB root-finding on surface moisture [m3/m3].
        moisture_picard: convergence tolerance for the mixed-form soil
            moisture Picard solve [m3/m3].
        moisture_bounds: admissible post-solve soil moisture overshoot before
            clipping or failure [m3/m3].
        coupling_temp: tolerance for soil temperature in coupling.
        coupling_mois: tolerance for soil moisture in coupling.
    """
    sfc_flux: float
    seb_root: float
    smb_flux: float
    moisture_picard: float
    moisture_bounds: float
    coupling_temp: float
    coupling_mois: float

@dataclass(frozen=True)
class NumericsConfig:
    """Numerical scheme parameters.

    Attributes:
        heat_diffusion_back_weight: Backward weighting factor for the
            soil heat diffusion solver (0.5 for Crank-Nicolson).
        iterations: a dataclass holding numerical iteration limits.
        tolerances: a dataclass holding numerical convergence criteria.
    """
    heat_diffusion_back_weight: float
    iterations: IterationsConfig
    tolerances: TolerancesConfig

@dataclass(frozen=True)
class TimeConfig:
    """Time-related parameters for the simulation.

    Attributes:
        utc_start: The starting time of the simulation in UTC seconds
            from midnight.
        utc_year: The year (UTC) at the start of the simulation.
        julian_day: The starting Julian day of the year.
    """
    utc_start: int
    utc_year: int
    julian_day: int

@dataclass(frozen=True)
class GridConfig:
    """Grid and spatial discretization parameters.

    Attributes:
        nx: Number of grid points in the x-direction.
        ny: Number of grid points in the y-direction.
        nz: Number of soil layers (grid points in the z-direction).
        z: Soil layer depths [m].
    """
    nx: int
    ny: int
    nz: int
    z : NDArray[np.float64]

@dataclass(frozen=True)
class SurfaceConfig:
    """Surface-related parameters.

    Attributes:
        z_o: Aerodynamic roughness length [m].
        z_t: Thermal roughness length [m].
        z_m: Measurement height for wind speed [m].
        z_s: Measurement height for temperature and humidity [m].
        albedo: Surface albedo (dimensionless).
        emissivity: Surface emissivity (dimensionless).
        model: Surface layer model selector ('most').
        psi_stable: Stable (z/L>=0) MOST integrated stability correction (ψ).
        zeta_max: Maximum |z/L| used to clamp Obukhov length for MOST.
        gustiness: Additional wind-speed magnitude [m/s] added in quadrature.
        gustiness_stable_only: Apply gustiness only when L>=0 if True.
    """
    z_o: float
    z_t: float
    z_m: float
    z_s: float
    albedo: float
    emissivity: float
    model: str
    psi_stable: str = "dyer-hicks"
    zeta_max: float = 5.0
    gustiness: float = 0.0
    gustiness_stable_only: bool = True

@dataclass(frozen=True)
class SoilConfig:
    """Soil model configuration.

    Attributes:
        properties: Name of soil property dataset (e.g., 'cosby-1984') or path
            to custom JSON file.
        model: Soil model selector ('brooks-corey', 'campbell', 'van-genuchten').
    """
    properties: str
    model: str

@dataclass(frozen=True)
class RadiationConfig:
    """Radiation model configuration.

    Attributes:
        model: Radiation model selector ('forcing', 'basic').
        latitude: Site latitude [degrees].
        longitude: Site longitude [degrees].
    """
    model: str
    latitude: float
    longitude: float

@dataclass(frozen=True)
class CanopyConfig:
    """Canopy / vegetation model configuration.

    Scalar parameters broadcast to every column when parsed; sequence
    parameters must have length ``ncol`` and retain their per-column
    values. ``model == 'none'`` skips canopy physics entirely
    (bare-soil fallback).

    Attributes:
        model: Canopy model selector ('none' or 'jarvis').
        lai: Leaf area index [m^2/m^2].
        veg_fraction: Vegetated surface fraction in [0, 1].
        rooting_depth: Depth over which roots integrate to unity [m].
        beta: Jackson-1996 root distribution parameter (dimensionless).
        rs_min: Minimum bulk stomatal resistance [s/m].
        rs_max: Maximum (cuticular) resistance [s/m].
        rg_half: Half-saturation net radiation for f1(R) [W/m^2]
            (Jarvis).
        vpd_coef: VPD sensitivity for f2(VPD) [1/Pa] (Jarvis).
        t_opt: Optimum air temperature for f3(T) [K] (Jarvis).
        t_coef: Width of f3(T) parabola [1/K^2] (Jarvis).
        r_ground: In-canopy aerodynamic resistance for heat transport
            from the radiative skin to the soil top [s/m]. Acts in
            series with the top-cell soil conductive resistance,
            scaled by ``veg_fraction``. ``0.0`` recovers the
            bare-skin behaviour where ``T_skin = T_soil_top``.
        water_capacity_lai: Canopy water holding capacity per LAI
            [kg/m^2 per LAI]. Total column capacity is
            ``veg_fraction * lai * water_capacity_lai``.
        wet_cooling_max: Maximum diagnostic nighttime wet-canopy cooling
            below the soil/radiative skin used for dewfall [K].
    """
    model: str = 'none'
    lai: FloatOrArrayLike = 0.0
    veg_fraction: FloatOrArrayLike = 0.0
    rooting_depth: FloatOrArrayLike = 0.0
    beta: FloatOrArrayLike = 0.965
    rs_min: FloatOrArrayLike = 40.0
    rs_max: FloatOrArrayLike = 5000.0
    rg_half: FloatOrArrayLike = 100.0
    vpd_coef: FloatOrArrayLike = 1.0e-4
    t_opt: FloatOrArrayLike = 298.0
    t_coef: FloatOrArrayLike = 1.6e-3
    r_ground: FloatOrArrayLike = 0.0
    water_capacity_lai: FloatOrArrayLike = 0.2
    wet_cooling_max: FloatOrArrayLike = 3.0

@dataclass(frozen=True)
class OutputConfig:
    """Output file configuration.

    Attributes:
        save: Boolean flag to enable or disable saving output.
        fields: A list of strings specifying which variables to save.
    """
    save: bool
    fields: list[str]
