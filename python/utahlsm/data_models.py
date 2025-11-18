#
# UtahLSM
#
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
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
from typing import List
import numpy as np
from numpy.typing import NDArray

# --- State Data Models ---

@dataclass
class AtmosphericState:
    """Holds the state of the atmosphere at a given time step.

    This data is provided by an external forcing file or a coupled
    atmospheric model, and represents the near-surface atmospheric
    conditions driving the land-surface model.

    Attributes:
        wind_speed: Wind speed [m/s].
        temperature: Air temperature [K].
        specific_humidity: Specific humidity [kg/kg].
        pressure: Atmospheric pressure [Pa].
        radiation_net: Net radiation [W/m^2].
    """
    wind_speed: float = 0.0
    temperature: float = 0.0
    specific_humidity: float = 0.0
    pressure: float = 0.0
    radiation_net: float = 0.0

@dataclass
class SoilState:
    """Holds the prognostic state of the soil column.

    This class represents the vertical profile of temperature and moisture
    within the soil, which evolves over time by the model's diffusion solvers.

    Attributes:
        temperature: Soil temperature profile [K].
        moisture: Soil moisture profile [m^3/m^3].
        type: Soil type name for each layer (string, e.g., 'clay', 'sand',
            'b11'). Names are lowercase and must match keys in the loaded
            soil properties dataset.
    """
    temperature: NDArray[np.float64] = field(
        default_factory=lambda: np.array([]))
    moisture: NDArray[np.float64] = field(
        default_factory=lambda: np.array([]))
    type: NDArray = field(
        default_factory=lambda: np.array([]))

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
        temperature: Surface temperature [K].
        moisture: Surface moisture content [kg/kg].
        specific_humidity: Surface-air specific humidity [kg/kg].
        fluxes: A dataclass containing all surface fluxes.
        turbulence: A dataclass containing turbulence scales.

    """
    temperature: float = 0.0
    moisture: float = 0.0
    specific_humidity: float = 0.0
    fluxes: SurfaceFluxes = field(default_factory=SurfaceFluxes)
    turbulence: TurbulenceScales = field(default_factory=TurbulenceScales)

@dataclass
class SolverState:
    """Holds intermediate variables used by the numerical solvers.

    This is a convenience class to store values that are calculated in one
    part of a solver and needed in another, avoiding recalculation and
    clarifying the data flow within complex numerical schemes.

    Attributes:
        conductivity_thermal_mid: Thermal conductivity at the midpoint between
            the top two soil layers [W/m/K].
    """
    conductivity_thermal_mid: float = 0.0

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
    atmos: List[AtmosphericState]

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
        smb_flux: iterations to solve soil moisture flux.
    """
    sfc_flux: int
    seb_bracket: int
    seb_root: int
    smb_flux: int

@dataclass(frozen=True)
class TolerancesConfig:
    """Convergence criteria for fields requiring an iterative solution.

    Several fields require an iterative approach to solve. This
    dataclass sets a tolerance needed to achieve convergence.

    Attributes:
        sfc_flux: tolerance for Obukhov length.
        seb_root: tolerance for seb root.
        smb_flux: tolerance for soil moisture flux.
    """
    sfc_flux: float
    seb_root: float
    smb_flux: float

@dataclass(frozen=True)
class NumericsConfig:
    """Numerical scheme parameters.

    Attributes:
        diffusion_back_weight: Backward weighting factor for the
            diffusion solver (0.5 for Crank-Nicolson).
        iterations: a dataclass holding numerical iteration limits.
        tolerances: a dataclass holding numerical convergence criteria.
    """
    diffusion_back_weight: float
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
        model: Integer ID for the surface layer model to use.
    """
    z_o: float
    z_t: float
    z_m: float
    z_s: float
    albedo: float
    emissivity: float
    model: int

@dataclass(frozen=True)
class SoilConfig:
    """Soil model configuration.

    Attributes:
        properties: Name of soil property dataset (e.g., 'cosby-1984') or path
            to custom JSON file.
        model: Integer ID for the soil physics model to use (1=BrooksCorey,
            2=Campbell, 3=VanGenuchten).
    """
    properties: str
    model: int

@dataclass(frozen=True)
class RadiationConfig:
    """Radiation model configuration.

    Attributes:
        model: Integer ID for the radiation model to use.
        latitude: Site latitude [degrees].
        longitude: Site longitude [degrees].
    """
    model: int
    latitude: float
    longitude: float

@dataclass(frozen=True)
class OutputConfig:
    """Output file configuration.

    Attributes:
        save: Boolean flag to enable or disable saving output.
        fields: A list of strings specifying which variables to save.
    """
    save: bool
    fields: List[str]
