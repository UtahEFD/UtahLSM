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
"""
data_models.py

This module defines the core, shared data structures (using dataclasses) for the UtahLSM model.
These classes provide a structured blueprint for configuration, initial conditions, and forcing data,
ensuring consistency across different parts of the model.
"""

from dataclasses import dataclass
import numpy as np
from typing import List

# --- Core Physical State Dataclasses ---

@dataclass
class AtmosphericState:
    """Represents the atmospheric conditions at a single point in time."""
    U: float      # Wind speed [m/s]
    T: float      # Air temperature [K]
    q: float      # Specific humidity [kg/kg]
    p: float      # Air pressure [Pa]
    R_net: float  # Net radiation [W/m^2]

@dataclass
class SoilState:
    """Represents the state of the soil column."""
    T: np.ndarray     # soil temperature [K]
    q: np.ndarray     # soil moisture [g/g]
    type: np.ndarray  # soil type [category]

@dataclass
class SurfaceState:
    """Represents the surface conditions at a single point in time."""
    Ts:  float      # surface skin temperature [K]
    qs:  float      # surface skin water content [g/g]
    qa:  float      # surface skin mixing ratio [g/g]
    ust: np.ndarray # friction velocity [m/s]
    obl: np.ndarray # obukhov length [m]
    wT:  np.ndarray # kinematic heat flux [K m/s]
    wq:  np.ndarray # kinematic moisture flux [m/s]
    shf: np.ndarray # sensible heat flux [W/m^2]
    lhf: np.ndarray # latent heat flux [W/m^2]
    ghf: np.ndarray # ground heat flux [W/m^2]

@dataclass
class SolverState:
    """Holds temporary variables for use in numerical solvers."""
    Kmid: float = 0.0  # soil thermal conductivity for SEB[W/(m*K)]

@dataclass(frozen=True)
class ForcingData:
    """Represents the entire time-series of meteorological forcing data."""
    ntime: int
    tstep: float
    atmos: List[AtmosphericState] # A list of atmospheric states, one for each timestep

# --- Configuration Dataclasses (from Namelist) ---

@dataclass(frozen=True)
class GeneralConfig:
    log_level: str
    
@dataclass(frozen=True)
class NumericsConfig:
    diffusion_back_weight: float

@dataclass(frozen=True)
class TimeConfig:
    step_seb: int
    step_dif: int
    utc_start: int
    julian_day: int

@dataclass(frozen=True)
class GridConfig:
    nx: int
    ny: int
    nz: int
    z : np.ndarray

@dataclass(frozen=True)
class SurfaceConfig:
    z_o: float
    z_t: float
    z_m: float
    z_s: float
    albedo: float
    emissivity: float
    model: int
    flux_iter_max: int
    flux_criteria: float
    temperature_reference: float

@dataclass(frozen=True)
class SoilConfig:
    param: int
    model: int

@dataclass(frozen=True)
class RadiationConfig:
    model: int
    latitude: float
    longitude: float

@dataclass(frozen=True)
class OutputConfig:
    save: bool
    fields: List[str]
