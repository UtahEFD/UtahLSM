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
from typing import List, Optional

# --- Core Physical State Dataclasses ---

@dataclass
class AtmosphericData:
    """Represents the atmospheric conditions at a single point in time."""
    U: float      # Wind speed [m/s]
    T: float      # Air temperature [K]
    q: float      # Specific humidity [kg/kg]
    p: float      # Air pressure [Pa]
    R_net: float  # Net radiation [W/m^2]

@dataclass
class SoilData:
    """Represents the state of the soil column."""
    z: np.ndarray     # vertical grid levels distance from surface [m]
    T: np.ndarray     # soil temperature [m/s]
    q: np.ndarray     # soil moisture [g/g]
    type: np.ndarray  # soil type [category]

@dataclass
class ForcingData:
    """Represents the entire time-series of meteorological forcing data."""
    ntime: int
    tstep: float
    atmos: List[AtmosphericData] # A list of atmospheric states, one for each timestep

# --- Configuration Dataclasses (from Namelist) ---

@dataclass
class TimeConfig:
    step_seb: int
    step_dif: int
    utc_start: int
    julian_day: int

@dataclass
class GridConfig:
    nx: int
    ny: int
    nz: int
    z : np.ndarray

@dataclass
class SurfaceConfig:
    z_o: float
    z_t: float
    z_m: float
    z_s: float
    albedo: float
    emissivity: float
    model: int

@dataclass
class SoilConfig:
    param: int
    model: int

@dataclass
class RadiationConfig:
    model: int
    latitude: float
    longitude: float

@dataclass
class OutputConfig:
    save: bool
    fields: List[str]
