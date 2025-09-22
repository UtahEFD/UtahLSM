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
from dataclasses import dataclass

@dataclass(frozen=True)
class _Thermodynamic:
    """Constants related to thermodynamics."""
    GAS_CONSTANT_DRY: float         = 287.0        # Gas constant for dry air [J/kg-K]
    GAS_CONSTANT_VAPOR: float       = 461.4        # Gas constant for water vapor [J/kg-K]
    EPSILON: float                  = 0.6220199393 # Ratio of dry/vapor gas constants
    SPECIFIC_HEAT: float            = 1004.0       # Specific heat of air [J/kg-K]
    LATENT_HEAT_VAPORIZATION: float = 2.45e6       # Latent heat of vaporization [J/kg]

@dataclass(frozen=True)
class _Physical:
    """Fundamental physical constants."""
    VON_KARMAN: float       = 0.41             # Von Karman constant []
    GRAVITY: float          = 9.81             # Gravitational acceleration [m/s^2]
    PI: float               = 3.14159265358979 # Pi
    STEFAN_BOLTZMANN: float = 5.6697e-8        # Stefan-Boltzmann constant [W/m^2-K^4]
    SOLAR_CONSTANT: float   = 1367.0           # Solar constant [W/m^2]

@dataclass(frozen=True)
class _Water:
    """Constants related to water properties."""
    DENSITY: float       = 1000.0  # Density of water [kg/m^3]
    SPECIFIC_HEAT: float = 4.184e6 # Volumetric heat capacity of water [J/m^3-K]

@dataclass(frozen=True)
class _Air:
    """Constants related to air properties."""
    DENSITY: float = 1.204 # Density of air [kg/m^3]

air: _Air                     = _Air()
physical: _Physical           = _Physical()
thermodynamic: _Thermodynamic = _Thermodynamic()
water: _Water                 = _Water()
