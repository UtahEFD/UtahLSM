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
"""Defines physical and thermodynamic constants for UtahLSM.

This module provides a centralized location for all physical constants
used throughout the land-surface model. Constants are grouped into logical
namespaces to provide a clean, organized interface.
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class _Thermodynamic:
    GAS_CONSTANT_DRY: float = 287.0        # Gas constant for dry air [J/kg-K]
    GAS_CONSTANT_VAPOR: float = 461.4      # Gas constant for water vapor [J/kg-K]
    EPSILON: float = 0.6220199393          # Ratio of dry/vapor gas constants
    EPSILON_VIRTUAL_TEMPERATURE: float = 0.608  # Virtual temperature correction
    SPECIFIC_HEAT: float = 1004.0          # Specific heat of air [J/kg-K]
    LATENT_HEAT_VAPORIZATION: float = 2.45e6   # Latent heat of vaporization [J/kg]
    TETENS_A: float = 17.269               # Tetens parameter A (dimensionless)
    TETENS_B: float = 35.86                # Tetens parameter B [K]
    ES_REF: float = 610.78                 # Reference vapor pressure at 0°C [Pa]


@dataclass(frozen=True)
class _Numerical:
    EPSILON: float = 1e-9                  # Small constant to prevent division by zero


@dataclass(frozen=True)
class _Physical:
    VON_KARMAN: float = 0.41              # Von Karman constant []
    GRAVITY: float = 9.81                 # Gravitational acceleration [m/s^2]
    PI: float = float(np.pi)              # Pi


@dataclass(frozen=True)
class _Water:
    DENSITY: float = 1000.0               # Density of water [kg/m^3]
    VOLUMETRIC_HEAT_CAPACITY: float = 4.184e6  # Volumetric heat capacity [J/m^3-K]


@dataclass(frozen=True)
class _Air:
    DENSITY_REF: float = 1.204            # Reference density of air [kg/m^3]
    TEMPERATURE_REF: float = 273.15       # Reference temperature [K]


@dataclass(frozen=True)
class _Radiation:
    DECLINATION_AMPLITUDE: float = 23.45  # Solar declination [degrees]
    SOLSTICE_DAY: int = 173               # Summer solstice reference day
    DAYS_PER_YEAR: float = 365.25         # Days per year (leap year adjusted)
    STEFAN_BOLTZMANN: float = 5.6697e-8   # Stefan-Boltzmann constant [W/m^2-K^4]
    SOLAR_CONSTANT: float = 1367.0        # Solar constant [W/m^2]


@dataclass(frozen=True)
class _Soil:
    CONDUCTIVITY_PF_THRESHOLD: float = 5.1    # Johansen (1975) pore fraction threshold
    CONDUCTIVITY_COEFF: float = 418.46        # Empirical coefficient [W/m-K]
    CONDUCTIVITY_EXP: float = 2.7             # Exponential parameter
    CONDUCTIVITY_MIN: float = 0.172           # Minimum conductivity [W/m-K]
    PSI_FIELD_CAPACITY: float = -3.3          # Matric potential at field capacity [m head]
    PSI_WILTING_POINT: float = -150.0         # Matric potential at permanent wilting [m head]


thermodynamic = _Thermodynamic()
numerical = _Numerical()
physical = _Physical()
water = _Water()
air = _Air()
radiation = _Radiation()
soil = _Soil()
