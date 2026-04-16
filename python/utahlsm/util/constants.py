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
"""Defines physical and thermodynamic constants for UtahLSM.

This module provides a centralized location for all physical constants
used throughout the land-surface model. Constants are grouped into logical
namespaces to provide a clean, organized interface.
"""
from types import SimpleNamespace

import numpy as np

# Thermodynamic constants
thermodynamic = SimpleNamespace(
    GAS_CONSTANT_DRY=287.0, # Gas constant for dry air [J/kg-K]
    GAS_CONSTANT_VAPOR=461.4, # Gas constant for water vapor [J/kg-K]
    EPSILON=0.6220199393, # Ratio of dry/vapor gas constants
    EPSILON_VIRTUAL_TEMPERATURE=0.608, # Ratio of dry/vapor gas constants
    SPECIFIC_HEAT=1004.0, # Specific heat of air [J/kg-K]
    LATENT_HEAT_VAPORIZATION=2.45e6, # Latent heat of vaporization [J/kg]
    TETENS_A=17.269, # Tetens parameter A (dimensionless)
    TETENS_B=35.86, # Tetens parameter B [K]
    ES_REF=610.78, # Reference vapor pressure at 0°C [Pa]
)

# Numerical constants
numerical = SimpleNamespace(
    EPSILON=1e-9, # Small constant to prevent division by zero
)

# Fundamental physical constants
physical = SimpleNamespace(
    VON_KARMAN=0.41, # Von Karman constant []
    GRAVITY=9.81, # Gravitational acceleration [m/s^2]
    PI=np.pi, # Pi
)

# Water properties
water = SimpleNamespace(
    DENSITY=1000.0, # Density of water [kg/m^3]
    VOLUMETRIC_HEAT_CAPACITY=4.184e6, # Volumetric heat capacity of water [J/m^3-K]
)

# Air properties
air = SimpleNamespace(
    DENSITY_REF=1.204, # Reference density of air [kg/m^3]
    TEMPERATURE_REF=273.15, # Reference temperature [K]
)

# Radiation model parameters (Spencer 1971)
radiation = SimpleNamespace(
    DECLINATION_AMPLITUDE=23.45, # Solar declination [degrees]
    SOLSTICE_DAY=173, # Summer solstice reference day
    DAYS_PER_YEAR=365.25, # Days per year (leap year adjusted)
    STEFAN_BOLTZMANN=5.6697e-8, # Stefan-Boltzmann constant [W/m^2-K^4]
    SOLAR_CONSTANT=1367.0, # Solar constant [W/m^2]
)

# Soil physics model parameters
soil = SimpleNamespace(
    CONDUCTIVITY_PF_THRESHOLD=5.1, # Johansen (1975) pore fraction threshold
    CONDUCTIVITY_COEFF=418.46, # Empirical coefficient [W/m-K]
    CONDUCTIVITY_EXP=2.7, # Exponential parameter
    CONDUCTIVITY_MIN=0.172, # Minimum conductivity [W/m-K]
    PSI_FIELD_CAPACITY=-3.3, # Matric potential at field capacity [m head]
    PSI_WILTING_POINT=-150.0, # Matric potential at permanent wilting [m head]
)
