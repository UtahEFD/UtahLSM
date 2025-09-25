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

"""Utah Land-Surface Model (UtahLSM).

A fast-response land-surface model designed for simulating the exchange
of heat, moisture, and momentum between the land surface and the atmosphere.

Developed collaboratively at the University of Utah and the NOAA National
Severe Storms Laboratory, UtahLSM is engineered for flexibility, allowing
for various physics parameterizations for the surface, soil, and radiation.

Main Features:
    - Modular physics components for easy extension and testing.
    - Support for different soil, radiation, and surface layer schemes.
    - NetCDF input/output for compatibility with standard atmospheric data formats.
"""

# --- Top-Level Imports ---
# These make classes and data models directly accessible when
# users import 'utahlsm'.

# The main land-surface model class
from .core import UtahLSM

# Core data structures for atmospheric, surface, soil, and solver states
#    as well as configuration settings from the namelist file
from .data_models import (
    AtmosphericState, SoilState, SurfaceState, SolverState, ForcingData,
    GeneralConfig, NumericsConfig, TimeConfig, GridConfig, SurfaceConfig,
    SoilConfig, RadiationConfig, OutputConfig
)

# Core I/O classes for handling model input and output
from .util.io import Input, Output
