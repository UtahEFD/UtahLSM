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

"""UtahLSM Physics Subpackage.

This package provides the core physics modules for the land-surface model.
It contains the abstract base classes for radiation, soil, and surface
processes, which define the common interface for different physics
parameterizations.

By importing the base classes here, they are made directly accessible under
the `utahlsm.physics` namespace.
"""

from .radiation import Radiation
from .soil import Soil
from .surface import Surface

__all__ = ['Radiation', 'Soil', 'Surface']
