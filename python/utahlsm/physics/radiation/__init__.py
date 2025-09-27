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
"""UtahLSM Radiation Physics Subpackage.

This module serves as the entry point for the radiation physics component
of the land-surface model. It imports and exposes the main `Radiation`
abstract base class, making it accessible to the rest of the model.
"""
from .radiation import Radiation

__all__ = ['Radiation']
