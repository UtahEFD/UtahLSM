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
"""UtahLSM Soil Physics Subpackage.

This module serves as the entry point for the soil physics component
of the land-surface model. It imports and exposes the main `Soil`
abstract base class, making it accessible to the rest of the model.
"""
from .soil import Soil

__all__ = ['Soil']
