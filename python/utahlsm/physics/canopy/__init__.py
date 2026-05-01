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
"""UtahLSM Canopy Physics Subpackage.

Exposes the `Canopy` abstract base class used by all vegetation
parameterizations (Jarvis and future variants). Concrete models are
constructed through `factory.get_canopy_model`.
"""
from .canopy import Canopy

__all__ = ['Canopy']
