#!/usr/bin/env python
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

"""Thermodynamic helpers shared across the physics modules.

Centralizes pure-function physical relationships (saturation vapor
pressure, saturation specific humidity) so that the surface, soil, and
canopy modules all reference a single implementation.
"""

from typing import Union

import numpy as np
from numpy.typing import NDArray

from ..util import constants as c

Number = Union[float, NDArray[np.float64]]


def saturation_vapor_pressure(T: Number) -> Number:
    """Saturation vapor pressure over liquid water via the Tetens formula.

    Args:
        T: Temperature [K].

    Returns:
        Saturation vapor pressure [Pa], same shape as input.
    """
    return c.thermodynamic.ES_REF * np.exp(
        c.thermodynamic.TETENS_A * (T - c.air.TEMPERATURE_REF)
        / (T - c.thermodynamic.TETENS_B)
    )


def saturation_specific_humidity(T: Number, p: Number) -> Number:
    """Saturation specific humidity q_sat [kg/kg].

    Args:
        T: Temperature [K].
        p: Air pressure [Pa].

    Returns:
        Saturation specific humidity [kg/kg], broadcast-shaped.
    """
    es = saturation_vapor_pressure(T)
    return c.thermodynamic.EPSILON * es / (p - 0.378 * es)
