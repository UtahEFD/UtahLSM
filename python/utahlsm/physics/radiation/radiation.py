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
"""Abstract base class for radiation models in UtahLSM.

This module defines the interface for all radiation physics parameterizations
used within the UtahLSM framework. It provides the `Radiation` abstract base
class, which ensures that any concrete radiation model implements the necessary
methods, and a factory function (`get_model`) for creating instances of those
models.
"""

from abc import ABC, abstractmethod
from typing import TypeVar

from ..._types import FloatOrArray
from ...data_models import AtmosphericState, SurfaceState
from ...util.io import logging_helper

RT = TypeVar('RT', bound='Radiation')
logger = logging_helper.get_logger('Radiation')

class Radiation(ABC):
    """Abstract base class for radiation models.

    Concrete subclasses must implement :meth:`compute_components`.
    """

    # Abstract methods ---
    @abstractmethod
    def compute_components(
        self,
        julian_day: int,
        time_utc: float,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
    ) -> tuple[FloatOrArray, FloatOrArray,
               FloatOrArray, FloatOrArray]:
        """Computes the four surface radiation components.

        Concrete subclasses return the downwelling and upwelling shortwave
        and longwave fluxes at the surface; the caller derives net
        radiation as ``sw_in - sw_out + lw_in - lw_out``.

        Args:
            julian_day: Current Julian day of the year.
            time_utc: Current time in UTC seconds from midnight.
            atm_state: Current state of the atmosphere.
            sfc_state: Current state of the surface.

        Returns:
            Tuple ``(sw_in, sw_out, lw_in, lw_out)`` in W/m^2, each
            shaped consistently with the column layout.
        """
        raise NotImplementedError
