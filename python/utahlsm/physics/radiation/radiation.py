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

import logging
from abc import ABC, abstractmethod
from typing import Tuple, TypeVar

import numpy as np
from numpy.typing import NDArray

from ...data_models import AtmosphericState, SurfaceState
from ...util.io import logging_helper

RT = TypeVar('RT', bound='Radiation')
logger = logging_helper.get_logger('RAD')

class Radiation(ABC):
    """Abstract base class for radiation models.

    This class defines the standard interface for radiation calculations
    and acts as a factory for creating specific radiation model instances.
    It cannot be instantiated directly.

    Attributes:
        logger: A logger for this class.
        latitude: The site latitude in degrees.
        longitude: The site longitude in degrees.
        albedo: The surface albedo (dimensionless).
        emissivity: The surface emissivity (dimensionless).
    """
    def __init__(self, latitude: float, longitude: float, albedo: float,
                 emissivity: float):
        """Initializes the Radiation base class.

        Args:
            latitude: The site latitude in degrees.
            longitude: The site longitude in degrees.
            albedo: The surface albedo (dimensionless).
            emissivity: The surface emissivity (dimensionless).
        """
        self.logger: logging.Logger = logging_helper.get_logger('RAD')
        # Store as radians for internal trig usage.
        self.latitude: float = float(np.deg2rad(latitude))
        self.longitude: float = float(np.deg2rad(longitude))
        self.albedo: float = albedo
        self.emissivity: float = emissivity

    # Abstract methods ---
    @abstractmethod
    def compute_components(
        self,
        julian_day: int,
        time_utc: float,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64],
               NDArray[np.float64], NDArray[np.float64]]:
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
