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
"""A basic radiation model for UtahLSM.

This module provides a simple, clear-sky radiation parameterization. It includes
methods for calculating incoming and outgoing shortwave and longwave radiation
based on fundamental physical principles and empirical relationships.
"""
import logging

import numpy as np

from ..._types import FloatOrArray
from ...data_models import AtmosphericState, SurfaceState
from ...util import constants as c
from ...util.io import logging_helper
from .radiation import Radiation


class RadBasic(Radiation):
    """A basic clear-sky radiation model.

    This class implements the `Radiation` interface and provides simple
    calculations for the four components of the surface radiation budget.
    """

    def __init__(
        self, latitude: float, longitude: float, albedo: float,
        emissivity: float
    ):
        """Initializes the RadBasic model.

        Args:
            latitude: The site latitude in degrees.
            longitude: The site longitude in degrees.
            albedo: The surface albedo (dimensionless).
            emissivity: The surface emissivity (dimensionless).
        """
        self.logger: logging.Logger = logging_helper.get_logger('RAD: Basic')
        self.logger.info('Using the basic model')
        super().__init__(latitude, longitude, albedo, emissivity)

    def compute_components(
        self,
        julian_day: int,
        time_utc: float,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
    ) -> tuple[FloatOrArray, FloatOrArray,
               FloatOrArray, FloatOrArray]:
        """Computes the four surface radiation components.

        Args:
            julian_day: Current Julian day of the year.
            time_utc: Current time in UTC seconds from midnight.
            atm_state: Current state of the atmosphere.
            sfc_state: Current state of the surface.

        Returns:
            Tuple ``(sw_in, sw_out, lw_in, lw_out)`` in W/m^2.
        """
        sw_in = self._shortwave_in(julian_day, time_utc)
        sw_out = self._shortwave_out(sw_in)
        lw_in = self._longwave_in(atm_state, sfc_state)
        lw_out = self._longwave_out(sfc_state)
        return sw_in, sw_out, lw_in, lw_out

    def _shortwave_in(self, julian_day: int, time_utc: float) -> FloatOrArray:
        """Computes downward shortwave radiation for clear-sky conditions.

        Args:
            julian_day: The current Julian day of the year.
            time_utc: The current time in UTC seconds from midnight.

        Returns:
            The incoming shortwave radiation in W/m^2, or 0 if the sun is down.
        """
        # local constants
        PI = c.physical.PI
        SC = c.radiation.SOLAR_CONSTANT

        declination = c.radiation.DECLINATION_AMPLITUDE * (PI/180.0) * \
            np.cos(2.0*PI*(julian_day-c.radiation.SOLSTICE_DAY) /
                   c.radiation.DAYS_PER_YEAR)
        sin_elevation = (np.sin(self.latitude)*np.sin(declination) -
                         np.cos(self.latitude)*np.cos(declination) *
                         np.cos((2*PI*time_utc/(24.0*3600.0))-
                                self.longitude))
        transmissivity = 0.6 + 0.2*sin_elevation
        sw_in = np.where(
            sin_elevation > 0,
            SC * transmissivity * sin_elevation,
            0.0
        )
        return sw_in

    def _shortwave_out(self, sw_in: FloatOrArray) -> FloatOrArray:
        """Computes upward shortwave radiation based on surface albedo.

        Args:
            sw_in: The incoming shortwave radiation in W/m^2.

        Returns:
            The outgoing shortwave radiation in W/m^2.
        """
        return self.albedo*sw_in

    def _longwave_in(
        self, atm_state: AtmosphericState, sfc_state: SurfaceState
    ) -> FloatOrArray:
        """Computes clear-sky downwelling longwave radiation.

        Uses the Brutsaert (1975) emissivity relation.

        Args:
            atm_state: The current state of the atmosphere.
            sfc_state: The current state of the surface.

        Returns:
            The incoming longwave radiation in W/m^2.
        """
        # local constants
        EPSILON = c.thermodynamic.EPSILON
        SB = c.radiation.STEFAN_BOLTZMANN

        # local references to atmospheric state
        pa = atm_state.pressure
        qa = atm_state.specific_humidity
        Ta = atm_state.temperature

        # vapor pressure and effective emissivity (Brutsaert 1975)
        vapor_pressure = (pa * qa) / (EPSILON + qa)
        emissivity_eff = 1.24 * (vapor_pressure / Ta) ** (1 / 7.0)

        return emissivity_eff * SB * (Ta ** 4)

    def _longwave_out(self, sfc_state: SurfaceState) -> FloatOrArray:
        """Computes upward longwave radiation using the Stefan-Boltzmann law.

        Args:
            sfc_state: The current state of the surface.

        Returns:
            The outgoing longwave radiation in W/m^2.
        """
        # local constants
        SB = c.radiation.STEFAN_BOLTZMANN

        # local references to surface state
        Ts = sfc_state.temperature

        return self.emissivity * SB * (Ts ** 4)
