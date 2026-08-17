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

This module provides a simple, clear-sky radiation parameterization for the
incoming components. The outgoing components are produced by the default
:meth:`Radiation.compute_outgoing` from the trial surface temperature, the
incoming fluxes, and the surface albedo/emissivity.
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

    Provides parameterized downwelling shortwave (clear-sky geometric
    optics) and longwave (Brutsaert 1975 effective emissivity).
    Outgoing components fall through to the base class.
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
        self.logger: logging.Logger = logging_helper.get_logger('Radiation')
        self.logger.info('Using the basic model')
        # Store as radians for internal trig usage.
        self.latitude: float = float(np.deg2rad(latitude))
        self.longitude: float = float(np.deg2rad(longitude))
        self.albedo: float = albedo
        self.emissivity: float = emissivity

    def compute_incoming(
        self,
        julian_day: int,
        time_utc: float,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
    ) -> tuple[FloatOrArray, FloatOrArray]:
        """Computes the downwelling shortwave and longwave components.

        Args:
            julian_day: Current Julian day of the year.
            time_utc: Current time in UTC seconds from midnight.
            atm_state: Current state of the atmosphere.
            sfc_state: Current state of the surface (unused).

        Returns:
            Tuple ``(sw_in, lw_in)`` in W/m^2.
        """
        del sfc_state
        sw_in = self._shortwave_in(julian_day, time_utc)
        lw_in = self._longwave_in(atm_state)
        return sw_in, lw_in

    def _shortwave_in(self, julian_day: int, time_utc: float) -> FloatOrArray:
        """Computes downward shortwave radiation for clear-sky conditions.

        Longitude is east-positive, so a site further east reaches solar
        noon earlier in UTC: ``t_noon = 43200 - lon * 86400 / (2 * pi)``
        seconds. That is produced by the ``+ self.longitude`` term inside
        the cosine below; the leading minus on the whole cosine term is
        what places the hour-angle origin at noon rather than midnight.
        Using ``- self.longitude`` instead shifts the diurnal cycle by
        ``2 * lon`` and is correct only on the prime meridian.

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
                         np.cos((2*PI*time_utc/(24.0*3600.0))+
                                self.longitude))
        transmissivity = 0.6 + 0.2*sin_elevation
        sw_in = np.where(
            sin_elevation > 0,
            SC * transmissivity * sin_elevation,
            0.0
        )
        return sw_in

    def _longwave_in(self, atm_state: AtmosphericState) -> FloatOrArray:
        """Computes clear-sky downwelling longwave radiation.

        Uses the Brutsaert (1975) effective emissivity relation
        ``eps = 1.24 * (e / T)^(1/7)``. The 1.24 coefficient is calibrated
        for vapor pressure in hPa (mb), so the vapor pressure is converted
        from Pa before the exponent is applied.

        Args:
            atm_state: The current state of the atmosphere.

        Returns:
            The incoming longwave radiation in W/m^2.
        """
        # local constants
        EPSILON = c.thermodynamic.EPSILON
        SB = c.radiation.STEFAN_BOLTZMANN
        PA_PER_HPA = 100.0

        # local references to atmospheric state
        pa = atm_state.pressure
        qa = atm_state.specific_humidity
        Ta = atm_state.temperature

        # Vapor pressure [Pa], then converted to hPa for Brutsaert (1975).
        # Applying the relation to Pa inflates the emissivity by a factor
        # 100^(1/7) = 1.93, which pushes it above unity and makes lw_in
        # exceed blackbody emission at the same temperature.
        vapor_pressure = (pa * qa) / (EPSILON + qa)
        emissivity_eff = 1.24 * (
            vapor_pressure / PA_PER_HPA / Ta
        ) ** (1 / 7.0)

        return emissivity_eff * SB * (Ta ** 4)
