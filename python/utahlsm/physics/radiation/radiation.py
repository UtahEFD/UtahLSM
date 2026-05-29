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

Radiation models are split into two halves:

* :meth:`compute_incoming` returns the downwelling shortwave and longwave
  fluxes. These are forcings on the surface and are evaluated once per
  timestep, either from a parameterization (clear-sky model) or from
  observed/coupled inputs.
* :meth:`compute_outgoing` returns the upwelling shortwave and longwave
  fluxes as a function of the *trial* surface temperature, surface
  optical properties, and the incoming components. It is evaluated every
  Brent iteration inside the surface energy budget so the
  4 * episilon * sigma * T_s^3 longwave restoring is correctly expressed.

The default :meth:`compute_outgoing` implements the standard one-source
model ``sw_out = alpha* sw_in`` and ``lw_out = epsilon * sigma * T_s^4 + (1 - epsilon) * lw_in``.
Subclasses need only override it for richer schemes (e.g. two-source
canopy/ground or two-stream multiple scattering).
"""

from abc import ABC, abstractmethod
from typing import TypeVar

import numpy as np
from numpy.typing import NDArray

from ..._types import FloatOrArray
from ...data_models import AtmosphericState, SurfaceState
from ...util import constants as c
from ...util.io import logging_helper

RT = TypeVar('RT', bound='Radiation')
logger = logging_helper.get_logger('Radiation')

class Radiation(ABC):
    """Abstract base class for radiation models.

    Concrete subclasses must implement :meth:`compute_incoming`. The
    default :meth:`compute_outgoing` covers the standard one-source case
    (single albedo and emissivity) and reads ``self.albedo`` and
    ``self.emissivity`` set on the subclass.
    """

    albedo: float
    emissivity: float

    # Abstract methods ---
    @abstractmethod
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
            sfc_state: Current state of the surface.

        Returns:
            Tuple ``(sw_in, lw_in)`` in W/m^2.
        """
        raise NotImplementedError

    def compute_outgoing(
        self,
        sfc_T: NDArray[np.float64],
        sw_in: NDArray[np.float64],
        lw_in: NDArray[np.float64],
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Computes the upwelling shortwave and longwave components.

        Default one-source implementation. Both incoming and outgoing
        components are treated as physical responses to ``sfc_T`` and
        the surface optical properties, so the SEB iterator can capture
        the longwave emission feedback during root-finding.

        ``sw_out = alpha * sw_in`` (single broadband albedo).

        ``lw_out = epsilon * sigma * T_s^4 + (1 - epsilon) * lw_in``
        (Stefan-Boltzmann emission plus reflection of the downwelling
        longwave). Keeping the reflected term makes the surface
        Kirchhoff-consistent: absorptivity equals emissivity, so the net
        longwave reduces to ``epsilon * (lw_in - sigma * T_s^4)`` rather
        than spuriously absorbing 100% of ``lw_in`` while emitting at
        ``epsilon < 1``.

        Args:
            sfc_T: Trial surface temperature [K].
            sw_in: Downwelling shortwave [W/m^2], same shape as sfc_T.
            lw_in: Downwelling longwave [W/m^2], same shape as sfc_T.
            atm_state: Current atmospheric state (unused by the default
                implementation; available for richer subclasses).
            sfc_state: Current surface state (unused by the default
                implementation; available for richer subclasses).

        Returns:
            Tuple ``(sw_out, lw_out)`` in W/m^2.
        """
        del atm_state, sfc_state  # unused in the one-source default
        SB = c.radiation.STEFAN_BOLTZMANN
        sw_out = self.albedo * np.asarray(sw_in)
        lw_out = (
            self.emissivity * SB * np.asarray(sfc_T) ** 4
            + (1.0 - self.emissivity) * np.asarray(lw_in)
        )
        return sw_out, lw_out
