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
"""Forcing-passthrough radiation model for UtahLSM."""

from ..._types import FloatOrArray
from ...data_models import AtmosphericState, SurfaceState
from ...util.io import logging_helper
from .radiation import Radiation


class RadForcing(Radiation):
    """Radiation model that passes incoming components through from forcing.

    Used when downwelling fluxes are provided externally (observation-based
    forcing or a coupled atmospheric model). Outgoing components are still
    treated as a response to surface state and are produced by the base
    class :meth:`compute_outgoing` from the trial T_s, the incoming
    fluxes, and the surface albedo/emissivity. Forcing files therefore
    only need to supply ``sw_in`` and ``lw_in``.
    """

    def __init__(self, albedo: float, emissivity: float) -> None:
        """Initializes the RadForcing model.

        Args:
            albedo: The surface albedo (dimensionless), used for sw_out.
            emissivity: The surface emissivity (dimensionless), used for
                lw_out.
        """
        logger = logging_helper.get_logger('Radiation')
        logger.info('Using forcing radiation data')
        self.albedo: float = albedo
        self.emissivity: float = emissivity

    def compute_incoming(
        self,
        julian_day: int,
        time_utc: float,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
    ) -> tuple[FloatOrArray, FloatOrArray]:
        """Returns the downwelling components from atmospheric forcing.

        Args:
            julian_day: Unused; present for interface compatibility.
            time_utc: Unused; present for interface compatibility.
            atm_state: Current atmospheric state carrying ``sw_in`` and
                ``lw_in`` from the forcing file.
            sfc_state: Unused; present for interface compatibility.

        Returns:
            Tuple ``(sw_in, lw_in)`` taken directly from ``atm_state``.
        """
        del julian_day, time_utc, sfc_state
        return atm_state.sw_in, atm_state.lw_in
