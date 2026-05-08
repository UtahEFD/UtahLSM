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
    """Radiation model that passes through components from atmospheric forcing.

    Used when radiation fields are provided externally (e.g. from a coupled
    model or observation-based forcing file) rather than computed online.
    ``compute_components`` simply returns the four component fields already
    loaded onto ``atm_state``.
    """

    def __init__(self) -> None:
        """Initializes the RadForcing model."""
        logger = logging_helper.get_logger('Radiation')
        logger.info('Using forcing radiation data')

    def compute_components(
        self,
        julian_day: int,
        time_utc: float,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
    ) -> tuple[FloatOrArray, FloatOrArray, FloatOrArray, FloatOrArray]:
        """Returns the four radiation components from atmospheric forcing.

        Args:
            julian_day: Unused; present for interface compatibility.
            time_utc: Unused; present for interface compatibility.
            atm_state: Current atmospheric state carrying sw_in, sw_out,
                lw_in, lw_out from the forcing file.
            sfc_state: Unused; present for interface compatibility.

        Returns:
            Tuple ``(sw_in, sw_out, lw_in, lw_out)`` taken directly from
            ``atm_state``.
        """
        return (
            atm_state.sw_in,
            atm_state.sw_out,
            atm_state.lw_in,
            atm_state.lw_out,
        )
