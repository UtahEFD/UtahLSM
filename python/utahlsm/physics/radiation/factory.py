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
"""Factory for creating radiation model instances."""

from ...data_models import RadiationConfig, SurfaceConfig
from ...exceptions import NamelistError
from ...util.io import logging_helper
from .rad_basic import RadBasic
from .rad_forcing import RadForcing
from .radiation import Radiation

logger = logging_helper.get_logger('RAD')


def get_radiation_model(
    radiation: RadiationConfig,
    surface: SurfaceConfig,
) -> Radiation:
    """Factory function to select and instantiate a radiation model.

    Args:
        radiation: Radiation configuration from the namelist.
        surface: Surface configuration supplying albedo and emissivity.

    Returns:
        An instance of a concrete :class:`Radiation` subclass.

    Raises:
        NamelistError: If ``radiation.model`` is not a recognised option.
    """
    if radiation.model == 'forcing':
        return RadForcing(
            surface.albedo,
            surface.emissivity,
        )

    if radiation.model == 'basic':
        return RadBasic(
            radiation.latitude,
            radiation.longitude,
            surface.albedo,
            surface.emissivity,
        )

    valid = ['forcing', 'basic']
    error_msg = f"'{radiation.model}' is an invalid radiation model."
    logger.error('x' * 62)
    logger.error('Namelist Error: %s', error_msg)
    logger.error('Valid options are: %s', ', '.join(valid))
    logger.error('x' * 62)
    raise NamelistError(error_msg)
