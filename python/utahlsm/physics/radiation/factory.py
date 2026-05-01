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

from ...exceptions import NamelistError
from ...util.io import logging_helper
from .rad_basic import RadBasic
from .radiation import Radiation

logger = logging_helper.get_logger('RAD')


def get_radiation_model(key: int, latitude: float, longitude: float,
                        albedo: float, emissivity: float) -> Radiation:
    """Factory function to select and instantiate a radiation model.

    Based on the integer key provided in the namelist, this function creates
    and returns an instance of the corresponding radiation model class.

    Args:
        key: An integer identifying the radiation model to use.
        latitude: The site latitude in degrees.
        longitude: The site longitude in degrees.
        albedo: The surface albedo (dimensionless).
        emissivity: The surface emissivity (dimensionless).

    Returns:
        An instance of a concrete `Radiation` subclass.

    Raises:
        NamelistError: If the provided `key` is not a valid model ID.
    """
    # Dictionary to map keys to classes
    rad_models = {
        1: RadBasic,
    }

    try:
        return rad_models[key](latitude, longitude, albedo, emissivity)
    except KeyError as e:
        error_msg = f'{key} is an invalid radiation model.'
        logger.error('x' * 62)
        logger.error('Namelist Error: %s', error_msg)
        logger.error('Valid options are:')
        for k, v in rad_models.items():
            logger.error('\t%d (%s)', k, v.__name__)
        logger.error('x' * 62)
        raise NamelistError(error_msg) from e
