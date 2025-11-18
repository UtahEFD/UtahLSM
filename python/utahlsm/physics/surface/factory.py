#
# UtahLSM
#
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#
"""Factory for creating surface model instances."""

from ...exceptions import NamelistError
from ...util.io import logging_helper
from .sfc import Surface
from .sfc_most import SurfaceMOST

logger = logging_helper.get_logger('SFC')


def get_surface_model(key: int) -> Surface:
    """Factory function to select and instantiate a surface layer model.

    Args:
        key: An integer identifying the surface layer model to use.

    Returns:
        An instance of a concrete `Surface` subclass.

    Raises:
        NamelistError: If the provided `key` is not a valid model ID.
    """
    # Dictionary to map keys to classes
    sfc_models = {
        1: SurfaceMOST,
    }

    try:
        return sfc_models[key]()
    except KeyError as e:
        error_msg = f'{key} is an invalid surface model.'
        logger.error('x' * 62)
        logger.error('Namelist Error: %s', error_msg)
        logger.error('Valid options are:')
        for k, v in sfc_models.items():
            logger.error('\t%d (%s)', k, v.__name__)
        logger.error('x' * 62)
        raise NamelistError(error_msg) from e
