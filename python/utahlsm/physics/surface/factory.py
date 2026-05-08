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
"""Factory for creating surface model instances."""

from ...data_models import SurfaceConfig
from ...exceptions import NamelistError
from ...util.io import logging_helper
from .sfc import Surface
from .sfc_most import SurfaceMOST

logger = logging_helper.get_logger('SFC')


def get_surface_model(surface: SurfaceConfig) -> Surface:
    """Factory function to select and instantiate a surface layer model.

    Args:
        surface: Surface layer configuration.

    Returns:
        An instance of a concrete `Surface` subclass.

    Raises:
        NamelistError: If the provided model name is not valid.
    """
    if surface.model == 'most':
        return SurfaceMOST(psi_stable=surface.psi_stable)

    valid = ['most']
    error_msg = f"'{surface.model}' is an invalid surface model."
    logger.error('x' * 62)
    logger.error('Namelist Error: %s', error_msg)
    logger.error('Valid options are: %s', ', '.join(valid))
    logger.error('x' * 62)
    raise NamelistError(error_msg)
