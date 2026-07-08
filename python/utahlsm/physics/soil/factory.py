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
"""Factory for creating soil model instances."""


from typing import Any

from ...exceptions import NamelistError
from ...util.io import logging_helper
from .soil import Soil
from .soil_brookscorey import BrooksCorey
from .soil_campbell import Campbell
from .soil_vangenuchten import VanGenuchten

logger = logging_helper.get_logger('Soil')

def get_soil_model(
    key: str,
    properties_dict: dict[str, dict[str, Any]],
    soil_type_names: list[str],
    dataset_name: str = 'custom',
    thermal_conductivity_model: str = 'mccumber-pielke',
) -> Soil:
    """Factory function to select and instantiate a soil model.

    Args:
        key: Model name ('brooks-corey', 'campbell', or 'van-genuchten').
        properties_dict: Dictionary mapping soil type names to property dicts.
        soil_type_names: List of soil type names for each layer.
        dataset_name: Human-readable name of the dataset being used.
        thermal_conductivity_model: Thermal conductivity parameterization
            ('mccumber-pielke' or 'johansen').

    Returns:
        An instance of a concrete `Soil` subclass.

    Raises:
        NamelistError: If the provided `key` is not a valid model name.
    """
    soil_models: dict[str, type[Soil]] = {
        'brooks-corey': BrooksCorey,
        'campbell': Campbell,
        'van-genuchten': VanGenuchten,
    }

    # Soil types defined with native van Genuchten retention parameters
    # (alpha/n) carry only a derived 'b' (= 1/(n-1)), not a measured
    # Campbell pore-size index, so they are valid only for the
    # van-genuchten model.
    if key in ('brooks-corey', 'campbell'):
        vg_native = sorted({
            name.lower() for name in soil_type_names
            if properties_dict.get(name.lower(), {}).get('parameterization')
            == 'van-genuchten'
        })
        if vg_native:
            error_msg = (
                f"Soil model '{key}' cannot use soil types defined with "
                f"native van Genuchten parameters (alpha/n): "
                f"{', '.join(vg_native)}. Use the 'van-genuchten' model or "
                f"a dataset providing measured 'b'/'psi_sat' for these types."
            )
            logger.error('x' * 62)
            logger.error('Namelist Error: %s', error_msg)
            logger.error('x' * 62)
            raise NamelistError(error_msg)

    try:
        return soil_models[key](
            properties_dict, soil_type_names, dataset_name,
            thermal_conductivity_model
        )
    except KeyError as e:
        error_msg = f"'{key}' is an invalid soil model."
        logger.error('x' * 62)
        logger.error('Namelist Error: %s', error_msg)
        logger.error('Valid options are: %s', ', '.join(soil_models))
        logger.error('x' * 62)
        raise NamelistError(error_msg) from e
