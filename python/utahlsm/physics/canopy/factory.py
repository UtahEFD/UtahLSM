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
"""Factory for creating canopy (vegetation) model instances."""

from typing import Optional

import numpy as np
from numpy.typing import NDArray

from ...data_models import CanopyConfig
from ...exceptions import NamelistError
from ...util.io import logging_helper
from .canopy import Canopy
from .canopy_jarvis import CanopyJarvis

logger = logging_helper.get_logger('CANOPY')


def get_canopy_model(
    config: CanopyConfig,
    z: NDArray[np.float64],
    ncol: int,
) -> Optional[Canopy]:
    """Selects and instantiates a canopy model from namelist config.

    Args:
        config: Parsed `CanopyConfig` dataclass from the namelist.
        z: Soil-layer node depths [m], shape (nz,).
        ncol: Number of horizontal columns.

    Returns:
        A concrete :class:`Canopy` subclass, or ``None`` when
        ``config.model == 'none'`` (bare-soil mode).

    Raises:
        NamelistError: If ``config.model`` is not a recognised option.
    """
    if config is None or config.model == 'none':
        return None

    # Helper to coerce a namelist scalar or sequence into (ncol,).
    def as_col(name: str, value) -> NDArray[np.float64]:
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            return np.full(ncol, float(arr))
        if arr.size == ncol:
            return arr.reshape(ncol)
        if arr.size == 1:
            return np.full(ncol, float(arr.flat[0]))
        raise NamelistError(
            f"Canopy parameter '{name}' has size {arr.size}; expected "
            f"1 or ncol={ncol}."
        )

    common = dict(
        lai=as_col('lai', config.lai),
        veg_fraction=as_col('veg_fraction', config.veg_fraction),
        rooting_depth=as_col('rooting_depth', config.rooting_depth),
        beta=as_col('beta', config.beta),
        rs_min=as_col('rs_min', config.rs_min),
        rs_max=as_col('rs_max', config.rs_max),
        z=z,
    )

    if config.model == 'jarvis':
        return CanopyJarvis(
            rg_half=as_col('rg_half', config.rg_half),
            vpd_coef=as_col('vpd_coef', config.vpd_coef),
            t_opt=as_col('t_opt', config.t_opt),
            t_coef=as_col('t_coef', config.t_coef),
            **common,
        )

    valid = ['none', 'jarvis']
    logger.error('x' * 62)
    logger.error("Namelist Error: '%s' is not a valid canopy model.",
                 config.model)
    logger.error('Valid options are: %s', ', '.join(valid))
    logger.error('x' * 62)
    raise NamelistError(f"Invalid canopy model '{config.model}'.")
