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

    def validate_col(
        name: str,
        values: NDArray[np.float64],
        *,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
        min_inclusive: bool = True,
        max_inclusive: bool = True,
    ) -> None:
        arr = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(arr)):
            raise NamelistError(
                f"Canopy parameter '{name}' must contain only finite values."
            )
        if min_value is not None:
            bad = arr < min_value if min_inclusive else arr <= min_value
            if np.any(bad):
                op = '>=' if min_inclusive else '>'
                bad_value = float(np.min(arr[bad]))
                raise NamelistError(
                    f"Canopy parameter '{name}' must be {op} {min_value}; "
                    f"got {bad_value}."
                )
        if max_value is not None:
            bad = arr > max_value if max_inclusive else arr >= max_value
            if np.any(bad):
                op = '<=' if max_inclusive else '<'
                bad_value = float(np.max(arr[bad]))
                raise NamelistError(
                    f"Canopy parameter '{name}' must be {op} {max_value}; "
                    f"got {bad_value}."
                )

    common = {
        'lai': as_col('lai', config.lai),
        'veg_fraction': as_col('veg_fraction', config.veg_fraction),
        'rooting_depth': as_col('rooting_depth', config.rooting_depth),
        'beta': as_col('beta', config.beta),
        'rs_min': as_col('rs_min', config.rs_min),
        'rs_max': as_col('rs_max', config.rs_max),
        'r_ground': as_col('r_ground', config.r_ground),
        'z': z,
    }
    validate_col('lai', common['lai'], min_value=0.0)
    validate_col(
        'veg_fraction', common['veg_fraction'], min_value=0.0, max_value=1.0
    )
    validate_col('rooting_depth', common['rooting_depth'], min_value=0.0)
    validate_col(
        'beta',
        common['beta'],
        min_value=0.0,
        max_value=1.0,
        min_inclusive=False,
        max_inclusive=False,
    )
    validate_col('rs_min', common['rs_min'], min_value=0.0, min_inclusive=False)
    validate_col('rs_max', common['rs_max'], min_value=0.0, min_inclusive=False)
    validate_col('r_ground', common['r_ground'], min_value=0.0)
    if np.any(common['rs_max'] < common['rs_min']):
        raise NamelistError(
            "Canopy parameter 'rs_max' must be >= 'rs_min' in every column."
        )

    if config.model == 'jarvis':
        rg_half = as_col('rg_half', config.rg_half)
        vpd_coef = as_col('vpd_coef', config.vpd_coef)
        t_opt = as_col('t_opt', config.t_opt)
        t_coef = as_col('t_coef', config.t_coef)
        validate_col('rg_half', rg_half, min_value=0.0, min_inclusive=False)
        validate_col('vpd_coef', vpd_coef, min_value=0.0)
        validate_col('t_opt', t_opt, min_value=0.0, min_inclusive=False)
        validate_col('t_coef', t_coef, min_value=0.0)
        return CanopyJarvis(
            rg_half=rg_half,
            vpd_coef=vpd_coef,
            t_opt=t_opt,
            t_coef=t_coef,
            **common,
        )

    valid = ['none', 'jarvis']
    logger.error('x' * 62)
    logger.error("Namelist Error: '%s' is not a valid canopy model.",
                 config.model)
    logger.error('Valid options are: %s', ', '.join(valid))
    logger.error('x' * 62)
    raise NamelistError(f"Invalid canopy model '{config.model}'.")
