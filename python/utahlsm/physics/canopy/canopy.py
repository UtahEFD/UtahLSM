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
"""Abstract base class for canopy / vegetation models in UtahLSM.

Defines the interface all canopy parameterizations must implement and
provides shared utilities (root-fraction profile construction,
root-weighted soil-state reductions).

Design axioms (see CLAUDE.md / design discussion):

* Single surface temperature — no leaf/ground split, no canopy heat
  storage.
* A simple prognostic wet-canopy water reservoir may be used by the
  orchestrator for dewfall/interception evaporation.
* r_s is refreshed once per outer SEB/SMB coupling iteration, not
  inside Brent evaluations.
* All per-column quantities carry shape ``(ncol,)``; per-layer-per-
  column quantities carry shape ``(nz, ncol)``.
"""

import logging
from abc import ABC, abstractmethod
from typing import TypeVar, cast

import numpy as np
from numpy.typing import NDArray

from ..._types import FloatOrArrayLike
from ...data_models import AtmosphericState, SoilState, SurfaceState
from ...util.io import logging_helper

CP = TypeVar('CP', bound='Canopy')
logger = logging_helper.get_logger('CANOPY')


class Canopy(ABC):
    """Abstract base class for canopy models.

    Subclasses must implement :meth:`compute_resistance`. The base class
    handles the root-distribution bookkeeping (per-column Jackson-1996
    profile truncated to ``rooting_depth``) and provides
    :meth:`root_zone_mean` for taking root-weighted averages of any
    layer field.

    Attributes:
        lai: Leaf area index per column [m^2/m^2], shape (ncol,).
        veg_fraction: Vegetated surface fraction per column, shape
            (ncol,). `(1 - veg_fraction)` is treated as bare soil.
        rooting_depth: Depth over which the root profile integrates to
            unity per column [m], shape (ncol,).
        beta: Jackson-1996 β parameter per column (dimensionless),
            shape (ncol,).
        rs_min: Minimum bulk stomatal resistance per column [s/m],
            shape (ncol,).
        rs_max: Maximum (cuticular) resistance per column [s/m],
            shape (ncol,).
        r_ground: In-canopy aerodynamic resistance for heat transport
            from the radiative skin to the soil top [s/m] per column,
            shape (ncol,). Used in series with the soil's top-cell
            conductive resistance to attenuate ground heat flux under
            vegetation.
        water_capacity_lai: Wet-canopy water holding capacity per LAI
            [kg/m^2 per LAI], shape (ncol,).
        wet_cooling_max: Maximum diagnostic nighttime wet-canopy cooling
            below the soil/radiative skin [K], shape (ncol,).
        root_fraction: Precomputed per-layer root fraction, shape
            (nz, ncol). Columns sum to 1.
    """

    def __init__(
        self,
        lai: NDArray[np.float64],
        veg_fraction: NDArray[np.float64],
        rooting_depth: NDArray[np.float64],
        beta: NDArray[np.float64],
        rs_min: NDArray[np.float64],
        rs_max: NDArray[np.float64],
        r_ground: NDArray[np.float64],
        water_capacity_lai: NDArray[np.float64],
        wet_cooling_max: NDArray[np.float64],
        z: NDArray[np.float64],
    ) -> None:
        """Initializes the canopy base.

        Args:
            lai: Leaf area index (ncol,).
            veg_fraction: Vegetated fraction in [0, 1] (ncol,).
            rooting_depth: Rooting depth [m], positive (ncol,).
            beta: Jackson-1996 β parameter (ncol,). Values near 0.96-0.97
                are typical for grasslands.
            rs_min: Minimum stomatal resistance [s/m] (ncol,).
            rs_max: Maximum resistance [s/m] (ncol,).
            r_ground: In-canopy aerodynamic resistance to ground heat
                transport [s/m] (ncol,).
            water_capacity_lai: Wet-canopy water holding capacity per LAI
                [kg/m^2 per LAI] (ncol,).
            wet_cooling_max: Maximum diagnostic nighttime wet-canopy cooling
                below the soil/radiative skin [K] (ncol,).
            z: Soil layer node depths [m], shape (nz,). Values are
                non-positive with `z[0] = 0` at the surface.
        """
        self.logger: logging.Logger = logging_helper.get_logger('CANOPY')
        ncol = self._infer_column_count(
            {
                'lai': lai,
                'veg_fraction': veg_fraction,
                'rooting_depth': rooting_depth,
                'beta': beta,
                'rs_min': rs_min,
                'rs_max': rs_max,
                'r_ground': r_ground,
                'water_capacity_lai': water_capacity_lai,
                'wet_cooling_max': wet_cooling_max,
            }
        )
        self.ncol = ncol
        self.lai = self._as_column_param('lai', lai, ncol)
        self.veg_fraction = self._as_column_param(
            'veg_fraction', veg_fraction, ncol
        )
        self.rooting_depth = self._as_column_param(
            'rooting_depth', rooting_depth, ncol
        )
        self.beta = self._as_column_param('beta', beta, ncol)
        self.rs_min = self._as_column_param('rs_min', rs_min, ncol)
        self.rs_max = self._as_column_param('rs_max', rs_max, ncol)
        self.r_ground = self._as_column_param('r_ground', r_ground, ncol)
        self.water_capacity_lai = self._as_column_param(
            'water_capacity_lai', water_capacity_lai, ncol
        )
        self.wet_cooling_max = self._as_column_param(
            'wet_cooling_max', wet_cooling_max, ncol
        )
        self.z = np.asarray(z, dtype=float)

        self.root_fraction: NDArray[np.float64] = self._build_root_profile()

    # --- Shared utilities ---

    @staticmethod
    def _infer_column_count(params: dict[str, FloatOrArrayLike]) -> int:
        """Infers the common column count from scalar-or-column parameters."""
        sizes: dict[str, int] = {}
        for name, value in params.items():
            arr = np.asarray(value, dtype=float)
            if arr.ndim > 0 and arr.size != 1:
                sizes[name] = arr.size
        unique_sizes = set(sizes.values())
        if len(unique_sizes) > 1:
            detail = ', '.join(
                f"{name}={size}" for name, size in sorted(sizes.items())
            )
            raise ValueError(
                "Canopy per-column parameter sizes disagree: "
                f"{detail}."
            )
        return unique_sizes.pop() if unique_sizes else 1

    @staticmethod
    def _as_column_param(
        name: str,
        value: FloatOrArrayLike,
        ncol: int,
    ) -> NDArray[np.float64]:
        """Returns a scalar-or-column canopy parameter as shape ``(ncol,)``."""
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            return np.full(ncol, float(arr))
        if arr.size == ncol:
            return cast(NDArray[np.float64], arr.reshape(ncol))
        if arr.size == 1:
            return np.full(ncol, float(arr.flat[0]))
        raise ValueError(
            f"Canopy parameter '{name}' has size {arr.size}; expected "
            f"1 or ncol={ncol}."
        )

    def _build_root_profile(self) -> NDArray[np.float64]:
        """Constructs the normalized per-layer root fraction.

        Uses the Jackson (1996) cumulative root distribution
        ``Y(d) = 1 - β^(100·d)`` with d in metres (β expressed per cm).
        Each layer's raw fraction is the difference of Y at its upper
        and lower interfaces, truncated at ``rooting_depth`` and
        renormalized so each column sums to 1.

        Returns:
            Root fraction array of shape (nz, ncol).
        """
        # Layer interfaces: z is negative-downward with z[0]=0. Build
        # per-layer top/bottom depths (positive, metres).
        nz = self.z.size
        ncol = self.ncol
        d_top = np.zeros(nz)
        d_bot = np.zeros(nz)
        # Treat each soil node as owning the half-cells above/below it.
        # Top of layer i:
        for i in range(nz):
            if i == 0:
                d_top[i] = 0.0
            else:
                d_top[i] = -0.5 * (self.z[i - 1] + self.z[i])
            if i == nz - 1:
                d_bot[i] = -self.z[i] + 0.5 * (self.z[i - 1] - self.z[i])
            else:
                d_bot[i] = -0.5 * (self.z[i] + self.z[i + 1])

        # Broadcast to (nz, ncol).
        d_top2 = np.broadcast_to(d_top[:, None], (nz, ncol))
        d_bot2 = np.broadcast_to(d_bot[:, None], (nz, ncol))
        rooting = self.rooting_depth[None, :]
        beta = self.beta[None, :]

        # Truncate each layer's lower bound at the column's rooting_depth.
        d_top_c = np.minimum(d_top2, rooting)
        d_bot_c = np.minimum(d_bot2, rooting)

        # Jackson cumulative in cm; convert m -> cm via ×100.
        y_top = 1.0 - np.power(beta, 100.0 * d_top_c)
        y_bot = 1.0 - np.power(beta, 100.0 * d_bot_c)
        raw = np.maximum(y_bot - y_top, 0.0)

        total = np.sum(raw, axis=0, keepdims=True)
        # Guard against degenerate config (e.g. rooting_depth = 0).
        total = np.where(total > 0.0, total, 1.0)
        return raw / total

    def root_zone_mean(
        self, field: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Returns the root-weighted column mean of a per-layer field.

        Args:
            field: Array with shape (nz,) or (nz, ncol).

        Returns:
            Per-column mean of shape (ncol,).
        """
        f = np.asarray(field)
        if f.ndim == 1:
            f = f[:, None]
        return cast(
            NDArray[np.float64], np.sum(self.root_fraction * f, axis=0))

    # --- Abstract interface ---

    @abstractmethod
    def compute_resistance(
        self,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
        soil_state: SoilState,
        theta_wilt: NDArray[np.float64],
        theta_fc: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Bulk stomatal + cuticular resistance [s/m] per column.

        Args:
            atm_state: Current atmospheric forcing.
            sfc_state: Current surface state (temperature, humidity,
                friction velocity). The caller guarantees these reflect
                the most recent SEB pass of the outer coupling loop.
            soil_state: Current soil profile (for root-zone moisture).
            theta_wilt: Per-layer permanent wilting point [m^3/m^3],
                shape (nz,) or (nz, ncol).
            theta_fc: Per-layer field capacity [m^3/m^3], shape (nz,)
                or (nz, ncol).

        Returns:
            Bulk resistance r_s [s/m] of shape (ncol,).
        """
        raise NotImplementedError
