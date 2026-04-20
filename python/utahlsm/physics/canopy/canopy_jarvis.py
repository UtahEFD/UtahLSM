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
"""Jarvis-style stomatal resistance parameterization.

Implements the multiplicative Jarvis (1976) / Noilhan-Planton (1989)
form

    r_s = r_s,min / (LAI · f1(R) · f2(VPD) · f3(T) · f4(θ_root))

with each stress factor in ``[0, 1]``. When the radiation, VPD,
temperature, or moisture factor collapses, ``r_s`` saturates at
``r_s,max`` (cuticular ceiling).
"""

import numpy as np
from numpy.typing import NDArray

from ...data_models import AtmosphericState, SoilState, SurfaceState
from ...physics import thermo
from ...util import constants as c
from .canopy import Canopy


class CanopyJarvis(Canopy):
    """Jarvis multiplicative-stress canopy model.

    Attributes:
        rg_half: Half-saturation radiation [W/m^2] for f1(R), per column.
        vpd_coef: VPD sensitivity parameter [1/Pa] for f2(VPD), per column.
        t_opt: Optimal leaf temperature [K] for f3(T), per column.
        t_coef: Width of the temperature optimum [1/K^2], per column.
    """

    def __init__(
        self,
        lai: NDArray[np.float64],
        veg_fraction: NDArray[np.float64],
        rooting_depth: NDArray[np.float64],
        beta: NDArray[np.float64],
        rs_min: NDArray[np.float64],
        rs_max: NDArray[np.float64],
        rg_half: NDArray[np.float64],
        vpd_coef: NDArray[np.float64],
        t_opt: NDArray[np.float64],
        t_coef: NDArray[np.float64],
        z: NDArray[np.float64],
    ) -> None:
        super().__init__(
            lai=lai,
            veg_fraction=veg_fraction,
            rooting_depth=rooting_depth,
            beta=beta,
            rs_min=rs_min,
            rs_max=rs_max,
            z=z,
        )
        self.logger.info('Using the Jarvis canopy model')
        self.rg_half = np.asarray(rg_half, dtype=float)
        self.vpd_coef = np.asarray(vpd_coef, dtype=float)
        self.t_opt = np.asarray(t_opt, dtype=float)
        self.t_coef = np.asarray(t_coef, dtype=float)

    def compute_resistance(
        self,
        atm_state: AtmosphericState,
        sfc_state: SurfaceState,
        soil_state: SoilState,
        theta_wilt: NDArray[np.float64],
        theta_fc: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        f1 = self._f_radiation(atm_state.radiation_net)
        f2 = self._f_vpd(atm_state)
        f3 = self._f_temperature(sfc_state.temperature)
        f4 = self._f_moisture(soil_state.moisture, theta_wilt, theta_fc)

        F = np.clip(f1 * f2 * f3 * f4, 1e-6, 1.0)
        lai_eff = np.maximum(self.lai, 1e-6)
        r_s = self.rs_min / (lai_eff * F)
        return np.minimum(r_s, self.rs_max)

    # --- Stress functions ---

    def _f_radiation(
        self, rad_net: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """f1(R) — radiation stress.

        Uses the Noilhan-Planton saturating form ``(1 + R/R_½) /
        (1 + LAI·R_½/R_½,min)``. We simplify to the widely-used
        variant ``f = R / (R + R_½)`` with R clipped at zero so
        stomata close at night.

        Args:
            rad_net: Net radiation [W/m^2] (ncol,).

        Returns:
            Radiation stress factor in [0, 1] (ncol,).
        """
        R = np.maximum(np.asarray(rad_net, dtype=float), 0.0)
        return R / (R + self.rg_half)

    def _f_vpd(self, atm_state: AtmosphericState) -> NDArray[np.float64]:
        """f2(VPD) — atmospheric dryness stress.

        ``f = 1 / (1 + VPD_coef · VPD)`` with VPD computed from air
        temperature, pressure, and specific humidity.
        """
        q_sat = thermo.saturation_specific_humidity(
            atm_state.temperature, atm_state.pressure
        )
        # Convert specific-humidity deficit to vapor-pressure deficit
        # (Pa) via e ≈ q·p/ε for small q.
        vpd = np.maximum(
            (q_sat - atm_state.specific_humidity)
            * atm_state.pressure / c.thermodynamic.EPSILON,
            0.0,
        )
        return 1.0 / (1.0 + self.vpd_coef * vpd)

    def _f_temperature(
        self, leaf_T: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """f3(T) — leaf temperature stress (parabolic about t_opt).

        In this big-leaf formulation the surface temperature doubles as
        the leaf temperature, so the caller passes ``sfc_state.temperature``
        rather than the atmospheric forcing. This preserves stress signal
        during high-insolation/low-wind conditions where T_leaf ≫ T_air.
        """
        stress = 1.0 - self.t_coef * (self.t_opt - leaf_T) ** 2
        return np.clip(stress, 0.0, 1.0)

    def _f_moisture(
        self,
        soil_moisture: NDArray[np.float64],
        theta_wilt: NDArray[np.float64],
        theta_fc: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """f4(θ_root) — root-zone moisture stress.

        Linear ramp from 0 at θ_wilt to 1 at θ_fc, computed from the
        root-fraction-weighted moisture, wilting point, and field
        capacity. Broadcasting handles both (nz,) and (nz, ncol)
        soil-property arrays.
        """
        theta_r = self.root_zone_mean(soil_moisture)
        theta_wilt_r = self.root_zone_mean(theta_wilt)
        theta_fc_r = self.root_zone_mean(theta_fc)
        denom = np.maximum(theta_fc_r - theta_wilt_r, 1e-6)
        return np.clip((theta_r - theta_wilt_r) / denom, 0.0, 1.0)
