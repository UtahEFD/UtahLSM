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
"""Surface layer parameterization based on Monin-Obukhov Similarity Theory.

This module provides an implementation of the `Surface` abstract base class
using standard Monin-Obukhov Similarity Theory (MOST) functions to describe
the stability and flux-profile relationships in the atmospheric surface layer.
"""
import logging
from typing import cast

import numpy as np

from ..._types import FloatOrArray
from ...util import constants as c
from ...util.io import logging_helper
from .sfc import Surface


class SurfaceMOST(Surface):
    """Implements surface layer physics using MOST.

    This class provides concrete implementations for the stability correction
    functions for momentum and heat.
    """
    def __init__(self, *, psi_stable: str = "dyer-hicks"):
        """Initializes the SurfaceMOST model.

        Args:
            psi_stable: Stable (z/L >= 0) integrated stability correction
                function (ψ) to use for momentum and heat. Options:
                - "dyer-hicks"
                - "beljaars-holtslag"
                - "cheng-brutsaert"
        """
        self.logger: logging.Logger = logging_helper.get_logger("SFC")
        psi_key = str(psi_stable).strip().lower()

        valid = {"dyer-hicks", "beljaars-holtslag", "cheng-brutsaert"}
        if psi_key not in valid:
            raise ValueError(
                f"Invalid psi_stable={psi_stable!r}. Valid options: "
                + ", ".join(sorted(valid))
            )
        self.psi_stable = psi_key

        self.logger.info("Using the MOST model (psi_stable=%s)", self.psi_stable)
        super().__init__()

    def _cap_obukhov_length(self, obukL: FloatOrArray,
                            min_val: float = 0.1) -> FloatOrArray:
        """Caps and preserves sign of Obukhov length.

        Ensures |obukL| >= min_val while preserving the original sign.
        This prevents division by zero and numerical instability in stability
        functions when Obukhov length is very small.

        Args:
            obukL: Original Obukhov length [m].
            min_val: Minimum magnitude threshold [m].

        Returns:
            Capped Obukhov length with original sign preserved.
        """
        obukL_arr = np.asarray(obukL, dtype=float)
        obukL_mag = np.maximum(np.abs(obukL_arr), min_val)
        return cast(FloatOrArray, np.copysign(obukL_mag, obukL_arr))

    def phim(self, z: float, obukL: FloatOrArray) -> FloatOrArray:
        """Computes the dimensionless stability function for momentum (phi_m).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of phi_m.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        stable = self.phim_stable(zeta)
        unstable = self.phim_unstable(zeta)
        return np.where(zeta >= 0, stable, unstable)

    def phim_stable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes momentum stability function for stable conditions.

        Computes phi_m using the standard Holtslag and De Bruin stability
        function for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Momentum stability function value [dimensionless].
        """
        return 1.0 + 5.0 * zeta

    def phim_unstable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes momentum stability function for unstable conditions.

        Computes phi_m using the Beljaars and Holtslag stability function
        for unstable stratification (zeta < 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Momentum stability function value [dimensionless].
        """
        # Clamp zeta <= 0 to ensure (1 - 16*zeta) > 0 for the power operation.
        # When called from vectorized code, stable zeta values are masked out
        # by np.where anyway.
        zeta = np.minimum(zeta, 0.0)
        return (1.0 - (16.0 * zeta))**(-0.25)

    def phih(self, z: float, obukL: FloatOrArray) -> FloatOrArray:
        """Computes the dimensionless stability function for heat (phi_h).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of phi_h.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        stable = self.phih_stable(zeta)
        unstable = self.phih_unstable(zeta)
        return np.where(zeta >= 0, stable, unstable)

    def phih_stable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes heat stability function for stable conditions.

        Computes phi_h using the standard Holtslag and De Bruin stability
        function for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Heat stability function value [dimensionless].
        """
        return 1.0 + 5.0 * zeta

    def phih_unstable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes heat stability function for unstable conditions.

        Computes phi_h using the Beljaars and Holtslag stability function
        for unstable stratification (zeta < 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Heat stability function value [dimensionless].
        """
        zeta = np.minimum(zeta, 0.0)
        return (1.0 - (16.0 * zeta))**(-0.50)

    def psim(self, z: float, obukL: FloatOrArray) -> FloatOrArray:
        """Computes the integrated stability function for momentum (psi_m).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of psi_m.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        stable = self.psim_stable(zeta)
        unstable = self.psim_unstable(zeta)
        return np.where(zeta >= 0, stable, unstable)

    def psim_stable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes integrated momentum stability function for stable conditions.

        Computes psi_m (the height-integrated stability function for momentum)
        for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Integrated momentum stability function value [dimensionless].
        """
        zeta = np.maximum(zeta, 0.0)
        if self.psi_stable == "dyer-hicks":
            return -5.0 * zeta
        if self.psi_stable == "beljaars-holtslag":
            # Holtslag & de Bruin (1988) / Beljaars & Holtslag (1991)
            a = 1.0
            b = 2.0 / 3.0
            c_ = 5.0
            d = 0.35
            c_over_d = c_ / d  # Pre-compute to avoid repeated division
            return -(
                a * zeta
                + b * (zeta - c_over_d) * np.exp(-d * zeta)
                + b * c_over_d
            )
        if self.psi_stable == "cheng-brutsaert":
            # Cheng & Brutsaert (2005)
            a = 6.1
            b = 2.5
            inner = zeta + (1.0 + zeta**b) ** (1.0 / b)
            return -a * np.log(inner)
        raise AssertionError("Unhandled psi_stable")  # pragma: no cover

    def psim_unstable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes integrated momentum stability function for unstable conditions.

        Computes psi_m (the height-integrated stability function for momentum)
        for unstable stratification (zeta < 0) using the Paulson formulation.

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Integrated momentum stability function value [dimensionless].
        """
        PI = c.physical.PI
        zeta = np.minimum(zeta, 0.0)
        x = (1.0 - (16.0 * zeta))**(0.25)
        return (
            2.0 * np.log((1.0 + x) / 2.0)
            + np.log((1.0 + x**2.0) / 2.0)
            - 2.0 * np.arctan(x)
            + PI / 2.0
        )

    def psih(self, z: float, obukL: FloatOrArray) -> FloatOrArray:
        """Computes the integrated stability function for heat (psi_h).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of psi_h.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        stable = self.psih_stable(zeta)
        unstable = self.psih_unstable(zeta)
        return np.where(zeta >= 0, stable, unstable)

    def psih_stable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes integrated heat stability function for stable conditions.

        Computes psi_h (the height-integrated stability function for heat)
        for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Integrated heat stability function value [dimensionless].
        """
        zeta = np.maximum(zeta, 0.0)
        if self.psi_stable == "dyer-hicks":
            return -5.0 * zeta
        if self.psi_stable == "beljaars-holtslag":
            a = 1.0
            b = 2.0 / 3.0
            c_ = 5.0
            d = 0.35
            c_over_d = c_ / d  # Pre-compute to avoid repeated division
            return -(
                a * zeta
                + b * (zeta - c_over_d) * np.exp(-d * zeta)
                + b * c_over_d
            )
        if self.psi_stable == "cheng-brutsaert":
            a = 5.3
            b = 1.1
            inner = zeta + (1.0 + zeta**b) ** (1.0 / b)
            return -a * np.log(inner)
        raise AssertionError("Unhandled psi_stable")  # pragma: no cover

    def psih_unstable(self, zeta: FloatOrArray) -> FloatOrArray:
        """Computes integrated heat stability function for unstable conditions.

        Computes psi_h (the height-integrated stability function for heat)
        for unstable stratification (zeta < 0) using the Paulson formulation.

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Integrated heat stability function value [dimensionless].
        """
        zeta = np.minimum(zeta, 0.0)
        x = (1.0 - (16.0 * zeta))**(0.50)
        return 2.0 * np.log((1.0 + x) / 2.0)

    def fm(self, z1: float, z0: float,
           obukhov_l: FloatOrArray) -> FloatOrArray:
        """Computes the log-law stability function for momentum.

        Args:
            z1: Upper height [m].
            z0: Lower height (roughness length) [m].
            obukhov_l: Obukhov length [m].

        Returns:
            The stability-corrected log-law function value.
        """
        VK = c.physical.VON_KARMAN
        return cast(FloatOrArray, VK / (
            np.log(z1 / z0) - self.psim(z1, obukhov_l) + self.psim(z0, obukhov_l)
        ))

    def fh(self, z1: float, z0h: float,
           obukhov_l: FloatOrArray) -> FloatOrArray:
        """Computes the log-law stability function for heat.

        Args:
            z1: Upper height [m].
            z0h: Lower height (thermal roughness length) [m].
            obukhov_l: Obukhov length [m].

        Returns:
            The stability-corrected log-law function value.
        """
        VK = c.physical.VON_KARMAN
        return cast(FloatOrArray, VK / (
            np.log(z1 / z0h) - self.psih(z1, obukhov_l) + self.psih(z0h, obukhov_l)
        ))
