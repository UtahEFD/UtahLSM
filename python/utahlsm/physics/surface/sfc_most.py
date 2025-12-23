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
"""Surface layer parameterization based on Monin-Obukhov Similarity Theory.

This module provides an implementation of the `Surface` abstract base class
using standard Monin-Obukhov Similarity Theory (MOST) functions to describe
the stability and flux-profile relationships in the atmospheric surface layer.
"""
import logging
import math

import numpy as np

from ...util import constants as c
from ...util.io import logging_helper
from .sfc import Surface


class SurfaceMOST(Surface):
    """Implements surface layer physics using MOST.

    This class provides concrete implementations for the stability correction
    functions for momentum and heat based on the widely used Businger-Dyer
    relations.
    """
    def __init__(self, *, psi_stable: str = "dyer-hicks"):
        """Initializes the SurfaceMOST model.

        Args:
            psi_stable: Stable (z/L >= 0) integrated stability correction
                function (ψ) to use for momentum and heat. Options:
                - "dyer-hicks" (legacy linear form)
                - "beljaars-holtslag"
                - "cheng-brutsaert"
        """
        self.logger: logging.Logger = logging_helper.get_logger("SFC")
        psi_key = str(psi_stable).strip().lower()
        if psi_key == "beljaars-holtslag-1991":
            psi_key = "beljaars-holtslag"
        if psi_key == "cheng-brutsaert-2005":
            psi_key = "cheng-brutsaert"

        valid = {"dyer-hicks", "beljaars-holtslag", "cheng-brutsaert"}
        if psi_key not in valid:
            raise ValueError(
                f"Invalid psi_stable={psi_stable!r}. Valid options: "
                + ", ".join(sorted(valid))
            )
        self.psi_stable = psi_key

        self.logger.info("Using the MOST model (psi_stable=%s)", self.psi_stable)
        super().__init__()

    def _cap_obukhov_length(self, obukL: float, min_val: float = 0.1) -> float:
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
        obukL_mag = max(abs(obukL), min_val)
        return np.copysign(obukL_mag, obukL)

    def phim(self, z: float, obukL: float) -> float:
        """Computes the dimensionless stability function for momentum (phi_m).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of phi_m.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        return self.phim_stable(zeta) if zeta >= 0 else self.phim_unstable(zeta)

    def phim_stable(self, zeta: float) -> float:
        """Computes momentum stability function for stable conditions.

        Computes phi_m using the standard Holtslag and De Bruin stability
        function for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Momentum stability function value [dimensionless].
        """
        return 1. + 5.*zeta

    def phim_unstable(self, zeta: float) -> float:
        """Computes momentum stability function for unstable conditions.

        Computes phi_m using the Beljaars and Holtslag stability function
        for unstable stratification (zeta < 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Momentum stability function value [dimensionless].
        """
        return (1.-(16.*zeta))**(-0.25)

    def phih(self,z: float, obukL: float) -> float:
        """Computes the dimensionless stability function for heat (phi_h).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of phi_h.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        return self.phih_stable(zeta) if zeta >= 0 else self.phih_unstable(zeta)

    def phih_stable(self, zeta: float) -> float:
        """Computes heat stability function for stable conditions.

        Computes phi_h using the standard Holtslag and De Bruin stability
        function for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Heat stability function value [dimensionless].
        """
        return 1. + 5.*zeta

    def phih_unstable(self, zeta: float) -> float:
        """Computes heat stability function for unstable conditions.

        Computes phi_h using the Beljaars and Holtslag stability function
        for unstable stratification (zeta < 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Heat stability function value [dimensionless].
        """
        return (1.-(16.*zeta))**(-0.50)

    def psim(self,z: float,obukL: float) -> float:
        """Computes the integrated stability function for momentum (psi_m).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of psi_m.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        return self.psim_stable(zeta) if zeta >= 0 else self.psim_unstable(zeta)

    def psim_stable(self, zeta: float) -> float:
        """Computes integrated momentum stability function for stable conditions.

        Computes psi_m (the height-integrated stability function for momentum)
        for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Integrated momentum stability function value [dimensionless].
        """
        zeta = max(float(zeta), 0.0)
        if self.psi_stable == "dyer-hicks":
            return -5.0 * zeta
        if self.psi_stable == "beljaars-holtslag":
            # Holtslag & de Bruin (1988) / Beljaars & Holtslag (1991)
            a = 1.0
            b = 2.0 / 3.0
            c_ = 5.0
            d = 0.35
            return -(
                a * zeta
                + b * (zeta - (c_ / d)) * math.exp(-d * zeta)
                + b * (c_ / d)
            )
        if self.psi_stable == "cheng-brutsaert":
            # Cheng & Brutsaert (2005)
            a = 6.1
            b = 2.5
            inner = zeta + (1.0 + zeta**b) ** (1.0 / b)
            return -a * math.log(inner)
        raise AssertionError("Unhandled psi_stable")  # pragma: no cover

    def psim_unstable(self, zeta: float) -> float:
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
        x = (1.-(16.*zeta))**(0.25)
        return (2.*np.log((1.+x)/2.) + np.log((1.+x**2.)/2.) -
                2.*math.atan2(1., self.phim_unstable(zeta)) + PI/2.)

    def psih(self,z: float,obukL: float) -> float:
        """Computes the integrated stability function for heat (psi_h).

        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].

        Returns:
            The value of psi_h.
        """
        obukL_cap = self._cap_obukhov_length(obukL)

        zeta = z / obukL_cap
        return self.psih_stable(zeta) if zeta >= 0 else self.psih_unstable(zeta)

    def psih_stable(self, zeta: float) -> float:
        """Computes integrated heat stability function for stable conditions.

        Computes psi_h (the height-integrated stability function for heat)
        for stable stratification (zeta >= 0).

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Integrated heat stability function value [dimensionless].
        """
        zeta = max(float(zeta), 0.0)
        if self.psi_stable == "dyer-hicks":
            return -5.0 * zeta
        if self.psi_stable == "beljaars-holtslag":
            a = 1.0
            b = 2.0 / 3.0
            c_ = 5.0
            d = 0.35
            return -(
                a * zeta
                + b * (zeta - (c_ / d)) * math.exp(-d * zeta)
                + b * (c_ / d)
            )
        if self.psi_stable == "cheng-brutsaert":
            a = 5.3
            b = 1.1
            inner = zeta + (1.0 + zeta**b) ** (1.0 / b)
            return -a * math.log(inner)
        raise AssertionError("Unhandled psi_stable")  # pragma: no cover

    def psih_unstable(self, zeta: float) -> float:
        """Computes integrated heat stability function for unstable conditions.

        Computes psi_h (the height-integrated stability function for heat)
        for unstable stratification (zeta < 0) using the Paulson formulation.

        Args:
            zeta: Dimensionless height parameter (z/L), where L is the
                Obukhov length [dimensionless].

        Returns:
            Integrated heat stability function value [dimensionless].
        """
        x = (1.-(16.*zeta))**(0.50)
        return 2.*np.log((1.+x)/2.)

    def fm(self, z1: float, z0: float, obukL: float) -> float:
        """Computes the log-law stability function for momentum.

        Args:
            z1: Upper height [m].
            z0: Lower height (roughness length) [m].
            obukL: Obukhov length [m].

        Returns:
            The stability-corrected log-law function value.
        """
        VK = c.physical.VON_KARMAN
        fm = VK / (np.log(z1/z0) - self.psim(z1,obukL) + self.psim(z0,obukL))
        return fm

    def fh(self, z1: float, z0h: float, obukL: float) -> float:
        """Computes the log-law stability function for heat.

        Args:
            z1: Upper height [m].
            z0h: Lower height (thermal roughness length) [m].
            obukL: Obukhov length [m].

        Returns:
            The stability-corrected log-law function value.
        """
        VK = c.physical.VON_KARMAN
        fh = VK / (np.log(z1/z0h) - self.psih(z1,obukL) + self.psih(z0h,obukL))
        return fh
