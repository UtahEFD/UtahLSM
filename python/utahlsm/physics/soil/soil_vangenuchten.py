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
"""Van Genuchten (1980) soil physics parameterization.

This module provides an implementation of the Soil abstract base class using
the hydraulic relationships described by van Genuchten (1980). This is a
widely used, flexible model for describing the soil water retention curve.
"""

import logging
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray

from ..._types import FloatOrArray
from ...util.io import logging_helper
from .soil import Soil


class VanGenuchten(Soil):
    """Implements the van Genuchten (1980) soil physics model.

    This class provides concrete implementations for calculating water
    potential,
    hydraulic conductivity, and diffusivity based on the van Genuchten model.
    """
    def __init__(self, properties_dict: dict[str, dict[str, Any]],
                 soil_type_names: list[str],
                 dataset_name: str = 'custom'):
        """Initializes the VanGenuchten soil model.

        Args:
            properties_dict: Dictionary mapping soil type names to properties.
            soil_type_names: List of soil type names for each layer.
            dataset_name: Human-readable name of the dataset being used.
        """
        self.logger: logging.Logger = logging_helper.get_logger('SOIL')
        self.logger.info('Using the Van Genuchten model')
        super().__init__(properties_dict, soil_type_names, dataset_name)

    def water_potential(
        self, soil_q: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
        """Computes soil water potential from soil moisture.

        Args:
            soil_q: Soil moisture content [m^3/m^3]. Can be a scalar for a
                single level or an array for the entire column.
            level: The specific soil layer index. Required if `soil_q` is a
                scalar, ignored if it is an array. Defaults to None.

        Returns:
            The soil water potential in meters [m].

        Raises:
            ValueError: If soil_q is out of valid bounds.
        """
        if level is not None:
            b        = self.properties.b[level]
            psi_sat  = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
        else:
            _2d = np.asarray(soil_q).ndim == 2
            b        = self.properties.b[:, None]        if _2d else self.properties.b
            psi_sat  = self.properties.psi_sat[:, None]  if _2d else self.properties.psi_sat
            porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
            residual = self.properties.residual[:, None] if _2d else self.properties.residual

        Se = (soil_q-residual)/(porosity-residual)
        m = 1 / (1+b)

        # The Van Genuchten equation involves (Se^(-1/m) - 1)^(1-m)
        # At saturation (Se=1): (1 - 1)^(1-m) = 0^(1-m) → 0 (correct physically)
        # At residual (Se=0): (inf - 1)^(1-m) → inf (correct physically)
        # However, computing 0^(1-m) with float exponent causes RuntimeWarning
        # We suppress the warning and handle edge cases explicitly
        with np.errstate(divide='ignore', invalid='ignore'):
            inner_term = (Se**(-1/m)) - 1
            psi = psi_sat * (inner_term**(1-m))

            # Handle edge cases that may produce inf/nan
            if np.isscalar(psi):
                if Se >= 0.9999:  # Essentially saturated, psi → 0
                    psi = 0.0
                elif not np.isfinite(psi):  # Handle any inf or nan
                    psi = -np.inf if psi_sat < 0 else np.inf
            else:
                # For arrays: set saturated values to 0, keep others as computed
                psi = np.where(Se >= 0.9999, 0.0, psi)

        return cast(FloatOrArray, psi)

    def water_content(
        self, psi: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
        """Computes soil moisture from water potential."""
        if level is not None:
            b        = self.properties.b[level]
            psi_sat  = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
        else:
            _2d = np.asarray(psi).ndim == 2
            b        = self.properties.b[:, None]        if _2d else self.properties.b
            psi_sat  = self.properties.psi_sat[:, None]  if _2d else self.properties.psi_sat
            porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
            residual = self.properties.residual[:, None] if _2d else self.properties.residual

        soil_e = porosity - residual
        m = 1.0 / (1.0 + b)
        psi_arr = np.asarray(psi, dtype=float)
        ratio = np.maximum(psi_arr / psi_sat, 0.0)
        Se = (1.0 + ratio ** (1.0 / (1.0 - m))) ** (-m)
        soil_q = residual + soil_e * Se
        soil_q = np.where(psi_arr >= 0.0, porosity, soil_q)
        soil_q = np.clip(soil_q, residual, porosity)
        return float(soil_q.item()) if np.isscalar(psi) else soil_q

    def moisture_capacity(
        self, psi: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
        """Computes specific moisture capacity dθ/dψ."""
        if level is not None:
            b        = self.properties.b[level]
            psi_sat  = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
        else:
            _2d = np.asarray(psi).ndim == 2
            b        = self.properties.b[:, None]        if _2d else self.properties.b
            psi_sat  = self.properties.psi_sat[:, None]  if _2d else self.properties.psi_sat
            porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
            residual = self.properties.residual[:, None] if _2d else self.properties.residual

        m = 1.0 / (1.0 + b)
        n = 1.0 / (1.0 - m)
        soil_e = porosity - residual
        psi_arr = np.asarray(psi, dtype=float)
        ratio = np.maximum(psi_arr / psi_sat, 0.0)

        with np.errstate(divide='ignore', invalid='ignore'):
            capacity = (
                soil_e
                * (-m * n / psi_sat)
                * ratio ** (n - 1.0)
                * (1.0 + ratio ** n) ** (-m - 1.0)
            )
        capacity = np.where(np.isfinite(capacity), capacity, 0.0)
        capacity = np.where(psi_arr >= 0.0, 0.0, capacity)
        return float(capacity) if np.isscalar(psi) else capacity

    def conductivity_moisture(
        self, soil_q: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
        """Computes soil hydraulic conductivity from soil moisture.

        Args:
            soil_q: Soil moisture content [m^3/m^3]. Can be a scalar or
                an array.
            level: The specific soil layer index if `soil_q` is a scalar.
                Defaults to None.

        Returns:
            The soil hydraulic conductivity [m/s].

        Raises:
            ValueError: If soil_q is out of valid bounds.
        """
        if level is not None:
            b        = self.properties.b[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
            K_sat    = self.properties.K_sat[level]
        else:
            _2d = np.asarray(soil_q).ndim == 2
            b        = self.properties.b[:, None]        if _2d else self.properties.b
            porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
            residual = self.properties.residual[:, None] if _2d else self.properties.residual
            K_sat    = self.properties.K_sat[:, None]    if _2d else self.properties.K_sat

        Se = (soil_q-residual)/(porosity-residual)
        m = 1 / (1+b)
        # Cache Se**(1/m) to avoid repeated exponentiation
        Se_pow_inv_m = Se**(1/m)
        inner = 1 - (1 - Se_pow_inv_m)**m
        conductivity = K_sat * np.sqrt(Se) * inner**2

        return conductivity

    def conductivity_gradient(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes the linearized dK/dθ for the Van Genuchten model.

        K'_lin = K(Se) / Se / (φ - θ_r), the secant from Se=0 converted
        to θ-space. Clamped to 0 where Se ≈ 0 to avoid divergence.

        Args:
            soil_q: Soil moisture content for all layers [m^3/m^3].

        Returns:
            Linearized dK/dθ for all layers [m/s].
        """
        _2d = soil_q.ndim == 2
        b        = self.properties.b[:, None]        if _2d else self.properties.b
        porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
        residual = self.properties.residual[:, None] if _2d else self.properties.residual
        K_sat    = self.properties.K_sat[:, None]    if _2d else self.properties.K_sat
        soil_e = porosity - residual
        Se = (soil_q - residual) / soil_e
        m = 1.0 / (1.0 + b)
        Se_pow = Se ** (1.0 / m)
        inner = 1.0 - (1.0 - Se_pow) ** m
        # K(Se)/Se = K_sat * Se^(-0.5) * inner^2; diverges as Se -> 0
        with np.errstate(divide='ignore', invalid='ignore'):
            gradient = K_sat * inner ** 2 / (np.sqrt(Se) * soil_e)
        gradient = np.where(Se > 1e-10, gradient, 0.0)
        return gradient

    def diffusivity_moisture(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes soil moisture diffusivity for the entire soil column.

        Args:
            soil_q: The soil moisture content for all layers [m^3/m^3].

        Returns:
            The soil moisture diffusivity for all layers [m^2/s].

        Raises:
            ValueError: If soil_q is out of valid bounds.
        """
        _2d = soil_q.ndim == 2
        b        = self.properties.b[:, None]        if _2d else self.properties.b
        psi_sat  = self.properties.psi_sat[:, None]  if _2d else self.properties.psi_sat
        porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
        residual = self.properties.residual[:, None] if _2d else self.properties.residual
        K_sat    = self.properties.K_sat[:, None]    if _2d else self.properties.K_sat
        Se = (soil_q-residual)/(porosity-residual)
        soil_e = porosity-residual
        m = 1 / (1+b)
        A = (1-m)*K_sat*psi_sat / (m*soil_e)
        # Cache Se**(1/m) to avoid repeated exponentiation
        Se_pow_inv_m = Se**(1/m)
        one_minus_Se_pow = 1 - Se_pow_inv_m
        C = Se**(0.5-(1/m)) * (one_minus_Se_pow**(-m) + one_minus_Se_pow**m - 2)
        diffusivity = A*C

        return diffusivity
