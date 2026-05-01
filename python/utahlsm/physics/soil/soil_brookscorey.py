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
"""Brooks and Corey (1964) soil physics parameterization.

This module provides an implementation of the Soil abstract base class using
the hydraulic relationships described by Brooks and Corey (1964). This is a
widely used model for describing soil hydraulic properties.
"""

import logging
from typing import Union

import numpy as np
from numpy.typing import NDArray

from ...util.io import logging_helper
from .soil import Soil


class BrooksCorey(Soil):
    """Implements the Brooks and Corey (1964) soil physics model.

    This class provides concrete implementations for calculating water
    potential,
    hydraulic conductivity, and diffusivity based on the Brooks-Corey model.
    """
    def __init__(self, properties_dict: dict, soil_type_names: list,
                 dataset_name: str = 'custom'):
        """Initializes the BrooksCorey soil model.

        Args:
            properties_dict: Dictionary mapping soil type names to properties.
            soil_type_names: List of soil type names for each layer.
            dataset_name: Human-readable name of the dataset being used.
        """
        self.logger: logging.Logger = logging_helper.get_logger('SOIL')
        self.logger.info('--- Using the Brooks-Corey model')
        super().__init__(properties_dict, soil_type_names, dataset_name)

    def water_potential(
        self, soil_q: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
        """Computes soil water potential from soil moisture.

        Args:
            soil_q: Soil moisture content [m^3/m^3]. Can be a scalar for a
                single level or a NumPy array for the entire column.
            level: The specific soil layer index. Required if `soil_q` is a
                scalar, ignored if it is an array. Defaults to None.

        Returns:
            The soil water potential in meters [m].

        Raises:
            ValueError: If soil_q is out of valid bounds.
        """
        if level is not None:
            b = self.properties.b[level]
            psi_sat = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
        else:
            soil_q_arr = np.asarray(soil_q)
            b = self._expand_profile_property(self.properties.b, soil_q_arr)
            psi_sat = self._expand_profile_property(
                self.properties.psi_sat, soil_q_arr)
            porosity = self._expand_profile_property(
                self.properties.porosity, soil_q_arr)
            residual = self._expand_profile_property(
                self.properties.residual, soil_q_arr)

        Se = np.maximum((soil_q-residual)/(porosity-residual), 1e-12)
        psi = psi_sat*( Se**(-b) )

        return psi

    def water_content(
        self, psi: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
        """Computes soil moisture from water potential."""
        if level is not None:
            b = self.properties.b[level]
            psi_sat = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
        else:
            psi_arr = np.asarray(psi)
            b = self._expand_profile_property(self.properties.b, psi_arr)
            psi_sat = self._expand_profile_property(
                self.properties.psi_sat, psi_arr)
            porosity = self._expand_profile_property(
                self.properties.porosity, psi_arr)
            residual = self._expand_profile_property(
                self.properties.residual, psi_arr)

        soil_e = porosity - residual
        psi_arr = np.asarray(psi, dtype=float)
        safe_psi = np.where(np.abs(psi_arr) < 1e-12, np.nan, psi_arr)
        with np.errstate(divide='ignore', invalid='ignore'):
            Se = np.abs(psi_sat / safe_psi) ** (1.0 / b)
        soil_q = residual + soil_e * Se
        soil_q = np.where(psi_arr >= psi_sat, porosity, soil_q)
        soil_q = np.clip(soil_q, residual, porosity)
        return soil_q.item() if np.isscalar(psi) else soil_q

    def moisture_capacity(
        self, psi: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
        """Computes specific moisture capacity dθ/dψ."""
        soil_q = self.water_content(psi, level=level)
        psi_arr = np.asarray(psi, dtype=float)
        soil_q_arr = np.asarray(soil_q, dtype=float)
        if level is not None:
            b = self.properties.b[level]
            residual = self.properties.residual[level]
        else:
            b = self._expand_profile_property(self.properties.b, psi_arr)
            residual = self._expand_profile_property(
                self.properties.residual, psi_arr)

        with np.errstate(divide='ignore', invalid='ignore'):
            capacity = -(soil_q_arr - residual) / (b * psi_arr)
        capacity = np.where(np.isfinite(capacity), capacity, 0.0)
        capacity = np.where(psi_arr >= 0.0, 0.0, capacity)
        return float(capacity) if np.isscalar(psi) else capacity

    def conductivity_moisture(
        self, soil_q: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
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
            b = self.properties.b[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
            K_sat = self.properties.K_sat[level]
        else:
            soil_q_arr = np.asarray(soil_q)
            b = self._expand_profile_property(self.properties.b, soil_q_arr)
            porosity = self._expand_profile_property(
                self.properties.porosity, soil_q_arr)
            residual = self._expand_profile_property(
                self.properties.residual, soil_q_arr)
            K_sat = self._expand_profile_property(
                self.properties.K_sat, soil_q_arr)

        Se = (soil_q-residual)/(porosity-residual)
        conductivity = K_sat*( Se**(2.*b+3.) )

        return conductivity

    def conductivity_gradient(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes the linearized dK/dθ for the Brooks-Corey model.

        K'_lin = K_sat · Se^(2b+2) / (φ - θ_r), the secant linearization
        using effective saturation.

        Args:
            soil_q: Soil moisture content for all layers [m^3/m^3].

        Returns:
            Linearized dK/dθ for all layers [m/s].
        """
        b = self._expand_profile_property(self.properties.b, soil_q)
        porosity = self._expand_profile_property(
            self.properties.porosity, soil_q)
        residual = self._expand_profile_property(
            self.properties.residual, soil_q)
        K_sat = self._expand_profile_property(self.properties.K_sat, soil_q)
        soil_e = porosity - residual
        Se = (soil_q - residual) / soil_e
        gradient = K_sat * Se ** (2.0 * b + 2.0) / soil_e
        return gradient

    def diffusivity_moisture(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes soil moisture diffusivity for the entire soil column.

        Args:
            soil_q: A NumPy array of soil moisture content for all
                layers [m^3/m^3].

        Returns:
            The soil moisture diffusivity for all layers [m^2/s].

        Raises:
            ValueError: If soil_q is out of valid bounds.
        """
        b = self._expand_profile_property(self.properties.b, soil_q)
        psi_sat = self._expand_profile_property(self.properties.psi_sat, soil_q)
        porosity = self._expand_profile_property(
            self.properties.porosity, soil_q)
        residual = self._expand_profile_property(
            self.properties.residual, soil_q)
        K_sat = self._expand_profile_property(self.properties.K_sat, soil_q)
        Se = (soil_q-residual)/(porosity-residual)
        diffusivity = -b*K_sat*psi_sat*( Se**(b+2.) ) / (porosity-residual)

        return diffusivity
