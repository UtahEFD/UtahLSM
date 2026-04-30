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
"""Campbell (1974) soil physics parameterization.

This module provides an implementation of the Soil abstract base class using
the hydraulic relationships described by Campbell (1974). This model is
often used for its simplicity and effectiveness in representing soil hydraulic
properties.
"""

import logging
from typing import Union

import numpy as np
from numpy.typing import NDArray

from ...util.io import logging_helper
from .soil import Soil


class Campbell(Soil):
    """Implements the Campbell (1974) soil physics model.

    This class provides concrete implementations for calculating water
    potential,
    hydraulic conductivity, and diffusivity based on the Campbell model.
    """

    def __init__(self, properties_dict: dict, soil_type_names: list,
                 dataset_name: str = 'custom'):
        """Initializes the Campbell soil model.

        Args:
            properties_dict: Dictionary mapping soil type names to properties.
            soil_type_names: List of soil type names for each layer.
            dataset_name: Human-readable name of the dataset being used.
        """
        self.logger: logging.Logger = logging_helper.get_logger('SOIL')
        self.logger.info('Using the Campbell model')
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
        else:
            soil_q_arr = np.asarray(soil_q)
            b = self._expand_profile_property(self.properties.b, soil_q_arr)
            psi_sat = self._expand_profile_property(
                self.properties.psi_sat, soil_q_arr)
            porosity = self._expand_profile_property(
                self.properties.porosity, soil_q_arr)

        ratio = np.maximum(soil_q/porosity, 1e-12)
        psi = psi_sat*(ratio**(-b))

        return psi

    def water_content(
        self, psi: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
        """Computes soil moisture from water potential."""
        if level is not None:
            b = self.properties.b[level]
            psi_sat = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
        else:
            psi_arr = np.asarray(psi)
            b = self._expand_profile_property(self.properties.b, psi_arr)
            psi_sat = self._expand_profile_property(
                self.properties.psi_sat, psi_arr)
            porosity = self._expand_profile_property(
                self.properties.porosity, psi_arr)

        psi_arr = np.asarray(psi, dtype=float)
        safe_psi = np.where(np.abs(psi_arr) < 1e-12, np.nan, psi_arr)
        with np.errstate(divide='ignore', invalid='ignore'):
            soil_q = porosity * np.abs(psi_sat / safe_psi) ** (1.0 / b)
        soil_q = np.where(psi_arr >= psi_sat, porosity, soil_q)
        soil_q = np.clip(soil_q, 0.0, porosity)
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
        else:
            b = self._expand_profile_property(
                self.properties.b, psi_arr)

        with np.errstate(divide='ignore', invalid='ignore'):
            capacity = -soil_q_arr / (b * psi_arr)
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
            K_sat = self.properties.K_sat[level]
        else:
            soil_q_arr = np.asarray(soil_q)
            b = self._expand_profile_property(self.properties.b, soil_q_arr)
            porosity = self._expand_profile_property(
                self.properties.porosity, soil_q_arr)
            K_sat = self._expand_profile_property(
                self.properties.K_sat, soil_q_arr)
        conductivity = K_sat * ((soil_q/porosity)**(2.*b+3.))

        return conductivity

    def conductivity_gradient(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes the linearized dK/dθ for the Campbell model.

        For Campbell (θ_r = 0): K'_lin = K_sat/φ · (θ/φ)^(2b+2) = K(θ)/θ.

        Args:
            soil_q: Soil moisture content for all layers [m^3/m^3].

        Returns:
            Linearized dK/dθ for all layers [m/s].
        """
        b = self._expand_profile_property(self.properties.b, soil_q)
        porosity = self._expand_profile_property(
            self.properties.porosity, soil_q)
        K_sat = self._expand_profile_property(self.properties.K_sat, soil_q)
        gradient = K_sat / porosity * (soil_q / porosity) ** (2.0 * b + 2.0)
        return gradient

    def diffusivity_moisture(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes soil moisture diffusivity for the entire soil column.

        Args:
            soil_q: Soil moisture content for all layers [m^3/m^3].

        Returns:
            The soil moisture diffusivity for all layers [m^2/s].

        Raises:
            ValueError: If soil_q is out of valid bounds.
        """
        b = self._expand_profile_property(self.properties.b, soil_q)
        psi_sat = self._expand_profile_property(self.properties.psi_sat, soil_q)
        porosity = self._expand_profile_property(
            self.properties.porosity, soil_q)
        K_sat = self._expand_profile_property(self.properties.K_sat, soil_q)
        diffusivity  = -b*K_sat*psi_sat*( (soil_q/porosity)**(b+2.) ) / porosity

        return diffusivity
