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
"""Van Genuchten (1980) soil physics parameterization.

This module provides an implementation of the Soil abstract base class using
the hydraulic relationships described by van Genuchten (1980). This is a
widely used, flexible model for describing the soil water retention curve.
"""

import logging
from typing import Union

import numpy as np
from numpy.typing import NDArray

from .soil import Soil
from ...util.io import logging_helper

class VanGenuchten(Soil):
    """Implements the van Genuchten (1980) soil physics model.

    This class provides concrete implementations for calculating water
    potential,
    hydraulic conductivity, and diffusivity based on the van Genuchten model.
    """
    def __init__(self, properties_dict: dict, soil_type_names: list,
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

    def surface_water_content(self, psi_sfc: float) -> float:
        """Computes surface soil water content from surface water potential.

        Args:
            psi_sfc: The soil water potential at the surface [m].

        Returns:
            The volumetric soil moisture content at the surface [m^3/m^3].
        """
        b = self.properties.b[0]
        psi_sat = self.properties.psi_sat[0]
        porosity = self.properties.porosity[0]
        residual = self.properties.residual[0]
        soil_e = porosity-residual
        m = 1 / (1+b)
        soil_q = residual + soil_e * (1/(1 + (psi_sfc/psi_sat)**(1/(1-m))))**(m)

        return soil_q

    def water_potential(
        self, soil_q: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
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
        self._validate_moisture_bounds(soil_q, level)

        if level is not None:
            b = self.properties.b[level]
            psi_sat = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
        else:
            b = self.properties.b
            psi_sat = self.properties.psi_sat
            porosity = self.properties.porosity
            residual = self.properties.residual

        Se = (soil_q-residual)/(porosity-residual)
        m = 1 / (1+b)
        psi = psi_sat*( ( (Se**(-1/m))-1 )**(1-m) )

        return psi

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
        self._validate_moisture_bounds(soil_q, level)

        if level is not None:
            b = self.properties.b[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
            K_sat = self.properties.K_sat[level]
        else:
            b = self.properties.b
            porosity = self.properties.porosity
            residual = self.properties.residual
            K_sat = self.properties.K_sat

        Se = (soil_q-residual)/(porosity-residual)
        m = 1 / (1+b)
        conductivity = (K_sat * np.sqrt(Se) *
                        ((1 - (1 - (Se**(1/m)) )**m )**2))

        return conductivity

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
        self._validate_moisture_bounds(soil_q)

        b = self.properties.b
        psi_sat = self.properties.psi_sat
        porosity = self.properties.porosity
        residual = self.properties.residual
        K_sat = self.properties.K_sat
        Se = (soil_q-residual)/(porosity-residual)
        soil_e = porosity-residual
        m = 1 / (1+b)
        A = (1-m)*K_sat*psi_sat / (m*soil_e)
        C = Se**(0.5-(1/m))*( (1 - Se**(1/m))**(-m) + (1- Se**(1/m))**m - 2 )
        diffusivity  = A*C

        return diffusivity
