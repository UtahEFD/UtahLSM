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
        # Guard against psi_sfc == 0 (saturated soil); return porosity
        psi_sfc = np.where(psi_sfc == 0, np.nan, psi_sfc)
        soil_q = residual+soil_e*( (psi_sat/psi_sfc)**(1./b) )
        soil_q = np.where(np.isnan(psi_sfc), porosity, soil_q)

        return soil_q

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
