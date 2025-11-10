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
import numpy as np
from numpy.typing import NDArray
from typing import Union
from .soil import Soil
from ...util import constants as c
from ...util.io import Input, logging_helper

class BrooksCorey(Soil):
    """Implements the Brooks and Corey (1964) soil physics model.
    
    This class provides concrete implementations for calculating water potential,
    hydraulic conductivity, and diffusivity based on the Brooks-Corey model.
    """
    def __init__(self, dataset_id: int, soil_type_array: NDArray[np.int_]):
        """Initializes the BrooksCorey soil model.
        
        Args:
            dataset_id: An integer ID for the soil parameter dataset to use.
            soil_type_array: A NumPy array of soil type IDs for each layer.
        """
        self.logger: logging.Logger = logging_helper.get_logger("SOIL")
        self.logger.info("--- Using the Brooks-Corey model")
        super().__init__(dataset_id, soil_type_array)
    
    def surface_water_content(self, psi: float) -> float:
        """Computes surface soil water content from surface water potential.
        
        Args:
            psi: The soil water potential at the surface [m].
        
        Returns:
            The volumetric soil moisture content at the surface [m^3/m^3].
        """
        b = self.properties.b[0]
        psi_sat = self.properties.psi_sat[0]
        porosity = self.properties.porosity[0]
        residual = self.properties.residual[0]
        soil_e = porosity-residual
        soil_q = residual+soil_e*( (psi_sat/psi)**(1./b) )
        
        return soil_q
    
    def water_potential(self, soil_q: Union[float,  NDArray[np.float64]], level: int = None) -> Union[float,  NDArray[np.float64]]:
        """Computes soil water potential from soil moisture.
        
        Args:
            soil_q: Soil moisture content [m^3/m^3]. Can be a scalar for a
                single level or a NumPy array for the entire column.
            level: The specific soil layer index. Required if `soil_q` is a
                scalar, ignored if it is an array. Defaults to None.
        
        Returns:
            The soil water potential in meters [m].
        """
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
        psi = psi_sat*( Se**(-b) )
        
        return psi
        
    def conductivity_moisture(self, soil_q: Union[float,  NDArray[np.float64]], level: int = None) -> Union[float,  NDArray[np.float64]]:
        """Computes soil hydraulic conductivity from soil moisture.
        
        Args:
            soil_q: Soil moisture content [m^3/m^3]. Can be a scalar or an array.
            level: The specific soil layer index if `soil_q` is a scalar.
                Defaults to None.
        
        Returns:
            The soil hydraulic conductivity [m/s].
        """
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
        conductivity = K_sat*( Se**(2.*b+3.) )
        
        return conductivity
    
    def diffusivity_moisture(self, soil_q:  NDArray[np.float64]) ->  NDArray[np.float64]:
        """Computes soil moisture diffusivity for the entire soil column.
        
        Args:
            soil_q: A NumPy array of soil moisture content for all layers [m^3/m^3].
        
        Returns:
            The soil moisture diffusivity for all layers [m^2/s].
        """
        b = self.properties.b
        psi_sat = self.properties.psi_sat
        porosity = self.properties.porosity
        residual = self.properties.residual
        K_sat = self.properties.K_sat
        Se = (soil_q-residual)/(porosity-residual)
        diffusivity = -b*K_sat*psi_sat*( Se**(b+2.) ) / (porosity-residual)
        
        return diffusivity