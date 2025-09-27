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
import numpy as np
from numpy.typing import NDArray
from typing import Union
from .soil import Soil
from ...util import constants as c
from ...util.io import Input, logging_helper

class Campbell(Soil):
    """Implements the Campbell (1974) soil physics model.
    
    This class provides concrete implementations for calculating water potential,
    hydraulic conductivity, and diffusivity based on the Campbell model.
    """

    def __init__(self,input: Input):
        """Initializes the Campbell soil model.
        
        Args:
            input: An `Input` object with the model's configuration settings.
        """
        self.logger = logging_helper.get_logger("SOIL")
        self.logger.info("Using the Campbell model")
        super().__init__(input)
    
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
        soil_q = porosity*(np.abs(psi_sat/psi)**(1./b))
                
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
        else:
            b = self.properties.b
            psi_sat = self.properties.psi_sat
            porosity = self.properties.porosity
        
        psi = psi_sat*((soil_q/porosity)**(-b))
        
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
            K_sat = self.properties.K_sat[level]
        else:
            b = self.properties.b
            porosity = self.properties.porosity
            K_sat = self.properties.K_sat
        conductivity = K_sat*( (soil_q/porosity)**(2.*b+3.) )
        
        return conductivity
    
    def diffusivity_moisture(self, soil_q:  NDArray[np.float64]) ->  NDArray[np.float64]:
        """Computes soil moisture diffusivity for the entire soil column.
        
        Args:
            soil_q: Soil moisture content for all layers [m^3/m^3].
        
        Returns:
            The soil moisture diffusivity for all layers [m^2/s].
        """
        b = self.properties.b
        psi_sat  = self.properties.psi_sat
        porosity = self.properties.porosity
        K_sat = self.properties.K_sat
        diffusivity  = -b*K_sat*psi_sat*( (soil_q/porosity)**(b+2.) ) / porosity
        
        return diffusivity