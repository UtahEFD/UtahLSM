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

import numpy as np
from typing import Union
from .soil import Soil
from ...util import constants as c
from ...util.io import Input, logging_helper

class BrooksCorey(Soil):

    # class initialization
    def __init__(self,input: Input):
        self.logger =  logging_helper.get_logger("SOIL")
        self.logger.info("--- Using the Brooks-Corey model")
        # initialize parent class
        super().__init__(input)
    
    # Compute soil surface moisture
    def surface_water_content(self, psi: float) -> float:
        b        = self.properties.b[0]
        psi_sat  = self.properties.psi_sat[0]
        porosity = self.properties.porosity[0]
        residual = self.properties.residual[0]
        soil_e   = porosity-residual
        soil_q   = residual+soil_e*( (psi_sat/psi)**(1./b) )
        
        return soil_q
    
    # Compute soil water potential (column)
    def water_potential(self, soil_q: Union[float, np.ndarray], level: int = None) -> Union[float, np.ndarray]:
        """
        Computes soil water potential.
        
        :param soil_q: Soil moisture content (scalar or array).
        :param level: The soil level for a scalar calculation (optional).
        :return: Soil water potential (scalar or array).
        """
        if level is not None:
            b        = self.properties.b[level]
            psi_sat  = self.properties.psi_sat[level]
            porosity = self.properties.porosity[level]
            residual = self.properties.residual[level]
        else:
            b        = self.properties.b
            psi_sat  = self.properties.psi_sat
            porosity = self.properties.porosity
            residual = self.properties.residual
        
        Se  = (soil_q-residual)/(porosity-residual)
        psi = psi_sat*( Se**(-b) )
        
        return psi
        
    # Computes soil moisture conductivity.
    def conductivity_moisture(self, soil_q: Union[float, np.ndarray], level: int = None) -> Union[float, np.ndarray]:
        """
        Computes soil moisture conductivity.
        
        :param soil_q: Soil moisture content (scalar or array).
        :param level: The soil level for a scalar calculation (optional).
        :return: Soil moisture conductivity (scalar or array).
        """
        if level is not None:
            b            = self.properties.b[level]
            porosity     = self.properties.porosity[level]
            residual     = self.properties.residual[level]
            K_sat        = self.properties.K_sat[level]
        else:
            b            = self.properties.b
            porosity     = self.properties.porosity
            residual     = self.properties.residual
            K_sat        = self.properties.K_sat
        
        Se           = (soil_q-residual)/(porosity-residual)
        conductivity = K_sat*( Se**(2.*b+3.) )
        
        return conductivity
    
    # Computes soil moisture diffusivity
    def diffusivity_moisture(self, soil_q: np.ndarray) -> np.ndarray:
        b            = self.properties.b
        psi_sat      = self.properties.psi_sat
        porosity     = self.properties.porosity
        residual     = self.properties.residual
        K_sat        = self.properties.K_sat
        Se           = (soil_q-residual)/(porosity-residual)
        diffusivity  = -b*K_sat*psi_sat*( Se**(b+2.) ) / (porosity-residual)
        
        return diffusivity