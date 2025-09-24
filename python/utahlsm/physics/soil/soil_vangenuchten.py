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

class VanGenuchten(Soil):

    # class initialization
    def __init__(self,input: Input):
        
        self.logger =  logging_helper.get_logger("SOIL")
        self.logger.info("[UtahLSM: Soil] \tUsing the Van Genuchten model")
        
        # initialize parent class
        super().__init__(input)
        
    # Compute soil surface moisture
    def surface_water_content(self, psi: float) -> float:
        b        = self.properties.b[0]
        psi_sat  = self.properties.psi_sat[0]
        porosity = self.properties.porosity[0]
        residual = self.properties.residual[0]
        soil_e   = porosity-residual
        m        = 1 / (1+b)
        soil_q   = residual + soil_e * (1 / ( 1 + (psi/psi_sat)**(1/(1-m)) ))**(m)
        
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
        
        Se       = (soil_q-residual)/(porosity-residual)
        m        = 1 / (1+b)
        psi      = psi_sat*( ( (Se**(-1/m))-1 )**(1-m) )
        
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
        m            = 1 / (1+b)
        conductivity = K_sat*np.sqrt(Se)*( (1 - (1 - (Se**(1/m)) )**m )**2 )
        
        return conductivity
    
    # Computes soil moisture diffusivity
    def diffusivity_moisture(self, soil_q: np.ndarray) -> np.ndarray:
        b            = self.properties.b
        psi_sat      = self.properties.psi_sat
        porosity     = self.properties.porosity
        residual     = self.properties.residual
        K_sat        = self.properties.K_sat
        Se           = (soil_q-residual)/(porosity-residual)
        soil_e       = porosity-residual
        m            = 1 / (1+b)
        A            = (1-m)*K_sat*psi_sat / (m*soil_e)
        C            = Se**(0.5-(1/m))*( (1 - Se**(1/m))**(-m) + (1- Se**(1/m))**m - 2 )
        diffusivity  = A*C
        
        return diffusivity