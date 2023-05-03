# 
# UtahLSM
# 
# Copyright (c) 2017–2023 Jeremy A. Gibbs
# Copyright (c) 2017–2023 Rob Stoll
# Copyright (c) 2017–2023 Eric Pardyjak
# Copyright (c) 2017–2023 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 
import sys
import numpy as np
from .soil import Soil
from util import constants as c

class Campbell(Soil):

    # class initialization
    def __init__(self,input):
        print("[UtahLSM: Soil] \t --- using the Campbell model with")
        # initialize parent class
        super().__init__(input)
    
    # Compute soil surface moisture
    def surface_water_content(self, psi):
        b        = self.properties[0].b
        psi_sat  = self.properties[0].psi_sat
        porosity = self.properties[0].porosity
        soil_q   = porosity*(np.abs(psi_sat/psi)**(1./b))
                
        return soil_q

    # Estimate soil surface moisture from surface mixing ratio
    def surface_water_content_estimate(self, sfc_T, sfc_q, atm_p):
        b        = self.properties[0].b
        psi_sat  = self.properties[0].psi_sat
        porosity = self.properties[0].porosity
        es       = 6.1078*np.exp(17.269*(sfc_T-273.15)/(sfc_T-35.86))
        qs       = 0.622*(es/(atm_p-0.378*es))
        ln       = np.log(sfc_q/qs)
        soil_q   = porosity*((c.Rv*sfc_T*ln/(c.grav*psi_sat))**(-1./b))
        return soil_q
    
    # Compute soil water potential (single level)
    def water_potential(self, soil_q, level):
        b        = self.properties[level].b
        psi_sat  = self.properties[level].psi_sat
        porosity = self.properties[level].porosity
        psi      = psi_sat*((soil_q/porosity)**(-b))
        if (psi>psi_sat): psi = psi_sat
        return psi
    
    # Computes soil moisture conductivity.
    def conductivity_moisture(self, soil_q, level):
        b            = self.properties[level].b
        porosity     = self.properties[level].porosity
        K_sat        = self.properties[level].K_sat
        conductivity = K_sat*( (soil_q/porosity)**(2.*b+3.) )
        
        return conductivity
    
    # Computes soil moisture diffusivity
    def diffusivity_moisture(self, soil_q, level):
        b            = self.properties[level].b
        psi_sat      = self.properties[level].psi_sat
        porosity     = self.properties[level].porosity
        K_sat        = self.properties[level].K_sat
        diffusivity  = -b*K_sat*psi_sat*( (soil_q/porosity)**(b+2.) ) / porosity
        
        return diffusivity