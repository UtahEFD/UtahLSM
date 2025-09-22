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
from .soil_type import SoilType
from ...util import constants as c
from ...util.io import logging_helper

class Soil(object):
    
    def __init__(self,input):
        
        self.logger =  logging_helper.get_logger("SOIL")
        self.input  = input    
        nz          = self.input.grid.nz
        dataset     = self.input.soil.param
        
        DATASET_NAMES = {
            1: "Clapp/Hornberger",
            2: "Cosby et al",
            3: "Rawls/Brakensiek"
        }
        self.logger.info(f"Using the {DATASET_NAMES[dataset]} dataset")
        
        # fill properties
        self.properties = np.empty(nz).astype(np.object_)
        for k in range(0,nz):
            self.properties[k] = SoilType.get_properties(dataset,self.input.initial.type[k])
    
    @staticmethod
    def get_model(key,input):
        
        # import soil sub-classes
        from .soil_brookscorey import BrooksCorey
        from .soil_campbell import Campbell
        from .soil_vangenuchten import VanGenuchten
        
        # dictionary to map keys to classes
        SOIL_MODELS = {
            1: BrooksCorey,
            2: Campbell,
            3: VanGenuchten,
        }
        
        # return class or throw error
        try:
            # look up model class from dictionary
            return SOIL_MODELS[key](input)
        except KeyError as e:
            self.logger.error(e)
            raise
    
    # Compute heat capacity
    def heat_capacity(self, soil_q, level):
        
        # Local constants
        CI_W = c.water.SPECIFIC_HEAT
        CP_A = c.thermodynamic.SPECIFIC_HEAT
        
        porosity = self.properties[level].porosity
        Ci       = self.properties[level].ci
        Ks       = (1.-porosity)*Ci + soil_q*CI_W + (porosity-soil_q)*CP_A
        
        return Ks
        
    # Compute surface mixing ratio
    def surface_mixing_ratio(self, sfc_T, sfc_q, atm_p):
        
        # Local constants
        G  = c.physical.GRAVITY
        RV = c.thermodynamic.GAS_CONSTANT_VAPOR
         
        psi      = self.water_potential(sfc_q, 0)
        h        = np.exp(G*psi/(RV*sfc_T))
        es       = 6.1078*np.exp(17.269*(sfc_T-273.15)/(sfc_T-35.86))
        hum_sat  = 0.622*(es/(atm_p-0.378*es))
        hum_spec = h*hum_sat
         
        return hum_spec
    
    # Compute soil thermal conductivity
    def conductivity_thermal(self, soil_q, level):
        psi = 100.*self.water_potential(soil_q,level)
        pf = np.log10(np.abs(psi))
        if (pf <= 5.1):
            conductivity = 418.46*np.exp(-(pf+2.7))
        else:
            conductivity = 0.172
        return conductivity
    
    # Compute soil thermal diffusivity
    def diffusivity_thermal(self, soil_q, level):
        heat_cap     = self.heat_capacity(soil_q,level)
        conductivity = self.conductivity_thermal(soil_q, level)
        diffusivity  = conductivity / heat_cap
        return diffusivity