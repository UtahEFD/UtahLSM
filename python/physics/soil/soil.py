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

import logging
import numpy as np

from util import constants as c
from .soil_type import SoilType

# local logger
logger = logging.getLogger("SOIL")

class Soil(object):
    
    def __init__(self,input):
        
        self.input = input    
        nz         = self.input.nsoil
        dataset    = self.input.soil_param
        
        if dataset==1:
            set_name = "Clapp/Hornberger"
        elif dataset==2:
            set_name = "Cosby et al"
        elif dataset==3:
            set_name="Rawls/Brakensiek"
        
        logger.info("--- the %s dataset"%set_name)
        
        # fill properties
        self.properties = np.empty(nz).astype(np.object_)
        for k in range(0,nz):
            self.properties[k] = SoilType.get_properties(dataset,self.input.soil_type[k])
    
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
        
        # look up model class from dictionary
        model_class = SOIL_MODELS.get(key)
        
        # return class or throw error
        try:
            # look up model class from dictionary
            return SOIL_MODELS[key](input)
        except KeyError as e:
            logger.error("x"*62)
            logger.error(f"Namelist Error: {e} is an invalid soil model.")
            logger.error(f"Valid options are:")
            for k,v in SOIL_MODELS.items():
                logger.error(f"\t{k} ({v.__name__})")
            logger.error("x"*62)
            raise SystemExit(1)
    
    # Compute heat capacity
    def heat_capacity(self, soil_q, level):
        porosity = self.properties[level].porosity
        Ci       = self.properties[level].ci
        Ks       = (1.-porosity)*Ci + soil_q*c.Ci_wat + (porosity-soil_q)*c.Cp_air
        
        return Ks

        
    # Compute surface mixing ratio
    def surface_mixing_ratio(self, sfc_T, sfc_q, atm_p):
        psi      = self.water_potential(sfc_q, 0)
        h        = np.exp(c.grav*psi/(c.Rv*sfc_T))
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