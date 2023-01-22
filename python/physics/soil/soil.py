import numpy as np

from util import constants as c
from .soil_type import SoilType

class SOIL(object):
    
    def __init__(self,inputLSM):
        
        # TODO: limit what is passed in
        # TODO: implement error check on models
        
        self.input = inputLSM    
        nz = self.input.nsoil
        
        # fill properties
        self.properties = np.empty(nz).astype(np.object_)
        for k in range(0,nz):
            self.properties[k] = SoilType.get_properties(self.input.param,self.input.soil_type[k])
    
    @staticmethod
    def get_model(key,inputLSM):
        if key == 1:
            from .soil_brookscorey import BrooksCorey
            return BrooksCorey(inputLSM)
        if key == 2:
            from .soil_campbell import Campbell
            return Campbell(inputLSM)
        if key == 3:
            from .soil_vangenuchten import VanGenuchten
            return VanGenuchten(inputLSM)
    
    # Compute heat capacity
    def heatCapacity(self, soil_q, level):
        porosity = self.properties[level].porosity
        Ci       = self.properties[level].ci
        return (1-porosity)*Ci + soil_q*c.Ci_wat + (porosity-soil_q)*c.Cp_air
    
    # Compute surface mixing ratio
    def surfaceMixingRatio(self, sfc_T, sfc_q, atm_p):
        psi      = self.waterPotential(sfc_q, 0)
        h        = np.exp(c.grav*psi/(c.Rv*sfc_T))
        es       = 6.1078*np.exp(17.269*(sfc_T-273.15)/(sfc_T-35.86))
        hum_sat  = 0.622*(es/(atm_p-0.378*es))
        hum_spec = h*hum_sat
        return hum_spec
    
    # Compute soil thermal conductivity
    def conductivityThermal(soil_q, level):
        psi = 100.*self.waterPotential(soil_q,level)
        pf = np.log10(np.abs(psi))
        if (pf <= 5.1):
            conductivity = 418.46*np.exp(-(pf+2.7))
        else:
            conductivity = 0.172
        return conductivity
    
    # Compute soil thermal diffusivity
    def diffusivityThermal(soil_q, level):
        heat_cap     = self.heatCapacity(soil_q,level)
        conductivity = self.conductivityThermal(soil_q, level)
        diffusivity  = conductivity / heat_cap
        return diffusivity