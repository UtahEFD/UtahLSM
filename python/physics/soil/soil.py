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
            from .soil_brookscorey import BROOKSCOREY
            return BROOKSCOREY(inputLSM)
        if key == 2:
            from .soil_campbell import CAMPBELL
            return CAMPBELL(inputLSM)
        if key == 3:
            from .soil_vangenuchten import VANGENUCHTEN
            return VANGENUCHTEN(inputLSM)
    
    # Soil::Soil(const std::vector<int>& soil_type, const int soil_param, const int soil_model, const int levels) {
    #     
    #     // Set up soil properties in the column
    #     properties.resize(levels);
    #     for (int i=0; i<levels; i++) {
    #         properties[i] = SoilType::getProperties(soil_type[i],soil_param);
    #     }
    # };
    
    
    # def heatCapacity(self, soil_q, level):
    #     porosity = self.properties[level].porosity
    #     Ci       = self.properties[level].ci
    #     return (1-porosity)*Ci + soil_q*c.Ci_wat + (porosity-soil_q)*c.Cp_air;
    
    # // Compute heat capacity
    # double Soil::heatCapacity(const double soil_q, const int level) {
    # 
    #     double porosity = properties[level]->porosity;
    #     double Ci = properties[level]->ci;
    # 
    #     return (1-porosity)*Ci + soil_q*c::Ci_wat + (porosity-soil_q)*c::Cp_air;
    # }
    
    
    # // Compute surface mixing ratio
    # double Soil::surfaceMixingRatio(const double sfc_T, const double sfc_q,
    #                                 const double atm_p) {
    #     
    #     double psi      = waterPotential(sfc_q, 0);
    #     double h        = std::exp(c::grav*psi/(c::Rv*sfc_T));
    #     double es       = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
    #     double hum_sat  = 0.622*(es/(atm_p-0.378*es));
    #     double hum_spec = h*hum_sat;
    #     //std::cout<<"IN surfaceMixingRatio"<<std::endl;
    #     return hum_spec;
    # }
    
    
    # // Compute soil thermal conductivity
    # double Soil::conductivityThermal(const double soil_q, const int level) {
    # 
    #     double conductivity;
    #     double psi = 100.*waterPotential(soil_q,level);
    #     double pf = std::log10(std::abs(psi));
    #     if (pf <= 5.1) {
    #         conductivity = 418.46*std::exp(-(pf+2.7));
    #     } else {
    #         conductivity = 0.172;
    #     }
    #     //std::cout<<"IN conductivityThermal"<<std::endl;
    #     return conductivity;
    # }
    
    
    
    # // Compute soil thermal diffusivity
    # double Soil::diffusivityThermal(const double soil_q, const int level) {
    # 
    #     double heat_cap     = heatCapacity(soil_q,level);
    #     double conductivity = conductivityThermal(soil_q, level);
    #     double diffusivity  = conductivity / heat_cap;
    #     //std::cout<<"IN diffusivityThermal"<<std::endl;
    #     return diffusivity;
    # }