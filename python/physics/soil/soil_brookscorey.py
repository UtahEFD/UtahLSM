import sys
from .soil import SOIL

class BrooksCorey(SOIL):

    # class initialization
    def __init__(self,inputLSM):
        
        print("[UtahLSM: SOIL] \t--- running with the Brooks-Corey model")
        # initialize parent class
        super().__init__(inputLSM)
    
    # Compute soil surface moisture
    def surfaceWaterContent(psi):
        b        = self.properties[0].b
        psi_sat  = self.properties[0].psi_sat
        porosity = self.properties[0].porosity
        residual = self.properties[0].residual
        soil_e   = porosity-residual
        soil_q   = residual+soil_e*( (psi_sat/psi)**(1./b) )
        return soil_q
    
    # Estimate soil surface moisture from surface mixing ratio
    def surfaceWaterContentEstimate(sfc_T, sfc_q, atm_p):
        b        = self.properties[0].b
        psi_sat  = self.properties[0].psi_sat
        porosity = self.properties[0].porosity
        residual = self.properties[0].residual
        es       = 6.1078*np.exp(17.269*(sfc_T-273.15)/(sfc_T-35.86))
        qs       = 0.622*(es/(atm_p-0.378*es))
        ln       = np.log(sfc_q/qs)
        soil_e   = porosity-residual
        soil_q   = residual+soil_e*( ( (c.Rv*sfc_T*ln) / (c.grav*psi_sat) )**(-1./b) )
        return soil_q
    
    # Compute soil water potential (single level)
    def waterPotential(soil_q, level):
        b        = self.properties[level].b
        psi_sat  = self.properties[level].psi_sat
        porosity = self.properties[level].porosity
        residual = self.properties[level].residual
        Se       = (soil_q-residual)/(porosity-residual)
        psi      = psi_sat*( Se**(-b) )
        if (psi>psi_sat):
            psi = psi_sat
        # TODO: Add checks on moisture
        # if (Se<0):
        #     sys.exit(0)
        # 
        # if (Se > 1):
        #     // std::cout<<"NOPE WALK THAT BACK"<<std::endl
        #     // std::cout<<Se<<std::endl
        #     // std::exit(0)
        return psi
    
    # Computes soil moisture conductivity.
    def conductivityMoisture(soil_q, level):
        b            = self.properties[level].b
        porosity     = self.properties[level].porosity
        residual     = self.properties[level].residual
        K_sat        = self.properties[level].K_sat
        Se           = (soil_q-residual)/(porosity-residual)
        conductivity = K_sat*( Se**(2.*b+3.) )
        return conductivity
    
    # Computes soil moisture diffusivity
    def diffusivityMoisture(soil_q, level):
        b            = self.properties[level].b
        psi_sat      = self.properties[level].psi_sat
        porosity     = self.properties[level].porosity
        residual     = self.properties[level].residual
        K_sat        = self.properties[level].K_sat
        Se           = (soil_q-residual)/(porosity-residual)
        diffusivity  = -b*K_sat*psi_sat*( Se**(b+2.) ) / (porosity-residual)
        return diffusivity