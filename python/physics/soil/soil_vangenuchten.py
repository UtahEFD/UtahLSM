from .soil import SOIL

class VanGenuchten(SOIL):

    # class initialization
    def __init__(self,inputLSM):
        
        # TODO: Check these computations, especially power terms
        
        print("[UtahLSM: SOIL] \t--- running with the Van Genuchten scheme")
        # initialize parent class
        super().__init__(inputLSM)
        
        # Compute soil surface moisture
        def surface_water_content(psi):
            b        = self.properties[0].b
            psi_sat  = self.properties[0].psi_sat
            porosity = self.properties[0].porosity
            residual = self.properties[0].residual
            soil_e   = porosity-residual
            m        = 1 / (1+b)
            soil_q   = residual + soil_e* ( ( (1 + ( (psi/psi_sat)**(1/(1-m)))  ))**(-m) )
            return soil_q
        
        # Estimate soil surface moisture from surface mixing ratio
        def surface_water_content_estimate(sfc_T, sfc_q, atm_p):
            b        = self.properties[0].b
            psi_sat  = self.properties[0].psi_sat
            porosity = self.properties[0].porosity
            residual = self.properties[0].residual
            es       = 6.1078*np.exp(17.269*(sfc_T-273.15)/(sfc_T-35.86))
            qs       = 0.622*(es/(atm_p-0.378*es))
            ln       = np.log(sfc_q/qs)
            soil_e   = porosity-residual
            m        = 1 / (1+b)
            soil_q   = residual+soil_e*( (1+ ( ((c.Rv*sfc_T*ln)/(c.grav*psi_sat))**(1/(1-m)) ))**(-m) )
            return soil_q
        
        # Compute soil water potential (single level)
        def water_potential(soil_q, level):
            b        = self.properties[level].b
            psi_sat  = self.properties[level].psi_sat
            porosity = self.properties[level].porosity
            residual = self.properties[level].residual
            Se       = (soil_q-residual)/(porosity-residual)
            m        = 1 / (1+b)
            psi      = psi_sat*( ( (Se**(-1/m))-1 )**(1-m) )
            if (psi>psi_sat): psi = psi_sat
            return psi
        
        # Computes soil moisture conductivity.
        def conductivity_moisture(soil_q, level):
            b            = self.properties[level].b
            porosity     = self.properties[level].porosity
            residual     = self.properties[level].residual
            K_sat        = self.properties[level].K_sat
            Se           = (soil_q-residual)/(porosity-residual)
            m            = 1 / (1+b)
            conductivity = K_sat*np.sqrt(Se)*( (1 - (1 - (Se**(1/m)) )**m )**2 )
            return conductivity
        
        # Computes soil moisture diffusivity
        def diffusivity_moisture(soil_q, level):
            b            = self.properties[level].b
            psi_sat      = self.properties[level].psi_sat
            porosity     = self.properties[level].porosity
            residual     = self.properties[level].residual
            K_sat        = self.properties[level].K_sat
            Se           = (soil_q-residual)/(porosity-residual)
            soil_e       = porosity-residual
            m            = 1 / (1+b)
            A            = (1-m)*K_sat*psi_sat / (m*soil_e)
            C            = (Se**(0.5-1/m)*( ( (1-( Se**(1/m))**(-m)) ) + ((1-( Se**(1/m)))**m)) - 2 )
            diffusivity  = A*C
            return diffusivity