import numpy as np
from util import constants as c
from .sfc import Surface

class SurfaceMOST(Surface):

    # class initialization
    def __init__(self,inputLSM, soil):
        
        print("[UtahLSM: Surface] \t\t --- running with the MOST scheme")
        # initialize parent class
        super().__init__(inputLSM)
        
    # main run loop
    def run(self, utc, atm_U, atm_T, atm_p, R_net, sfc_T, sfc_q, soil_T, soil_q, soil_z):
        
        # config variables
        max_iterations = 200
        converged = False
        last_L = 0.1
        criteria = 0.1
        ref_T = 300
        z_m = self.input.z_m
        z_s = self.input.z_s
        z_o = self.input.z_o
        z_t = self.input.z_t
        
        # Compute surface mixing ratio
        gnd_q  = self.soil.surface_mixing_ratio(sfc_T,sfc_q,atm_p)
        
    #     # loop to solve for sensible flux, latent flux, ustar, and L
    #     for _ in range (0,max_iterations):
    #         
    #         # Compute ground flux
    #         # First time through we estimate based on Santanello and Friedl (2003)
    #         if ( (first)):
    #             if (soil_q[0]>=0.4):
    #                 A = 0.31
    #                 B = 74000
    #             elif (soil_q[0]<0.4 and soil_q[0] >= 0.25):
    #                 A = 0.33
    #                 B = 85000
    #             else:
    #                 A = 0.35
    #                 B = 100000
    #             flux_gr = R_net*A*np.cos((2*c.pi*(utc)+10800)/B)
    #         else:
    #             K0 = soil.conductivity_thermal(soil_q[0],0)
    #             K1 = soil.conductivity_thermal(soil_q[1],1)
    #             Kmid = 0.5*(K0 + K1)
    #             flux_gr = Kmid*(sfc_T - soil_T[1])/(soil_z[0]-soil_z[1]);
    #         
    #         # Compute friction velocity
    #         ustar = atm_U*self.fm(z_m, z_o,zeta_m,zeta_o)
    #         
    #         # Compute heat flux
    #         flux_wT = (sfc_T-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t)
    #         
    #     dz    = self.input.dz
    #     z0    = self.input.z0
    #     z0h   = self.input.z0h
    #     z1    = dz / 2
    #     T_0   = self.input.Tref
    #     g     = constants.grav
    #     vk    = constants.vonk
    #     zf    = self.input.zf
    #     
    #     # surface information
    #     Ms = np.sqrt(U**2 + V**2)
    #     
    #     # loop to solve for L and ustar
    #     for _ in range (0,4):
    # 
    #         # compute fluxes and Obukhov length
    #         self.ustar = Ms*self.fm(z1,z0)
    #         self.hflux = (self.Ts[tidx]-T)*self.ustar*self.fh(z1,z0h)
    #         
    #         # protect against purely neutral conditions
    #         if self.hflux == 0:
    #             self.obukL = 0
    #         else:
    #             self.obukL = -self.ustar**3*T_0/(vk*g*self.hflux)
    #         
    #     self.tau13 = -(self.ustar**2)*(U/Ms)
    #     self.tau23 = -(self.ustar**2)*(V/Ms)
    #     if self.obukL == 0: 
    #         self.zeta = 0
    #     else:
    #        self.zeta = zf/self.obukL
    #       
    #     