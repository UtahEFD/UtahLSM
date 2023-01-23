import numpy as np
from util import constants as c
from .radiation import Radiation

class RadBasic(Radiation):

    # class initialization
    def __init__(self,inputLSM):
        
        print("[UtahLSM: RAD] \t\t --- running with the basic model")
        # initialize parent class
        super().__init__(inputLSM)
        
    # Computes the surface net radiation
    def compute_net(self, julian_day, time_utc, sfc_T):
        
        # Compute incoming shortwave radiation
        sw_in = self.shortwave_in(julian_day, time_utc, self.latitude, self.longitude)
        
        # Compute outgoing shortwave radiation
        sw_out = self.shortwave_out(self.albedo, sw_in)
        
        # Compute outgoing longwave radiation
        lw_out = self.longwave_out(self.emissivity, sfc_T)
        
        # Compute incoming longwave radiation (current hack is net longwave of -50)
        lw_in = lw_out - 50.0
        
        # Compute net radiation
        return sw_in - sw_out + lw_in - lw_out
    
    # Computes the downward longwave radiation at the surface
    def longwave_in(self):
        return 0
    
    # Computes the upward longwave radiation at the surface
    def longwave_out(self, emissivity, sfc_T):
        return emissivity * c.sb * (sfc_T**4)
    
    # Computes the downward shortwave radiation at the surface
    def shortwave_in(self, julian_day, time_utc, latitude, longitude):
        sw_in = 0
        declination = 23.45*(c.pi/180.0)*np.cos(2.0*c.pi*(julian_day-173)/365.25)
        sin_elevation = np.sin(latitude)*np.sin(declination) - np.cos(latitude)*np.cos(declination)* np.cos( (2*c.pi*time_utc/(24.0*3600.0)) - longitude )
        if (sin_elevation > 0):
            transmissivity = (0.6 + 0.2*sin_elevation)
            sw_in = c.sc * transmissivity * sin_elevation
        return sw_in
        
    # Computes the upward shortwave radiation at the surface
    def shortwave_out(self, albedo, sw_in):
        return albedo * sw_in