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
from data_models import AtmosphericState, SurfaceState
from .radiation import Radiation
from ...util import constants as c
from ...util.io import logging_helper

class RadBasic(Radiation):

    # class initialization
    def __init__(self, input):
        
        # local logger
        self.logger = logging_helper.get_logger("RAD: Basic")
        self.logger.info("Using the basic model")
        # initialize parent class
        super().__init__(input)
        
    # Computes the surface net radiation
    def compute_net(self, julian_day:int, time_utc:int, atm_state:AtmosphericState, sfc_state:SurfaceState):
        
        # local copies of constants
        latitude    = self.input.radiation.latitude
        longitude   = self.input.radiation.longitude
        albedo      = self.input.surface.albedo
        emissivity  = self.input.surface.emissivity
        
        # Compute incoming shortwave radiation
        sw_in = self.shortwave_in(julian_day, time_utc, latitude, longitude)
        
        # Compute outgoing shortwave radiation
        sw_out = self.shortwave_out(albedo, sw_in)
        
        # Compute outgoing longwave radiation
        lw_out = self.longwave_out(emissivity, sfc_state)
        
        # Compute incoming longwave radiation (current hack is net longwave of -50)
        lw_in = self.longwave_in(atm_state,sfc_state)

        # Compute net radiation
        return sw_in - sw_out + lw_in - lw_out
    
    # Computes the downward longwave radiation at the surface
    def longwave_in(self, atm_state, sfc_state):
        """simple clear-sky downwelling lw computation from Brutsaert (1975)"""
        
        # local references to atmospheric state
        pa = atm_state.p
        qa = sfc_state.qa
        Ts = sfc_state.Ts
        
        # vapor pressure
        ea = (pa*qa) / (c.epsilon + qa)
        
        # effective emissivity
        emissivity = 1.24*(ea/Ts)**(1/7.)
        
        # downward longwave
        lw_in = emissivity * c.sb * (Ts**4)
        
        return lw_in
            
    # Computes the upward longwave radiation at the surface
    def longwave_out(self, emissivity, sfc_state):
        Ts = sfc_state.Ts
        return emissivity * c.sb * (Ts**4)
    
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