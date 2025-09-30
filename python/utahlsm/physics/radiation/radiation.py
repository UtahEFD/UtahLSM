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
"""Abstract base class for radiation models in UtahLSM.

This module defines the interface for all radiation physics parameterizations
used within the UtahLSM framework. It provides the `Radiation` abstract base
class, which ensures that any concrete radiation model implements the necessary
methods, and a factory function (`get_model`) for creating instances of those
models.
"""
from abc import ABC, abstractmethod
from ...data_models import AtmosphericState, SurfaceState
from ...util.io import logging_helper

class Radiation(ABC):
    """Abstract base class for radiation models.
    
    This class defines the standard interface for radiation calculations
    and acts as a factory for creating specific radiation model instances.
    It cannot be instantiated directly.
    
    Attributes:
        logger: A logger for this class.
        latitude: The site latitude in degrees.
        longitude: The site longitude in degrees.
        albedo: The surface albedo (dimensionless).
        emissivity: The surface emissivity (dimensionless).
    """
    def __init__(self, latitude: float, longitude: float, albedo: float, emissivity: float):
        """Initializes the Radiation base class.
        
        Args:
            latitude: The site latitude in degrees.
            longitude: The site longitude in degrees.
            albedo: The surface albedo (dimensionless).
            emissivity: The surface emissivity (dimensionless).
        """
        self.logger = logging_helper.get_logger("RAD")
        self.latitude = latitude
        self.longitude = longitude
        self.albedo = albedo
        self.emissivity = emissivity
        
    @staticmethod
    def get_model(key: int, latitude: float, longitude: float, albedo: float, emissivity: float):
        """Factory method to select and instantiate a radiation model.
        
        Based on the integer key provided in the namelist, this method imports
        and returns an instance of the corresponding radiation model class.
        
        Args:
            key: An integer identifying the radiation model to use.
            latitude: The site latitude in degrees.
            longitude: The site longitude in degrees.
            albedo: The surface albedo (dimensionless).
            emissivity: The surface emissivity (dimensionless).
        
        Returns:
            An instance of a concrete `Radiation` subclass.
        
        Raises:
            SystemExit: If the provided `key` is not a valid model ID.
        """
        # import radiation sub-classes
        from .rad_basic import RadBasic
        
        # dictionary to map keys to classes
        RAD_MODELS = {
            1: RadBasic,
        }
            
        # return class or throw error
        try:
            # look up model class from dictionary
            return RAD_MODELS[key](latitude, longitude, albedo, emissivity)
        except KeyError as e:
            self.logger.error("x"*62)
            self.logger.error(f"Namelist Error: {key} is an invalid radiation model.")
            self.logger.error(f"Valid options are:")
            for k,v in RAD_MODELS.items():
                self.logger.error(f"\t{k} ({v.__name__})")
            self.logger.error("x"*62)
            raise SystemExit(1)
    
    # Abstract methods ---
    @abstractmethod
    def compute_net(self, julian_day: int, utc: float, atm_state: AtmosphericState, sfc_state: SurfaceState) -> float:
        """Computes the net radiation at the surface.
        
        This is an abstract method that must be implemented by any concrete
        subclass. It calculates the net radiation flux (shortwave and longwave)
        at the land surface.
        
        Args:
            julian_day: The current Julian day of the year.
            utc: The current time in UTC seconds from midnight.
            atm_state: The current state of the atmosphere.
            sfc_state: The current state of the surface.
        
        Returns:
            The net radiation in W/m^2.
        """
        raise NotImplementedError
