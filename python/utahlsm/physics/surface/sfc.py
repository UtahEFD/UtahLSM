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
"""Abstract base class for surface layer models in UtahLSM.

This module defines the interface for all surface layer physics schemes
used within the UtahLSM framework. It provides the `Surface` abstract base
class, which ensures that any concrete surface layer model implements the
necessary stability functions, and a factory function (`get_model`) for
creating instances of those models.
"""
from abc import ABC, abstractmethod
from ...util.io import logging_helper

class Surface(ABC):
    """Abstract base class for surface layer models.
    
    This class defines the standard interface for surface layer calculations
    and acts as a factory for creating specific model instances. It cannot
    be instantiated directly.
    
    Attributes:
        logger: A logger for this class.
    """
    def __init__(self):
        """Initializes the Surface base class."""
        self.logger =  logging_helper.get_logger("SFC")
   
    @staticmethod
    def get_model(key):
        """Factory method to select and instantiate a surface layer model.
        
        Args:
            key: An integer identifying the surface layer model to use.
        
        Returns:
            An instance of a concrete `Surface` subclass.
        
        Raises:
            SystemExit: If the provided `key` is not a valid model ID.
        """
        # import surface sub-classes
        from .sfc_most import SurfaceMOST
       
        # dictionary to map keys to classes
        SFC_MODELS = {
            1: SurfaceMOST,
        }
       
        # return class or throw error
        try:
            # look up model class from dictionary
            return SFC_MODELS[key]()
        except KeyError as e:
            self.logger.error("x"*62)
            self.logger.error(f"Namelist Error: {e} is an invalid surface model.")
            self.logger.error(f"Valid options are:")
            for k,v in SFC_MODELS.items():
                self.logger.error(f"\t{k} ({v.__name__})")
            self.logger.error("x"*62)
            raise SystemExit(1)
   
    #--- Abstract methods ---   
    @abstractmethod
    def fm(self, z_m: float, z_o: float, obl: float) -> float:
        """Computes the stability function for momentum.
        
        This is an abstract method that must be implemented by any concrete
        subclass.
        
        Args:
            z_m: Measurement height for wind speed [m].
            z_o: Aerodynamic roughness length [m].
            obl: Obukhov length [m].
        
        Returns:
            The dimensionless stability correction factor for momentum.
        """
        raise NotImplementedError
   
    @abstractmethod
    def fh(self, z_s: float, z_t: float, obl: float) -> float:
        """Computes the stability function for heat.
        
        This is an abstract method that must be implemented by any concrete
        subclass.
        
        Args:
            z_s: Measurement height for temperature and humidity [m].
            z_t: Thermal roughness length [m].
            obl: Obukhov length [m].
        
        Returns:
            The dimensionless stability correction factor for heat and scalars.
        """
        raise NotImplementedError
