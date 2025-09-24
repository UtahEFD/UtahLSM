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
from abc import ABC, abstractmethod
from ...util.io import logging_helper

class Surface(ABC):
    
    def __init__(self):
        self.logger =  logging_helper.get_logger("SFC")
    
    @staticmethod
    def get_model(key):
        
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
    
    # Abstract methods    
    @abstractmethod
    def fm(self, z_m: float, z_o: float, obl: float) -> float:
        """Computes the stability function for momentum."""
        raise NotImplementedError
    
    @abstractmethod
    def fh(self, z_s: float, z_t: float, obl: float) -> float:
        """Computes the stability function for heat."""
        raise NotImplementedError
