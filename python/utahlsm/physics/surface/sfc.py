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
import logging

# local logger
logger = logging.getLogger("SFC")

class Surface(object):
    
    def __init__(self):pass
    
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
            logger.error("x"*62)
            logger.error(f"Namelist Error: {e} is an invalid surface model.")
            logger.error(f"Valid options are:")
            for k,v in SFC_MODELS.items():
                logger.error(f"\t{k} ({v.__name__})")
            logger.error("x"*62)
            raise SystemExit(1)
