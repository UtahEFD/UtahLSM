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
logger = logging.getLogger("RAD")

class Radiation(object):
    
    def __init__(self, input):
    
        self.input = input
        
    @staticmethod
    def get_model(key, input):
        
        # import radiation sub-classes
        from .rad_basic import RadBasic
        
        # dictionary to map keys to classes
        RAD_MODELS = {
            1: RadBasic,
        }
            
        # return class or throw error
        try:
            # look up model class from dictionary
            return RAD_MODELS[key](input)
        except KeyError as e:
            logger.error("x"*62)
            logger.error(f"Namelist Error: {key} is an invalid radiation model.")
            logger.error(f"Valid options are:")
            for k,v in RAD_MODELS.items():
                logger.error(f"\t{k} ({v.__name__})")
            logger.error("x"*62)
            raise SystemExit(1)
