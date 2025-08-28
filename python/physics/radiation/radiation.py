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
import sys

class Radiation(object):
    
    def __init__(self, input):
    
        self.latitude   = input.latitude
        self.longitude  = input.longitude
        self.albedo     = input.albedo
        self.emissivity = input.emissivity
        
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
            print("x"*62)
            print(f"Namelist Error: {key} is an invalid radiation model.")
            print(f"Valid options are:")
            for k,v in RAD_MODELS.items():
                print(f"\t{k} ({v.__name__})")
            print("x"*62)
            
            sys.exit(1)
