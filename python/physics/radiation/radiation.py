# 
# UtahLSM
# 
# Copyright (c) 2017–2023 Jeremy A. Gibbs
# Copyright (c) 2017–2023 Rob Stoll
# Copyright (c) 2017–2023 Eric Pardyjak
# Copyright (c) 2017–2023 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

class Radiation(object):
    
    def __init__(self, input):
    
        self.latitude   = input.latitude
        self.longitude  = input.longitude
        self.albedo     = input.albedo
        self.emissivity = input.emissivity
        
    @staticmethod
    def get_model(key, input):
        if key == 1:
            from .rad_basic import RadBasic
            return RadBasic(input)