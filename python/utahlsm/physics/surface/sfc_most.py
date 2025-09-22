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

import math
import numpy as np
from .sfc import Surface
from ...util import constants as c
from ...util.io import logging_helper

class SurfaceMOST(Surface):

    # class initialization
    def __init__(self):
        
        self.logger = logging_helper.get_logger("SFC")
        self.logger.info("Using the MOST model")
        # initialize parent class
        super().__init__()
        
    # gradient momentum function
    def phim(self,z, obukL):
        
        # z/L
        zeta = 0 if obukL==0 else z/obukL
        
        # compute gradient momentum function
        return self.phim_stable(zeta) if zeta >= 0 else self.phim_unstable(zeta)
    
    # gradient momentum function under neutral/stable conditions
    def phim_stable(self,zeta):
        return 1. + 5.*zeta
    
    # gradient momentum function under unstable conditions
    def phim_unstable(self,zeta):
        return (1.-(16.*zeta))**(-0.25)
    
    # gradient scalar function
    def phih(self,z, obukL):
        
        # z/L
        zeta = 0 if obukL==0 else z/obukL
        
        # compute gradient scalar function
        return self.phih_stable(zeta) if zeta >= 0 else self.phih_unstable(zeta)
        
    # gradient scalar function under neutral/stable conditions
    def phih_stable(self,zeta):
        return 1. + 5.*zeta
    
    # gradient scalar function under unstable conditions
    def phih_unstable(self,zeta):
        return (1.-(16.*zeta))**(-0.50)
    
    # integral stability correction for momentum
    def psim(self,z,obukL):
        
        # z/L
        zeta = 0 if obukL==0 else z/obukL
        
        # compute integral stability correction
        return self.psim_stable(zeta) if zeta >= 0 else self.psim_unstable(zeta)
    
    # integral stability correction for momentum under neutral/stable conditions
    def psim_stable(self,zeta):
        return -5.*zeta
    
    # integral stability correction for momentum under unstable conditions
    def psim_unstable(self,zeta):
        
        # local constants
        PI = c.physical.PI
        
        x = (1.-(16.*zeta))**(0.25)
        
        return 2.*np.log((1.+x)/2.)+np.log((1.+x**2.)/2.)-2.*math.atan2(1.,self.phim_unstable(zeta))+PI/2.
    
    # integral stability correction for scalars
    def psih(self,z,obukL):
        # z/L
        zeta = 0 if obukL==0 else z/obukL
        
        # compute integral stability correction
        return self.psih_stable(zeta) if zeta >= 0 else self.psih_unstable(zeta)
        
    # integral stability correction for scalars under neutral/stable conditions
    def psih_stable(self,zeta):
        return -5.*zeta
    
    # integral stability correction for scalars under unstable conditions
    def psih_unstable(self,zeta):
        x = (1.-(16.*zeta))**(0.50)
        
        return 2.*np.log((1.+x)/2.)
    
    # common log-law function for momentum
    def fm(self, z1, z0, obukL):
        
        # local constants
        VK = c.physical.VON_KARMAN
        
        fm = VK / (np.log(z1/z0) - self.psim(z1,obukL) + self.psim(z0,obukL))
        return fm
    
    # common log-law function for heat
    def fh(self, z1, z0h, obukL):
        
        # local constants
        VK = c.physical.VON_KARMAN
        
        fh = VK / (np.log(z1/z0h) - self.psih(z1,obukL) + self.psih(z0h,obukL))
        return fh
