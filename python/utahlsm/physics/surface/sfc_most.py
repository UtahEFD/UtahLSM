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
"""Surface layer parameterization based on Monin-Obukhov Similarity Theory.

This module provides an implementation of the `Surface` abstract base class
using standard Monin-Obukhov Similarity Theory (MOST) functions to describe
the stability and flux-profile relationships in the atmospheric surface layer.
"""
import math
import logging
import numpy as np
from .sfc import Surface
from ...util import constants as c
from ...util.io import logging_helper

class SurfaceMOST(Surface):
    """Implements surface layer physics using MOST.
    
    This class provides concrete implementations for the stability correction
    functions for momentum and heat based on the widely used Businger-Dyer
    relations.
    """
    def __init__(self):
        """Initializes the SurfaceMOST model."""
        self.logger: logging.Logger = logging_helper.get_logger("SFC")
        self.logger.info("Using the MOST model")
        super().__init__()
        
    def phim(self, z: float, obukL: float) -> float:
        """Computes the dimensionless stability function for momentum (phi_m).
        
        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].
        
        Returns:
            The value of phi_m.
        """
        obukL_min = 0.1
        obukL_mag = max(abs(obukL), obukL_min)
        obukL_cap = np.copysign(obukL_mag, obukL)
        
        zeta = z / (obukL_cap)
        return self.phim_stable(zeta) if zeta >= 0 else self.phim_unstable(zeta)
    
    def phim_stable(self,zeta: float) -> float:
        """Computes phi_m for stable conditions (zeta >= 0)."""
        return 1. + 5.*zeta
    
    def phim_unstable(self,zeta: float) -> float:
        """Computes phi_m for unstable conditions (zeta < 0)."""
        return (1.-(16.*zeta))**(-0.25)
    
    def phih(self,z: float, obukL: float) -> float:
        """Computes the dimensionless stability function for heat (phi_h).
        
        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].
        
        Returns:
            The value of phi_h.
        """
        obukL_min = 0.1
        obukL_mag = max(abs(obukL), obukL_min)
        obukL_cap = np.copysign(obukL_mag, obukL)
        
        zeta = z / (obukL_cap)
        return self.phih_stable(zeta) if zeta >= 0 else self.phih_unstable(zeta)
        
    def phih_stable(self,zeta: float) -> float:
        """Computes phi_h for stable conditions (zeta >= 0)."""
        return 1. + 5.*zeta
    
    def phih_unstable(self,zeta:float) -> float:
        """Computes phi_h for unstable conditions (zeta < 0)."""
        return (1.-(16.*zeta))**(-0.50)
    
    def psim(self,z: float,obukL: float) -> float:
        """Computes the integrated stability function for momentum (psi_m).
        
        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].
        
        Returns:
            The value of psi_m.
        """
        obukL_min = 0.1
        obukL_mag = max(abs(obukL), obukL_min)
        obukL_cap = np.copysign(obukL_mag, obukL)
        
        zeta = z / (obukL_cap)
        return self.psim_stable(zeta) if zeta >= 0 else self.psim_unstable(zeta)
    
    def psim_stable(self,zeta: float) -> float:
        """Computes psi_m for stable conditions (zeta >= 0)."""
        return -5.*zeta
    
    def psim_unstable(self,zeta: float) -> float:
        """Computes psi_m for unstable conditions (zeta < 0)."""
        PI = c.physical.PI
        x = (1.-(16.*zeta))**(0.25)
        return 2.*np.log((1.+x)/2.)+np.log((1.+x**2.)/2.)-2.*math.atan2(1.,self.phim_unstable(zeta))+PI/2.
    
    def psih(self,z: float,obukL: float) -> float:
        """Computes the integrated stability function for heat (psi_h).
        
        Args:
            z: Height above the surface [m].
            obukL: Obukhov length [m].
        
        Returns:
            The value of psi_h.
        """
        obukL_min = 0.1
        obukL_mag = max(abs(obukL), obukL_min)
        obukL_cap = np.copysign(obukL_mag, obukL)
        
        zeta = z / (obukL_cap)
        return self.psih_stable(zeta) if zeta >= 0 else self.psih_unstable(zeta)
        
    def psih_stable(self,zeta: float) -> float:
        """Computes psi_h for stable conditions (zeta >= 0)."""
        return -5.*zeta
    
    def psih_unstable(self,zeta: float) -> float:
        """Computes psi_h for unstable conditions (zeta < 0)."""
        x = (1.-(16.*zeta))**(0.50)
        return 2.*np.log((1.+x)/2.)
    
    def fm(self, z1: float, z0: float, obukL: float) -> float:
        """Computes the log-law stability function for momentum.
        
        Args:
            z1: Upper height [m].
            z0: Lower height (roughness length) [m].
            obukL: Obukhov length [m].
        
        Returns:
            The stability-corrected log-law function value.
        """
        VK = c.physical.VON_KARMAN
        fm = VK / (np.log(z1/z0) - self.psim(z1,obukL) + self.psim(z0,obukL))
        return fm
    
    def fh(self, z1: float, z0h: float, obukL: float) -> float:
        """Computes the log-law stability function for heat.
        
        Args:
            z1: Upper height [m].
            z0h: Lower height (thermal roughness length) [m].
            obukL: Obukhov length [m].
        
        Returns:
            The stability-corrected log-law function value.
        """
        VK = c.physical.VON_KARMAN
        fh = VK / (np.log(z1/z0h) - self.psih(z1,obukL) + self.psih(z0h,obukL))
        return fh
