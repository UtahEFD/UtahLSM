class Surface(object):
    
    def __init__(self,inputLSM, soil):
        self.input = inputLSM
        self.ustar = 0
        self.obukL = 0
        self.soil  = soil
        
    @staticmethod
    def get_model(key,inputLSM):
        if key == 1:
            from .sfc_most import SurfaceMOST
            return SurfaceMOST(inputLSM)
    
    # function to return the momentum and scalar fluxes
    def get_fluxes(self):
        return self.tau13, self.tau23, self.ustar, self.hflux
        
    # gradient momentum function
    def phim(self,z):
        zeta = np.zeros(self.obukLT.shape)
        mask = self.obukLT!=0
        if any(mask): zeta[mask] = z/self.obukLT[mask]
        for _ in zeta:
            if _ > 0:
                return self.phim_stable(_)
            else:
                return self.phim_unstable(_)
    
    # gradient momentum function under neutral/stable conditions
    def phim_stable(self,zeta):
        return 1 + 5*zeta
    
    # gradient momentum function under unstable conditions
    def phim_unstable(self,zeta):
        return (1-(16*zeta))**(-0.25)
    
    # gradient scalar function
    def phih(self,z):
        zeta = np.zeros(self.obukLT.shape)
        mask = self.obukLT!=0
        if any(mask): zeta[mask] = z/self.obukLT[mask]
        for _ in zeta:
            if _ > 0:
                return self.phih_stable(_)
            else:
                return self.phih_unstable(_)
        
    # gradient scalar function under neutral/stable conditions
    def phih_stable(self,zeta):
        return 1 + 5*zeta
    
    # gradient scalar function under unstable conditions
    def phih_unstable(self,zeta):
        return (1-(16*zeta))**(-0.50)
    
    # integral stability correction for momentum
    def psim(self,z,obukLT):
        zeta = np.zeros(self.obukLT.shape)
        psim = np.zeros(self.obukLT.shape)
        
        mask = self.obukLT!=0
        if any(mask): zeta[mask] = z/self.obukLT[mask]
        
        psim = np.where(zeta > 0, self.psim_stable(zeta), self.psim_unstable(zeta))
        return psim
    
    # integral stability correction for momentum under neutral/stable conditions
    def psim_stable(self,zeta):
        return -5*zeta
    
    # integral stability correction for momentum under unstable conditions
    def psim_unstable(self,zeta):
        x = (self.phim_unstable(zeta))**(-1)
        return 2*np.log((1+x)/2)+np.log((1+x**2)/2)-2*np.arctan(x)+np.pi/2
    
    # integral stability correction for scalars
    def psih(self,z,obukLT):
        zeta = np.zeros(self.obukLT.shape)
        psih = np.zeros(obukLT.shape)
        
        mask = self.obukLT!=0
        if any(mask): zeta[mask] = z/self.obukLT[mask]
        
        psih = np.where(zeta > 0, self.psih_stable(zeta), self.psih_unstable(zeta))
        
        return psih
        
    # integral stability correction for scalars under neutral/stable conditions
    def psih_stable(self,zeta):
        return -5*zeta
    
    # Integral stability correction for scalars under unstable conditions
    def psih_unstable(self,zeta):
        x = (self.phih_unstable(zeta))**(-1)
        return 2*np.log((1+x)/2)
    
    # common log-law function for momentum
    def fm(self, z1, z0, obukLT):
        vk = constants.vonk
        fm = np.zeros(obukLT.shape)
        fm = vk / (np.log(z1/z0) - self.psim(z1,obukLT) + self.psim(z0,obukLT))
        return fm
    
    # common log-law function for heat
    def fh(self, z1, z0h, obukLT):
        vk = constants.vonk
        fh = np.zeros(obukLT.shape)
        fh = vk / (np.log(z1/z0h) - self.psih(z1,obukLT) + self.psih(z0h,obukLT))
        
        return fh