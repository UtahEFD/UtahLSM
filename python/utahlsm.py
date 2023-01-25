#!/usr/bin/env python
import argparse
import numpy as np
import os
import sys
import time
from physics import Radiation, Soil, Surface
from util import constants as c, io, matrix

# custom error message for user case entry
class InvalidCase(Exception):
    pass

# land-surface model class
class UtahLSM:
    """This is the main UtahLSM class
    
    :param inputLSM: A handle to the :class:`util.io.Input`
    :param outputLSM: A handle to the :class:`util.io.Output`
    """
    
    # model class initialization
    def __init__(self,inputLSM, outputLSM, ustar, flux_wT, flux_wq):
        """Constructor method
        """
        # TODO: finish writing main init function
        # set the input and output fields
        self.input  = inputLSM
        self.output = outputLSM
        
        print("[UtahLSM: Setup] \t Reading input settings")
        # Input time section
        self.dt_seb = inputLSM.step_seb
        self.dt_dif = inputLSM.step_dif
        
        # Input grid section
        self.nx = inputLSM.nx
        self.ny = inputLSM.ny
        
        # Input length scale section
        self.z_o = self.input.z_o
        self.z_t = self.input.z_t
        self.z_m = self.input.z_m
        self.z_s = self.input.z_s
        
        # Input soil section
        self.nz             = self.input.nsoil
        self.soil_param     = self.input.param
        self.soil_model     = self.input.model
        
        print("[UtahLSM: Setup] \t Reading input data")
        self.soil_z    = self.input.soil_z
        self.soil_T    = self.input.soil_T
        self.soil_q    = self.input.soil_q
        self.soil_type = self.input.soil_type
        
        # Initialize new surface values for first run
        self.sfc_T_new = self.soil_T[0]
        self.sfc_q_new = self.soil_q[0]
        
        # Initialize history arrays for first run
        soil_T_last = self.soil_T
        soil_q_last = self.soil_q
        
        # Modify soil levels to be negative
        self.soil_z = -1*self.soil_z
        
        # Input radiation section
        self.comp_rad = self.input.comp_rad
        if (self.comp_rad==1):
            print("[UtahLSM: Setup] \t Creating radiation model")
            self.albedo     = self.input.albedo
            self.emissivity = self.input.emissivity
            self.latitude   = self.input.latitude
            self.longitude  = self.input.longitude
            self.julian_day = self.input.julian_day
            
            # convert latitude and longitude into radians
            self.latitude  = self.latitude * constants.pi / 180.0
            self.longitude = self.longitude * constants.pi / 180.0
            
            # Create radiation model
            self.rad = Radiation.get_model(1,self.input)
        else:
            print("[UtahLSM: Radiation] \t Using offline data, no model")
        
        # Create soil model
        print("[UtahLSM: Setup] \t Creating soil model")
        self.soil = Soil.get_model(self.soil_model,self.input)
        
        print("[UtahLSM: Setup] \t Creating surface model")
        # choose surface model
        self.sfc = Surface.get_model(1)

        print("[UtahLSM: Setup] \t Creating additional fields")
        # initialize flux arrays
        self.ust     = np.array([ustar])
        self.flux_wT = np.array([flux_wT])
        self.flux_wq = np.array([flux_wq])
        self.obl     = np.zeros(1)
        self.shf     = np.zeros(1)
        self.lhf     = np.zeros(1)
        self.ghf     = np.zeros(1)
        
        # Local atmospheric data
        self.atm_U = 0
        self.atm_T = 0
        self.atm_q = 0
        self.atm_p = 0
        self.R_net = 0
        
        # Local time data
        self.first      = True # flag whether first time step or not
        self.step_count = 0    # number of times the LSM has been called
        self.tstep      = 0    # current time step
        self.runtime    = 0    # current elapsed time
        self.utc        = 0    # current time in UTC

        print("[UtahLSM: Setup] \t Creating output file")    
        # set reference to output dimensions
        self.output_dims = {
            't':0,
            'z':self.nz
        }
        self.output.set_dims(self.output_dims)

        # set reference to output fields
        self.output_fields = {
            'soil_z':self.soil_z,
            'soil_T':self.soil_T,
            'soil_q':self.soil_q,
            'soil_type':self.soil_type,
            'ust':self.ust,
            'obl':self.obl,
            'shf':self.shf,
            'lhf':self.lhf,
            'ghf':self.ghf
        }
        self.output.set_fields(self.output_fields)

        # write initial data
        self.output.save(self.output_fields,0,0,initial=True)
        
    # Update atmospheric quantities prior to solving
    def update_fields(self,dt, u, T, q, p, rad=0):
        self.tstep    = dt
        self.atm_U    = u
        self.atm_T    = T
        self.atm_q    = q
        self.atm_p    = p
        self.runtime += tstep
        
        # Run radiation model and update time/date if needed
        if (self.comp_rad==1):
            utc         = np.fmod(self.runtime,86400)
            julian_day += int(runtime/86400);
            self.R_net  = self.rad.computeNet(julian_day,utc,self.soil_T[0])
        else:
            self.R_net = rad
        
        # Keep winds from being exactly zero
        if (self.atm_U==0): self.atm_U = 0.1
            
    # Run the model
    def run(self):
                
        # Set initial new temp and moisture
        self.sfc_T_new = self.soil_T[0];
        self.sfc_q_new = self.soil_q[0];
        
        # Check if time to re-compute balances
        if ( (self.step_count % self.dt_seb)==0 ):
            pass
            #self.solveSEB()
            #self.solveSMB()
        else:
            # just return new fluxes
            self.compute_fluxes(self.soil_T[0],self.soil_q[0])
        
        # Save current temperature and moisture
        self.soil_T_last = self.soil_T
        self.soil_q_last = self.soil_q
        
        # check if time to compute diffusion
        if ( (self.step_count % self.dt_dif)==0 ):
            
            # Solve heat diffusion
            self.solve_diffusion_heat()
        
            # solve moisture diffusion
            self.solve_diffusion_mois()
        
        # Change flag of whether initial time
        if self.first: self.first = False
        
        # Increment step counter
        self.step_count += 1
        
    # Save output fields
    # TODO: write save function
    def save(): pass
    
    # Compute fluxes using similarity theory
    # TODO: Clean compute_fluxes up
    def compute_fluxes(self, sfc_T, sfc_q):
        
        # Local variables
        max_iterations = 200
        converged      = False
        last_L         = 0.1
        criteria       = 0.1
        ref_T          = 300
        
        # Compute surface mixing ratio
        gnd_q  = self.soil.surface_mixing_ratio(sfc_T,sfc_q,self.atm_p)
        
        # Sensible flux, latent flux, ustar, and L
        for i in range(0,max_iterations):
            # First time through we estimate based on Santanello and Friedl (2003)
            if self.first:
                if (self.soil_q[0]>=0.4):
                    A = 0.31
                    B = 74000
                elif (self.soil_q[0]<0.4 and self.soil_q[0] >= 0.25):
                    A = 0.33
                    B = 85000
                else:
                    A = 0.35
                    B = 100000
                self.ghf = self.R_net*A*np.cos((2*c.pi*(utc)+10800)/B)
            else:
                K0       = self.soil.conductivity_thermal(self.soil_q[0],0)
                K1       = self.soil.conductivity_thermal(self.soil_q[1],1)
                Kmid     = 0.5*(K0 + K1)
                self.ghf[:] = Kmid*(sfc_T - self.soil_T[1])/(self.soil_z[0]-self.soil_z[1])
            
            # Compute friction velocity
            self.ust[:] = self.atm_U*self.sfc.fm(self.z_m, self.z_o, self.obl[:])
            
            # Compute heat flux
            self.flux_wT[:] = (sfc_T-self.atm_T)*self.ust[:]*self.sfc.fh(self.z_s, self.z_t, self.obl[:])
            
            # Compute latent flux
            if ( (self.first) and (i == 0)):
                self.flux_wq[:] = (self.R_net - self.ghf[:] - self.flux_wT[:]*c.rho_air*c.Cp_air)/(c.rho_air*c.Lv)
                gnd_q = self.atm_q + self.flux_wq[:] / (self.ust[:]*self.sfc.fh(self.z_s,self.z_t,self.obl[:]))
                self.soil_q[0] = self.soil.surface_water_content_estimate(self.soil_T[0],gnd_q, self.atm_p)
                self.sfc_q_new = self.soil_q[0]
            else:
                self.flux_wq[:] = (gnd_q-self.atm_q)*self.ust[:]*self.sfc.fh(self.z_s,self.z_t,self.obl[:])
            
            # Compute virtual heat flux
            flux_wTv = self.flux_wT[:] + ref_T*0.61*self.flux_wq[:]
            
            # Compute L
            last_L = self.obl[:]
            self.obl[:] = -(self.ust[:]**3)*ref_T/(c.vonk*c.grav*flux_wTv)
            
            # Bounds check on L
            if (self.z_m/self.obl[:] > 5.):  self.obl[:] = 5.
            if (self.z_m/self.obl[:] < -5.): self.obl[:] = 5.

            # Check for convergence
            converged = np.abs(last_L-self.obl[:]) <= criteria
            if (converged):
                self.shf[:] = c.rho_air*c.Cp_air*self.flux_wT[:]
                self.lhf[:] = c.rho_air*c.Lv*self.flux_wq[:]
                break
    
    # Compute the surface energy budget
    def compute_seb(self, sfc_T):
        
        # Compute fluxes using passed in values
        self.compute_fluxes(sfc_T,self.sfc_q_new);
        
        # Write sensible and latent heat fluxes in [W/m^2]
        Qh = c.rho_air*c.Cp_air*self.flux_wT
        Ql = c.rho_air*c.Lv*self.flux_wq
        Qg = self.ghf
        
        # Compute surface energy balance
        SEB = self.R_net - Qg - Qh - Ql
        return SEB
    
    # Compute the derivative of the surface energy budget
    def compute_dseb(self, sfc_T):
        
        # Compute derivative of SEB wrt temperature
        heat_cap = self.soil.heat_capacity(self.sfc_q_new,0)
        dSEB_dT  = 4.*self.emissivity*c.sb*(sfc_T**3)
        + c.rho_air*c.Cp_air*self.ust*self.sfc.fh(z_s,z_t,self.obl[:])
        + heat_cap*(self.soil_z[0]-self.soil_z[1])/(2*self.tstep)
        return dSEB_dT
    
    # Solve the surface energy budget
    def solve_seb(self):
        
        # Local variables
        max_iter_temp = 200
        max_iter_flux = 200
        temp_1        = self.soil_T[0] - 1
        temp_2        = self.soil_T[0] + 1
        temp_criteria = 0.001
        flux_criteria = 0.001
        
        # Compute SEB using current bracketed temperatures
        SEB_l = self.compute_seb(temp_1)
        SEB_h = self.compute_seb(temp_2)
        
        # Dynamic bracket adjustments
        out_of_bracket = (SEB_l > 0.0 and SEB_h > 0.0) or (SEB_l < 0.0 and SEB_h < 0.0)
        
        while (out_of_bracket):
            
            # Expand brackets by 1 K
            temp_1 -= 1
            temp_2 += 1
            
            # Recompute SEB at brackets
            SEB_l = self.compute_seb(temp_1)
            SEB_h = self.compute_seb(temp_2)
            
            # Check for proper brackets
            out_of_bracket = (SEB_l > 0.0 and SEB_h > 0.0) or (SEB_l < 0.0 and SEB_h < 0.0)
        
        if ((SEB_l > 0.0 and SEB_h > 0.0) or (SEB_l < 0.0 and SEB_h < 0.0)):
            throw("Please adjust brackets for Ts")
        
        # If SEB from low bracket Ts = 0, then that value of Ts is solution
        if (SEB_l == 0.0): self.sfc_T_new = temp_1
        
        # If SEB from high bracket Ts = 0, then that value of Ts is solution
        if (SEB_h == 0.0): self.sfc_T_new = temp_2
        
        # Orient the solutions such that SEB(temp_l) < 0;
        if (SEB_l < 0.0):
            temp_l = temp_1
            temp_h = temp_2
        else:
            temp_l = temp_2
            temp_h = temp_1
        
        # Prepare for convergence looping
        dTs = np.abs(temp_h-temp_l)
        
        # Convergence loop for flux
        for ff in range(0,max_iter_flux):
                    
            # Convergence loop for temperature
            for tt in range(0,max_iter_temp):
                
                # Compute SEB and dSEB_dTs
                SEB     = self.compute_seb(self.sfc_T_new);
                dSEB_dT = self.compute_dseb(self.sfc_T_new)
                
                # Update brackets
                if (SEB<0.): temp_l = self.sfc_T_new
                if (SEB>0.): temp_h = self.sfc_T_new
                
                # Bracket and bisect temperature if Newton out of range
                if ((((self.sfc_T_new-temp_h)*dSEB_dT-SEB)*((self.sfc_T_new-temp_l)*dSEB_dT-SEB)>0.0)
                    or (np.abs(2*SEB) > np.abs(dTs_old*dSEB_dT))):
                    dTs_old   = dTs
                    dTs       = 0.5*(temp_h-temp_l)
                    last_T    = self.sfc_T_new
                    self.sfc_T_new = temp_l + dTs
                    if (temp_l == self.sfc_T_new): break
                else:
                    dTs_old   = dTs
                    dTs       = SEB / dSEB_dT
                    last_T    = self.sfc_T_new
                    self.sfc_T_new = self.sfc_T_new - dTs
                    if (last_T == self.sfc_T_new): break
                
                # Check for convergence
                if (np.abs( (self.sfc_T_new-last_T)/last_T) <= temp_criteria): break
                
                # If convergence fails, recompute flux
                # computeFluxes(soil_T[0],soil_q[0]);
            
            # Save current flux for convergence criteria
            last_F = self.flux_wT
            
            # Recompute heat flux using new temperature
            self.compute_fluxes(self.sfc_T_new,self.sfc_q_new)
            
            # Check for convergence
            if (np.abs(self.flux_wT-last_F) <= flux_criteria):
                Qh = c.rho_air*c.Cp_air*self.flux_wT
                Ql = c.rho_air*c.Lv*self.flux_wq
                Qg = self.ghf
                break
            
            # If flux fails to converge, split temperature difference
            self.sfc_T_new = 0.5*(self.sfc_T_new + last_T)
    
    # Solve the surface moisture budget
    def solve_smb(self):
        
        # Local variables
        max_iter_flux = 200
        delta         = 0.5 
        flux_criteria = .001
        
        # Moisture potential at first two levels below ground
        psi0 = self.soil.water_potential(self.soil_q[0], 0)
        psi1 = self.soil.water_potential(self.soil_q[1], 1)
        
        # Compute initial soil moisture flux
        K0    = self.soil.conductivity_moisture(self.soil_q[0],0)
        K1    = self.soil.conductivity_moisture(self.soil_q[1],1)
        K_avg = 0.5*(K0+K1)
        D0    = self.soil.diffusivity_moisture(self.soil_q[0],0)
        D1    = self.soil.diffusivity_moisture(self.soil_q[1],1)
        D_avg = 0.5*(D0+D1)
        flux_sm  = c.rho_wat*K_avg*((psi0 - psi1)/(self.soil_z[0]-self.soil_z[1]) + 1)
        flux_sm = c.rho_wat*D_avg*(self.soil_q[0]-self.soil_q[1])/(self.soil_z[0]-self.soil_z[1]) + c.rho_wat*K_avg
        
        # Compute evaporation
        E = c.rho_air*self.flux_wq
        
        # Convergence loop for moisture flux
        for ff in range(0,max_iter_flux):
            
            # Save soil moisture flux for convergence test
            flux_sm_last = flux_sm
            
            # Compute new weighted soil moisture flux
            flux_sm = delta*flux_sm_last - (1.-delta)*E
            
            # Re-compute moisture potential
            psi0    = psi1 + (self.soil_z[0]-self.soil_z[1])*((flux_sm/(c.rho_wat*K_avg))-1)
            if (psi0 > self.soil.properties[0].psi_sat):
                psi0 = self.soil.properties[0].psi_sat
            
            # Update soil moisture
            sfc_q_new = self.soil.surface_water_content(psi0)
            gnd_q     = self.soil.surface_mixing_ratio(self.sfc_T_new,self.sfc_q_new,self.atm_p)
            E = c.rho_air*(gnd_q-self.atm_q)*self.ust*self.sfc.fh(z_s, z_t,self.obl[:])
            
            # Update soil moisture transfer
            K0    = self.soil.conductivity_moisture(self.sfc_q_new,0)
            K1    = self.soil.conductivity_moisture(self.sfc_q_new,1)
            K_avg = 0.5*(K0+K1)
            
            # Check for convergence
            converged = np.abs((E + flux_sm)/E) <=flux_criteria
        
            if (converged): break
    
    # Solve the diffusion equation for soil heat
    def solve_diffusion_heat(self):
        
        # Local variables
        AB  = 1.0
        AF  = 1.0-AB
        dz  = self.soil_z[0] - self.soil_z[1]
        dz2 = dz**2
        
        K     = np.zeros(self.nz)
        K_mid = np.zeros(self.nz-1)
        z_mid = np.zeros(self.nz-1)
        r     = np.zeros(self.nz-1)
        e     = np.zeros(self.nz-1)
        f     = np.zeros(self.nz-1)
        g     = np.zeros(self.nz-1)
        
        for i in range(0,self.nz-1):
            K[i]     = self.soil.diffusivity_thermal(self.soil_q[i],i)
            K[i+1]   = self.soil.diffusivity_thermal(self.soil_q[i+1],i+1)
            K_mid[i] = 0.5*(K[i]+K[i+1])
            z_mid[i] = 0.5*(self.soil_z[i]+self.soil_z[i+1])
        
        # Get the time step restriction
        dt_T = 1
        
        # Loop through diffusion by substep
        t = 0
        while (t<=self.tstep):
                
            # Set up and solve a tridiagonal matrix
            # AT(n+1) = r(n), where n denotes the time level
            # e, f, g the components of A matrix
            # T(n+1)  the soil temperature vector at t=n+1
            # r(n)    the soil temperature vector at t=n multiplied by coefficients
        
            # Matrix coefficients for first level below surface
            Cp  = self.dt_dif * dt_T * K_mid[0] / dz2
            Cm  = self.dt_dif * dt_T * K_mid[1] / dz2
            CBp = -AB * Cp
            CBm = -AB * Cm
            CB  = 1 - CBp - CBm
            CFp = AF * Cp
            CFm = AF * Cm
            CF  = 1 - CFp - CFm
        
            e[0] = 0
            f[0] = CB
            g[0] = CBm
            r[0] = CFp * self.soil_T[0] + CF * self.soil_T[1] + CFm * self.soil_T[2] - CBp * self.sfc_T_new
            
            # Matrix coefficients for the interior levels
            for i in range(1,self.nz-2):
        
                # for soil_T in this loop:
                # i   -> j+1 level
                # i+1 -> j   level
                # i+2 -> j-1 level
                Cp  = self.dt_dif * dt_T * K_mid[i] / dz2
                Cm  = self.dt_dif * dt_T * K_mid[i+1] / dz2
                CBp = -AB * Cp
                CBm = -AB * Cm
                CB  = 1 - CBp - CBm
                CFp = AF * Cp
                CFm = AF * Cm
                CF  = 1 - CFp - CFm
        
                e[i] = CBp
                f[i] = CB
                g[i] = CBm
                r[i] = CFp * self.soil_T[i] + CF * self.soil_T[i+1] + CFm * self.soil_T[i+2]
        
            # Matrix coefficients for bottom level
            j = self.nz-2
        
            Cp  = self.dt_dif * dt_T * K_mid[j] / dz2
            Cm  = self.dt_dif * dt_T * K_mid[j] / dz2
            CBp = -AB * Cp
            CBm = -AB * Cm
            CB  = 1 - CBp - CBm
            CFp = AF * Cp
            CFm = AF * Cm
            CF  = 1 - CFp - CFm
        
            e[j] = (CBp - CBm)
            f[j] = (CB + 2 * CBm)
            g[j] = 0
            r[j] = (CFp - CFm) * self.soil_T[j] + (CF + 2* CFm) * self.soil_T[j+1]
                    
            # now we can add new sfc T to column array
            self.soil_T[0] = self.sfc_T_new
        
            # Solve the tridiagonal system
            try:
                # we only need to send the layers below surface
                matrix.tridiagonal(e,f,g,r,self.soil_T[1::])
            except:
                sys.exit(0)
        
            # solve K for this step to get a new dt
            for i in range(0, self.nz-1):
                K[i]     = self.soil.diffusivity_thermal(self.soil_q[i],i)
                K[i+1]   = self.soil.diffusivity_thermal(self.soil_q[i+1],i+1)
                K_mid[i] = 0.5*(K[i]+K[i+1])
                z_mid[i] = 0.5*(self.soil_z[i]+self.soil_z[i+1])
        
            # Get the time step restriction
            Kmax = np.max(K)
            dt_T = dz2 / (2*Kmax)
            t+=dt_T
    
    # Solve the diffusion equation for soil moisture
    def solve_diffusion_mois(self):
        # Local variables
        AB  = 1.0
        AF  = 1.0-AB
        dz  = self.soil_z[0] - self.soil_z[1]
        dz2 = dz**2
        
        K_lin = np.zeros(self.nz)
        D     = np.zeros(self.nz)
        D_mid = np.zeros(self.nz-1)
        z_mid = np.zeros(self.nz-1)
        r     = np.zeros(self.nz-1)
        e     = np.zeros(self.nz-1)
        f     = np.zeros(self.nz-1)
        g     = np.zeros(self.nz-1)
        
        # Get the time step restriction
        dt_q = 1
        
        # Loop through diffusion by substep
        t = 0
        while (t<=self.tstep):
        
            # Set up and solve a tridiagonal matrix
            # AT(n+1) = r(n), where n denotes the time level
            # e, f, g the components of A matrix
            # T(n+1)  the soil temperature vector at t=n+1
            # r(n)    the soil temperature vector at t=n multiplied by coefficients
        
            # first soil level below the surface
            # common coefficients
            Cpd  = self.dt_dif * dt_q * D_mid[0] / dz2
            Cmd  = self.dt_dif * dt_q * D_mid[1] / dz2
            Cpk  = self.dt_dif * dt_q * K_lin[0] / (2*dz)
            Cmk  = self.dt_dif * dt_q * K_lin[2] / (2*dz)
        
            # coefficients for backward scheme
            CBpd = -AB * Cpd
            CBmd = -AB * Cmd
            CBpk = -AB * Cpk
            CBmk = -AB * Cmk
            CB   = (1 - CBpd - CBmd)
            CBp  = CBpd + CBpk
            CBm  = CBmd - CBmk
        
            # coefficients for forward scheme
            CFpd = AF * Cpd
            CFmd = AF * Cmd
            CFpk = AF * Cpk
            CFmk = AF * Cmk
            CF   = (1 - CFpd - CFmd)
            CFp  = CFpd + CFpk
            CFm  = CFmd - CFmk
        
            # matrix components
            e[0] = 0
            f[0] = CB
            g[0] = CBm
            r[0] = CFp * self.soil_q[0] + CF * self.soil_q[1] + CFm * self.soil_q[2] - CBp * self.sfc_q_new
        
            # interior soil levels
            for i in range(1,self.nz-2):
        
                # for soil_T in this loop:
                # i   -> j+1 level
                # i+1 -> j   level
                # i+2 -> j-1 level# 
                
                # common coefficients
                Cpd  = self.dt_dif * dt_q * D_mid[i] / dz2
                Cmd  = self.dt_dif * dt_q * D_mid[i+1] / dz2
                Cpk  = self.dt_dif * dt_q * K_lin[i] / (2*dz)
                Cmk  = self.dt_dif * dt_q * K_lin[i+2] / (2*dz)
        
                # coefficients for backward scheme
                CBpd = -AB * Cpd
                CBmd = -AB * Cmd
                CBpk = -AB * Cpk
                CBmk = -AB * Cmk
                CB   = (1 - CBpd - CBmd)
                CBp  = CBpd + CBpk
                CBm  = CBmd - CBmk
        
                # coefficients for forward scheme
                CFpd = AF * Cpd
                CFmd = AF * Cmd
                CFpk = AF * Cpk
                CFmk = AF * Cmk
                CF   = (1 - CFpd - CFmd)
                CFp  = CFpd + CFpk
                CFm  = CFmd - CFmk
        
                # matrix components
                e[i] = CBp
                f[i] = CB
                g[i] = CBm
                r[i] = CFp * self.soil_q[i] + CF * self.soil_q[i+1] + CFm * self.soil_q[i+2]
        
            # Matrix coefficients for bottom level
            j = self.nz-2
        
            # common coefficients
            Cpd  = self.dt_dif * dt_q * D_mid[j] / dz2
            Cmd  = self.dt_dif * dt_q * D_mid[j] / dz2
            Cpk  = self.dt_dif * dt_q * K_lin[j] / (2*dz)
            Cmk  = self.dt_dif * dt_q * K_lin[j] / (2*dz)
        
            # coefficients for backward scheme
            CBpd = -AB * Cpd
            CBmd = -AB * Cmd
            CBpk = -AB * Cpk
            CBmk = -AB * Cmk
            CB   = (1 - CBpd - CBmd)
            CBp  = CBpd + CBpk
            CBm  = CBmd - CBmk
        
            # coefficients for forward scheme
            CFpd = AF * Cpd
            CFmd = AF * Cmd
            CFpk = AF * Cpk
            CFmk = AF * Cmk
            CF   = (1 - CFpd - CFmd)
            CFp  = CFpd + CFpk
            CFm  = CFmd - CFmk
        
            # matrix components
            e[j] = (CBp - CBm)
            f[j] = (CB + 2 * CBm)
            g[j] = 0
            r[j] = (CFp - CFm) * self.soil_q[j] + (CF + 2 * CFm) * self.soil_q[j+1]
        
            # now we can add new sfc q to column array
            self.soil_q[0] = self.sfc_q_new
        
            # Solve the tridiagonal system
            try:
                # we only need the layers below the surface
                matrix.tridiagonal(e,f,g,r,soil_q[1::])
            except:
                sys.exit(0)
        
            # solve for D to get new dt
            for i in range(0,self.nz-1):
                D[i]     = self.soil.diffusivity_moisture(self.soil_q[i],i)
                D[i+1]   = self.soil.diffusivity_moisture(self.soil_q[i+1],i+1)
                D_mid[i] = 0.5*(D[i]+D[i+1])
                z_mid[i] = 0.5*(self.soil_z[i]+self.soil_z[i+1])
        
                # linearized K
                K_lin[i] = self.soil.conductivity_moisture(self.soil_q[i],i)/self.soil_q[i]
                if (i==self.nz-2):
                    K_lin[i+1] = soil.conductivity_moisture(self.soil_q[i+1],i+1)/self.soil_q[i+1]
        
            # get new dt
            Dmax = np.max(D)
            dt_q = dz2 / (2*Dmax)
            t+=dt_q

# main program to run the LSM
if __name__ == "__main__":
    
    # let's time this thing
    t1 = time.time()

    # a nice welcome message
    print("##############################################################")
    print("#                                                            #")
    print("#                     Welcome to UtahLSM                     #")
    print("#   A land surface model created at the University of Utah   #")
    print("#       and the NOAA National Severe Storms Laboratory       #")
    print("#                                                            #")
    print("##############################################################")
    
    # get case from user
    parser = argparse.ArgumentParser(description="Run a case with UtahLSM")
    parser.add_argument("-c", "--case", dest='case', required=True,
                        action='store', type=str, help="Case name")
    parser.add_argument("-o", "--output", dest='outfile', 
                        action='store', type=str, help="Output file name")
    args = parser.parse_args()
    case = args.case
    outf = args.outfile
    
    # create Input instance
    try:
        if os.path.exists('../cases/%s/'%case):
            namelist    = '../cases/%s/lsm_namelist.json'%case
            initfile    = '../cases/%s/lsm_init.nc'%case 
            offlinefile = '../cases/%s/lsm_offline.nc'%case  
            inputLSM    = io.Input(namelist,initfile,offlinefile)
        else:
            raise InvalidCase('Error: The folder ../cases/%s does not exist.'%case)
    except InvalidCase as e:
        print(e)
        sys.exit(1)
    
    print("[UtahLSM: Info] \t Running offline for the %s case"%case)
    
    # create Output instance
    if not outf:
        outf='lsm_%s_py.nc'%case
    outputLSM = io.Output(outf)
    
    # grid information
    nx = inputLSM.nx
    ny = inputLSM.ny
    
    # get offline input data
    ntime = inputLSM.ntime
    tstep = inputLSM.tstep
    atm_U = inputLSM.atm_U
    atm_T = inputLSM.atm_T
    atm_q = inputLSM.atm_q
    atm_p = inputLSM.atm_p
    R_net = inputLSM.r_net
    
    # local fluxes to be modified by lsm
    ustar   = 0.0
    flux_wq = 0.0
    flux_wT = 0.0
    lsm     = UtahLSM(inputLSM,outputLSM,ustar,flux_wq,flux_wT)
    
    # Loop through each time
    utc = 0
    for t in range(0,ntime):
        utc += tstep
        sys.stdout.write("\r[UtahLSM: Run] \t\t Running for time %05.2f of %05.2f"%(utc,ntime*tstep))
        sys.stdout.flush()
        
        # update user-specified fields
        lsm.update_fields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t])
        lsm.run()
        # loop through each lsm instance
        #for j in range(0,ny):
        #    for i in range(0,nx):
                # update user-specified fields
                #global_utah_lsm[k].update_fields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t])
                
                # run the model
                #global_utah_lsm[k].run()
                
                # save output
                #global_utah_lsm[k].save(outputLSM)
    
    # time info
    t2 = time.time()
    tt = t2 - t1
    print("\n[UtahLSM: Run] \t\t Done! Completed in %0.2f seconds"%tt)
    print("##############################################################")