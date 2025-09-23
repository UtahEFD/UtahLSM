#!/usr/bin/env python
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

from dataclasses import replace
import logging
import numpy as np

from .data_models import AtmosphericState, SurfaceState, SolverState
from .physics import Radiation, Soil, Surface
from .util import constants as c, solvers
from .util.io import Input, Output, logging_helper

# land-surface model class
class UtahLSM:
    """This is the main UtahLSM class
    
    :param input_lsm: A handle to the :class:`util.io.Input`
    :param output_lsm: A handle to the :class:`util.io.Output`
    """
    
    # model class initialization
    def __init__(self,input_lsm, output_lsm):
        """constructor method
        """
        
        # local logger
        self.logger = logging_helper.get_logger("UtahLSM")
        
        # set the input and output fields
        self.input  = input_lsm
        self.output = output_lsm
        
        # copy mutable data 
        self.soil_state  = replace(self.input.initial)
        
        # radiation model
        if (self.input.radiation.model):            
            self.rad = Radiation.get_model(self.input.radiation.model,self.input)
        else:
            self.logger.info("Using radiation forcing data")
        
        # soil model
        self.soil = Soil.get_model(self.input.soil.model,self.input)
        
        # surface model
        self.sfc = Surface.get_model(self.input.surface.model) 
        
        # initialize surface state
        self.sfc_state = SurfaceState(Ts=0,
                                      qs=0,
                                      qa=0,
                                      ust = np.zeros(1),
                                      obl = np.zeros(1),
                                      wT  = np.zeros(1),
                                      wq  = np.zeros(1),
                                      shf = np.zeros(1),
                                      lhf = np.zeros(1),
                                      ghf = np.zeros(1)
                                     )
        
        # initialize local atmospheric data
        self.atm_state: AtmosphericState = None
        
        # initialize solver state
        self.solver_state: SolverState = SolverState()
        
        # initialize local time data
        self.tstep: float = 0    # current time step

        # set reference to output dimensions
        self.output_dims = {
            't':0,
            'z':self.input.grid.nz
        }
        self.output.set_dims(self.output_dims)

        # set reference to output fields
        self.output_fields = {
            'ust'   :self.sfc_state.ust,
            'obl'   :self.sfc_state.obl,
            'shf'   :self.sfc_state.shf,
            'lhf'   :self.sfc_state.lhf,
            'ghf'   :self.sfc_state.ghf,
            'soil_z':self.input.grid.z,
            'soil_T':self.soil_state.T,
            'soil_q':self.soil_state.q,
        }
        self.output.set_fields(self.output_fields)

        # write initial data
        self.output.save(self.output_fields,0,0,initial=True)
        
    # update atmospheric quantities prior to solving
    def update(self, dt: float, runtime: float, atm_state: AtmosphericState):
        
        self.logger.info(f"[time = {runtime:7.1f}]")
        self.logger.info("Updating atmospheric state")
        
        # update model state
        self.tstep     = dt
        self.atm_state = atm_state
        
        # update surface state
        pa = self.atm_state.p
        Ts = self.soil_state.T[0]
        qs = self.soil_state.q[0]
        
        qa = self.soil.surface_mixing_ratio(Ts,qs,pa)
            
        self.sfc_state.Ts = Ts
        self.sfc_state.qs = qs
        self.sfc_state.qa = qa
        
        # run radiation model and update time/date if needed
        if (self.input.radiation.model):
            utc        = np.fmod((self.input.time.utc_start+runtime),86400)
            julian_day = self.input.time.julian_day + int(utc/86400)
            self.atm_state.R_net = self.rad.compute_net(julian_day,utc,self.atm_state,self.sfc_state)
        
        # Keep winds from being exactly zero
        if (self.atm_state.U==0): self.atm_state.U = 1E-4
        
    # Run the model
    def run(self, step_count: int, runtime: float) -> SurfaceState:
        
        self.logger.info(f"Solving soil state")
        
        # Set initial new temp and moisture
        self.sfc_state.Ts = self.soil_state.T[0]
        self.sfc_state.qs = self.soil_state.q[0]
        
        # Check if time to re-compute balances
        if ( (step_count % self.input.time.step_seb)==0 ):
            self.solve_seb()
            self.solve_smb()
        else:
            # just return new fluxes
            self.compute_fluxes(self.soil_state.T[0],self.soil_state.q[0])
        
        # check if time to compute diffusion
        if ( (step_count % self.input.time.step_dif)==0 ):
            
            # Solve heat diffusion
            self.solve_diffusion_heat()
            
            # solve moisture diffusion
            self.solve_diffusion_mois()
        
    # Save output fields
    def save(self, step_count: int, runtime: float):
        # write output
        self.logger.info(f"Saving data to file\n{'-'*19}")
        self.output.save(self.output_fields,step_count,runtime)
    
    # Compute fluxes using similarity theory
    def compute_fluxes(self, sfc_T, sfc_q):
        
        # Local constants
        VK  = c.physical.VON_KARMAN
        G   = c.physical.GRAVITY
        RHO = c.air.DENSITY
        CP  = c.thermodynamic.SPECIFIC_HEAT
        LV  = c.thermodynamic.LATENT_HEAT_VAPORIZATION
        
        # Local variables
        converged = False
        iter_max  = self.input.surface.flux_iter_max
        criteria  = self.input.surface.flux_criteria
        ref_T     = self.atm_state.T
        
        # Compute surface mixing ratio
        gnd_q  = self.soil.surface_mixing_ratio(sfc_T,sfc_q,self.atm_state.p)
        
        # Compute ground flux
        Kmid = self.solver_state.Kmid
        self.sfc_state.ghf[0] = Kmid*(sfc_T - self.soil_state.T[1])/(self.input.grid.z[0]-self.input.grid.z[1])
        
        # Sensible flux, latent flux, ustar, and L
        for i in range(0,iter_max):
            
            # Compute stability functions
            fm = self.sfc.fm(self.input.surface.z_m, self.input.surface.z_o, self.sfc_state.obl[0])
            fh = self.sfc.fh(self.input.surface.z_s, self.input.surface.z_t, self.sfc_state.obl[0])
            
            # Compute friction velocity
            self.sfc_state.ust[0] = self.atm_state.U*fm
            
            # Compute heat flux
            self.sfc_state.wT[0] = (sfc_T-self.atm_state.T)*self.sfc_state.ust[0]*fh
            
            # Compute latent flux
            self.sfc_state.wq[0] = (gnd_q-self.atm_state.q)*self.sfc_state.ust[0]*fh
            
            # Compute virtual heat flux
            flux_wTv = self.sfc_state.wT[0] + ref_T*0.61*self.sfc_state.wq[0]
            
            # Compute L
            last_L = self.sfc_state.obl[0]
            self.sfc_state.obl[0] = -(self.sfc_state.ust[0]**3)*ref_T/(VK*G*flux_wTv)
            
            # Bounds check on L
            if (self.input.surface.z_m/self.sfc_state.obl[0] > 5.): 
                self.sfc_state.obl[0] =  self.input.surface.z_m/5.
            if (self.input.surface.z_m/self.sfc_state.obl[0] < -5.): 
                self.sfc_state.obl[0] = -self.input.surface.z_m/5.
            
            # Check for convergence
            converged = np.abs(last_L-self.sfc_state.obl[0]) <= criteria
            if (converged):
                self.sfc_state.shf[0] = RHO*CP*self.sfc_state.wT[0]
                self.sfc_state.lhf[0] = RHO*LV*self.sfc_state.wq[0]
                break
    
    # Solve the surface energy budget using a custom implementation of Brent's Method
    def solve_seb(self):
        """
        Solves the Surface Energy Budget (SEB) to find the surface temperature.
        This function first establishes a valid temperature bracket [a, b] where the
        SEB function changes sign, then uses a robust root-finding algorithm
        (solvers.root_brent) to find the precise temperature.
        """
        
        # Calculate and store Kmid in the solver_state object
        # Calculate thermal conductivity for the entire soil column at once
        K_all = self.soil.conductivity_thermal(self.soil_state.q)
        
        # Store the mid-point conductivity between the top two layers
        self.solver_state.Kmid = 0.5 * (K_all[0] + K_all[1])
        
        # Objective function for the root finder. The root is found when SEB is zero.
        def seb_function(sfc_T):
            return self.compute_seb(sfc_T)
        
        # 1. Establish an initial temperature bracket
        temp_a = self.soil_state.T[0] - 1.0
        temp_b = self.soil_state.T[0] + 1.0
        seb_a  = seb_function(temp_a)
        seb_b  = seb_function(temp_b)
        
        # 2. Aggressively expand the bracket if the root is not contained within it.
        #    This loop ensures f(a) and f(b) have opposite signs.
        max_bracket_iter = 50
        iter_count = 0
        while seb_a * seb_b > 0 and iter_count < max_bracket_iter:
            if abs(seb_a) < abs(seb_b):
                temp_a -= 5.0  # Expand bracket by a larger, fixed step
                seb_a = seb_function(temp_a)
            else:
                temp_b += 5.0
                seb_b = seb_function(temp_b)
            iter_count += 1
        
        if iter_count >= max_bracket_iter:
            self.logger.error("Failed to find a valid bracket for solve_seb after %d iterations.", max_bracket_iter)
            raise SystemExit(1)
        
        # 3. Call the custom root-finder to get the surface temperature
        try:
            temp_root, converged = solvers.root_brent(seb_function, temp_a, temp_b)
            if not converged:
                self.logger.warning("SEB root-finder did not converge within the maximum iterations.")
            
            self.sfc_state.Ts = temp_root
            
            # Final flux calculation with the converged temperature
            self.compute_fluxes(self.sfc_state.Ts, self.sfc_state.qs)
            self.logger.debug(f"SEB converged to T_sfc = {self.sfc_state.Ts:.3f} K")
        
        except Exception as e:
            self.logger.error(f"An exception occurred during SEB root finding: {e}")
            raise SystemExit(1)

    # Compute the surface energy budget
    def compute_seb(self, sfc_T):

        # Compute fluxes using passed in values
        self.compute_fluxes(sfc_T,self.sfc_state.qs);
        
        # Compute surface energy balance
        SEB = self.atm_state.R_net - self.sfc_state.ghf[0] - self.sfc_state.shf[0] - self.sfc_state.lhf[0]
        
        return SEB
    
    # Solve the surface moisture budget
    def solve_smb(self):
        
        # Local constants
        RHO_W = c.water.DENSITY
        RHO_A = c.air.DENSITY
        
        # Local variables
        max_iter_flux = 200
        delta         = 0.5 
        flux_criteria = .001
        
        psi_all = self.soil.water_potential(self.soil_state.q)
        D_all   = self.soil.diffusivity_moisture(self.soil_state.q)
        K_all   = self.soil.conductivity_moisture(self.soil_state.q)
        
        # Moisture potential at the top two levels
        psi0 = psi_all[0]
        psi1 = psi_all[1]
        
        # Average diffusivity and conductivity between the top two levels
        Dmid = 0.5 * (D_all[0] + D_all[1])
        Kmid = 0.5 * (K_all[0] + K_all[1])
        
        #D0    = self.soil.diffusivity_moisture(self.soil_state.q[0],0)
        #D1    = self.soil.diffusivity_moisture(self.soil_state.q[1],1) 
        #D_avg = 0.5*(D0+D1)
        
        flux_sm  = RHO_W*Kmid*((psi0 - psi1)/(self.input.grid.z[0]-self.input.grid.z[1]) + 1.0)
        #flux_sm  = c.rho_wat*D_avg*(self.soil_state.q[0]-self.soil_state.q[1])/(self.input.grid.z[0]-self.input.grid.z[1]) + c.rho_wat*Kmid
        
        # Compute evaporation
        E = RHO_A*self.sfc_state.wq[0]
        
        # Convergence loop for moisture flux
        for ff in range(0,max_iter_flux):
            
            # Save soil moisture flux for convergence test
            flux_sm_last = flux_sm
            
            # Compute new weighted soil moisture flux
            flux_sm = delta*flux_sm_last - (1.0-delta)*E
            
            # Re-compute moisture potential
            
            psi0 = psi1 + (self.input.grid.z[0]-self.input.grid.z[1])*((flux_sm/(RHO_W*Kmid))-1.0)
            
            psi_sat = self.soil.properties.psi_sat[0]
            if (psi0 > psi_sat):
                psi0 = psi_sat
            
            # Update soil moisture
            self.sfc_state.qs = self.soil.surface_water_content(psi0)
            
            gnd_q = self.soil.surface_mixing_ratio(self.sfc_state.Ts,self.sfc_state.qs,self.atm_state.p)
            E     = RHO_A*(gnd_q-self.atm_state.q)*self.sfc_state.ust[0]*self.sfc.fh(self.input.surface.z_s,self.input.surface.z_t,self.sfc_state.obl[0])
            
            # Update soil moisture transfer
            K0    = self.soil.conductivity_moisture_scalar(self.sfc_state.qs,0)
            K1    = K_all[1]
            Kmid = 0.5*(K0+K1)
            
            # Check for convergence
            converged = np.abs((E + flux_sm)/E) <=flux_criteria
            
            if (converged):
                self.solver_state.Kmid = Kmid
                break

    # Solve the diffusion equation for soil heat
    def solve_diffusion_heat(self):
        
        # Local variables
        AB       = self.input.numerics.diffusion_back_weight
        AF       = 1.0-AB
        nz       = self.input.grid.nz
        dz       = self.input.grid.z[0] - self.input.grid.z[1]
        dz2      = dz**2
        step_dif = self.input.time.step_dif
        
        r, e, f, g = [np.zeros(self.input.grid.nz - 1) for _ in range(4)]
        
        K_all = self.soil.diffusivity_thermal(self.soil_state.q)
        
        K_mid = K_mid = 0.5 * (K_all[:-1] + K_all[1:]) 
        
        # Get the time step restriction
        dt_T = self.tstep
        
        # loop through diffusion by sub-step
        t = 0
        while (t<=self.tstep):
            
            # Get the time step restriction
            #K_max = np.max(K_all)
            dt_T  = self.tstep#dz2 / (2.0 * K_max)
            
            # Ensure the last sub-step lands exactly on the main time step
            # if t + dt_T > self.tstep:
            #     dt_T = self.tstep - t
            
            # Set up and solve a tridiagonal matrix
            # AT(n+1) = r(n), where n denotes the time level
            # e, f, g the components of A matrix
            # T(n+1)  the soil temperature vector at t=n+1
            # r(n)    the soil temperature vector at t=n multiplied by coefficients
        
            # Matrix coefficients for first level below surface
            Cp  = float(step_dif) * dt_T * K_mid[0] / dz2
            Cm  = float(step_dif) * dt_T * K_mid[1] / dz2
            CBp = -AB * Cp
            CBm = -AB * Cm
            CB  = 1.0 - CBp - CBm
            CFp = AF * Cp
            CFm = AF * Cm
            CF  = 1.0 - CFp - CFm
        
            e[0] = 0
            f[0] = CB
            g[0] = CBm
            r[0] = CFp * self.soil_state.T[0] + CF * self.soil_state.T[1] + CFm * self.soil_state.T[2] - CBp * self.sfc_state.Ts
            
            # Matrix coefficients for the interior levels
            for i in range(1,self.input.grid.nz-2):
        
                # for soil_T in this loop:
                # i   -> j+1 level
                # i+1 -> j   level
                # i+2 -> j-1 level
                Cp  = float(step_dif) * dt_T * K_mid[i] / dz2
                Cm  = float(step_dif) * dt_T * K_mid[i+1] / dz2
                CBp = -AB * Cp
                CBm = -AB * Cm
                CB  = 1.0 - CBp - CBm
                CFp = AF * Cp
                CFm = AF * Cm
                CF  = 1.0 - CFp - CFm
        
                e[i] = CBp
                f[i] = CB
                g[i] = CBm
                r[i] = CFp * self.soil_state.T[i] + CF * self.soil_state.T[i+1] + CFm * self.soil_state.T[i+2]
        
            # Matrix coefficients for bottom level
            j = self.input.grid.nz-2
        
            Cp  = float(step_dif) * dt_T * K_mid[j] / dz2
            Cm  = float(step_dif) * dt_T * K_mid[j] / dz2
            CBp = -AB * Cp
            CBm = -AB * Cm
            CB  = 1.0 - CBp - CBm
            CFp = AF * Cp
            CFm = AF * Cm
            CF  = 1.0 - CFp - CFm
        
            e[j] = (CBp - CBm)
            f[j] = (CB + 2.0 * CBm)
            g[j] = 0
            r[j] = (CFp - CFm) * self.soil_state.T[j] + (CF + 2.0* CFm) * self.soil_state.T[j+1]
            
            # now we can add new sfc T to column array
            self.soil_state.T[0] = self.sfc_state.Ts
        
            # Solve the tridiagonal system
            # we only need to send the layers below surface
            self.soil_state.T[1::] = solvers.tridiagonal(e,f,g,r)
            
            # update time
            t+=dt_T
    
    # Solve the diffusion equation for soil moisture
    def solve_diffusion_mois(self):
        
        # Local variables
        AB  = self.input.numerics.diffusion_back_weight
        AF  = 1.0-AB
        dz  = self.input.grid.z[0] - self.input.grid.z[1]
        dz2 = dz**2
        
        r, e, f, g = [np.zeros(self.input.grid.nz - 1) for _ in range(4)]
        
        # loop through diffusion by sub-step
        t = 0
        while (t<=self.tstep):
            
            D_all = self.soil.diffusivity_moisture(self.soil_state.q)
            K_all = self.soil.conductivity_moisture(self.soil_state.q)
            D_mid = 0.5 * (D_all[:-1] + D_all[1:])
            K_lin = K_all/self.soil_state.q
            
            #D_max = np.max(D_all)
            dt_q = self.tstep#dz2 / (2.0 * D_max)
            
            # if t + dt_q > self.tstep:
            #     dt_q = self.tstep - t
            
            # Set up and solve a tridiagonal matrix
            # AT(n+1) = r(n), where n denotes the time level
            # e, f, g the components of A matrix
            # T(n+1)  the soil temperature vector at t=n+1
            # r(n)    the soil temperature vector at t=n multiplied by coefficients
            
            # first soil level below the surface
            # common coefficients
            Cpd  = float(self.input.time.step_dif) * dt_q * D_mid[0] / dz2
            Cmd  = float(self.input.time.step_dif) * dt_q * D_mid[1] / dz2
            Cpk  = float(self.input.time.step_dif) * dt_q * K_lin[0] / (2*dz)
            Cmk  = float(self.input.time.step_dif) * dt_q * K_lin[2] / (2*dz)
            
            # coefficients for backward scheme
            CBpd = -AB * Cpd
            CBmd = -AB * Cmd
            CBpk = -AB * Cpk
            CBmk = -AB * Cmk
            CB   = (1.0 - CBpd - CBmd)
            CBp  = CBpd + CBpk
            CBm  = CBmd - CBmk
            
            # coefficients for forward scheme
            CFpd = AF * Cpd
            CFmd = AF * Cmd
            CFpk = AF * Cpk
            CFmk = AF * Cmk
            CF   = (1.0 - CFpd - CFmd)
            CFp  = CFpd + CFpk
            CFm  = CFmd - CFmk
            
            # matrix components
            e[0] = 0
            f[0] = CB
            g[0] = CBm
            r[0] = CFp * self.soil_state.q[0] + CF * self.soil_state.q[1] + CFm * self.soil_state.q[2] - CBp * self.sfc_state.qs
            
            # interior soil levels
            for i in range(1,self.input.grid.nz-2):
                # for soil_T in this loop:
                # i   -> j+1 level
                # i+1 -> j   level
                # i+2 -> j-1 level
                
                # common coefficients
                Cpd  = float(self.input.time.step_dif) * dt_q * D_mid[i] / dz2
                Cmd  = float(self.input.time.step_dif) * dt_q * D_mid[i+1] / dz2
                Cpk  = float(self.input.time.step_dif) * dt_q * K_lin[i] / (2*dz)
                Cmk  = float(self.input.time.step_dif) * dt_q * K_lin[i+2] / (2*dz)
                
                # coefficients for backward scheme
                CBpd = -AB * Cpd
                CBmd = -AB * Cmd
                CBpk = -AB * Cpk
                CBmk = -AB * Cmk
                CB   = (1.0 - CBpd - CBmd)
                CBp  = CBpd + CBpk
                CBm  = CBmd - CBmk
                
                # coefficients for forward scheme
                CFpd = AF * Cpd
                CFmd = AF * Cmd
                CFpk = AF * Cpk
                CFmk = AF * Cmk
                CF   = (1.0 - CFpd - CFmd)
                CFp  = CFpd + CFpk
                CFm  = CFmd - CFmk
                
                # matrix components
                e[i] = CBp
                f[i] = CB
                g[i] = CBm
                r[i] = CFp * self.soil_state.q[i] + CF * self.soil_state.q[i+1] + CFm * self.soil_state.q[i+2]
            
            # Matrix coefficients for bottom level
            j = self.input.grid.nz-2
            
            # common coefficients
            Cpd  = float(self.input.time.step_dif) * dt_q * D_mid[j] / dz2
            Cmd  = float(self.input.time.step_dif) * dt_q * D_mid[j] / dz2
            Cpk  = float(self.input.time.step_dif) * dt_q * K_lin[j] / (2*dz)
            Cmk  = float(self.input.time.step_dif) * dt_q * K_lin[j] / (2*dz)
            
            # coefficients for backward scheme
            CBpd = -AB * Cpd
            CBmd = -AB * Cmd
            CBpk = -AB * Cpk
            CBmk = -AB * Cmk
            CB   = (1.0 - CBpd - CBmd)
            CBp  = CBpd + CBpk
            CBm  = CBmd - CBmk
            
            # coefficients for forward scheme
            CFpd = AF * Cpd
            CFmd = AF * Cmd
            CFpk = AF * Cpk
            CFmk = AF * Cmk
            CF   = (1.0 - CFpd - CFmd)
            CFp  = CFpd + CFpk
            CFm  = CFmd - CFmk
            
            # matrix components
            e[j] = (CBp - CBm)
            f[j] = (CB + 2.0 * CBm)
            g[j] = 0
            r[j] = (CFp - CFm) * self.soil_state.q[j] + (CF + 2.0 * CFm) * self.soil_state.q[j+1]
            
            # now we can add new sfc q to column array
            self.soil_state.q[0] = self.sfc_state.qs
            
            # solve the tridiagonal system
            # we only need the layers below the surface
            self.soil_state.q[1::] = solvers.tridiagonal(e,f,g,r)
            
            # update time
            t+=dt_q
