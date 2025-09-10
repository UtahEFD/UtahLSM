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

import argparse
from dataclasses import replace
import logging
import numpy as np
import os
import sys
import time
from data_models import AtmosphericData, Fluxes
from physics import Radiation, Soil, Surface
from util import constants as c, solvers
from util.io import Input, Output

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
        # set the input and output fields
        self.input  = input_lsm
        self.output = output_lsm
        
        # copy mutable data
        logger.info("Copying input data")    
        self.soil_state  = replace(self.input.initial)
        
        # initialize new surface values for first run
        self.sfc_T_new = self.soil_state.T.item(0)
        self.sfc_q_new = self.soil_state.q.item(0)
        
        # radiation model
        logger.info("Creating radiation model")
        if (self.input.radiation.model):            
            self.rad = Radiation.get_model(self.rad_model,self.input)
        else:
            logger.info("--- using offline data, no model")
        
        # soil model
        logger.info("Creating soil model")
        self.soil = Soil.get_model(self.input.soil.model,self.input)
        
        # surface model
        logger.info("Creating surface model")
        self.sfc = Surface.get_model(self.input.surface.model)
        
        # data for output file
        logger.info("Creating output file")    
        
        # initialize fluxes
        self.fluxes = Fluxes(ust = np.zeros(1),
                             obl = np.zeros(1),
                             wT  = np.zeros(1),
                             wq  = np.zeros(1),
                             shf = np.zeros(1),
                             lhf = np.zeros(1),
                             ghf = np.zeros(1)
                            )
        
        # initialize local atmospheric data
        self.atm_state: AtmosphericData = None
        
        # initialize local time data
        self.tstep = 0    # current time step

        # set reference to output dimensions
        self.output_dims = {
            't':0,
            'z':self.input.grid.nz
        }
        self.output.set_dims(self.output_dims)

        # set reference to output fields
        self.output_fields = {
            'ust':self.fluxes.ust,
            'obl':self.fluxes.obl,
            'shf':self.fluxes.shf,
            'lhf':self.fluxes.lhf,
            'ghf':self.fluxes.ghf,
            'soil_z':self.input.grid.z,
            'soil_T':self.soil_state.T,
            'soil_q':self.soil_state.q,
        }
        self.output.set_fields(self.output_fields)

        # write initial data
        self.output.save(self.output_fields,0,0,initial=True)
        
    # update atmospheric quantities prior to solving
    def update(self, dt: float, runtime: float, atm_state: AtmosphericData):
        
        # update model state
        self.tstep     = dt
        self.atm_state = atm_state
        
        # run radiation model and update time/date if needed
        if (self.input.radiation.model):
            utc        = np.fmod((self.input.time.utc_start+runtime),86400)
            julian_day = self.input.time.julian_day + int(utc_total/86400)
            self.R_net = self.rad.compute_net(julian_day,utc,self.soil_state.T[0])
        else:
            self.R_net = self.atm_state.R_net
        
        # Keep winds from being exactly zero
        if (self.atm_state.U==0): self.atm_state.U = 1E-4
        
    # Run the model
    def run(self, step_count: int) -> Fluxes:
                
        # Set initial new temp and moisture
        self.sfc_T_new = self.soil_state.T[0]
        self.sfc_q_new = self.soil_state.q[0]
        
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
        self.output.save(self.output_fields,step_count,runtime)
    
    # Compute fluxes using similarity theory
    def compute_fluxes(self, sfc_T, sfc_q):
        
        # Local variables
        max_iterations = 200
        converged      = False
        last_L         = 1000.0
        criteria       = 0.1
        ref_T          = 300.0
        
        # Compute surface mixing ratio
        gnd_q  = self.soil.surface_mixing_ratio(sfc_T,sfc_q,self.atm_state.p)
        
        # Compute ground flux
        K0          = self.soil.conductivity_thermal(self.soil_state.q[0],0)
        K1          = self.soil.conductivity_thermal(self.soil_state.q[1],1)
        Kmid        = 0.5*(K0 + K1)
        self.fluxes.ghf[0] = Kmid*(sfc_T - self.soil_state.T[1])/(self.input.grid.z[0]-self.input.grid.z[1])
        
        # Sensible flux, latent flux, ustar, and L
        for i in range(0,max_iterations):
            
            # Compute stability functions
            fm = self.sfc.fm(self.input.surface.z_m, self.input.surface.z_o, self.fluxes.obl[0])
            fh = self.sfc.fh(self.input.surface.z_s, self.input.surface.z_t, self.fluxes.obl[0])
            
            # Compute friction velocity
            self.fluxes.ust[0] = self.atm_state.U*fm
            
            # Compute heat flux
            self.fluxes.wT[0] = (sfc_T-self.atm_state.T)*self.fluxes.ust[0]*fh
            
            # Compute latent flux
            self.fluxes.wq[0] = (gnd_q-self.atm_state.q)*self.fluxes.ust[0]*fh
                
            # Compute virtual heat flux
            flux_wTv = self.fluxes.wT[0] + ref_T*0.61*self.fluxes.wq[0]
            
            # Compute L
            last_L = self.fluxes.obl[0]
            self.fluxes.obl[0] = -(self.fluxes.ust[0]**3)*ref_T/(c.vonk*c.grav*flux_wTv)
            
            # Bounds check on L
            if (self.input.surface.z_m/self.fluxes.obl[0] > 5.):  
                self.fluxes.obl[0] =  self.input.surface.z_m/5.
            if (self.input.surface.z_m/self.fluxes.obl[0] < -5.): 
                self.fluxes.obl[0] = -self.input.surface.z_m/5.
            
            # Check for convergence
            converged = np.abs(last_L-self.fluxes.obl[0]) <= criteria
            if (converged):
                self.fluxes.shf[0] = c.rho_air*c.Cp_air*self.fluxes.wT[0]
                self.fluxes.lhf[0] = c.rho_air*c.Lv*self.fluxes.wq[0]
                break
    
    # Solve the surface energy budget using a custom implementation of Brent's Method
    def solve_seb(self):
        """
        Solves the Surface Energy Budget (SEB) to find the surface temperature.
        This function first establishes a valid temperature bracket [a, b] where the
        SEB function changes sign, then uses a robust root-finding algorithm
        (solvers.root_brent) to find the precise temperature.
        """
        # Objective function for the root finder. The root is found when SEB is zero.
        def seb_function(sfc_T):
            return self.compute_seb(sfc_T)
        
        # 1. Establish an initial temperature bracket
        temp_a = self.soil_state.T[0] - 1.0
        temp_b = self.soil_state.T[0] + 1.0
        seb_a = seb_function(temp_a)
        seb_b = seb_function(temp_b)
        
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
            logger.error("Failed to find a valid bracket for solve_seb after %d iterations.", max_bracket_iter)
            sys.exit(1)
        
        # 3. Call the custom root-finder to get the surface temperature
        try:
            temp_root, converged = solvers.root_brent(seb_function, temp_a, temp_b)
            if not converged:
                logger.warning("SEB root-finder did not converge within the maximum iterations.")
            
            self.sfc_T_new = temp_root
            
            # Final flux calculation with the converged temperature
            self.compute_fluxes(self.sfc_T_new, self.sfc_q_new)
            logger.debug(f"SEB converged to T_sfc = {self.sfc_T_new:.3f} K")
        
        except Exception as e:
            logger.error(f"An exception occurred during SEB root finding: {e}")
            sys.exit(1)

    # Compute the surface energy budget
    def compute_seb(self, sfc_T):

        # Compute fluxes using passed in values
        self.compute_fluxes(sfc_T,self.sfc_q_new);
        
        # Compute surface energy balance
        SEB = self.atm_state.R_net - self.fluxes.ghf[0] - self.fluxes.shf[0] - self.fluxes.lhf[0]
        
        return SEB
    
    # Solve the surface moisture budget
    def solve_smb(self):
        
        # Local variables
        max_iter_flux = 200
        delta         = 0.5 
        flux_criteria = .001
        
        # Moisture potential at first two levels below ground
        psi0 = self.soil.water_potential(self.soil_state.q[0], 0)
        psi1 = self.soil.water_potential(self.soil_state.q[1], 1)
        
        # Compute initial soil moisture flux
        K0    = self.soil.conductivity_moisture(self.soil_state.q[0],0)
        K1    = self.soil.conductivity_moisture(self.soil_state.q[1],1)
        K_avg = 0.5*(K0+K1)
        
        D0    = self.soil.diffusivity_moisture(self.soil_state.q[0],0)
        D1    = self.soil.diffusivity_moisture(self.soil_state.q[1],1) 
        D_avg = 0.5*(D0+D1)
        
        flux_sm  = c.rho_wat*K_avg*((psi0 - psi1)/(self.input.grid.z[0]-self.input.grid.z[1]) + 1.0)
        #flux_sm  = c.rho_wat*D_avg*(self.soil_state.q[0]-self.soil_state.q[1])/(self.input.grid.z[0]-self.input.grid.z[1]) + c.rho_wat*K_avg
        
        # Compute evaporation
        E = c.rho_air*self.fluxes.wq[0]
        
        # Convergence loop for moisture flux
        for ff in range(0,max_iter_flux):
            
            # Save soil moisture flux for convergence test
            flux_sm_last = flux_sm
            
            # Compute new weighted soil moisture flux
            flux_sm = delta*flux_sm_last - (1.0-delta)*E
            
            # Re-compute moisture potential
            
            psi0 = psi1 + (self.input.grid.z[0]-self.input.grid.z[1])*((flux_sm/(c.rho_wat*K_avg))-1.0)
            
            if (psi0 > self.soil.properties[0].psi_sat):
                psi0 = self.soil.properties[0].psi_sat
            
            # Update soil moisture
            self.sfc_q_new = self.soil.surface_water_content(psi0)
            
            gnd_q = self.soil.surface_mixing_ratio(self.sfc_T_new,self.sfc_q_new,self.atm_state.p)
            E     = c.rho_air*(gnd_q-self.atm_state.q)*self.fluxes.ust[0]*self.sfc.fh(self.input.surface.z_s,self.input.surface.z_t,self.fluxes.obl[0])
            
            # Update soil moisture transfer
            K0    = self.soil.conductivity_moisture(self.sfc_q_new,0)
            K1    = self.soil.conductivity_moisture(self.soil_state.q[1],1)
            K_avg = 0.5*(K0+K1)
            
            # Check for convergence
            converged = np.abs((E + flux_sm)/E) <=flux_criteria
            
            if (converged): 
                break

    # Solve the diffusion equation for soil heat
    def solve_diffusion_heat(self):
        
        # Local variables
        AB  = 1.0
        AF  = 1.0-AB
        dz  = self.input.grid.z[0] - self.input.grid.z[1]
        dz2 = dz**2
        
        K     = np.zeros(self.input.grid.nz)
        K_mid = np.zeros(self.input.grid.nz-1)
        z_mid = np.zeros(self.input.grid.nz-1)
        r     = np.zeros(self.input.grid.nz-1)
        e     = np.zeros(self.input.grid.nz-1)
        f     = np.zeros(self.input.grid.nz-1)
        g     = np.zeros(self.input.grid.nz-1)
        
        for i in range(0,self.input.grid.nz-1):
            K[i]     = self.soil.diffusivity_thermal(self.soil_state.q[i],i)
            K[i+1]   = self.soil.diffusivity_thermal(self.soil_state.q[i+1],i+1)
            K_mid[i] = 0.5*(K[i]+K[i+1])
            z_mid[i] = 0.5*(self.input.grid.z[i]+self.input.grid.z[i+1])
        
        # Get the time step restriction
        dt_T = 1.0
        
        # loop through diffusion by sub-step
        t = 0
        while (t<=self.tstep):

            # Set up and solve a tridiagonal matrix
            # AT(n+1) = r(n), where n denotes the time level
            # e, f, g the components of A matrix
            # T(n+1)  the soil temperature vector at t=n+1
            # r(n)    the soil temperature vector at t=n multiplied by coefficients
        
            # Matrix coefficients for first level below surface
            Cp  = float(self.input.time.step_dif) * dt_T * K_mid[0] / dz2
            Cm  = float(self.input.time.step_dif) * dt_T * K_mid[1] / dz2
            CBp = -AB * Cp
            CBm = -AB * Cm
            CB  = 1.0 - CBp - CBm
            CFp = AF * Cp
            CFm = AF * Cm
            CF  = 1.0 - CFp - CFm
        
            e[0] = 0
            f[0] = CB
            g[0] = CBm
            r[0] = CFp * self.soil_state.T[0] + CF * self.soil_state.T[1] + CFm * self.soil_state.T[2] - CBp * self.sfc_T_new
            
            # Matrix coefficients for the interior levels
            for i in range(1,self.input.grid.nz-2):
        
                # for soil_T in this loop:
                # i   -> j+1 level
                # i+1 -> j   level
                # i+2 -> j-1 level
                Cp  = float(self.input.time.step_dif) * dt_T * K_mid[i] / dz2
                Cm  = float(self.input.time.step_dif) * dt_T * K_mid[i+1] / dz2
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
        
            Cp  = float(self.input.time.step_dif) * dt_T * K_mid[j] / dz2
            Cm  = float(self.input.time.step_dif) * dt_T * K_mid[j] / dz2
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
            self.soil_state.T[0] = self.sfc_T_new
        
            # Solve the tridiagonal system
            # we only need to send the layers below surface
            self.soil_state.T[1::] = solvers.tridiagonal(e,f,g,r)
            
            # update conductivities for sub-step
            for i in range(0, self.input.grid.nz-1):
                K[i]     = self.soil.diffusivity_thermal(self.soil_state.q[i],i)
                K[i+1]   = self.soil.diffusivity_thermal(self.soil_state.q[i+1],i+1)
                K_mid[i] = 0.5*(K[i]+K[i+1])
                z_mid[i] = 0.5*(self.input.grid.z[i]+self.input.grid.z[i+1])
            
            # adjust time step if not at final time
            if (t!=self.tstep):
                
                # compute new diffusion time step
                Kmax = np.max(K)
                dt_T = dz2 / (2.0*Kmax)
                
                # check if we need to relax dt to meet end time exactly
                if (t+dt_T>self.tstep):
                    dt_T = self.tstep - t
            
            # update time
            t+=dt_T
    
    # Solve the diffusion equation for soil moisture
    def solve_diffusion_mois(self):
        
        # Local variables
        AB  = 1.0
        AF  = 1.0-AB
        dz  = self.input.grid.z[0] - self.input.grid.z[1]
        dz2 = dz**2
        
        K_lin = np.zeros(self.input.grid.nz)
        D     = np.zeros(self.input.grid.nz)
        D_mid = np.zeros(self.input.grid.nz-1)
        z_mid = np.zeros(self.input.grid.nz-1)
        r     = np.zeros(self.input.grid.nz-1)
        e     = np.zeros(self.input.grid.nz-1)
        f     = np.zeros(self.input.grid.nz-1)
        g     = np.zeros(self.input.grid.nz-1)
        
        # Get the time step restriction
        dt_q = 1.0
        
        # loop through diffusion by sub-step
        t = 0
        while (t<=self.tstep):
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
            r[0] = CFp * self.soil_state.q[0] + CF * self.soil_state.q[1] + CFm * self.soil_state.q[2] - CBp * self.sfc_q_new
            
            # interior soil levels
            for i in range(1,self.input.grid.nz-2):
                # for soil_T in this loop:
                # i   -> j+1 level
                # i+1 -> j   level
                # i+2 -> j-1 level# 
                
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
            self.soil_state.q[0] = self.sfc_q_new
            
            # solve the tridiagonal system
            # we only need the layers below the surface
            self.soil_state.q[1::] = solvers.tridiagonal(e,f,g,r)
                
            # update diffusivities and conductivities for sub-step
            for i in range(0,self.input.grid.nz-1):
                D[i]     = self.soil.diffusivity_moisture(self.soil_state.q[i],i)
                D[i+1]   = self.soil.diffusivity_moisture(self.soil_state.q[i+1],i+1)
                D_mid[i] = 0.5*(D[i]+D[i+1])
                z_mid[i] = 0.5*(self.input.grid.z[i]+self.input.grid.z[i+1])
                
                # linearized K
                K_lin[i] = self.soil.conductivity_moisture(self.soil_state.q[i],i)/self.soil_state.q[i]
                if (i==self.input.grid.nz-2):
                    K_lin[i+1] = self.soil.conductivity_moisture(self.soil_state.q[i+1],i+1)/self.soil_state.q[i+1]
            
            # adjust time step if not at final time
            if (t!=self.tstep):
                
                # compute new diffusion time step
                Dmax = np.max(D)
                dt_q = dz2 / (2.0*Dmax)
                
                # check if we need to relax dt to meet end time exactly
                if (t+dt_q>self.tstep):
                    dt_q = self.tstep - t
            
            # update time
            t+=dt_q

# main program to run the LSM
if __name__ == "__main__":
    
    # let's time this thing
    t1 = time.time()
    
    # configure logging
    log_format = '{asctime} [{levelname:^8s}] {name:^20s} {message}'
    logging.basicConfig(level=logging.DEBUG,
                        format=log_format,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        style='{',
                        filename='utahlsm.log',
                        filemode='w')    
    
    # create a console handler for printing to the screen
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(log_format, "%Y-%m-%d %H:%M:%S",style='{',))
    logging.getLogger('').addHandler(console_handler)
    
    # local logger
    logger = logging.getLogger("UtahLSM")
    
    # a nice welcome message
    logger.info("##############################################################")
    logger.info("#                                                            #")
    logger.info("#                     Welcome to UtahLSM                     #")
    logger.info("#   A land surface model created at the University of Utah   #")
    logger.info("#       and the NOAA National Severe Storms Laboratory       #")
    logger.info("#                                                            #")
    logger.info("##############################################################")
    
    # get case from user
    parser = argparse.ArgumentParser(description="Run a case with UtahLSM")
    parser.add_argument("-c", "--case", dest='case', required=True,
                        action='store', type=str, help="Case name")
    parser.add_argument("-o", "--output", dest='outfile', 
                        action='store', type=str, help="Output file name")
    args = parser.parse_args()
    case = args.case
    outf = args.outfile
    
    # file paths to input data
    namelist    = '../cases/%s/lsm_namelist.json'%case
    initfile    = '../cases/%s/lsm_init.nc'%case 
    offlinefile = '../cases/%s/lsm_offline.nc'%case
    
    # create Input instance
    try:
        input_lsm = Input(namelist,initfile,offlinefile)
    except:
        raise SystemExit(1)

    logger.info("Running offline for the %s case"%case)
    
    # create Output instance
    if not outf:
        outf='lsm_%s_py.nc'%case
    output_lsm = Output(outf)
    
    # grid information (not used yet)
    nx = input_lsm.grid.nx
    ny = input_lsm.grid.ny
    
    # --- Main Time Loop ---
    logger.info("Starting simulation time loop...")
    
    # local fluxes to be modified by lsm
    lsm = UtahLSM(input_lsm,output_lsm)
    
    # Loop through each time
    runtime = 0
    tstep = input_lsm.forcing.tstep
    for step_count, atm_state in enumerate(input_lsm.forcing.atmos):
        runtime += tstep
        logger.info(f"Running for time: {runtime:8.2f} of {input_lsm.forcing.ntime*tstep:8.2f}")
        
        # update user-specified fields
        lsm.update(tstep, runtime, atm_state)
        lsm.run(runtime)
        lsm.save(step_count,runtime)
    
    # time info
    t2 = time.time()
    tt = t2 - t1
    logger.info("Done! Completed in %0.4f seconds"%tt)
    logger.info("##############################################################")