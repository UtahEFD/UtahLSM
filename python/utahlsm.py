#!/usr/bin/env python
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
        # set the input and output fields
        self.input  = inputLSM
        self.output = outputLSM
        
        print("[UtahLSM: Setup] \t Reading input settings")
        # Input time section
        self.dt_seb     = self.input.step_seb
        self.dt_dif     = self.input.step_dif
        self.utc        = self.input.utc_start
        self.julian_day = self.input.julian_day
        
        # Input grid section
        self.nx = self.input.nx
        self.ny = self.input.ny
        
        # Input length scale section
        self.z_o = self.input.z_o
        self.z_t = self.input.z_t
        self.z_m = self.input.z_m
        self.z_s = self.input.z_s
        
        # Input soil section
        self.nz         = self.input.nsoil
        self.soil_param = self.input.param
        self.soil_model = self.input.model
        
        print("[UtahLSM: Setup] \t Reading input data")
        self.soil_z    = self.input.soil_z
        self.soil_T    = self.input.soil_T
        self.soil_q    = self.input.soil_q
        self.soil_type = self.input.soil_type
        
        # Initialize new surface values for first run
        self.sfc_T_new = self.soil_T.item(0)
        self.sfc_q_new = self.soil_q.item(0)
        
        # Initialize history arrays for first run
        self.soil_T_last = self.soil_T
        self.soil_q_last = self.soil_q
        
        # Modify soil levels to be negative away from surface
        self.soil_z = -1*self.soil_z
        
        # Input radiation section
        self.comp_rad   = self.input.comp_rad
        self.albedo     = self.input.albedo
        self.emissivity = self.input.emissivity
        self.latitude   = self.input.latitude
        self.longitude  = self.input.longitude
        if (self.comp_rad==1):
            print("[UtahLSM: Setup] \t Creating radiation model")
                        
            # convert latitude and longitude into radians
            self.latitude  = self.latitude * constants.pi / 180.0
            self.longitude = self.longitude * constants.pi / 180.0
            
            # Create radiation model
            self.rad = Radiation.get_model(1,self.input)
        else:
            print("[UtahLSM: Radiation] \t --- using offline data, no model")
        
        # Create soil model
        print("[UtahLSM: Setup] \t Creating soil model")
        self.soil = Soil.get_model(self.soil_model,self.input)
        
        print("[UtahLSM: Setup] \t Creating surface model")
        # choose surface model
        self.sfc = Surface.get_model(1)

        print("[UtahLSM: Setup] \t Creating output file")    
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

        # set reference to output dimensions
        self.output_dims = {
            't':0,
            'z':self.nz
        }
        self.output.set_dims(self.output_dims)

        # set reference to output fields
        self.output_fields = {
            'ust':self.ust,
            'obl':self.obl,
            'shf':self.shf,
            'lhf':self.lhf,
            'ghf':self.ghf,
            'soil_z':self.soil_z,
            'soil_T':self.soil_T,
            'soil_q':self.soil_q,
        }
        self.output.set_fields(self.output_fields)

        # write initial data
        self.output.save(self.output_fields,0,0,initial=True)
        
    # Update atmospheric quantities prior to solving
    def update_fields(self, dt, u, T, q, p, rad=0):
        self.tstep    = dt
        self.atm_U    = u
        self.atm_T    = T
        self.atm_q    = q
        self.atm_p    = p
        self.runtime += tstep
        self.utc      = np.fmod(self.runtime,86400)
        
        # Run radiation model and update time/date if needed
        if (self.comp_rad==1):
            julian_day += int(runtime/86400);
            self.R_net  = self.rad.computeNet(julian_day,self.utc,self.soil_T[0])
        else:
            self.R_net = rad
        
        # Keep winds from being exactly zero
        if (self.atm_U==0): self.atm_U = 0.1
        
        # debugging
        if (self.runtime==29400.0):
            print('\r')
            print("--------------")
            io.Logger.print_hex(self.tstep, "update_fields\t\t", "tstep")
            io.Logger.print_hex(self.utc,   "update_fields\t\t", "t utc")
            io.Logger.print_hex(self.atm_U, "update_fields\t\t", "atm_U")
            io.Logger.print_hex(self.atm_T, "update_fields\t\t", "atm_T")
            io.Logger.print_hex(self.atm_q, "update_fields\t\t", "atm_q")
            io.Logger.print_hex(self.atm_p, "update_fields\t\t", "atm_p")
            io.Logger.print_hex(self.R_net, "update_fields\t\t", "R_net")
            print("--------------")
        
    # Run the model
    def run(self):
                
        # Set initial new temp and moisture
        self.sfc_T_new = self.soil_T[0]
        self.sfc_q_new = self.soil_q[0]
        
        if (self.runtime==29400.0):
            print("--------------")
            io.Logger.print_hex(self.sfc_T_new, "run\t\t\t\t\t", 'sfc_T_new')
            io.Logger.print_hex(self.sfc_q_new, "run\t\t\t\t\t", 'sfc_q_new')
            print("--------------")
        
        # Check if time to re-compute balances
        if ( (self.step_count % self.dt_seb)==0 ):
            self.solve_seb()
            self.solve_smb()
        else:
            # just return new fluxes
            self.compute_fluxes(self.soil_T[0],self.soil_q[0])
        
        # Save current temperature and moisture
        self.soil_T_last = self.soil_T.copy()
        self.soil_q_last = self.soil_q.copy()
        
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
    def save(self):
        # write output
        self.output.save(self.output_fields,self.step_count,self.runtime,initial=False)
        # 
        # # close output file       
        #self.output.close()
    
    # Compute fluxes using similarity theory
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
            
            # Compute ground flux
            K0          = self.soil.conductivity_thermal(self.soil_q[0],0)
            K1          = self.soil.conductivity_thermal(self.soil_q[1],1)
            Kmid        = 0.5*(K0 + K1)
            self.ghf[0] = Kmid*(sfc_T - self.soil_T[1])/(self.soil_z[0]-self.soil_z[1])
            
            # Compute friction velocity
            self.ust[0] = self.atm_U*self.sfc.fm(self.z_m, self.z_o, self.obl[0])
            
            # Compute heat flux
            self.flux_wT[0] = (sfc_T-self.atm_T)*self.ust[0]*self.sfc.fh(self.z_s, self.z_t, self.obl[0])
            
            # Compute latent flux
            if ( (self.first) and (i == 0)):
                self.flux_wq[0] = (self.R_net - self.ghf[0] - self.flux_wT[0]*c.rho_air*c.Cp_air)/(c.rho_air*c.Lv)
                gnd_q = self.atm_q + self.flux_wq.item(0) / (self.ust[0]*self.sfc.fh(self.z_s,self.z_t,self.obl[0]))
                self.sfc_q_new = self.soil.surface_water_content_estimate(self.soil_T[0],gnd_q, self.atm_p)
                self.soil_q[0] = self.sfc_q_new
            else:
                self.flux_wq[0] = (gnd_q-self.atm_q)*self.ust[0]*self.sfc.fh(self.z_s,self.z_t,self.obl[0])
                
            # Compute virtual heat flux
            flux_wTv = self.flux_wT[0] + ref_T*0.61*self.flux_wq[0]
                
            # Compute L
            last_L = self.obl[0]
            self.obl[0] = -(self.ust[0]**3)*ref_T/(c.vonk*c.grav*flux_wTv)
            
            if (self.runtime==29400.0):
                print("--------------")
                io.Logger.print_hex(gnd_q,           "compute_fluxes\t\t", 'gnd_q')
                io.Logger.print_hex(self.ghf[0],     "compute_fluxes\t\t", 'ghf')
                io.Logger.print_hex(self.ust[0],     "compute_fluxes\t\t", 'ust')
                io.Logger.print_hex(self.flux_wT[0], "compute_fluxes\t\t", 'wT')
                io.Logger.print_hex(self.flux_wq[0], "compute_fluxes\t\t", 'wq')
                io.Logger.print_hex(flux_wTv,        "compute_fluxes\t\t", 'wTv')
                io.Logger.print_hex(self.obl[0],     "compute_fluxes\t\t", 'obl')
                print("--------------")
            
            # Bounds check on L
            if (self.z_m/self.obl[0] > 5.):  self.obl[0] = self.z_m/5.
            if (self.z_m/self.obl[0] < -5.): self.obl[0] = self.z_m/5.
            
            # Check for convergence
            converged = np.abs(last_L-self.obl[0]) <= criteria
            if (converged):
                self.shf[0] = c.rho_air*c.Cp_air*self.flux_wT[0]
                self.lhf[0] = c.rho_air*c.Lv*self.flux_wq[0]
                break

        # Exit if L convergence fails
        if (not converged):
            print("[UtahLSM: Fluxes] \t Converge failed")
            sys.exit()
    
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
        dTs     = np.abs(temp_h-temp_l)
        dTs_old = dTs
        
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
                    or (np.abs(2.0*SEB) > np.abs(dTs_old*dSEB_dT))):
                    dTs_old        = dTs
                    dTs            = 0.5*(temp_h-temp_l)
                    last_T         = self.sfc_T_new
                    self.sfc_T_new = temp_l + dTs
                    if (temp_l == self.sfc_T_new): break
                else:
                    dTs_old        = dTs
                    dTs            = SEB / dSEB_dT
                    last_T         = self.sfc_T_new
                    self.sfc_T_new = self.sfc_T_new - dTs
                    if (last_T == self.sfc_T_new): break
                
                # Check for convergence
                if (np.abs( (self.sfc_T_new-last_T)/last_T) <= temp_criteria): break
                
                # If convergence fails, recompute flux
                # computeFluxes(soil_T[0],soil_q[0]);
            
            # Save current flux for convergence criteria
            last_F = self.flux_wT[0]
            
            # Recompute heat flux using new temperature
            self.compute_fluxes(self.sfc_T_new,self.sfc_q_new)
            
            # Check for convergence
            if (np.abs(self.flux_wT-last_F) <= flux_criteria):
                Qh = c.rho_air*c.Cp_air*self.flux_wT[0]
                Ql = c.rho_air*c.Lv*self.flux_wq[0]
                Qg = self.ghf[0]
                if (self.runtime==29400.0):
                    print("--------------")
                    io.Logger.print_hex(Qh, "compute_seb\t\t\t",'Qh')
                    io.Logger.print_hex(Ql, "compute_seb\t\t\t",'Ql')
                    io.Logger.print_hex(Qg, "compute_seb\t\t\t",'Qg')
                    print("--------------")
                break
            
            # If flux fails to converge, split temperature difference
            self.sfc_T_new = 0.5*(self.sfc_T_new + last_T)
    
    # Compute the surface energy budget
    def compute_seb(self, sfc_T):

        # Compute fluxes using passed in values
        self.compute_fluxes(sfc_T,self.sfc_q_new);
        
        # Write sensible and latent heat fluxes in [W/m^2]
        Qh = c.rho_air*c.Cp_air*self.flux_wT[0]
        Ql = c.rho_air*c.Lv*self.flux_wq[0]
        Qg = self.ghf[0]
        
        # Compute surface energy balance
        SEB = self.R_net - Qg - Qh - Ql
        
        if (self.runtime==29400.0):
            print("--------------")
            io.Logger.print_hex(Qh,         "compute_seb\t\t\t", 'Qh')
            io.Logger.print_hex(Ql,         "compute_seb\t\t\t", 'Ql')
            io.Logger.print_hex(Qg,         "compute_seb\t\t\t", 'Qg')
            io.Logger.print_hex(self.R_net, "compute_seb\t\t\t", 'Rn')
            io.Logger.print_hex(SEB,        "compute_seb\t\t\t", 'SEB')
            print("--------------")
        
        return SEB
    
    # Compute the derivative of the surface energy budget
    def compute_dseb(self, sfc_T):
        
        # Compute derivative of SEB wrt temperature
        heat_cap = self.soil.heat_capacity(self.sfc_q_new,0)
        dSEB_dT  = 4.0*self.emissivity*c.sb*(sfc_T**3) \
        + c.rho_air*c.Cp_air*self.ust[0]*self.sfc.fh(self.z_s,self.z_t,self.obl[0]) \
        + heat_cap*(self.soil_z[0]-self.soil_z[1])/(2.0*self.tstep)
        
        if (self.runtime==29400.0):
            print("--------------")
            io.Logger.print_hex(heat_cap, "compute_dseb\t\t", 'heat_cap')
            io.Logger.print_hex(dSEB_dT,  "compute_dseb\t\t", 'dSEB_dT')
            print("--------------")
        
        return dSEB_dT
    
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
        
        flux_sm  = c.rho_wat*K_avg*((psi0 - psi1)/(self.soil_z[0]-self.soil_z[1]) + 1.0)
        #flux_sm  = c.rho_wat*D_avg*(self.soil_q[0]-self.soil_q[1])/(self.soil_z[0]-self.soil_z[1]) + c.rho_wat*K_avg
        
        # Compute evaporation
        E = c.rho_air*self.flux_wq[0]
        
        # Convergence loop for moisture flux
        for ff in range(0,max_iter_flux):
            
            # Save soil moisture flux for convergence test
            flux_sm_last = flux_sm
            
            # Compute new weighted soil moisture flux
            flux_sm = delta*flux_sm_last - (1.0-delta)*E
            
            # Re-compute moisture potential
            psi0 = psi1 + (self.soil_z[0]-self.soil_z[1])*((flux_sm/(c.rho_wat*K_avg))-1.0)
            if (psi0 > self.soil.properties[0].psi_sat):
                psi0 = self.soil.properties[0].psi_sat
            
            # Update soil moisture
            self.sfc_q_new = self.soil.surface_water_content(psi0)
            
            gnd_q = self.soil.surface_mixing_ratio(self.sfc_T_new,self.sfc_q_new,self.atm_p)
            E     = c.rho_air*(gnd_q-self.atm_q)*self.ust[0]*self.sfc.fh(self.z_s,self.z_t,self.obl[0])
            
            # Update soil moisture transfer
            K0    = self.soil.conductivity_moisture(self.sfc_q_new,0)
            K1    = self.soil.conductivity_moisture(self.sfc_q_new,1)
            K_avg = 0.5*(K0+K1)
            
            # Check for convergence
            converged = np.abs((E + flux_sm)/E) <=flux_criteria
            if (self.runtime==29400.0): 
                io.Logger.print_hex(E,       "solve_smb\t\t\t", 'E')
                io.Logger.print_hex(flux_sm, "solve_smb\t\t\t", 'flux_sm')
            if (converged): 
                # if (self.runtime==29400.0): 
                #     io.Logger.print_hex(E,       'E')
                #     io.Logger.print_hex(flux_sm, 'flux_sm')                
                break

    # Solve the diffusion equation for soil heat
    def solve_diffusion_heat(self):
        
        if (self.runtime==29400.0): 
            print("----BEFORET---")
            io.Logger.print_hex(self.sfc_T_new,"diffusion_heat\t\t","sfc_T_new")
            for ii in range(self.nz):
                io.Logger.print_hex(self.soil_T[ii],"diffusion_heat\t\t","soil_T (%02d)"%ii) 
            print("--------------")
        
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
            Cp  = float(self.dt_dif) * dt_T * K_mid[0] / dz2
            Cm  = float(self.dt_dif) * dt_T * K_mid[1] / dz2
            CBp = -AB * Cp
            CBm = -AB * Cm
            CB  = 1.0 - CBp - CBm
            CFp = AF * Cp
            CFm = AF * Cm
            CF  = 1.0 - CFp - CFm
        
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
                Cp  = float(self.dt_dif) * dt_T * K_mid[i] / dz2
                Cm  = float(self.dt_dif) * dt_T * K_mid[i+1] / dz2
                CBp = -AB * Cp
                CBm = -AB * Cm
                CB  = 1.0 - CBp - CBm
                CFp = AF * Cp
                CFm = AF * Cm
                CF  = 1.0 - CFp - CFm
        
                e[i] = CBp
                f[i] = CB
                g[i] = CBm
                r[i] = CFp * self.soil_T[i] + CF * self.soil_T[i+1] + CFm * self.soil_T[i+2]
        
            # Matrix coefficients for bottom level
            j = self.nz-2
        
            Cp  = float(self.dt_dif) * dt_T * K_mid[j] / dz2
            Cm  = float(self.dt_dif) * dt_T * K_mid[j] / dz2
            CBp = -AB * Cp
            CBm = -AB * Cm
            CB  = 1.0 - CBp - CBm
            CFp = AF * Cp
            CFm = AF * Cm
            CF  = 1.0 - CFp - CFm
        
            e[j] = (CBp - CBm)
            f[j] = (CB + 2.0 * CBm)
            g[j] = 0
            r[j] = (CFp - CFm) * self.soil_T[j] + (CF + 2.0* CFm) * self.soil_T[j+1]
                    
            # now we can add new sfc T to column array
            self.soil_T[0] = self.sfc_T_new
        
            # Solve the tridiagonal system
            try:
                # we only need to send the layers below surface
                matrix.tridiagonal(e,f,g,r,self.soil_T[1::])
            except:
                sys.exit(0)
            
            # update conductivities for sub-step
            for i in range(0, self.nz-1):
                K[i]     = self.soil.diffusivity_thermal(self.soil_q[i],i)
                K[i+1]   = self.soil.diffusivity_thermal(self.soil_q[i+1],i+1)
                K_mid[i] = 0.5*(K[i]+K[i+1])
                z_mid[i] = 0.5*(self.soil_z[i]+self.soil_z[i+1])
            
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
        if (self.runtime==29400.0): 
            print("----AFTERT----")
            for ii in range(self.nz):
                io.Logger.print_hex(self.soil_T[ii],"diffusion_heat\t\t","soil_T (%02d)"%ii)
            print("--------------")
    
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
        dt_q = 1.0
        
        if (self.runtime==29400.0): 
            print("----BEFOREQ---")
            io.Logger.print_hex(self.sfc_q_new, "diffusion_mois\t\t","sfc_q_new")
            for ii in range(self.nz):
                io.Logger.print_hex(self.soil_q[ii],"diffusion_mois\t\t","soil_q (%02d)"%ii)
            print("--------------")
        
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
            Cpd  = float(self.dt_dif) * dt_q * D_mid[0] / dz2
            Cmd  = float(self.dt_dif) * dt_q * D_mid[1] / dz2
            Cpk  = float(self.dt_dif) * dt_q * K_lin[0] / (2*dz)
            Cmk  = float(self.dt_dif) * dt_q * K_lin[2] / (2*dz)
            
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
            r[0] = CFp * self.soil_q[0] + CF * self.soil_q[1] + CFm * self.soil_q[2] - CBp * self.sfc_q_new
            
            # interior soil levels
            for i in range(1,self.nz-2):
                # for soil_T in this loop:
                # i   -> j+1 level
                # i+1 -> j   level
                # i+2 -> j-1 level# 
                
                # common coefficients
                Cpd  = float(self.dt_dif) * dt_q * D_mid[i] / dz2
                Cmd  = float(self.dt_dif) * dt_q * D_mid[i+1] / dz2
                Cpk  = float(self.dt_dif) * dt_q * K_lin[i] / (2*dz)
                Cmk  = float(self.dt_dif) * dt_q * K_lin[i+2] / (2*dz)
                
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
                r[i] = CFp * self.soil_q[i] + CF * self.soil_q[i+1] + CFm * self.soil_q[i+2]
            
            # Matrix coefficients for bottom level
            j = self.nz-2
            
            # common coefficients
            Cpd  = float(self.dt_dif) * dt_q * D_mid[j] / dz2
            Cmd  = float(self.dt_dif) * dt_q * D_mid[j] / dz2
            Cpk  = float(self.dt_dif) * dt_q * K_lin[j] / (2*dz)
            Cmk  = float(self.dt_dif) * dt_q * K_lin[j] / (2*dz)
            
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
            r[j] = (CFp - CFm) * self.soil_q[j] + (CF + 2.0 * CFm) * self.soil_q[j+1]
            
            # now we can add new sfc q to column array
            self.soil_q[0] = self.sfc_q_new
            
            # solve the tridiagonal system
            try:
                # we only need the layers below the surface
                matrix.tridiagonal(e,f,g,r,self.soil_q[1::])
            except:
                sys.exit(0)
                
            # update diffusivities and conductivities for sub-step
            for i in range(0,self.nz-1):
                D[i]     = self.soil.diffusivity_moisture(self.soil_q[i],i)
                D[i+1]   = self.soil.diffusivity_moisture(self.soil_q[i+1],i+1)
                D_mid[i] = 0.5*(D[i]+D[i+1])
                z_mid[i] = 0.5*(self.soil_z[i]+self.soil_z[i+1])
                
                # linearized K
                K_lin[i] = self.soil.conductivity_moisture(self.soil_q[i],i)/self.soil_q[i]
                if (i==self.nz-2):
                    K_lin[i+1] = self.soil.conductivity_moisture(self.soil_q[i+1],i+1)/self.soil_q[i+1]
            
            
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
        
        if (self.runtime==29400.0): 
            print("----AFTERQ----")
            for ii in range(self.nz):
                io.Logger.print_hex(self.soil_q[ii],"diffusion_mois\t\t","soil_q (%02d)"%ii)
            print("--------------")

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
        sys.stdout.write("\r[UtahLSM: Run] \t\t Running for time: %05.2f of %05.2f"%(utc,ntime*tstep))
        sys.stdout.flush()
        
        # update user-specified fields
        lsm.update_fields(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t])
        lsm.run()
        lsm.save()
    
    # time info
    t2 = time.time()
    tt = t2 - t1
    print("\n[UtahLSM: Run] \t\t Done! Completed in %0.4f seconds"%tt)
    print("##############################################################")