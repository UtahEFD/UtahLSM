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
import logging
import numpy as np
import os
import sys
import time
from physics import Radiation, Soil, Surface
from util import constants as c, matrix
from util.io import Input, Output

# custom error message for user case entry
class InvalidCase(Exception):
    pass

# local logger
#logger = logging.getLogger(__name__)

# land-surface model class
class UtahLSM:
    """This is the main UtahLSM class
    
    :param inputLSM: A handle to the :class:`io.Input`
    :param outputLSM: A handle to the :class:`io.Output`
    """
    
    # model class initialization
    def __init__(self,inputLSM, outputLSM, ustar, flux_wT, flux_wq):
        """Constructor method
        """
        # set the input and output fields
        self.input  = inputLSM
        self.output = outputLSM
        
        logger.info("Reading input settings")
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
        
        # Input surface section
        self.sfc_model = self.input.sfc_model
        
        # Input soil section
        self.nz         = self.input.nsoil
        self.soil_param = self.input.soil_param
        self.soil_model = self.input.soil_model
        
        logger.info("Reading input data")
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
        self.rad_model  = self.input.rad_model
        self.albedo     = self.input.albedo
        self.emissivity = self.input.emissivity
        self.latitude   = self.input.latitude
        self.longitude  = self.input.longitude
        
        logger.info("Creating radiation model")
        if (self.rad_model):
                        
            # convert latitude and longitude into radians
            self.latitude  = self.latitude * c.pi / 180.0
            self.longitude = self.longitude * c.pi / 180.0
            
            # Create radiation model
            self.rad = Radiation.get_model(self.rad_model,self.input)
        else:
            logger.info("--- using offline data, no model")
        
        # Create soil model
        logger.info("Creating soil model")
        self.soil = Soil.get_model(self.soil_model,self.input)
        
        logger.info("Creating surface model")
        # choose surface model
        self.sfc = Surface.get_model(self.sfc_model)

        logger.info("Creating output file")    
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
    def update(self, dt, u, T, q, p, rad=0):
        self.tstep    = dt
        self.atm_U    = u
        self.atm_T    = T
        self.atm_q    = q
        self.atm_p    = p
        self.runtime += tstep
        self.utc      = np.fmod(self.runtime,86400)
        
        # Run radiation model and update time/date if needed
        if (self.rad_model):
            self.julian_day += int(self.runtime/86400);
            self.R_net  = self.rad.compute_net(self.julian_day,self.utc,self.soil_T[0])
        else:
            self.R_net = rad
        
        # Keep winds from being exactly zero
        if (self.atm_U==0): self.atm_U = 1E-4
        
        # debugging
        # if (self.runtime<1E7):
        #     print('\r')
        #     print("--------------")
        #     Logger.print_double(self.tstep, "update\t\t", "tstep")
        #     Logger.print_double(self.utc,   "update\t\t", "t utc")
        #     Logger.print_double(self.atm_U, "update\t\t", "atm_U")
        #     Logger.print_double(self.atm_T, "update\t\t", "atm_T")
        #     Logger.print_double(self.atm_q, "update\t\t", "atm_q")
        #     Logger.print_double(self.atm_p, "update\t\t", "atm_p")
        #     Logger.print_double(self.R_net, "update\t\t", "R_net")
        #     print("--------------")
        
    # Run the model
    def run(self):
                
        # Set initial new temp and moisture
        self.sfc_T_new = self.soil_T[0]
        self.sfc_q_new = self.soil_q[0]
        
        # if (self.runtime<1E7):
        #     print("--------------")
        #     Logger.print_double(self.sfc_T_new, "run\t\t\t\t\t", 'sfc_T_new')
        #     Logger.print_double(self.sfc_q_new, "run\t\t\t\t\t", 'sfc_q_new')
        #     print("--------------")
        
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
        last_L         = 1000.0
        criteria       = 0.1
        ref_T          = 300.0
        
        # Compute surface mixing ratio
        gnd_q  = self.soil.surface_mixing_ratio(sfc_T,sfc_q,self.atm_p)
        
        # Compute ground flux
        K0          = self.soil.conductivity_thermal(self.soil_q[0],0)
        K1          = self.soil.conductivity_thermal(self.soil_q[1],1)
        Kmid        = 0.5*(K0 + K1)
        self.ghf[0] = Kmid*(sfc_T - self.soil_T[1])/(self.soil_z[0]-self.soil_z[1])
        
        # Sensible flux, latent flux, ustar, and L
        for i in range(0,max_iterations):
            
            # Compute stability functions
            fm = self.sfc.fm(self.z_m, self.z_o, self.obl[0])
            fh = self.sfc.fh(self.z_s, self.z_t, self.obl[0])
            
            # Compute friction velocity
            self.ust[0] = self.atm_U*fm
            
            # Compute heat flux
            self.flux_wT[0] = (sfc_T-self.atm_T)*self.ust[0]*fh
            
            # Compute latent flux
            self.flux_wq[0] = (gnd_q-self.atm_q)*self.ust[0]*fh
                
            # Compute virtual heat flux
            flux_wTv = self.flux_wT[0] + ref_T*0.61*self.flux_wq[0]
            
            # Compute L
            last_L = self.obl[0]
            self.obl[0] = -(self.ust[0]**3)*ref_T/(c.vonk*c.grav*flux_wTv)
            
            # Bounds check on L
            if (self.z_m/self.obl[0] > 5.):  self.obl[0] = self.z_m/5.
            if (self.z_m/self.obl[0] < -5.): self.obl[0] = -self.z_m/5.
                
            # if (self.runtime==25800):
            #     print("--------------")
            #     Logger.print_double(gnd_q,           "compute_fluxes\t\t", 'gnd_q')
            #     Logger.print_double(self.ghf[0],     "compute_fluxes\t\t", 'ghf')
            #     Logger.print_double(self.ust[0],     "compute_fluxes\t\t", 'ust')
            #     Logger.print_double(self.atm_U,      "compute_fluxes\t\t", 'atm_U')
            #     Logger.print_double(fm,              "compute_fluxes\t\t", 'fm')
            #     Logger.print_double(fh,              "compute_fluxes\t\t", 'fh')
            #     Logger.print_double(self.flux_wT[0], "compute_fluxes\t\t", 'wT')
            #     Logger.print_double(self.flux_wq[0], "compute_fluxes\t\t", 'wq')
            #     Logger.print_double(flux_wTv,        "compute_fluxes\t\t", 'wTv')
            #     Logger.print_double(self.obl[0],     "compute_fluxes\t\t", 'obl')
            #     Logger.print_double(last_L,          "compute_fluxes\t\t", 'obl_old')
            #     Logger.print_double(np.abs(last_L-self.obl[0]), "compute_fluxes\t\t", 'obl_diff')
            #     Logger.print_double(criteria,     "compute_fluxes\t\t", 'criteria')
            #     print("--------------")
            
            # Check for convergence
            converged = np.abs(last_L-self.obl[0]) <= criteria
            if (converged):
                self.shf[0] = c.rho_air*c.Cp_air*self.flux_wT[0]
                self.lhf[0] = c.rho_air*c.Lv*self.flux_wq[0]
                break
        
        # Exit if L convergence fails
        # if (not converged):
        #     print("[UtahLSM: Fluxes] \t Converge failed")
        #     sys.exit()
    
    # Solve the surface energy budget using a custom implementation of Brent's Method
    def solve_seb(self):
        """
        Solves the Surface Energy Budget (SEB) to find the surface temperature.
        This function first establishes a valid temperature bracket [a, b] where the
        SEB function changes sign, then uses a robust root-finding algorithm
        (_solve_root_brent) to find the precise temperature.
        """
        # Objective function for the root finder. The root is found when SEB is zero.
        def seb_function(sfc_T):
            return self.compute_seb(sfc_T)
        
        # 1. Establish an initial temperature bracket
        temp_a = self.soil_T[0] - 1.0
        temp_b = self.soil_T[0] + 1.0
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
            temp_root, converged = self._solve_root_brent(seb_function, temp_a, temp_b)
            if not converged:
                logger.warning("SEB root-finder did not converge within the maximum iterations.")
            
            self.sfc_T_new = temp_root
            
            # Final flux calculation with the converged temperature
            self.compute_fluxes(self.sfc_T_new, self.sfc_q_new)
            logger.debug(f"SEB converged to T_sfc = {self.sfc_T_new:.3f} K")
        
        except Exception as e:
            logger.error(f"An exception occurred during SEB root finding: {e}")
            sys.exit(1)

    def _solve_root_brent(self, f, a, b, tol=1e-3, max_iter=100):
        """
        Custom implementation of Brent's method for finding the root of a function.
        It combines bisection, secant, and inverse quadratic interpolation methods.
        
        :param f: The function for which to find a root, f(x) = 0.
        :param a: The lower bound of the bracket.
        :param b: The upper bound of the bracket.
        :param tol: The desired tolerance for the root.
        :param max_iter: The maximum number of iterations to perform.
        :return: A tuple (root, converged_status).
        """
        fa = f(a)
        fb = f(b)
        
        if fa * fb >= 0:
            raise ValueError("Root not bracketed in _solve_root_brent (f(a) * f(b) >= 0).")
        
        # Ensure 'a' is the best current guess
        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa
        
        c = a  # c is the previous best approximation
        fc = fa
        mflag = True # Flag to indicate whether to use bisection
        s = 0 # Current iterate
        
        for i in range(max_iter):
            # Inverse Quadratic Interpolation
            if fa != fc and fb != fc:
                s = (a * fb * fc / ((fa - fb) * (fa - fc)) +
                     b * fa * fc / ((fb - fa) * (fb - fc)) +
                     c * fa * fb / ((fc - fa) * (fc - fb)))
            # Secant Method
            else:
                s = b - fb * (b - a) / (fb - fa)
                 
            # Check if the interpolated point is within bounds and efficient
            cond1 = (s < (3 * a + b) / 4 or s > b)
            cond2 = mflag and (abs(s - b) >= abs(b - c) / 2)
            cond3 = not mflag and (abs(s - b) >= abs(c - d) / 2)
            cond4 = mflag and (abs(b - c) < tol)
            cond5 = not mflag and (abs(c - d) < tol)
            
            if cond1 or cond2 or cond3 or cond4 or cond5:
                # If conditions are not met, fall back to bisection
                s = (a + b) / 2
                mflag = True
            else:
                mflag = False
                
            fs = f(s)
            d, c = c, b  # Update previous values
            
            # Update points for next iteration
            if fa * fs < 0:
                b = s
                fb = fs
            else:
                a = s
                fa = fs
                 
            # Ensure 'a' is the best current guess
            if abs(fa) < abs(fb):
                a, b = b, a
                fa, fb = fb, fa
                 
            # Check for convergence
            if abs(b - a) < tol:
                return b, True
                 
        return b, False
    
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
        
        # if (self.runtime<1E7):
        #     print("--------------")
        #     Logger.print_double(Qh,         "compute_seb\t\t\t", 'Qh')
        #     Logger.print_double(Ql,         "compute_seb\t\t\t", 'Ql')
        #     Logger.print_double(Qg,         "compute_seb\t\t\t", 'Qg')
        #     Logger.print_double(self.R_net, "compute_seb\t\t\t", 'Rn')
        #     Logger.print_double(SEB,        "compute_seb\t\t\t", 'SEB')
        #     print("--------------")
        
        return SEB
    
    # Compute the derivative of the surface energy budget
    def compute_dseb(self, sfc_T):
        
        # Compute derivative of SEB wrt temperature
        heat_cap = self.soil.heat_capacity(self.sfc_q_new,0)
        dSEB_dT  = 4.0*self.emissivity*c.sb*(sfc_T**3) \
        + c.rho_air*c.Cp_air*self.ust[0]*self.sfc.fh(self.z_s,self.z_t,self.obl[0]) \
        + heat_cap/(self.soil_z[0]-self.soil_z[1])
        
        # if (self.runtime<1E7):
        #     print("--------------")
        #     Logger.print_double(heat_cap, "compute_dseb\t\t", 'heat_cap')
        #     Logger.print_double(dSEB_dT,  "compute_dseb\t\t", 'dSEB_dT')
        #     print("--------------")
        
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
            K1    = self.soil.conductivity_moisture(self.soil_q[1],1)
            K_avg = 0.5*(K0+K1)
            
            # Check for convergence
            converged = np.abs((E + flux_sm)/E) <=flux_criteria
            
            # if (self.runtime<1E7): 
            #     Logger.print_double(E,       "solve_smb\t\t\t", 'E')
            #     Logger.print_double(flux_sm, "solve_smb\t\t\t", 'flux_sm')
            if (converged): 
                # if (self.runtime<1E7): 
                #     Logger.print_double(E,       'E')
                #     Logger.print_double(flux_sm, 'flux_sm')                
                break

    # Solve the diffusion equation for soil heat
    def solve_diffusion_heat(self):
        
        # if (self.runtime<1E7): 
        #     print("----BEFORET---")
        #     Logger.print_double(self.sfc_T_new,"diffusion_heat\t\t","sfc_T_new")
        #     for ii in range(self.nz):
        #         Logger.print_double(self.soil_T[ii],"diffusion_heat\t\t","soil_T (%02d)"%ii) 
        #     print("--------------")
        
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
            # we only need to send the layers below surface
            matrix.tridiagonal(e,f,g,r,self.soil_T[1::])
            
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
        # if (self.runtime<1E7): 
        #     print("----AFTERT----")
        #     for ii in range(self.nz):
        #         Logger.print_double(self.soil_T[ii],"diffusion_heat\t\t","soil_T (%02d)"%ii)
        #     print("--------------")
    
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
        
        # if (self.runtime<1E7): 
        #     print("----BEFOREQ---")
        #     Logger.print_double(self.sfc_q_new, "diffusion_mois\t\t","sfc_q_new")
        #     for ii in range(self.nz):
        #         Logger.print_double(self.soil_q[ii],"diffusion_mois\t\t","soil_q (%02d)"%ii)
        #     print("--------------")
        
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
            # we only need the layers below the surface
            matrix.tridiagonal(e,f,g,r,self.soil_q[1::])
                
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
        
        # if (self.runtime<1E7): 
        #     print("----AFTERQ----")
        #     for ii in range(self.nz):
        #         Logger.print_double(self.soil_q[ii],"diffusion_mois\t\t","soil_q (%02d)"%ii)
        #     print("--------------")

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
    
    # create Input instance
    try:
        if os.path.exists('../cases/%s/'%case):
            namelist    = '../cases/%s/lsm_namelist.json'%case
            initfile    = '../cases/%s/lsm_init.nc'%case 
            offlinefile = '../cases/%s/lsm_offline.nc'%case  
            inputLSM    = Input(namelist,initfile,offlinefile)
        else:
            raise InvalidCase('Error: The folder ../cases/%s does not exist.'%case)
    except InvalidCase as e:
        logger.error(e)
        raise SystemExit(1)

    logger.info("Running offline for the %s case"%case)
    
    # create Output instance
    if not outf:
        outf='lsm_%s_py.nc'%case
    outputLSM = Output(outf)
    
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
        logger.info(f"Running for time: {utc:8.2f} of {ntime*tstep:8.2f}")
        
        # update user-specified fields
        lsm.update(tstep,atm_U[t],atm_T[t],atm_q[t],atm_p[t],R_net[t])
        lsm.run()
        lsm.save()
    
    # time info
    t2 = time.time()
    tt = t2 - t1
    logger.info("Done! Completed in %0.4f seconds"%tt)
    logger.info("##############################################################")