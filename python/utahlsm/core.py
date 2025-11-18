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

"""The core UtahLSM model orchestrator.

This module contains the main `UtahLSM` class, which serves as the central
controller for the land-surface model. It initializes the various physics
components (radiation, soil, surface), manages the model's state, and steps
through the simulation in time.
"""

import logging
from dataclasses import replace
from typing import Callable, Optional

import numpy as np

from .data_models import (
    AtmosphericState, SoilState, SolverState, SurfaceState
)
from .exceptions import NamelistError, UtahLSMError
from .physics import Radiation, Soil, Surface
from .physics.radiation.factory import get_radiation_model
from .physics.soil.factory import get_soil_model
from .physics.surface.factory import get_surface_model
from .util import constants as c, solvers
from .util.io import Input, Output, logging_helper


class UtahLSM:
    """The main Utah Land-Surface Model class.

    This class orchestrates the entire simulation. It holds the model's
    state, calls the physics modules in the correct sequence, and manages
    the time-stepping and I/O operations.

    Attributes:
        input: An instance of the Input class containing all configuration.
        output: An instance of the Output class for writing simulation results.
        soil_state: The current state of the soil column (temperature/moisture).
        sfc_state: The current state of the surface diagnostics (fluxes, etc.).
        atm_state: The current near-surface atmospheric forcing conditions.
        solver_state: Intermediate variables for the numerical solvers.
        rad: The selected radiation model instance.
        soil: The selected soil model instance.
        sfc: The selected surface layer model instance.
        tstep: The current model time step [s].
        logger: A logger instance for this class.
        output_dims: A dictionary of dimensions for the output file.
        output_fields: A dictionary of fields to be written to the output file.
    """
    def __init__(self, input_lsm: Input, output_lsm: Output) -> None:
        """Initializes the UtahLSM model.

        Args:
            input_lsm: An `Input` object containing the model configuration.
            output_lsm: An `Output` object for handling data output.
        """
        self.logger: logging.Logger = logging_helper.get_logger('UtahLSM')
        self.input: Input = input_lsm
        self.output: Output = output_lsm
        self.tstep: float = 0.0

        self._setup_states()
        self._setup_physics()
        self._setup_output()

    #--- Public Methods ---

    def update(self, dt: float, runtime: float,
               atm_state: AtmosphericState) -> None:
        """Updates the model with new atmospheric forcing data for the
        current step.

        Args:
            dt: The time step duration [s].
            runtime: The total elapsed simulation time [s].
            atm_state: The atmospheric state for the current time step.
        """
        self.logger.info('[time = %7.1f]', runtime)
        self.logger.info('Updating atmospheric state')

        self.tstep = dt
        self.atm_state = atm_state  # pylint: disable=attribute-defined-outside-init

        # Update surface state from the top soil layer
        sfc_T = self.soil_state.temperature[0]
        sfc_q = self.soil_state.moisture[0]
        sfc_r = self.soil.surface_mixing_ratio(
            sfc_T, sfc_q, self.atm_state.pressure)
        self.sfc_state.temperature = sfc_T
        self.sfc_state.moisture = sfc_q
        self.sfc_state.specific_humidity = sfc_r

        # Run radiation model if configured
        if self.input.radiation.model:
            utc = np.fmod((self.input.time.utc_start+runtime), 86400)
            # Wrap julian day to stay in valid range [1, 365/366]
            # accounting for leap years
            days_per_year = (366 if self._is_leap_year(
                self.input.time.utc_year) else 365)
            julian_day = ((self.input.time.julian_day + int(utc/86400) - 1)
                          % days_per_year) + 1
            self.atm_state.radiation_net = self.rad.compute_net(
                julian_day, utc, self.atm_state, self.sfc_state)

    def run(self) -> None:
        """Runs the core model physics for a single time step.

        This includes solving the surface energy and moisture budgets and
        updating the soil profiles via diffusion solvers.
        """
        self.logger.info('Solving soil state')

        # Set initial guesses for new surface temp and moisture
        self.sfc_state.temperature = self.soil_state.temperature[0]
        self.sfc_state.moisture = self.soil_state.moisture[0]

        # Solve surface energy and moisture budgets
        self._solve_seb()
        self._solve_smb()

        # Solve diffusion equations for heat and moisture
        self._solve_diffusion_heat()
        self._solve_diffusion_mois()

    def save(self, step_count: int, runtime: float) -> None:
        """Saves the model's current state to the output file.

        Args:
            step_count: The current time step number.
            runtime: The total elapsed simulation time [s].
        """
        self.logger.info('Saving data to file\n-------------------')
        self.output.save(self.output_fields,step_count,runtime)

    # --- Internal Methods ---

    def _setup_states(self) -> None:
        """Initializes all state containers for the model."""
        self.logger.info('Setting up initial states')
        self.soil_state: SoilState = replace(self.input.initial)
        self.sfc_state: SurfaceState = SurfaceState()
        self.atm_state: AtmosphericState = AtmosphericState()
        self.solver_state: SolverState = SolverState()

    def _setup_physics(self) -> None:
        """Initializes the physics modules based on user configuration."""
        self.logger.info('Initializing physics modules')
        try:
            if self.input.radiation.model:
                self.rad: Radiation = get_radiation_model(
                    self.input.radiation.model,
                    self.input.radiation.latitude,
                    self.input.radiation.longitude,
                    self.input.surface.albedo,
                    self.input.surface.emissivity
                )
            else:
                self.rad: Radiation = None  # type: ignore
                self.logger.info('Using radiation forcing data')
            self.soil: Soil = get_soil_model(
                self.input.soil.model,
                self.input.soil_properties,
                self.input.soil_type_names,
                self.input.soil_properties_name
            )
            self.sfc: Surface = get_surface_model(self.input.surface.model)
        except NamelistError as e:
            self.logger.error('Failed to initialize physics modules: %s.', e)
            raise

    def _setup_output(self) -> None:
        """Sets up the output file dimensions and fields."""
        self.output_dims: dict = {
            't': 0,
            'z': self.input.grid.nz
        }
        self.output.set_dims(self.output_dims)

        self.output_fields: dict = {
            'ust': self.sfc_state.turbulence.friction_velocity,
            'obl': self.sfc_state.turbulence.obukhov_length,
            'shf': self.sfc_state.fluxes.sensible_heat,
            'lhf': self.sfc_state.fluxes.latent_heat,
            'ghf': self.sfc_state.fluxes.ground_heat,
            'soil_z': self.input.grid.z,
            'soil_T': self.soil_state.temperature,
            'soil_q': self.soil_state.moisture,
        }
        self.output.set_fields(self.output_fields)
        self.output.save(self.output_fields, 0, 0, initial=True)

    @staticmethod
    def _is_leap_year(year: int) -> bool:
        """Determines if a year is a leap year.

        A year is a leap year if:
        - It is divisible by 4 AND not divisible by 100, OR
        - It is divisible by 400

        Args:
            year: The year to check.

        Returns:
            True if the year is a leap year, False otherwise.
        """
        return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

    def _solve_seb(self) -> None:
        """Solves the Surface Energy Budget (SEB) to find surface temperature.

        This function first establishes a valid temperature bracket [a, b]
        where the SEB function changes sign, then uses a robust root-finding
        algorithm (solvers.root_brent) to find the precise temperature.
        """
        # Calculate thermal conductivity for the entire soil column
        K_all = self.soil.conductivity_thermal(self.soil_state.moisture)
        self.solver_state.K_mid = 0.5 * (K_all[0] + K_all[1])

        # Establish an initial temperature bracket
        temp_a = self.soil_state.temperature[0] - 1.0
        temp_b = self.soil_state.temperature[0] + 1.0
        seb_a  = self._compute_seb(temp_a)
        seb_b  = self._compute_seb(temp_b)

        # Expand the bracket if the root is not contained within it
        ITER_MAX = self.input.numerics.iterations.seb_bracket
        iter_count = 0
        while seb_a * seb_b > 0 and iter_count < ITER_MAX:
            if abs(seb_a) < abs(seb_b):
                temp_a -= 5.0
                seb_a = self._compute_seb(temp_a)
            else:
                temp_b += 5.0
                seb_b = self._compute_seb(temp_b)
            iter_count += 1

        if iter_count >= ITER_MAX:
            self.logger.error(
                'Failed to find a valid bracket for _solve_seb.')
            raise UtahLSMError(
                f'Failed to find a valid bracket for surface energy '
                f'balance after {ITER_MAX} iterations.'
            )

        # Find the root (surface temperature)
        try:
            ITER_MAX = self.input.numerics.iterations.seb_root
            TOLERANCE = self.input.numerics.tolerances.seb_root
            temp, converged = solvers.root_brent(
                self._compute_seb, temp_a, temp_b, ITER_MAX, TOLERANCE)
            if not converged:
                self.logger.warning('SEB root-finder did not converge.')
            self.sfc_state.temperature = temp
            self._compute_fluxes(self.sfc_state.temperature,
                                 self.sfc_state.moisture)
            self.logger.debug(
                'SEB converged to T_sfc = %.3f K',
                self.sfc_state.temperature)
        except Exception as e:
            self.logger.error('Error during SEB root finding: %s', e)
            raise

    # Compute the surface energy budget
    def _compute_seb(self, sfc_T: float) -> float:
        """Computes the surface energy budget residual for a given surface
            temperature.

        Args:
            sfc_T: The surface temperature [K] to test.

        Returns:
            The residual of the surface energy budget [W/m^2].
        """
        self._compute_fluxes(sfc_T, self.sfc_state.moisture)
        SEB = (self.atm_state.radiation_net
               - self.sfc_state.fluxes.ground_heat[0]
               - self.sfc_state.fluxes.sensible_heat[0]
               - self.sfc_state.fluxes.latent_heat[0])

        return SEB

    def _solve_smb(self) -> None:
        """Solves the Surface Moisture Budget (SMB).

        This function finds the soil moisture flux and surface evaporation.
        It then iteratively blends the two in time until convergence.
        """
        # Local constants and variables
        delta = 0.5

        RHO_W = c.water.DENSITY
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        TOL = self.input.numerics.tolerances.smb_flux
        ITER_MAX = self.input.numerics.iterations.smb_flux

        z_s = self.input.surface.z_s
        z_t = self.input.surface.z_t
        dz = self.input.grid.z[0] - self.input.grid.z[1]

        atm_T = self.atm_state.temperature
        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        sfc_T = self.sfc_state.temperature
        L = self.sfc_state.turbulence.obukhov_length[0]
        ust = self.sfc_state.turbulence.friction_velocity[0]
        rho_a = atm_p / (RD * atm_T)
        fh = self.sfc.fh(z_s, z_t, L)

        psi_sat = self.soil.properties.psi_sat[0]
        psi_all = self.soil.water_potential(self.soil_state.moisture)
        K_hydraulic = self.soil.conductivity_moisture(self.soil_state.moisture)
        psi0 = psi_all[0]
        psi1 = psi_all[1]
        K0 = K_hydraulic[0]
        K1 = K_hydraulic[1]
        K_mid = 0.5 * (K0 + K1)

        # Soil moisture flux and evaporation
        flux_sm  = RHO_W*K_mid*((psi0 - psi1)/dz + 1.0)
        E = rho_a*self.sfc_state.fluxes.kinematic_moisture[0]

        # Iteratively solve for moisture flux
        for _ in range(0,ITER_MAX):

            # New blended soil moisture flux
            flux_sm_last = flux_sm
            flux_sm = delta*flux_sm_last - (1.0-delta)*E

            # New evaporation
            psi0 = psi1 + dz*((flux_sm/(RHO_W*K_mid))-1.0)
            psi0 = min(psi0, psi_sat)
            self.sfc_state.moisture = self.soil.surface_water_content(
                psi0)
            gnd_q = self.soil.surface_mixing_ratio(
                sfc_T, self.sfc_state.moisture, atm_p)
            E = rho_a*(gnd_q-atm_q)*ust*fh

            K0 = self.soil.conductivity_moisture(
                self.sfc_state.moisture, level=0)
            K_mid = 0.5*(K0+K1)  # K_mid is hydraulic conductivity at midpoint

            if abs((E + flux_sm) / E) <= TOL:
                break

    def _compute_fluxes(self, sfc_T: float, sfc_q: float) -> None:
        """Computes surface fluxes using Monin-Obukhov Similarity Theory.

        This is an iterative process to find the friction velocity (ustar)
        and Obukhov length (L) that are consistent with the calculated
        sensible and latent heat fluxes.

        Args:
            sfc_T: Surface temperature [K].
            sfc_q: Surface moisture [m^3/m^3].

        Raises:
            Exception: If the flux computation encounters an unrecoverable error
                in the Monin-Obukhov iteration process.
        """
        # Local constants and variables
        converged = False

        VK = c.physical.VON_KARMAN
        G = c.physical.GRAVITY
        CP = c.thermodynamic.SPECIFIC_HEAT
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        LV = c.thermodynamic.LATENT_HEAT_VAPORIZATION
        TOL = self.input.numerics.tolerances.sfc_flux
        ITER_MAX = self.input.numerics.iterations.sfc_flux

        atm_T = self.atm_state.temperature
        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        atm_ws = self.atm_state.wind_speed
        ref_T = atm_T
        rho = atm_p / (RD * atm_T)

        z_m = self.input.surface.z_m
        z_o = self.input.surface.z_o
        z_s = self.input.surface.z_s
        z_t = self.input.surface.z_t
        dz = self.input.grid.z[0] - self.input.grid.z[1]

        # Compute surface-air specific humidity
        gnd_q  = self.soil.surface_mixing_ratio(sfc_T, sfc_q, atm_p)

        # Compute ground heat flux
        K_mid = self.solver_state.K_mid
        self.sfc_state.fluxes.ground_heat[0] = (
            K_mid * (sfc_T - self.soil_state.temperature[1]) / dz)

        # Iteratively solve for fluxes and stability
        L = self.sfc_state.turbulence.obukhov_length[0]
        for _ in range(0,ITER_MAX):

            # Stability functions
            fm = self.sfc.fm(z_m, z_o, L)
            fh = self.sfc.fh(z_s, z_t, L)

            # Friction velocity
            self.sfc_state.turbulence.friction_velocity[0] = atm_ws*fm
            ustar = self.sfc_state.turbulence.friction_velocity[0]

            # Kinematic fluxes
            flux_wT = (sfc_T-atm_T)*ustar*fh
            flux_wq = (gnd_q-atm_q)*ustar*fh
            self.sfc_state.fluxes.kinematic_heat[0] = flux_wT
            self.sfc_state.fluxes.kinematic_moisture[0] = flux_wq
            flux_wTv = flux_wT + ref_T*0.61*flux_wq

            # Obukhov length
            last_L = L
            if flux_wTv != 0:
                L = -(ustar**3) * ref_T / (VK * G * flux_wTv)
            else:
                L = 1e6 # Large positive for neutral

            # Bound L to prevent extreme instability/stability
            if z_m/L > 5.0:
                L = z_m/5.0
            elif z_m/L < -5.0:
                L = -z_m/5.0

            # Check for convergence
            if abs(last_L - L) <= TOL:
                self.sfc_state.turbulence.obukhov_length[0] = L
                self.sfc_state.fluxes.sensible_heat[0] = (
                    rho * CP * flux_wT)
                self.sfc_state.fluxes.latent_heat[0] = (
                    rho * LV * flux_wq)
                converged = True
                break

        # Set final Obukhov length if loop completed without converging
        if not converged:
            self.sfc_state.turbulence.obukhov_length[0] = L
            self.logger.warning(
                'Obukhov length did not converge. Final value = %f', L)

    def _solve_diffusion(
        self,
        state_field: np.ndarray,
        get_diffusivity: Callable,
        get_conductivity: Optional[Callable],
        sfc_boundary: float,
        field_name: str = 'field'
    ) -> None:
        """Solves a generic 1D diffusion equation using a theta scheme.

        This is a parameterized diffusion solver that handles both heat and
        moisture diffusion. The key difference is that moisture diffusion
        includes an additional hydraulic conductivity gradient term.

        Args:
            state_field: Reference to the field to update (temperature or
                moisture) [NDArray[np.float64]].
            get_diffusivity: Callable that computes diffusivity profile from
                soil moisture [Callable[[NDArray[np.float64]], NDArray[np.float64]]].
            get_conductivity: Optional callable that computes conductivity
                profile from soil moisture. None for heat diffusion,
                function for moisture diffusion
                [Optional[Callable[[NDArray[np.float64]], NDArray[np.float64]]]].
            sfc_boundary: Surface boundary value for Dirichlet BC [float].
            field_name: Name of the field for logging/documentation [str].

        Physics:
            - Diffusivity always depends on soil moisture (not the state
              being solved)
            - Conductivity (if present) also depends on soil moisture
            - Dirichlet BC at top (surface): uses sfc_boundary
            - Neumann BC at bottom: assumes zero gradient
            - Theta scheme parameterization:
              theta = 0.0 -> FTCS (explicit)
              theta = 0.5 -> Crank-Nicolson
              theta = 1.0 -> BTCS (implicit)
        """
        self.logger.debug('Solving %s diffusion', field_name)
        theta_b = self.input.numerics.diffusion_back_weight
        theta_f = 1.0 - theta_b
        nz = self.input.grid.nz
        dz = self.input.grid.z[0] - self.input.grid.z[1]
        dz2 = dz**2
        dt = self.tstep
        e, f, g, r = [np.zeros(nz - 1) for _ in range(4)]

        # Compute diffusivity using soil moisture (always, for both heat
        # and moisture)
        D = get_diffusivity(self.soil_state.moisture)
        D_mid = 0.5 * (D[:-1] + D[1:])

        # Compute conductivity terms if provided (moisture case)
        K_lin = None
        if get_conductivity is not None:
            K_hydraulic = get_conductivity(state_field)
            # Avoid division by zero in dry conditions
            K_lin = np.where(state_field > 1e-9, K_hydraulic / state_field, 0.0)

        # === First soil level below surface (i=0) ===
        Cp = dt * D_mid[0] / dz2
        Cm = dt * D_mid[1] / dz2

        # Backward (implicit) coefficients
        CBp = -theta_b * Cp
        CBm = -theta_b * Cm
        CB = 1.0 - CBp - CBm

        # Forward (explicit) coefficients
        CFp = theta_f * Cp
        CFm = theta_f * Cm
        CF = 1.0 - CFp - CFm

        # Add conductivity terms if applicable (moisture case)
        if K_lin is not None:
            Cpk = dt * K_lin[0] / (2 * dz)
            Cmk = dt * K_lin[2] / (2 * dz)
            CBpk = -theta_b * Cpk
            CBmk = -theta_b * Cmk
            CBp += CBpk
            CBm -= CBmk
            CFpk = theta_f * Cpk
            CFmk = theta_f * Cmk
            CFp += CFpk
            CFm -= CFmk

        f[0] = CB
        g[0] = CBm
        r[0] = (CFp * state_field[0] + CF * state_field[1] +
                CFm * state_field[2] - CBp * sfc_boundary)

        # === Interior soil levels ===
        for i in range(1, nz - 2):
            Cp = dt * D_mid[i] / dz2
            Cm = dt * D_mid[i + 1] / dz2

            CBp = -theta_b * Cp
            CBm = -theta_b * Cm
            CB = 1.0 - CBp - CBm

            CFp = theta_f * Cp
            CFm = theta_f * Cm
            CF = 1.0 - CFp - CFm

            # Add conductivity terms if applicable
            if K_lin is not None:
                Cpk = dt * K_lin[i] / (2 * dz)
                Cmk = dt * K_lin[i + 2] / (2 * dz)
                CBpk = -theta_b * Cpk
                CBmk = -theta_b * Cmk
                CBp += CBpk
                CBm -= CBmk
                CFpk = theta_f * Cpk
                CFmk = theta_f * Cmk
                CFp += CFpk
                CFm -= CFmk

            e[i] = CBp
            f[i] = CB
            g[i] = CBm
            r[i] = (CFp * state_field[i] + CF * state_field[i + 1] +
                    CFm * state_field[i + 2])

        # === Bottom level (Neumann BC: zero gradient) ===
        j = nz - 2
        Cp = dt * D_mid[j] / dz2
        Cm = dt * D_mid[j] / dz2

        CBp = -theta_b * Cp
        CBm = -theta_b * Cm
        CB = 1.0 - CBp - CBm

        CFp = theta_f * Cp
        CFm = theta_f * Cm
        CF = 1.0 - CFp - CFm

        # Add conductivity terms if applicable
        if K_lin is not None:
            Cpk = dt * K_lin[j] / (2 * dz)
            Cmk = dt * K_lin[j] / (2 * dz)
            CBpk = -theta_b * Cpk
            CBmk = -theta_b * Cmk
            CBp += CBpk
            CBm -= CBmk
            CFpk = theta_f * Cpk
            CFmk = theta_f * Cmk
            CFp += CFpk
            CFm -= CFmk

        e[j] = CBp - CBm
        f[j] = CB + 2.0 * CBm
        r[j] = ((CFp - CFm) * state_field[j] +
                (CF + 2.0 * CFm) * state_field[j + 1])

        # Solve and update
        state_field[0] = sfc_boundary
        state_field[1:] = solvers.tridiagonal(e, f, g, r)

    def _solve_diffusion_heat(self) -> None:
        """Solves the soil heat diffusion equation using a theta scheme."""
        self._solve_diffusion(
            state_field=self.soil_state.temperature,
            get_diffusivity=self.soil.diffusivity_thermal,
            get_conductivity=None,
            sfc_boundary=self.sfc_state.temperature,
            field_name='temperature'
        )

    def _solve_diffusion_mois(self) -> None:
        """Solves the soil moisture diffusion equation using a theta scheme."""
        self._solve_diffusion(
            state_field=self.soil_state.moisture,
            get_diffusivity=self.soil.diffusivity_moisture,
            get_conductivity=self.soil.conductivity_moisture,
            sfc_boundary=self.sfc_state.moisture,
            field_name='moisture'
        )
