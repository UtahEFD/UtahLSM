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

from .data_models import AtmosphericState, SoilState, SolverState, SurfaceState
from .exceptions import NamelistError, SolverError
from .physics import Radiation, Soil, Surface
from .physics.radiation.factory import get_radiation_model
from .physics.soil.factory import get_soil_model
from .physics.surface.factory import get_surface_model
from .util import constants as c
from .util import solvers
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
        self._load_atm_state(atm_state)

        # Update surface state from the top soil layer
        sfc_T = np.array(self.soil_state.temperature[0], copy=True)
        sfc_q = np.array(self.soil_state.moisture[0], copy=True)
        sfc_r = self.soil.surface_mixing_ratio(
            sfc_T, sfc_q, self.atm_state.pressure)
        self.sfc_state.temperature = sfc_T
        self.sfc_state.moisture = sfc_q
        self.sfc_state.specific_humidity = sfc_r

        # Run radiation model if configured
        if self.input.radiation.model:
            total_seconds = self.input.time.utc_start + runtime
            days_passed = int(total_seconds // 86400)
            current_utc = total_seconds % 86400
            days_per_year = (366 if self._is_leap_year(
                self.input.time.utc_year) else 365)
            julian_day = ((self.input.time.julian_day
                + days_passed - 1) % days_per_year) + 1
            self.atm_state.radiation_net = self._as_column_vector(
                self.rad.compute_net(
                    julian_day, current_utc, self.atm_state, self.sfc_state
                ),
                "radiation_net",
            )

    def run(self) -> None:
        """Runs the core model physics for a single time step.

        This includes solving the surface energy and moisture budgets and
        updating the soil profiles via diffusion solvers.
        """
        self.logger.info('Solving soil state')

        if (getattr(self.input, "numerics", None) is not None
            and getattr(self.input.numerics, "warm_start_turbulence", False)
            and not getattr(self, "_did_warm_start_turbulence", False)
        ):
            self._warm_start_turbulence()
            self._did_warm_start_turbulence = True

        # Set initial guesses for new surface temp and moisture
        self.sfc_state.temperature = np.array(
            self.soil_state.temperature[0], copy=True)
        self.sfc_state.moisture = np.array(
            self.soil_state.moisture[0], copy=True)

        # Solve surface energy and moisture budgets
        self._solve_surface_coupling()
        
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
        nx = getattr(self.input.grid, "nx", 1)
        ny = getattr(self.input.grid, "ny", 1)
        self.ncol = nx * ny
        self.soil_state: SoilState = replace(self.input.initial)
        self.soil_state.temperature = self._ensure_column_field(
            self.soil_state.temperature, "soil temperature")
        self.soil_state.moisture = self._ensure_column_field(
            self.soil_state.moisture, "soil moisture")

        self.sfc_state: SurfaceState = SurfaceState()
        self.sfc_state.temperature = np.zeros(self.ncol)
        self.sfc_state.moisture = np.zeros(self.ncol)
        self.sfc_state.specific_humidity = np.zeros(self.ncol)
        self.sfc_state.fluxes.kinematic_heat = np.zeros(self.ncol)
        self.sfc_state.fluxes.kinematic_moisture = np.zeros(self.ncol)
        self.sfc_state.fluxes.sensible_heat = np.zeros(self.ncol)
        self.sfc_state.fluxes.latent_heat = np.zeros(self.ncol)
        self.sfc_state.fluxes.ground_heat = np.zeros(self.ncol)
        self.sfc_state.turbulence.friction_velocity = np.zeros(self.ncol)
        self.sfc_state.turbulence.obukhov_length = np.zeros(self.ncol)

        self.atm_state: AtmosphericState = AtmosphericState()
        self.atm_state.wind_speed = np.zeros(self.ncol)
        self.atm_state.temperature = np.zeros(self.ncol)
        self.atm_state.specific_humidity = np.zeros(self.ncol)
        self.atm_state.pressure = np.zeros(self.ncol)
        self.atm_state.radiation_net = np.zeros(self.ncol)

        self.solver_state: SolverState = SolverState()
        self.solver_state.conductivity_thermal_mid = np.zeros(self.ncol)
        self._did_warm_start_turbulence: bool = False

    def _ensure_column_field(
        self, field: np.ndarray, name: str
    ) -> np.ndarray:
        """Ensures soil fields are shaped as (nz, ncol)."""
        data = np.asarray(field, dtype=float)
        nz = self.input.grid.nz
        ny = getattr(self.input.grid, "ny", 1)
        nx = getattr(self.input.grid, "nx", 1)
        ncol = getattr(self, "ncol", ny * nx)

        if data.ndim == 1:
            if data.shape[0] != nz:
                raise ValueError(
                    f"{name} length {data.shape[0]} does not match nz={nz}."
                )
            return data[:, None] if ncol > 1 else data[:, None]
        if data.ndim == 2:
            if data.shape == (nz, ncol):
                return data
            if data.shape == (nz, 1):
                return np.repeat(data, ncol, axis=1) if ncol > 1 else data
            raise ValueError(
                f"{name} shape {data.shape} does not match (nz, ncol)=({nz}, "
                f"{ncol})."
            )
        if data.ndim == 3:
            if data.shape == (nz, ny, nx):
                return data.reshape(nz, ncol)
            raise ValueError(
                f"{name} shape {data.shape} does not match (nz, ny, nx)=({nz}, "
                f"{ny}, {nx})."
            )
        raise ValueError(
            f"{name} has unsupported dimensions: {data.ndim}."
        )

    def _as_column_vector(self, value: object, name: str) -> np.ndarray:
        """Coerces scalars or (y, x) fields into (ncol,) arrays."""
        data = np.asarray(value, dtype=float)
        ny = getattr(self.input.grid, "ny", 1)
        nx = getattr(self.input.grid, "nx", 1)
        ncol = self.ncol

        if data.ndim == 0:
            return np.full(ncol, float(data))
        if data.ndim == 1:
            if data.size == ncol:
                return data
            if data.size == 1:
                return np.full(ncol, float(data[0]))
            raise ValueError(
                f"{name} length {data.size} does not match ncol={ncol}."
            )
        if data.ndim == 2:
            if data.shape == (ny, nx):
                return data.reshape(ncol)
            if data.shape == (1, 1):
                return np.full(ncol, float(data[0, 0]))
            raise ValueError(
                f"{name} shape {data.shape} does not match (ny, nx)=({ny}, "
                f"{nx})."
            )
        raise ValueError(
            f"{name} has unsupported dimensions: {data.ndim}."
        )

    def _copy_atm_state(self) -> AtmosphericState:
        """Returns a shallow copy of the current atmospheric state arrays."""
        return AtmosphericState(
            wind_speed=np.array(self.atm_state.wind_speed, copy=True),
            temperature=np.array(self.atm_state.temperature, copy=True),
            specific_humidity=np.array(self.atm_state.specific_humidity, copy=True),
            pressure=np.array(self.atm_state.pressure, copy=True),
            radiation_net=np.array(self.atm_state.radiation_net, copy=True),
        )

    def _load_atm_state(self, atm_state: AtmosphericState) -> None:
        """Loads atmospheric state data into column vectors."""
        self.atm_state.wind_speed = self._as_column_vector(
            atm_state.wind_speed, "wind_speed")
        self.atm_state.temperature = self._as_column_vector(
            atm_state.temperature, "temperature")
        self.atm_state.specific_humidity = self._as_column_vector(
            atm_state.specific_humidity, "specific_humidity")
        self.atm_state.pressure = self._as_column_vector(
            atm_state.pressure, "pressure")
        self.atm_state.radiation_net = self._as_column_vector(
            atm_state.radiation_net, "radiation_net")

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
                self.rad: Radiation = None
                self.logger.info('Using radiation forcing data')
            self.soil: Soil = get_soil_model(
                self.input.soil.model,
                self.input.soil_properties,
                self.input.soil_type_names,
                self.input.soil_properties_name
            )
            self.sfc: Surface = get_surface_model(self.input.surface)
        except NamelistError as e:
            self.logger.error('Failed to initialize physics modules: %s.', e)
            raise

    def _setup_output(self) -> None:
        """Sets up the output file dimensions and fields."""
        nx = getattr(self.input.grid, "nx", 1)
        ny = getattr(self.input.grid, "ny", 1)
        self.output_dims: dict = {
            't': 0,
            'z': self.input.grid.nz
        }
        if nx > 1 or ny > 1:
            self.output_dims.update({
                'y': ny,
                'x': nx
            })
        self.output.set_dims(self.output_dims)

        self.sfc_state.temperature = np.array(
            self.soil_state.temperature[0], copy=True)
        self.sfc_state.moisture = np.array(
            self.soil_state.moisture[0], copy=True)

        forcing0 = None
        if getattr(self.input, "forcing", None) is not None:
            atmos = getattr(self.input.forcing, "atmos", [])
            forcing0 = atmos[0] if atmos else None

        if forcing0 is not None and getattr(
            self.input.numerics, "initialize_surface_temperature_from_seb",False
        ):
            saved_atm_state = self._copy_atm_state()
            saved_tstep = getattr(self, "tstep", 0.0)
            self._load_atm_state(forcing0)
            self.tstep = float(getattr(self.input.forcing, "tstep", 0.0))

            self._solve_seb()
            self.soil_state.temperature[0] = self.sfc_state.temperature

            if hasattr(self.output, "outfile") and hasattr(
                self.output.outfile, "setncattr"
            ):
                self.output.outfile.setncattr(
                    "initial_surface_temperature",
                    "initialized from SEB using forcing[0]",
                )

            self._load_atm_state(saved_atm_state)
            self.tstep = saved_tstep

        if forcing0 is not None and getattr(
            self.input.numerics, "warm_start_turbulence", False
        ):
            self._warm_start_turbulence()
            self._did_warm_start_turbulence = True
            if hasattr(self.output, "outfile") and hasattr(
                self.output.outfile, "setncattr"
            ):
                self.output.outfile.setncattr(
                    "initial_diagnostics",
                    "warm_start_turbulence using forcing[0]",
                )

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

    def _warm_start_turbulence(self) -> None:
        """Warm-start MOST diagnostics using forcing[0] (offline mode only)."""
        if getattr(self.input, "forcing", None) is None:
            return
        atmos = getattr(self.input.forcing, "atmos", [])
        if not atmos:
            return
        if not hasattr(self, "_compute_fluxes"):
            return

        forcing0 = atmos[0]
        saved_atm_state = self._copy_atm_state()
        saved_tstep = getattr(self, "tstep", 0.0)

        self._load_atm_state(forcing0)
        self.tstep = float(getattr(self.input.forcing, "tstep", 0.0))

        sfc_T = np.array(self.soil_state.temperature[0], copy=True)
        sfc_q = np.array(self.soil_state.moisture[0], copy=True)
        self.sfc_state.temperature = sfc_T
        self.sfc_state.moisture = sfc_q
        self._compute_fluxes(sfc_T, sfc_q)

        self._load_atm_state(saved_atm_state)
        self.tstep = saved_tstep

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

        Uses a vectorized Brent's method to solve all columns simultaneously.
        This is significantly faster than column-by-column iteration.
        """
        # Calculate thermal conductivity for the entire soil column
        K_all = self.soil.conductivity_thermal(self.soil_state.moisture)
        self.solver_state.conductivity_thermal_mid = 0.5 * (K_all[0] + K_all[1])

        iter_max_bracket = self.input.numerics.iterations.seb_bracket
        iter_max_root = self.input.numerics.iterations.seb_root
        tol_root = self.input.numerics.tolerances.seb_root

        # Keep the initial Obukhov length fixed during root-finding
        initial_L = np.array(self.sfc_state.turbulence.obukhov_length, copy=True)

        # Create vectorized brackets around current temperatures
        current_T = self.sfc_state.temperature.copy()
        bracket_width = 1.0
        temp_a = current_T - bracket_width
        temp_b = current_T + bracket_width

        seb_a = self._compute_seb_vec(temp_a, initial_L)
        seb_b = self._compute_seb_vec(temp_b, initial_L)

        # Vectorized bracket expansion
        for _ in range(iter_max_bracket):
            needs_expansion = seb_a * seb_b > 0
            if not np.any(needs_expansion):
                break
            step = 5.0
            # Expand in the direction with smaller residual magnitude
            expand_left = needs_expansion & (np.abs(seb_a) < np.abs(seb_b))
            expand_right = needs_expansion & ~expand_left

            temp_a = np.where(expand_left, temp_a - step, temp_a)
            temp_b = np.where(expand_right, temp_b + step, temp_b)

            # Only recompute SEB for columns that changed
            if np.any(expand_left):
                seb_a_new = self._compute_seb_vec(temp_a, initial_L)
                seb_a = np.where(expand_left, seb_a_new, seb_a)
            if np.any(expand_right):
                seb_b_new = self._compute_seb_vec(temp_b, initial_L)
                seb_b = np.where(expand_right, seb_b_new, seb_b)
        else:
            # Check if any brackets failed after all iterations
            failed = seb_a * seb_b > 0
            if np.any(failed):
                bad_cols = np.where(failed)[0]
                raise SolverError(
                    f"SEB Bracket failed at cols={bad_cols.tolist()}. "
                    f"Residuals a: {seb_a[failed]}, b: {seb_b[failed]}"
                )

        # Vectorized Brent root-finding
        def seb_func(T: np.ndarray) -> np.ndarray:
            return self._compute_seb_vec(T, initial_L)

        self.sfc_state.temperature, converged = solvers.root_brent_vec(
            seb_func, temp_a, temp_b,
            iter_max=iter_max_root, tol=tol_root
        )

        if not converged.all():
            n_failed = int(np.sum(~converged))
            self.logger.warning(
                'SEB root-finding did not converge for %d of %d columns.',
                n_failed, self.ncol
            )

        # Final flux update with resolved temperatures (vectorized, once)
        self._compute_fluxes(self.sfc_state.temperature, self.sfc_state.moisture)

    def _compute_seb_vec(
        self,
        sfc_T: np.ndarray,
        initial_L: np.ndarray
    ) -> np.ndarray:
        """Computes the SEB residual for all columns without mutating state.

        This is a pure function used during root-finding iterations.

        Args:
            sfc_T: Surface temperature array [K] (ncol,).
            initial_L: Fixed Obukhov length array [m] (ncol,).

        Returns:
            Array of SEB residuals [W/m^2] (ncol,).
        """
        CP = c.thermodynamic.SPECIFIC_HEAT
        LV = c.thermodynamic.LATENT_HEAT_VAPORIZATION
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        VK = c.physical.VON_KARMAN
        G = c.physical.GRAVITY
        EVT = c.thermodynamic.EPSILON_VIRTUAL_TEMPERATURE
        TOL = self.input.numerics.tolerances.sfc_flux
        ITER_MAX = self.input.numerics.iterations.sfc_flux

        sfc_T = np.asarray(sfc_T)
        sfc_q = self.sfc_state.moisture
        atm_T = self.atm_state.temperature
        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        atm_ws = self.atm_state.wind_speed
        K_mid = self.solver_state.conductivity_thermal_mid
        soil_T1 = self.soil_state.temperature[1]

        z_m = self.input.surface.z_m
        z_o = self.input.surface.z_o
        z_s = self.input.surface.z_s
        z_t = self.input.surface.z_t
        zeta_max = self.input.surface.zeta_max
        gustiness = self.input.surface.gustiness
        gustiness_stable_only = self.input.surface.gustiness_stable_only
        dz = self.input.grid.z[0] - self.input.grid.z[1]

        ref_T = atm_T
        rho = atm_p / (RD * atm_T)

        # Surface-air specific humidity
        gnd_q = self.soil.surface_mixing_ratio(sfc_T, sfc_q, atm_p)

        # Ground heat flux
        ground_heat = K_mid * (sfc_T - soil_T1) / dz

        # Iterate for Obukhov length with fixed initial L
        L = np.array(initial_L, copy=True)
        converged = np.zeros_like(L, dtype=bool)

        for _ in range(ITER_MAX):
            fm = self.sfc.fm(z_m, z_o, L)
            fh = self.sfc.fh(z_s, z_t, L)

            wind_eff = atm_ws
            if gustiness > 0.0:
                gust = np.hypot(atm_ws, gustiness)
                if gustiness_stable_only:
                    wind_eff = np.where(L >= 0.0, gust, atm_ws)
                else:
                    wind_eff = gust

            ustar = wind_eff * fm
            flux_wT = (sfc_T - atm_T) * ustar * fh
            flux_wq = (gnd_q - atm_q) * ustar * fh
            flux_wTv = flux_wT + EVT * ref_T * flux_wq

            L_new = np.where(
                flux_wTv != 0.0,
                -(ustar**3) * ref_T / (VK * G * flux_wTv),
                1e6,
            )

            zeta = z_m / L_new
            L_new = np.where(zeta > zeta_max, z_m / zeta_max, L_new)
            L_new = np.where(zeta < -zeta_max, -z_m / zeta_max, L_new)

            diff = np.abs(L_new - L)
            converged |= diff <= TOL
            L = np.where(converged, L, L_new)

            if converged.all():
                break

        # Compute fluxes
        sensible = rho * CP * flux_wT
        latent = rho * LV * flux_wq

        # SEB residual
        seb = (
            self.atm_state.radiation_net
            - ground_heat
            - sensible
            - latent
        )

        return seb

    def _solve_smb(self) -> None:
        """Solves the Surface Moisture Budget (SMB).

        This function finds the soil moisture flux and surface evaporation.
        It then iteratively blends the two in time until convergence.
        """
        # Local constants and variables
        delta = 0.5#self.input.numerics.coupling_relaxation

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
        L = self.sfc_state.turbulence.obukhov_length
        ust = self.sfc_state.turbulence.friction_velocity
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
        flux_sm = RHO_W * K_mid * ((psi0 - psi1) / dz + 1.0)
        E = rho_a * self.sfc_state.fluxes.kinematic_moisture

        # Iteratively solve for moisture flux
        converged = np.zeros(self.ncol, dtype=bool)
        for _ in range(0, ITER_MAX):

            # New blended soil moisture flux
            flux_sm_last = flux_sm
            flux_sm_new = (1.0 - delta) * flux_sm_last - delta * E

            # New evaporation
            psi0_new = psi1 + dz * ((flux_sm_new / (RHO_W * K_mid)) - 1.0)
            psi0_new = np.minimum(psi0_new, psi_sat)
            sfc_moisture_new = self.soil.surface_water_content(psi0_new)
            gnd_q = self.soil.surface_mixing_ratio(
                sfc_T, sfc_moisture_new, atm_p)
            E_new = rho_a * (gnd_q - atm_q) * ust * fh

            K0_new = self.soil.conductivity_moisture(
                sfc_moisture_new, level=0)
            K_mid_new = 0.5 * (K0_new + K1)

            err = np.abs(E_new + flux_sm_new)
            newly_converged = err <= TOL
            converged |= newly_converged

            if converged.all():
                flux_sm = flux_sm_new
                psi0 = psi0_new
                self.sfc_state.moisture = sfc_moisture_new
                E = E_new
                K_mid = K_mid_new
                break

            mask = ~converged
            flux_sm = np.where(mask, flux_sm_new, flux_sm)
            psi0 = np.where(mask, psi0_new, psi0)
            self.sfc_state.moisture = np.where(
                mask, sfc_moisture_new, self.sfc_state.moisture)
            E = np.where(mask, E_new, E)
            K_mid = np.where(mask, K_mid_new, K_mid)

    def _compute_fluxes(self, sfc_T: np.ndarray, sfc_q: np.ndarray) -> None:
        """Computes surface fluxes using Monin-Obukhov Similarity Theory.

        This is an iterative process to find the friction velocity (ustar)
        and Obukhov length (L) that are consistent with the calculated
        sensible and latent heat fluxes. All columns are processed in parallel.

        Args:
            sfc_T: Surface temperature array [K] (ncol,).
            sfc_q: Surface moisture array [m^3/m^3] (ncol,).
        """
        VK = c.physical.VON_KARMAN
        G = c.physical.GRAVITY
        CP = c.thermodynamic.SPECIFIC_HEAT
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        LV = c.thermodynamic.LATENT_HEAT_VAPORIZATION
        EVT = c.thermodynamic.EPSILON_VIRTUAL_TEMPERATURE
        TOL = self.input.numerics.tolerances.sfc_flux
        ITER_MAX = self.input.numerics.iterations.sfc_flux

        sfc_T_vec = self._as_column_vector(sfc_T, "sfc_T")
        sfc_q_vec = self._as_column_vector(sfc_q, "sfc_q")
        atm_T = self.atm_state.temperature
        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        atm_ws = self.atm_state.wind_speed
        L = np.array(self.sfc_state.turbulence.obukhov_length, copy=True)
        K_mid = self.solver_state.conductivity_thermal_mid
        soil_T1 = self.soil_state.temperature[1]

        ref_T = atm_T
        rho = atm_p / (RD * atm_T)

        z_m = self.input.surface.z_m
        z_o = self.input.surface.z_o
        z_s = self.input.surface.z_s
        z_t = self.input.surface.z_t
        zeta_max = self.input.surface.zeta_max
        gustiness = self.input.surface.gustiness
        gustiness_stable_only = self.input.surface.gustiness_stable_only
        dz = self.input.grid.z[0] - self.input.grid.z[1]

        # Compute surface-air specific humidity
        gnd_q = self.soil.surface_mixing_ratio(sfc_T_vec, sfc_q_vec, atm_p)

        # Compute ground heat flux
        ground_heat = K_mid * (sfc_T_vec - soil_T1) / dz

        # Iteratively solve for fluxes and stability
        converged = np.zeros_like(L, dtype=bool)
        flux_wT = np.zeros_like(L)
        flux_wq = np.zeros_like(L)
        ustar = np.zeros_like(L)

        for _ in range(ITER_MAX):
            fm = self.sfc.fm(z_m, z_o, L)
            fh = self.sfc.fh(z_s, z_t, L)

            wind_eff = atm_ws
            if gustiness > 0.0:
                gust = np.hypot(atm_ws, gustiness)
                if gustiness_stable_only:
                    wind_eff = np.where(L >= 0.0, gust, atm_ws)
                else:
                    wind_eff = gust

            ustar = wind_eff * fm
            flux_wT = (sfc_T_vec - atm_T) * ustar * fh
            flux_wq = (gnd_q - atm_q) * ustar * fh
            flux_wTv = flux_wT + EVT * ref_T * flux_wq

            L_new = np.where(
                flux_wTv != 0.0,
                -(ustar**3) * ref_T / (VK * G * flux_wTv),
                1e6,
            )

            zeta = z_m / L_new
            L_new = np.where(zeta > zeta_max, z_m / zeta_max, L_new)
            L_new = np.where(zeta < -zeta_max, -z_m / zeta_max, L_new)

            diff = np.abs(L_new - L)
            converged |= diff <= TOL
            L = np.where(converged, L, L_new)

            if converged.all():
                break

        if not converged.all():
            self.logger.warning(
                'Obukhov length did not converge for %d columns.',
                int(np.sum(~converged)),
            )

        sensible = rho * CP * flux_wT
        latent = rho * LV * flux_wq

        self.sfc_state.turbulence.obukhov_length[:] = L
        self.sfc_state.turbulence.friction_velocity[:] = ustar
        self.sfc_state.fluxes.ground_heat[:] = ground_heat
        self.sfc_state.fluxes.kinematic_heat[:] = flux_wT
        self.sfc_state.fluxes.kinematic_moisture[:] = flux_wq
        self.sfc_state.fluxes.sensible_heat[:] = sensible
        self.sfc_state.fluxes.latent_heat[:] = latent

    def _solve_surface_coupling(self) -> None:
        """Iteratively solves the coupled SEB and SMB.

        This method performs a Picard iteration, alternating between solving
        the Surface Energy Budget (SEB) for temperature and the Surface
        Moisture Budget (SMB) for moisture until both state variables converge.
        """
        max_outer_iter = self.input.numerics.iterations.coupling
        tol_temp = self.input.numerics.tolerances.coupling_temp
        tol_mois = self.input.numerics.tolerances.coupling_mois

        for i in range(max_outer_iter):
            # Store previous states to check convergence
            prev_T = np.array(self.sfc_state.temperature, copy=True)
            prev_q = np.array(self.sfc_state.moisture, copy=True)

            # Solve SEB for temperature
            self._solve_seb()

            # Solve SMB for moisture
            self._solve_smb()

            # Check convergence
            diff_T = np.abs(self.sfc_state.temperature - prev_T)
            diff_q = np.abs(self.sfc_state.moisture - prev_q)

            if np.all(diff_T < tol_temp) and np.all(diff_q < tol_mois):
                self.logger.debug(
                    'Surface coupling converged in %d iterations.', i + 1)
                self._compute_fluxes(
                    self.sfc_state.temperature, self.sfc_state.moisture
                )
                return

        self.logger.warning(
            'Surface coupling did not converge after %d iterations. '
            'dT: %.4f, dq: %.4e', max_outer_iter,
            float(np.max(diff_T)), float(np.max(diff_q)))
        self._compute_fluxes(self.sfc_state.temperature, self.sfc_state.moisture)
    
    def _solve_diffusion(self,state_field: np.ndarray,get_diffusivity: Callable,
                         get_conductivity: Optional[Callable],
                         sfc_boundary: float,field_name: str = 'field') -> None:
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
        field = np.asarray(state_field, dtype=float)
        if field.ndim == 1:
            field = field[:, None]
        if field.ndim != 2 or field.shape[0] != nz:
            raise ValueError(
                f"{field_name} has unexpected shape {field.shape}."
            )
        ncol = field.shape[1]
        sfc_boundary_vec = self._as_column_vector(
            sfc_boundary, f"{field_name} boundary")
        if sfc_boundary_vec.size != ncol:
            raise ValueError(
                f"{field_name} boundary size {sfc_boundary_vec.size} does not "
                f"match ncol={ncol}."
            )

        e, f, g, r = [np.zeros((nz - 1, ncol)) for _ in range(4)]

        # Compute diffusivity using soil moisture (always, for both heat
        # and moisture)
        D = get_diffusivity(self.soil_state.moisture)
        D_mid = 0.5 * (D[:-1] + D[1:])

        # Compute conductivity terms if provided (moisture case)
        K_lin = None
        if get_conductivity is not None:
            K_hydraulic = get_conductivity(field)
            # Avoid division by zero in dry conditions
            K_lin = np.where(field > 1e-9, K_hydraulic / field, 0.0)

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
        r[0] = (CFp * field[0] + CF * field[1] +
                CFm * field[2] - CBp * sfc_boundary_vec)

        # === Interior soil levels (Vectorized) ===
        # Define slices to represent indices i, i+1, and i+2
        # Original loop: for i in range(1, nz - 2)
        # indices: 1, 2, ..., nz-3
        idx     = slice(1, nz - 2)  # corresponds to i
        idx_p1  = slice(2, nz - 1)  # corresponds to i+1
        idx_p2  = slice(3, nz)      # corresponds to i+2

        # Compute diffusion coefficients for all interior points
        # D_mid is size (nz-1), so we slice up to nz-2
        Cp = dt * D_mid[idx] / dz2
        Cm = dt * D_mid[idx_p1] / dz2

        # Backward (implicit) coefficients
        CBp = -theta_b * Cp
        CBm = -theta_b * Cm
        CB  = 1.0 - CBp - CBm

        # Forward (explicit) coefficients
        CFp = theta_f * Cp
        CFm = theta_f * Cm
        CF  = 1.0 - CFp - CFm

        # Add conductivity terms if applicable (moisture case)
        if K_lin is not None:
            # K_lin is size (nz)
            Cpk = dt * K_lin[idx] / (2 * dz)
            Cmk = dt * K_lin[idx_p2] / (2 * dz)
            CBpk = -theta_b * Cpk
            CBmk = -theta_b * Cmk
            CBp += CBpk
            CBm -= CBmk
            CFpk = theta_f * Cpk
            CFmk = theta_f * Cmk
            CFp += CFpk
            CFm -= CFmk

        # Assign coefficients to tridiagonal matrix arrays
        e[idx] = CBp
        f[idx] = CB
        g[idx] = CBm

        # Compute the Right Hand Side (RHS) vector r
        # state_field is size (nz)
        r[idx] = (CFp * field[idx] +
                  CF  * field[idx_p1] +
                  CFm * field[idx_p2])

        # === Bottom level (Neumann BC: zero gradient) ===
        j = nz - 2
        Cp = dt * D_mid[j] / dz2
        Cm = dt * D_mid[j] / dz2

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

        # Assign coefficients to tridiagonal matrix arrays
        e[j] = CBp - CBm
        f[j] = CB + 2.0 * CBm

        # Compute the Right Hand Side (RHS) vector r
        # state_field is size (nz)
        r[j] = ((CFp - CFm) * field[j] +
                (CF + 2.0 * CFm) * field[j + 1])

        # Solve and update
        field[0] = sfc_boundary_vec
        field[1:] = solvers.tridiagonal(e, f, g, r)

        if state_field.ndim == 1:
            state_field[:] = field[:, 0]
        else:
            state_field[:] = field

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
