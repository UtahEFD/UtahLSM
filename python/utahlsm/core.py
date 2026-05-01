#!/usr/bin/env python
#
# UtahLSM
#
# Copyright (c) 2017–2026 Jeremy A. Gibbs
# Copyright (c) 2017–2026 Rob Stoll
# Copyright (c) 2017–2026 Eric Pardyjak
# Copyright (c) 2017–2026 Pete Willemsen
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
    AtmosphericState,
    CanopyState,
    SoilState,
    SolverState,
    SurfaceState,
)
from .exceptions import NamelistError, SolverError
from .physics import Canopy, Radiation, Soil, Surface, thermo
from .physics.canopy.factory import get_canopy_model
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
        canopy: The selected canopy model instance, or None in bare-soil mode.
        canopy_state: Canopy / vegetation diagnostics updated each step.
        tstep: The current model time step [s].
        logger: A logger instance for this class.
        output_dims: A dictionary of dimensions for the output file.
        full_output_fields: A dictionary of all available output fields for the
            current model configuration.
        output_fields: A dictionary of fields selected for normal output.
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
        """Updates the model with new atmospheric forcing data.

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
        sfc_theta = np.array(self.soil_state.moisture[0], copy=True)
        sfc_q = self.soil.surface_specific_humidity(
            sfc_T, sfc_theta, self.atm_state.pressure)
        self.sfc_state.temperature = sfc_T
        self.sfc_state.moisture = sfc_theta
        self.sfc_state.specific_humidity = sfc_q

        # Run radiation model if configured. The model populates the four
        # component fields (sw_in, sw_out, lw_in, lw_out) on atm_state;
        # radiation_net is derived from them so the SEB and any consumers
        # of the components remain mutually consistent.
        if self.input.radiation.model:
            total_seconds = self.input.time.utc_start + runtime
            days_passed = int(total_seconds // 86400)
            current_utc = total_seconds % 86400
            days_per_year = (366 if self._is_leap_year(
                self.input.time.utc_year) else 365)
            julian_day = ((self.input.time.julian_day
                + days_passed - 1) % days_per_year) + 1
            sw_in, sw_out, lw_in, lw_out = self.rad.compute_components(
                julian_day, current_utc, self.atm_state, self.sfc_state
            )
            self.atm_state.sw_in = self._as_column_vector(sw_in, "sw_in")
            self.atm_state.sw_out = self._as_column_vector(sw_out, "sw_out")
            self.atm_state.lw_in = self._as_column_vector(lw_in, "lw_in")
            self.atm_state.lw_out = self._as_column_vector(lw_out, "lw_out")
            self.atm_state.radiation_net = (
                self.atm_state.sw_in - self.atm_state.sw_out
                + self.atm_state.lw_in - self.atm_state.lw_out
            )

    def run(self) -> None:
        """Runs the core model physics for a single time step.

        This includes solving the surface energy and moisture budgets and
        updating the soil profiles via diffusion solvers.
        """
        self.logger.info('Solving soil state')
        if hasattr(self, 'soil'):
            self.soil._validate_moisture_bounds(self.soil_state.moisture)

        if (self.input.numerics.warm_start_turbulence
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
        self.output.save(self.output_fields, step_count, runtime)

    # --- Internal Methods ---

    def _setup_states(self) -> None:
        """Initializes all state containers for the model."""
        self.logger.info('Setting up initial states')
        nx = self.input.grid.nx
        ny = self.input.grid.ny
        self.ncol = nx * ny
        self.soil_state: SoilState = replace(self.input.initial)
        self.soil_state.temperature = self._ensure_column_field(
            self.soil_state.temperature, "soil temperature")
        self.soil_state.moisture = self._ensure_column_field(
            self.soil_state.moisture, "soil moisture")

        self.sfc_state: SurfaceState = SurfaceState()
        self.sfc_state.temperature = np.zeros(self.ncol)
        self.sfc_state.soil_top_temperature = np.zeros(self.ncol)
        self.sfc_state.moisture = np.zeros(self.ncol)
        self.sfc_state.specific_humidity = np.zeros(self.ncol)
        self.sfc_state.fluxes.kinematic_heat = np.zeros(self.ncol)
        self.sfc_state.fluxes.kinematic_moisture = np.zeros(self.ncol)
        self.sfc_state.fluxes.sensible_heat = np.zeros(self.ncol)
        self.sfc_state.fluxes.latent_heat = np.zeros(self.ncol)
        self.sfc_state.fluxes.ground_heat = np.zeros(self.ncol)
        self.sfc_state.turbulence.friction_velocity = np.zeros(self.ncol)
        self.sfc_state.turbulence.obukhov_length = np.full(self.ncol, 1e6)

        self.atm_state: AtmosphericState = AtmosphericState()
        self.atm_state.wind_speed = np.zeros(self.ncol)
        self.atm_state.temperature = np.zeros(self.ncol)
        self.atm_state.specific_humidity = np.zeros(self.ncol)
        self.atm_state.pressure = np.zeros(self.ncol)
        self.atm_state.sw_in = np.zeros(self.ncol)
        self.atm_state.sw_out = np.zeros(self.ncol)
        self.atm_state.lw_in = np.zeros(self.ncol)
        self.atm_state.lw_out = np.zeros(self.ncol)
        self.atm_state.radiation_net = np.zeros(self.ncol)

        self.solver_state: SolverState = SolverState()
        self.solver_state.conductivity_thermal_mid = np.zeros(self.ncol)
        nz_diff = self.input.grid.nz - 1
        self.solver_state.diffusion_e = np.zeros((nz_diff, self.ncol))
        self.solver_state.diffusion_f = np.zeros((nz_diff, self.ncol))
        self.solver_state.diffusion_g = np.zeros((nz_diff, self.ncol))
        self.solver_state.diffusion_r = np.zeros((nz_diff, self.ncol))
        self._did_warm_start_turbulence: bool = False

        # Canopy defaults to bare-soil until _setup_physics instantiates
        # a concrete model from the namelist. Stubbing it here keeps
        # partial-initialization test paths that skip _setup_physics safe.
        self.canopy: Optional[Canopy] = None
        nz = self.input.grid.nz
        self.canopy_state: CanopyState = CanopyState(
            resistance=np.full(self.ncol, np.inf),
            theta_root=np.zeros(self.ncol),
            transpiration=np.zeros(self.ncol),
            evap_soil=np.zeros(self.ncol),
            latent_veg=np.zeros(self.ncol),
            latent_soil=np.zeros(self.ncol),
            root_uptake=np.zeros((nz, self.ncol)),
        )

    def _ensure_column_field(
        self, field: np.ndarray, name: str
    ) -> np.ndarray:
        """Ensures soil fields are shaped as (nz, ncol)."""
        data = np.asarray(field, dtype=float)
        nz = self.input.grid.nz
        ny = self.input.grid.ny
        nx = self.input.grid.nx
        ncol = self.ncol

        if data.ndim == 1:
            if data.shape[0] != nz:
                raise ValueError(
                    f"{name} length {data.shape[0]} does not match nz={nz}."
                )
            return data[:, None]
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
        ny = self.input.grid.ny
        nx = self.input.grid.nx
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
            sw_in=np.array(self.atm_state.sw_in, copy=True),
            sw_out=np.array(self.atm_state.sw_out, copy=True),
            lw_in=np.array(self.atm_state.lw_in, copy=True),
            lw_out=np.array(self.atm_state.lw_out, copy=True),
            radiation_net=np.array(self.atm_state.radiation_net, copy=True),
        )

    def _load_atm_state(self, atm_state: AtmosphericState) -> None:
        """Loads atmospheric state data into column vectors.

        Forcing-driven runs supply the four radiation components; the
        net is derived here so SEB residuals stay consistent with
        component-level diagnostics (e.g. Jarvis f1 reading sw_in).
        Built-in radiation runs overwrite all five fields in
        :meth:`update` after this method returns.
        """
        self.atm_state.wind_speed = self._as_column_vector(
            atm_state.wind_speed, "wind_speed")
        self.atm_state.temperature = self._as_column_vector(
            atm_state.temperature, "temperature")
        self.atm_state.specific_humidity = self._as_column_vector(
            atm_state.specific_humidity, "specific_humidity")
        self.atm_state.pressure = self._as_column_vector(
            atm_state.pressure, "pressure")
        self.atm_state.sw_in = self._as_column_vector(
            atm_state.sw_in, "sw_in")
        self.atm_state.sw_out = self._as_column_vector(
            atm_state.sw_out, "sw_out")
        self.atm_state.lw_in = self._as_column_vector(
            atm_state.lw_in, "lw_in")
        self.atm_state.lw_out = self._as_column_vector(
            atm_state.lw_out, "lw_out")
        self.atm_state.radiation_net = (
            self.atm_state.sw_in - self.atm_state.sw_out
            + self.atm_state.lw_in - self.atm_state.lw_out
        )

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
            self.canopy: Optional[Canopy] = get_canopy_model(
                self.input.canopy, self.input.grid.z, self.ncol
            )
        except NamelistError as e:
            self.logger.error('Failed to initialize physics modules: %s.', e)
            raise

    def _setup_output(self) -> None:
        """Sets up the output file dimensions and fields."""
        nx = self.input.grid.nx
        ny = self.input.grid.ny
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
        self.sfc_state.soil_top_temperature = np.array(
            self.soil_state.temperature[0], copy=True)
        self.sfc_state.moisture = np.array(
            self.soil_state.moisture[0], copy=True)

        forcing0 = None
        if self.input.forcing is not None:
            atmos = self.input.forcing.atmos
            forcing0 = atmos[0] if atmos else None

        if forcing0 is not None and (
            self.input.numerics.initialize_surface_temperature_from_seb
        ):
            saved_atm_state = self._copy_atm_state()
            saved_tstep = getattr(self, "tstep", 0.0)
            self._load_atm_state(forcing0)
            self.tstep = float(self.input.forcing.tstep)

            self._refresh_canopy_diagnostics()
            self._solve_seb()
            # The top soil cell holds the soil-top temperature, which
            # equals the radiative skin only when r_canopy_thermal = 0.
            self.soil_state.temperature[0] = self.sfc_state.soil_top_temperature

            if hasattr(self.output, "outfile") and hasattr(
                self.output.outfile, "setncattr"
            ):
                self.output.outfile.setncattr(
                    "initial_surface_temperature",
                    "initialized from SEB using forcing[0]",
                )

            self._load_atm_state(saved_atm_state)
            self.tstep = saved_tstep

        if forcing0 is not None and self.input.numerics.warm_start_turbulence:
            self._warm_start_turbulence()
            self._did_warm_start_turbulence = True
            if hasattr(self.output, "outfile") and hasattr(
                self.output.outfile, "setncattr"
            ):
                self.output.outfile.setncattr(
                    "initial_diagnostics",
                    "warm_start_turbulence using forcing[0]",
                )

        self.full_output_fields = {
            'ust': self.sfc_state.turbulence.friction_velocity,
            'obl': self.sfc_state.turbulence.obukhov_length,
            'shf': self.sfc_state.fluxes.sensible_heat,
            'lhf': self.sfc_state.fluxes.latent_heat,
            'ghf': self.sfc_state.fluxes.ground_heat,
            'soil_z': self.input.grid.z,
            'soil_T': self.soil_state.temperature,
            'soil_q': self.soil_state.moisture,
        }
        if getattr(self, 'canopy', None) is not None:
            self.full_output_fields.update({
                'r_s': self.canopy_state.resistance,
                'theta_root': self.canopy_state.theta_root,
                'lhf_soil': self.canopy_state.latent_soil,
                'lhf_veg': self.canopy_state.latent_veg,
            })
        self.output_fields = self._select_output_fields(self.full_output_fields)
        self.output.set_fields(self.output_fields)
        self.output.save(self.output_fields, 0, 0, initial=True)

    def _select_output_fields(
            self, available_fields: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """Selects the configured subset of output fields for this run.

        Args:
            available_fields: All fields available from the current model
                configuration.

        Returns:
            The subset of fields to write during normal output.

        Raises:
            ValueError: If the requested field list mixes ``all`` with explicit
                names, contains unknown fields, or requests fields that are not
                available in the current configuration.
        """
        if not self.input.output.save:
            return {}

        requested = list(self.input.output.fields)
        if 'all' in requested:
            if len(requested) != 1:
                raise ValueError(
                    "output.fields must be ['all'] or an explicit field list."
                )
            return dict(available_fields)

        supported_fields = Output.supported_fields()
        unknown_fields = [
            field for field in requested if field not in supported_fields
        ]
        if unknown_fields:
            raise ValueError(
                f'Unknown output field(s) requested: {unknown_fields}. '
                f'Supported fields: {sorted(supported_fields)}.'
            )

        unavailable_fields = [
            field for field in requested if field not in available_fields
        ]
        if unavailable_fields:
            raise ValueError(
                f'Output field(s) not available for this configuration: '
                f'{unavailable_fields}. Available fields: '
                f'{sorted(available_fields)}.'
            )

        return {field: available_fields[field] for field in requested}

    def _warm_start_turbulence(self) -> None:
        """Warm-start MOST diagnostics using forcing[0] (offline mode only)."""
        if self.input.forcing is None:
            return
        atmos = self.input.forcing.atmos
        if not atmos:
            return
        forcing0 = atmos[0]
        saved_atm_state = self._copy_atm_state()
        saved_tstep = getattr(self, "tstep", 0.0)

        self._load_atm_state(forcing0)
        self.tstep = float(self.input.forcing.tstep)

        sfc_T = np.array(self.soil_state.temperature[0], copy=True)
        sfc_q = np.array(self.soil_state.moisture[0], copy=True)
        self.sfc_state.temperature = sfc_T
        self.sfc_state.moisture = sfc_q
        self._refresh_canopy_diagnostics()
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

    def _refresh_canopy_diagnostics(self) -> None:
        """Refreshes r_s and θ_root from the current outer-iteration state.

        Called once per outer SEB/SMB coupling pass so that subsequent SEB
        root-finds and SMB Brent solves see a fixed stomatal resistance
        (per the fast-response design axiom - no per-Brent-evaluation
        re-evaluation of f1-f4).
        """
        canopy = getattr(self, 'canopy', None)
        if canopy is None:
            return
        theta_wilt = self.soil.theta_wilt
        theta_fc = self.soil.theta_fc
        r_s = canopy.compute_resistance(
            self.atm_state, self.sfc_state, self.soil_state,
            theta_wilt, theta_fc,
        )
        self.canopy_state.resistance[:] = r_s
        self.canopy_state.theta_root[:] = canopy.root_zone_mean(
            self.soil_state.moisture
        )

    def _partition_flux_wq(
        self,
        sfc_T: np.ndarray,
        gnd_q: np.ndarray,
        atm_q: np.ndarray,
        atm_p: np.ndarray,
        ustar: np.ndarray,
        fh: np.ndarray,
        cols: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Kinematic moisture flux with optional canopy partition.

        Bare-soil mode: single-source form
        ``flux_wq = (gnd_q - atm_q) · u* · f_h``.

        Canopy-active mode: two-source (Noilhan-Planton) combining bare-
        soil evaporation through the aerodynamic resistance with
        transpiration through ``r_a + r_s``. Both sources use the same
        surface temperature (big-leaf, single-T axiom).
        """
        canopy = getattr(self, 'canopy', None)
        if canopy is None:
            return (gnd_q - atm_q) * ustar * fh

        if cols is not None:
            f_veg = canopy.veg_fraction[cols]
            r_s = self.canopy_state.resistance[cols]
        else:
            f_veg = canopy.veg_fraction
            r_s = self.canopy_state.resistance

        q_sat = thermo.saturation_specific_humidity(sfc_T, atm_p)
        u_fh = ustar * fh
        e_soil = (1.0 - f_veg) * (gnd_q - atm_q) * u_fh
        # 1/(ra + rs) normalized: flux = f_veg·(q_sat - q_a)·u_fh/(1 + u_fh·r_s)
        denom = 1.0 + u_fh * r_s
        t_veg = f_veg * (q_sat - atm_q) * u_fh / denom
        return e_soil + t_veg

    def _finalize_canopy_partition(self) -> None:
        """Splits the converged total LH into soil and canopy components.

        Also records per-layer root uptake rate [m^3/m^3/s] used as a
        sink term in the soil moisture RHS. Layer 0's root fraction
        is folded into layer 1 so the surface BC is unaffected by
        transpiration (the SMB already balances the top layer).
        """
        LV = c.thermodynamic.LATENT_HEAT_VAPORIZATION
        RHO_W = c.water.DENSITY
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        EVT = c.thermodynamic.EPSILON_VIRTUAL_TEMPERATURE

        canopy = getattr(self, 'canopy', None)
        if canopy is None:
            if hasattr(self, 'canopy_state'):
                self.canopy_state.evap_soil[:] = 0.0
                self.canopy_state.transpiration[:] = 0.0
                self.canopy_state.latent_soil[:] = (
                    self.sfc_state.fluxes.latent_heat)
                self.canopy_state.latent_veg[:] = 0.0
                self.canopy_state.root_uptake[:] = 0.0
            return

        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        atm_T = self.atm_state.temperature
        sfc_T = self.sfc_state.temperature
        sfc_q = self.sfc_state.moisture
        rho_a = atm_p / (RD * atm_T * (1.0 + EVT * atm_q))

        ust = self.sfc_state.turbulence.friction_velocity
        fh = self.sfc.fh(
            self.input.surface.z_s, self.input.surface.z_t,
            self.sfc_state.turbulence.obukhov_length,
        )
        gnd_q = self.soil.surface_specific_humidity(sfc_T, sfc_q, atm_p)
        q_sat = thermo.saturation_specific_humidity(sfc_T, atm_p)
        f_veg = canopy.veg_fraction
        r_s = self.canopy_state.resistance
        u_fh = ust * fh

        E_soil_kin = (1.0 - f_veg) * (gnd_q - atm_q) * u_fh
        T_veg_kin = f_veg * (q_sat - atm_q) * u_fh / (1.0 + u_fh * r_s)

        E_soil_mass = rho_a * E_soil_kin          # [kg/m^2/s]
        T_veg_mass = rho_a * T_veg_kin            # [kg/m^2/s]

        self.canopy_state.evap_soil[:] = E_soil_mass
        self.canopy_state.transpiration[:] = T_veg_mass
        self.canopy_state.latent_soil[:] = LV * E_soil_mass
        self.canopy_state.latent_veg[:] = LV * T_veg_mass

        # Per-layer uptake rate [m^3/m^3/s] folded into layers 1..nz-1.
        nz = self.input.grid.nz
        root_frac = canopy.root_fraction.copy()  # (nz, ncol)
        if nz > 1:
            root_frac[1] += root_frac[0]
            root_frac[0] = 0.0
        dz = self.input.grid.z[0] - self.input.grid.z[1]
        # T_veg_mass broadcast over layers times root fraction / (rho_w dz).
        self.canopy_state.root_uptake[:] = (
            root_frac * T_veg_mass[None, :] / (RHO_W * dz)
        )

    def _solve_seb(self) -> None:
        """Solves the Surface Energy Budget (SEB) to find surface temperature.

        Uses a vectorized Brent's method to solve all columns simultaneously.
        This is significantly faster than column-by-column iteration.
        """
        # Calculate thermal conductivity for the entire soil column.
        # Harmonic mean treats the two half-layers as series resistors,
        # which is the physically consistent effective conductivity for
        # 1D Fourier conduction between layer centers.
        K_all = self.soil.conductivity_thermal(self.soil_state.moisture)
        K0, K1 = K_all[0], K_all[1]
        self.solver_state.conductivity_thermal_mid = np.where(
            (K0 + K1) > 0.0, 2.0 * K0 * K1 / (K0 + K1), 0.0
        )

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
                idx = np.where(expand_left)[0]
                seb_a[idx] = self._compute_seb_vec(
                    temp_a[idx], initial_L[idx], cols=idx)
            if np.any(expand_right):
                idx = np.where(expand_right)[0]
                seb_b[idx] = self._compute_seb_vec(
                    temp_b[idx], initial_L[idx], cols=idx)
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

    def _solve_most(
        self,
        sfc_T: np.ndarray,
        sfc_q: np.ndarray,
        L_init: np.ndarray,
        max_iter: int,
        cols: Optional[np.ndarray] = None,
    ) -> tuple:
        """Computes surface fluxes via Monin-Obukhov Similarity Theory.

        This is the shared MOST solver used by both the SEB residual
        computation and the final flux update. All columns are processed
        in parallel.

        Args:
            sfc_T: Surface temperature array [K] (ncol,) or (len(cols),).
            sfc_q: Surface moisture array [m^3/m^3], same shape as sfc_T.
            L_init: Initial Obukhov length array [m], same shape as sfc_T.
            max_iter: Maximum number of MOST iterations. Use 1 for a
                single-pass evaluation with fixed L.
            cols: Optional array of column indices to evaluate. When
                provided, only those columns are read from the model
                state arrays. sfc_T, sfc_q, and L_init must already be
                sliced to match len(cols).

        Returns:
            Tuple of (ustar, flux_wT, flux_wq, ground_heat,
                      soil_top_T, L, sensible, latent, converged).
        """
        CP = c.thermodynamic.SPECIFIC_HEAT
        LV = c.thermodynamic.LATENT_HEAT_VAPORIZATION
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        VK = c.physical.VON_KARMAN
        G = c.physical.GRAVITY
        EVT = c.thermodynamic.EPSILON_VIRTUAL_TEMPERATURE
        TOL = self.input.numerics.tolerances.sfc_flux

        canopy = getattr(self, 'canopy', None)
        if cols is not None:
            atm_T = self.atm_state.temperature[cols]
            atm_p = self.atm_state.pressure[cols]
            atm_q = self.atm_state.specific_humidity[cols]
            atm_ws = self.atm_state.wind_speed[cols]
            K_mid = self.solver_state.conductivity_thermal_mid[cols]
            soil_T1 = self.soil_state.temperature[1, cols]
            if canopy is not None:
                f_veg = canopy.veg_fraction[cols]
                r_g_aero = canopy.r_ground[cols]
            else:
                f_veg = np.zeros_like(soil_T1)
                r_g_aero = np.zeros_like(soil_T1)
        else:
            atm_T = self.atm_state.temperature
            atm_p = self.atm_state.pressure
            atm_q = self.atm_state.specific_humidity
            atm_ws = self.atm_state.wind_speed
            K_mid = self.solver_state.conductivity_thermal_mid
            soil_T1 = self.soil_state.temperature[1]
            if canopy is not None:
                f_veg = canopy.veg_fraction
                r_g_aero = canopy.r_ground
            else:
                f_veg = np.zeros_like(soil_T1)
                r_g_aero = np.zeros_like(soil_T1)

        z_m = self.input.surface.z_m
        z_o = self.input.surface.z_o
        z_s = self.input.surface.z_s
        z_t = self.input.surface.z_t
        zeta_max = self.input.surface.zeta_max
        gustiness = self.input.surface.gustiness
        gustiness_stable_only = self.input.surface.gustiness_stable_only
        dz = self.input.grid.z[0] - self.input.grid.z[1]

        ref_T = atm_T
        # Moist-air density uses virtual temperature (Tv = T(1 + 0.608 q)).
        # Tv > T for humid air, so dry-T gives ~0.5–1% high LE in practice.
        Tv = atm_T * (1.0 + EVT * atm_q)
        rho = atm_p / (RD * Tv)

        gnd_q = self.soil.surface_specific_humidity(sfc_T, sfc_q, atm_p)

        # Ground heat flux: skin → soil-top conduction with an in-canopy
        # aerodynamic resistance acting in series. Converting the
        # aerodynamic resistance r_ground [s/m] to a thermal resistance
        # uses the volumetric heat capacity of air ρ·Cp [J/m^3/K], so
        # r_canopy_thermal [K·m^2/W] = r_ground / (ρ·Cp). The resistance
        # is scaled by veg_fraction so a bare patch (f_veg=0) recovers
        # the original direct-conduction limit.
        r_soil = dz / np.maximum(K_mid, 1e-12)
        r_canopy_thermal = f_veg * r_g_aero / (rho * CP)
        r_total = r_soil + r_canopy_thermal
        ground_heat = (sfc_T - soil_T1) / r_total
        # Soil-top temperature: what the conduction equation alone (no
        # canopy) would need at the soil surface to carry G into layer
        # 1. Equals sfc_T when r_canopy_thermal = 0; otherwise cooler.
        soil_top_T = sfc_T - ground_heat * r_canopy_thermal

        L = np.array(L_init, copy=True)
        converged = np.zeros_like(L, dtype=bool)
        ustar = np.zeros_like(L)
        flux_wT = np.zeros_like(L)
        flux_wq = np.zeros_like(L)

        for _ in range(max_iter):
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
            flux_wq = self._partition_flux_wq(
                sfc_T, gnd_q, atm_q, atm_p, ustar, fh, cols=cols
            )
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

        sensible = rho * CP * flux_wT
        latent = rho * LV * flux_wq

        return (ustar, flux_wT, flux_wq, ground_heat, soil_top_T,
                L, sensible, latent, converged)

    def _compute_seb_vec(
        self,
        sfc_T: np.ndarray,
        initial_L: np.ndarray,
        cols: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Computes the SEB residual without mutating state.

        Uses a single-pass MOST evaluation with fixed Obukhov length for
        efficiency during root-finding iterations. The final consistent L
        is resolved by _compute_fluxes after the root is found.

        Args:
            sfc_T: Surface temperature array [K] (ncol,) or (len(cols),).
            initial_L: Fixed Obukhov length array [m], same shape as sfc_T.
            cols: Optional array of column indices to evaluate. When
                provided, only those columns are computed.

        Returns:
            Array of SEB residuals [W/m^2].
        """
        sfc_T = np.asarray(sfc_T)
        if cols is not None:
            sfc_q = self.sfc_state.moisture[cols]
            rad_net = self.atm_state.radiation_net[cols]
        else:
            sfc_q = self.sfc_state.moisture
            rad_net = self.atm_state.radiation_net

        _, _, _, ground_heat, _, _, sensible, latent, _ = self._solve_most(
            sfc_T, sfc_q, initial_L, max_iter=1, cols=cols
        )

        return rad_net - ground_heat - sensible - latent

    def _solve_smb(self) -> None:
        """Solves the Surface Moisture Budget (SMB) for surface moisture.

        Finds θ_sfc that balances the evaporative demand against the
        Darcy flux from the subsurface via bracketed Brent root-finding
        on θ_sfc ∈ [θ_residual, θ_porosity]. The residual is

            R(θ_sfc) = E(θ_sfc, T_sfc)
                       + rho_w K_mid (θ_sfc) [ (ψ(θ_sfc) - ψ_1) / Δz + 1 ]

        At θ_sfc = θ_residual: K → 0, ψ → -∞, gnd_q → 0, so
        R ≈ -rho_a atm_q u* f_h ≤ 0.
        At θ_sfc = θ_porosity: K → K_sat, ψ → 0, gnd_q is saturated, so
        R > 0. The monotonic sign change guarantees a bracketed root,
        which Brent's method finds robustly - including the dry-soil
        regime where the previous ψ-inversion approach was ill-conditioned.
        """
        RHO_W = c.water.DENSITY
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        EVT = c.thermodynamic.EPSILON_VIRTUAL_TEMPERATURE

        z_s = self.input.surface.z_s
        z_t = self.input.surface.z_t
        dz = self.input.grid.z[0] - self.input.grid.z[1]

        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        sfc_T = self.sfc_state.temperature
        L = self.sfc_state.turbulence.obukhov_length
        ust = self.sfc_state.turbulence.friction_velocity
        Tv = self.atm_state.temperature * (1.0 + EVT * atm_q)
        rho_a = atm_p / (RD * Tv)
        fh = self.sfc.fh(z_s, z_t, L)

        residual_q = float(self.soil.properties.residual[0])
        porosity = float(self.soil.properties.porosity[0])

        # Subsurface properties (fixed during SMB solve)
        psi1 = self.soil.water_potential(self.soil_state.moisture)[1]
        K1 = self.soil.conductivity_moisture(self.soil_state.moisture)[1]

        canopy = getattr(self, 'canopy', None)
        if canopy is not None:
            f_veg = canopy.veg_fraction
        else:
            f_veg = np.zeros_like(sfc_T)

        def smb_residual(theta: np.ndarray) -> np.ndarray:
            theta_c = np.clip(theta, residual_q, porosity)
            psi0 = self.soil.water_potential(theta_c, level=0)
            K0 = self.soil.conductivity_moisture(theta_c, level=0)
            # Geometric mean: K_h spans many orders of magnitude with
            # moisture, so arithmetic/harmonic means are dominated by the
            # wetter node. Geometric is the standard LSM pragmatic choice.
            K_mid = np.maximum(np.sqrt(K0 * K1), 1e-14)
            gnd_q = self.soil.surface_specific_humidity(sfc_T, theta_c, atm_p)
            # SMB is a bare-soil top-layer balance; transpiration is
            # removed from the root-zone moisture budget (diffusion RHS),
            # not the surface flux residual.
            E_soil = (1.0 - f_veg) * rho_a * (gnd_q - atm_q) * ust * fh
            return E_soil + RHO_W * K_mid * ((psi0 - psi1) / dz + 1.0)

        # Shrink the bracket slightly off the physical bounds to keep
        # ψ(θ) and K(θ) finite and well-defined at the endpoints.
        eps = 1e-6
        span = porosity - residual_q
        a = np.full_like(sfc_T, residual_q + eps * span)
        b = np.full_like(sfc_T, porosity - eps * span)

        iter_max = self.input.numerics.iterations.smb_flux
        tol = self.input.numerics.tolerances.smb_flux

        theta_sfc, converged = solvers.root_brent_vec(
            smb_residual, a, b, iter_max=iter_max, tol=tol
        )

        if not converged.all():
            self.logger.warning(
                'SMB root-finding did not converge for %d of %d columns.',
                int(np.sum(~converged)), self.ncol,
            )

        self.sfc_state.moisture = np.clip(theta_sfc, residual_q, porosity)

    def _compute_fluxes(self, sfc_T: np.ndarray, sfc_q: np.ndarray) -> None:
        """Computes surface fluxes using Monin-Obukhov Similarity Theory.

        This is an iterative process to find the friction velocity (ustar)
        and Obukhov length (L) that are consistent with the calculated
        sensible and latent heat fluxes. All columns are processed in parallel.

        Args:
            sfc_T: Surface temperature array [K] (ncol,).
            sfc_q: Surface moisture array [m^3/m^3] (ncol,).
        """
        ITER_MAX = self.input.numerics.iterations.sfc_flux

        sfc_T_vec = self._as_column_vector(sfc_T, "sfc_T")
        sfc_q_vec = self._as_column_vector(sfc_q, "sfc_q")
        L_init = np.array(self.sfc_state.turbulence.obukhov_length, copy=True)

        (ustar, flux_wT, flux_wq, ground_heat, soil_top_T,
         L, sensible, latent, converged) = self._solve_most(
            sfc_T_vec, sfc_q_vec, L_init, max_iter=ITER_MAX
        )

        if not converged.all():
            self.logger.warning(
                'Obukhov length did not converge for %d columns.',
                int(np.sum(~converged)),
            )

        self.sfc_state.turbulence.obukhov_length[:] = L
        self.sfc_state.turbulence.friction_velocity[:] = ustar
        self.sfc_state.fluxes.ground_heat[:] = ground_heat
        self.sfc_state.fluxes.kinematic_heat[:] = flux_wT
        self.sfc_state.fluxes.kinematic_moisture[:] = flux_wq
        self.sfc_state.fluxes.sensible_heat[:] = sensible
        self.sfc_state.fluxes.latent_heat[:] = latent
        self.sfc_state.soil_top_temperature = np.asarray(soil_top_T)

    def _solve_surface_coupling(self) -> None:
        """Iteratively solves the coupled SEB and SMB.

        Performs a Picard iteration, alternating between solving the
        Surface Energy Budget (SEB) for temperature and the Surface
        Moisture Budget (SMB) for moisture. Each inner solve is a robust
        bracketed Brent root-find, so no under-relaxation is needed to
        damp oscillations — convergence typically occurs in a few outer
        iterations.
        """
        max_outer_iter = self.input.numerics.iterations.coupling
        tol_temp = self.input.numerics.tolerances.coupling_temp
        tol_mois = self.input.numerics.tolerances.coupling_mois

        for i in range(max_outer_iter):
            prev_T = np.array(self.sfc_state.temperature, copy=True)
            prev_q = np.array(self.sfc_state.moisture, copy=True)

            # Refresh r_s and θ_root once per outer iteration so SEB/SMB
            # inner solves see a fixed stomatal resistance.
            self._refresh_canopy_diagnostics()

            self._solve_seb()
            self._solve_smb()

            diff_T = np.abs(self.sfc_state.temperature - prev_T)
            diff_q = np.abs(self.sfc_state.moisture - prev_q)

            if np.all(diff_T < tol_temp) and np.all(diff_q < tol_mois):
                self.logger.debug(
                    'Surface coupling converged in %d iterations.', i + 1)
                # Refresh once more so canopy diagnostics and the
                # transpiration partition reflect the final converged
                # surface state, not the state from the start of the
                # last Picard iteration.
                self._refresh_canopy_diagnostics()
                self._compute_fluxes(
                    self.sfc_state.temperature, self.sfc_state.moisture
                )
                self._finalize_canopy_partition()
                return

        self.logger.warning(
            'Surface coupling did not converge after %d iterations. '
            'dT: %.4f, dq: %.4e', max_outer_iter,
            float(np.max(diff_T)), float(np.max(diff_q)))
        self._refresh_canopy_diagnostics()
        self._compute_fluxes(self.sfc_state.temperature, self.sfc_state.moisture)
        self._finalize_canopy_partition()

    def _solve_diffusion(self,state_field: np.ndarray,get_diffusivity: Callable,
                         get_conductivity: Optional[Callable],
                         sfc_boundary: float,field_name: str = 'field',
                         source_term: Optional[np.ndarray] = None,
                         avg_diffusivity: str = 'arithmetic') -> None:
        """Solves a generic 1D diffusion equation using a theta scheme.

        This helper now serves the soil heat solve. The moisture equation
        uses a dedicated mixed-form Richards implementation.

        Args:
            state_field: Reference to the field to update (temperature or
                moisture) [NDArray[np.float64]].
            get_diffusivity: Callable that computes diffusivity profile from
                soil moisture [Callable[[NDArray[np.float64]], NDArray[np.float64]]].
            get_conductivity: Optional callable that computes an auxiliary
                conductivity-like profile used by the tridiagonal assembly.
                Retained for backwards compatibility within this helper
                [Optional[Callable[[NDArray[np.float64]], NDArray[np.float64]]]].
            sfc_boundary: Surface boundary value for Dirichlet BC [float].
            field_name: Name of the field for logging/documentation [str].
            source_term: Optional per-layer source (or sink, if negative)
                with units of the state field per second, shape
                (nz, ncol). Applied as ``dt · source`` to the RHS of the
                tridiagonal system for layers 1..nz-1.
            avg_diffusivity: Averaging rule for interface diffusivity. Use
                ``"arithmetic"`` for simple means or ``"geometric"`` for
                multiplicative averaging.

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
        theta_b = self.input.numerics.heat_diffusion_back_weight
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

        e = self.solver_state.diffusion_e
        f = self.solver_state.diffusion_f
        g = self.solver_state.diffusion_g
        r = self.solver_state.diffusion_r
        if e.shape != (nz - 1, ncol):
            e = self.solver_state.diffusion_e = np.zeros((nz - 1, ncol))
            f = self.solver_state.diffusion_f = np.zeros((nz - 1, ncol))
            g = self.solver_state.diffusion_g = np.zeros((nz - 1, ncol))
            r = self.solver_state.diffusion_r = np.zeros((nz - 1, ncol))
        else:
            e.fill(0.0)
            f.fill(0.0)
            g.fill(0.0)
            r.fill(0.0)

        # Compute diffusivity using soil moisture (always, for both heat
        # and moisture). Moisture diffusivity D_θ spans orders of
        # magnitude with θ, so arithmetic face-averaging is dominated by
        # the wetter node; geometric mean (Haverkamp & Vauclin, 1979) is
        # the standard choice there. Thermal diffusivity α varies far
        # less, so arithmetic is fine.
        D = get_diffusivity(self.soil_state.moisture)
        if avg_diffusivity == 'geometric':
            D_mid = np.sqrt(np.maximum(D[:-1] * D[1:], 0.0))
        else:
            D_mid = 0.5 * (D[:-1] + D[1:])

        # Compute conductivity terms if provided (moisture case)
        K_lin = None
        if get_conductivity is not None:
            K_lin = get_conductivity(field)

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

        # Apply optional per-layer source (e.g. root-uptake sink) for
        # layers 1..nz-1. source_term has shape (nz, ncol); r is
        # (nz-1, ncol), indexed so r[j] ↔ layer (j+1).
        if source_term is not None:
            src = np.asarray(source_term, dtype=float)
            if src.ndim == 1:
                src = src[:, None]
            if src.shape != (nz, ncol):
                raise ValueError(
                    f"source_term shape {src.shape} does not match "
                    f"(nz, ncol)=({nz}, {ncol})."
                )
            r += dt * src[1:]

        # Solve and update
        field[0] = sfc_boundary_vec
        field[1:] = solvers.tridiagonal(e, f, g, r)

        if state_field.ndim == 1:
            state_field[:] = field[:, 0]
        else:
            state_field[:] = field

    def _solve_diffusion_heat(self) -> None:
        """Solves the soil heat diffusion equation using a theta scheme.

        The Dirichlet BC is the soil-top temperature, not the radiative
        skin: when an in-canopy thermal resistance is configured the
        two differ by ``G * r_canopy_thermal``. With no canopy
        resistance the two are identical, recovering the original
        bare-skin behaviour.
        """
        self._solve_diffusion(
            state_field=self.soil_state.temperature,
            get_diffusivity=self.soil.diffusivity_thermal,
            get_conductivity=None,
            sfc_boundary=self.sfc_state.soil_top_temperature,
            field_name='temperature'
        )

    def _solve_mixed_moisture(
        self,
        source_term: Optional[np.ndarray] = None,
    ) -> None:
        """Solves soil moisture with a mixed-form Richards Picard iteration.

        The nonlinear solve iterates in pressure head while the prognostic
        state remains volumetric moisture content. Darcy fluxes are assembled
        at faces using a positive-downward depth coordinate and a lagged
        conductivity from the current Picard iterate.
        """
        nz = self.input.grid.nz
        dt = self.tstep
        if nz < 2:
            raise ValueError('Mixed-form moisture solve requires nz >= 2.')

        moisture = np.asarray(self.soil_state.moisture, dtype=float)
        squeeze = False
        if moisture.ndim == 1:
            moisture = moisture[:, None]
            squeeze = True
        elif moisture.ndim != 2 or moisture.shape[0] != nz:
            raise ValueError(
                f'moisture has unexpected shape {moisture.shape}.')

        ncol = moisture.shape[1]
        theta_sfc = self._as_column_vector(
            self.sfc_state.moisture, 'moisture boundary')
        if theta_sfc.size != ncol:
            raise ValueError(
                f'moisture boundary size {theta_sfc.size} does not match '
                f'ncol={ncol}.'
            )

        dx = abs(self.input.grid.z[1] - self.input.grid.z[0])
        if dx <= 0.0:
            raise ValueError(
                f'Expected non-zero soil spacing, got dz={dx}.'
            )
        dx2 = dx ** 2

        e = self.solver_state.diffusion_e
        f = self.solver_state.diffusion_f
        g = self.solver_state.diffusion_g
        r = self.solver_state.diffusion_r
        if e.shape != (nz - 1, ncol):
            e = self.solver_state.diffusion_e = np.zeros((nz - 1, ncol))
            f = self.solver_state.diffusion_f = np.zeros((nz - 1, ncol))
            g = self.solver_state.diffusion_g = np.zeros((nz - 1, ncol))
            r = self.solver_state.diffusion_r = np.zeros((nz - 1, ncol))
        else:
            e.fill(0.0)
            f.fill(0.0)
            g.fill(0.0)
            r.fill(0.0)

        residual = self.soil._expand_profile_property(
            self.soil.properties.residual, moisture)
        porosity = self.soil._expand_profile_property(
            self.soil.properties.porosity, moisture)
        span = np.maximum(porosity - residual, 1e-12)
        theta_for_head = np.clip(
            moisture, residual + 1e-6 * span, porosity
        )

        residual0 = float(self.soil.properties.residual[0])
        porosity0 = float(self.soil.properties.porosity[0])
        span0 = max(porosity0 - residual0, 1e-12)
        theta_sfc_head = np.clip(theta_sfc, residual0 + 1e-6 * span0, porosity0)
        psi_sfc = np.asarray(
            self.soil.water_potential(theta_sfc_head, level=0), dtype=float
        )

        psi_iter = np.asarray(
            self.soil.water_potential(theta_for_head), dtype=float
        )
        psi_iter[0] = psi_sfc
        theta_iter = np.array(theta_for_head, copy=True)
        theta_iter[0] = theta_sfc
        theta_old = np.array(moisture, copy=True)

        src = None
        if source_term is not None:
            src = np.asarray(source_term, dtype=float)
            if src.ndim == 1:
                src = src[:, None]
            if src.shape != (nz, ncol):
                raise ValueError(
                    f"source_term shape {src.shape} does not match "
                    f"(nz, ncol)=({nz}, {ncol})."
                )

        iter_max = self.input.numerics.iterations.moisture_picard
        tol = max(float(self.input.numerics.tolerances.moisture_picard), 1e-8)
        converged = np.zeros(ncol, dtype=bool)

        for i in range(iter_max):
            theta_iter[:] = np.asarray(
                self.soil.water_content(psi_iter), dtype=float
            )
            theta_iter[0] = theta_sfc

            capacity = np.asarray(
                self.soil.moisture_capacity(psi_iter), dtype=float
            )
            capacity = np.maximum(capacity, 0.0)
            K_node = np.asarray(
                self.soil.conductivity_moisture(theta_iter), dtype=float
            )
            K_face = np.sqrt(np.maximum(K_node[:-1] * K_node[1:], 0.0))

            q_up = K_face * (1.0 - (psi_iter[1:] - psi_iter[:-1]) / dx)
            q_dn = np.zeros_like(q_up)
            q_dn[:-1] = (
                K_face[1:]
                * (1.0 - (psi_iter[2:] - psi_iter[1:-1]) / dx)
            )
            # Zero gradient of pressure head at the lower boundary implies
            # free drainage: q = K.
            q_dn[-1] = K_node[-1]

            K_up = K_face
            K_dn = np.zeros_like(K_up)
            K_dn[:-1] = K_face[1:]

            e.fill(0.0)
            f[:] = capacity[1:] / dt + (K_up + K_dn) / dx2
            g.fill(0.0)
            e[1:] = -K_up[1:] / dx2
            g[:-1] = -K_dn[:-1] / dx2

            r[:] = (
                -(theta_iter[1:] - theta_old[1:]) / dt
                - (q_dn - q_up) / dx
            )
            if src is not None:
                r += src[1:]

            delta_psi = solvers.tridiagonal(e, f, g, r)
            psi_iter[1:] += delta_psi
            psi_iter[0] = psi_sfc

            theta_next = np.asarray(
                self.soil.water_content(psi_iter), dtype=float
            )
            theta_next[0] = theta_sfc
            delta_theta = np.max(
                np.abs(theta_next[1:] - theta_iter[1:]), axis=0
            )
            theta_iter[:] = theta_next
            converged = delta_theta < tol
            if converged.all():
                self.logger.debug(
                    'Mixed-form moisture solve converged in %d Picard '
                    'iterations.',
                    i + 1,
                )
                break

        if not converged.all():
            self.logger.warning(
                'Mixed-form moisture solve did not converge for %d of %d '
                'columns after %d iterations.',
                int(np.sum(~converged)),
                ncol,
                iter_max,
            )

        if squeeze:
            self.soil_state.moisture[:] = theta_iter[:, 0]
        else:
            self.soil_state.moisture[:] = theta_iter

    def _solve_diffusion_mois(self) -> None:
        """Solves soil moisture with a mixed-form Richards equation.

        When a canopy is active, the per-layer root extraction rate (a
        volumetric sink in [m^3/m^3/s]) is injected explicitly into the
        moisture tendency so transpiration removes water distributively
        from the root zone rather than from the bare-soil top-layer
        balance.
        """
        source = None
        if getattr(self, 'canopy', None) is not None:
            source = -self.canopy_state.root_uptake
        self._solve_mixed_moisture(source_term=source)
        self._enforce_soil_moisture_bounds()

    def _enforce_soil_moisture_bounds(self) -> None:
        """Clips tiny moisture overshoots and raises on material violations."""
        moisture = np.asarray(self.soil_state.moisture, dtype=float)
        residual = self.soil._expand_profile_property(
            self.soil.properties.residual, moisture
        )
        porosity = self.soil._expand_profile_property(
            self.soil.properties.porosity, moisture
        )
        tol = max(float(self.input.numerics.tolerances.moisture_bounds), 1e-8)

        below_hard = moisture < (residual - tol)
        above_hard = moisture > (porosity + tol)
        hard_mask = below_hard | above_hard
        if np.any(hard_mask):
            min_delta = float(np.min(moisture - residual))
            max_delta = float(np.max(moisture - porosity))
            bad_idx = np.argwhere(hard_mask)
            examples = ", ".join(
                f"{tuple(int(i) for i in idx)}={float(moisture[tuple(idx)]):.6f}"
                for idx in bad_idx[:5]
            )
            raise SolverError(
                'Soil moisture left physical bounds after the mixed moisture '
                'solve: '
                f'{int(np.sum(hard_mask))} cells outside [residual, porosity] '
                f'by more than tol={tol:.1e}. '
                f'Min(theta-residual)={min_delta:.3e}, '
                f'Max(theta-porosity)={max_delta:.3e}. '
                f'Examples: {examples}'
            )

        clip_mask = (moisture < residual) | (moisture > porosity)
        if np.any(clip_mask):
            self.logger.warning(
                'Clipping %d soil moisture values to [residual, porosity] '
                'after the mixed moisture solve (tol=%.1e).',
                int(np.sum(clip_mask)),
                tol,
            )
            np.clip(moisture, residual, porosity, out=moisture)
