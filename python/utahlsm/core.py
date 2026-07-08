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
from typing import Callable, Literal, Optional, cast

import numpy as np

from ._types import FloatOrArray
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

    # Class-level annotation for attributes set inside helper methods
    # called from __init__; lets mypy resolve the type when this attribute
    # is read in update()/save() before mypy has traced the helper.
    output_fields: dict[str, np.ndarray | list[str]]

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

        # Compute the downwelling radiation components. Outgoing
        # components are evaluated as a function of trial T_s inside the
        # SEB iterator (see _compute_seb_vec) so the longwave emission
        # feedback is captured during root-finding. After SEB converges,
        # _solve_seb writes the final sw_out, lw_out, and radiation_net
        # back to atm_state for downstream consumers and output.
        total_seconds = self.input.time.utc_start + runtime
        days_passed = int(total_seconds // 86400)
        current_utc = total_seconds % 86400
        days_per_year = (366 if self._is_leap_year(
            self.input.time.utc_year) else 365)
        julian_day = ((self.input.time.julian_day
            + days_passed - 1) % days_per_year) + 1
        sw_in, lw_in = self.rad.compute_incoming(
            julian_day, current_utc, self.atm_state, self.sfc_state
        )
        self.atm_state.sw_in = self._as_column_vector(sw_in, "sw_in")
        self.atm_state.lw_in = self._as_column_vector(lw_in, "lw_in")
        # Provisional outgoing diagnostics from current T_s; refreshed
        # after the SEB solve.
        self._refresh_radiation_diagnostics(self.sfc_state.temperature)

    def run(self) -> None:
        """Runs the core model physics for a single time step.

        This includes solving the surface energy and moisture budgets and
        updating the soil profiles via diffusion solvers.
        """
        self.logger.info('Solving soil state')

        # Set initial guesses for new surface temp and moisture
        self.sfc_state.temperature = np.array(
            self.soil_state.temperature[0], copy=True)
        self.sfc_state.moisture = np.array(
            self.soil_state.moisture[0], copy=True)

        self._prepare_precipitation_for_timestep()

        # Solve surface energy and moisture budgets
        self._solve_surface_coupling()

        # Post-coupling precipitation bookkeeping: advance the Green-Ampt
        # event state and build the macropore deposition profile before
        # the heat and moisture solves consume them.
        self._update_infiltration_history()
        self._distribute_bypass()

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
        self._refresh_reassigned_outputs()
        self.output.save(self.output_fields, step_count, runtime)

    def _refresh_reassigned_outputs(self) -> None:
        """Re-point output_fields at fields that are reassigned each step.

        Most outputs are updated in place, but the radiative skin temperature
        and outgoing longwave are reassigned (``self.x = np.array(...)``) by
        the SEB solve, so a reference captured once would go stale.
        """
        for name, value in (
            ('T_skin', self.sfc_state.temperature),
            ('lw_out', self.atm_state.lw_out),
        ):
            if name in self.output_fields:
                self.output_fields[name] = np.asarray(value)

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
        self.sfc_state.seb_residual = np.zeros(self.ncol)

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
        self.atm_state.seb_storage = np.zeros(self.ncol)
        self.atm_state.precipitation = np.zeros(self.ncol)
        self.sfc_state.fluxes.runoff = np.zeros(self.ncol)
        self._soil_precipitation = np.zeros(self.ncol)
        self._precipitation_prepared = False
        # Infiltrated water flux [kg/m^2/s] (post-interception, post-bypass,
        # post-runoff), set by the SMB and consumed by the rain
        # heat-advection source.
        self._infiltration_flux = np.zeros(self.ncol)
        # Macropore bypass: throughfall captured by cracks [kg/m^2/s] and
        # its per-layer deposition profile [kg/m^2/s], rebuilt each step.
        self._bypass_flux = np.zeros(self.ncol)
        self._bypass_deposit: Optional[np.ndarray] = None
        # Green-Ampt event state: cumulative matrix-infiltrated depth [m]
        # since the start of the current rain event (decays between events).
        self._ga_cum_infil = np.zeros(self.ncol)

        self.solver_state: SolverState = SolverState()
        self.solver_state.conductivity_thermal_mid = np.zeros(self.ncol)
        nz_diff = self.input.grid.nz - 1
        self.solver_state.diffusion_e = np.zeros((nz_diff, self.ncol))
        self.solver_state.diffusion_f = np.zeros((nz_diff, self.ncol))
        self.solver_state.diffusion_g = np.zeros((nz_diff, self.ncol))
        self.solver_state.diffusion_r = np.zeros((nz_diff, self.ncol))

        # Canopy defaults to bare-soil until _setup_physics instantiates
        # a concrete model from the namelist. Stubbing it here keeps
        # partial-initialization test paths that skip _setup_physics safe.
        self.canopy: Optional[Canopy] = None
        nz = self.input.grid.nz
        self.canopy_state: CanopyState = CanopyState(
            resistance=np.full(self.ncol, np.inf),
            theta_root=np.zeros(self.ncol),
            transpiration=np.zeros(self.ncol),
            wet_evaporation=np.zeros(self.ncol),
            evap_soil=np.zeros(self.ncol),
            water_storage=np.zeros(self.ncol),
            water_capacity=np.zeros(self.ncol),
            latent_veg=np.zeros(self.ncol),
            latent_wet=np.zeros(self.ncol),
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
                return np.broadcast_to(data, (nz, ncol)) if ncol > 1 else data
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

    def _ensure_writable_soil_field(
        self, field_name: Literal["temperature", "moisture"]
    ) -> np.ndarray:
        """Returns a writable soil state field, copying broadcast views lazily."""
        if field_name == "temperature":
            field = self.soil_state.temperature
            output_name = "soil_T"
        else:
            field = self.soil_state.moisture
            output_name = "soil_q"
        if not field.flags.writeable:
            field = np.array(field, copy=True)
            if field_name == "temperature":
                self.soil_state.temperature = field
            else:
                self.soil_state.moisture = field
            if hasattr(self, "output_fields"):
                self.output_fields[output_name] = field
        return field

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

    def _load_atm_state(self, atm_state: AtmosphericState) -> None:
        """Loads atmospheric state data into column vectors.

        Forcing-driven runs supply the downwelling radiation components
        (``sw_in``, ``lw_in``); the upwelling components and net
        radiation are derived as functions of T_s by the radiation
        model and refreshed by the SEB solver, so they are not read
        from forcing here.
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
        self.atm_state.lw_in = self._as_column_vector(
            atm_state.lw_in, "lw_in")
        seb_storage = self._as_column_vector(
            atm_state.seb_storage, "seb_storage")
        if isinstance(self.atm_state.seb_storage, np.ndarray):
            self.atm_state.seb_storage[:] = seb_storage
        else:
            self.atm_state.seb_storage = seb_storage
        precipitation = self._as_column_vector(
            atm_state.precipitation, "precipitation")
        if isinstance(self.atm_state.precipitation, np.ndarray):
            self.atm_state.precipitation[:] = precipitation
        else:
            self.atm_state.precipitation = precipitation
        self._precipitation_prepared = False
        RD = c.thermodynamic.GAS_CONSTANT_DRY
        EVT = c.thermodynamic.EPSILON_VIRTUAL_TEMPERATURE
        Tv = self.atm_state.temperature * (
            1.0 + EVT * self.atm_state.specific_humidity)
        self.sfc_state.air_density = self.atm_state.pressure / (RD * Tv)

    def _prepare_precipitation_for_timestep(self) -> None:
        """Routes liquid precipitation through canopy storage once per step.

        The raw forcing remains on ``atm_state.precipitation`` for output.
        The cached ``_soil_precipitation`` is throughfall plus bare-ground
        precipitation, and is what the SMB consumes during Picard coupling.
        """
        P_total = np.asarray(self.atm_state.precipitation, dtype=float)
        soil_precip = np.array(P_total, copy=True)

        canopy = getattr(self, 'canopy', None)
        dt = float(getattr(self, 'tstep', 0.0))
        if canopy is not None and dt > 0.0:
            storage = np.asarray(self.canopy_state.water_storage, dtype=float)
            capacity = np.asarray(self.canopy_state.water_capacity, dtype=float)
            room = np.maximum(capacity - storage, 0.0)
            interception_potential = np.maximum(P_total, 0.0) * canopy.veg_fraction
            interception = np.minimum(interception_potential, room / dt)
            self.canopy_state.water_storage[:] = np.clip(
                storage + interception * dt,
                0.0,
                capacity,
            )
            soil_precip = P_total - interception

        cached_precip = getattr(self, '_soil_precipitation', None)
        if cached_precip is None or np.asarray(cached_precip).shape != soil_precip.shape:
            self._soil_precipitation = np.array(soil_precip, copy=True)
        else:
            self._soil_precipitation[:] = soil_precip
        self._precipitation_prepared = True

    def _refresh_radiation_diagnostics(self, sfc_T: np.ndarray) -> None:
        """Updates outgoing radiation and net radiation on atm_state.

        Outgoing components are responses to T_s, so they are derived
        from the radiation model whenever T_s changes (after forcing
        load and after each SEB solve). Stored on ``atm_state`` for
        downstream consumers (output, diagnostics).

        Args:
            sfc_T: Surface temperature [K], shape (ncol,).
        """
        sw_in = np.asarray(self.atm_state.sw_in)
        lw_in = np.asarray(self.atm_state.lw_in)
        sw_out, lw_out = self.rad.compute_outgoing(
            np.asarray(sfc_T), sw_in, lw_in,
            self.atm_state, self.sfc_state,
        )
        self.atm_state.sw_out = np.asarray(sw_out)
        self.atm_state.lw_out = np.asarray(lw_out)
        self.atm_state.radiation_net = (
            sw_in - self.atm_state.sw_out
            + lw_in - self.atm_state.lw_out
        )

    def _setup_physics(self) -> None:
        """Initializes the physics modules based on user configuration."""
        self.logger.info('Initializing physics modules')
        try:
            self.rad: Radiation = get_radiation_model(
                self.input.radiation,
                self.input.surface,
            )
            self.soil: Soil = get_soil_model(
                self.input.soil.model,
                self.input.soil_properties,
                self.input.soil_type_names,
                self.input.soil_properties_name,
                self.input.soil.thermal_conductivity_model,
            )
            self.sfc: Surface = get_surface_model(self.input.surface)
            self.canopy = get_canopy_model(
                self.input.canopy, self.input.grid.z, self.ncol
            )
            if self.canopy is not None:
                self.canopy_state.water_capacity[:] = (
                    self.canopy.veg_fraction
                    * self.canopy.lai
                    * self.canopy.water_capacity_lai
                )
        except NamelistError as e:
            self.logger.error('Failed to initialize physics modules: %s.', e)
            raise

    def _setup_output(self) -> None:
        """Sets up the output file dimensions and fields."""
        nx = self.input.grid.nx
        ny = self.input.grid.ny
        self.output_dims: dict[str, int] = {
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

        if forcing0 is not None:
            # Drive the full SEB+SMB Picard coupling at forcing[0] so that
            # T_s, θ_sfc, and L are mutually consistent at t=0. Without
            # this, step 1 sees SMB shifting θ_sfc away from the initial
            # soil-top moisture, which redistributes the energy balance
            # (LE↔SHF) and produces a large, non-physical jump in u*/L.
            assert self.input.forcing is not None
            self._load_atm_state(forcing0)
            self.tstep = float(self.input.forcing.tstep)

            self._solve_surface_coupling()

            # Top soil cell carries soil-top temperature; equals the
            # radiative skin only when r_canopy_thermal = 0.
            self._ensure_writable_soil_field("temperature")
            self.soil_state.temperature[0] = self.sfc_state.soil_top_temperature

            if self.output.outfile is not None:
                self.output.outfile.setncattr(
                    "initial_state",
                    "SEB+SMB coupling using forcing[0]",
                )

        self.full_output_fields: dict[str, np.ndarray | list[str]] = {
            'ust': self.sfc_state.turbulence.friction_velocity,
            'obl': self.sfc_state.turbulence.obukhov_length,
            'shf': self.sfc_state.fluxes.sensible_heat,
            'lhf': self.sfc_state.fluxes.latent_heat,
            'ghf': self.sfc_state.fluxes.ground_heat,
            # Radiative skin temperature and emitted longwave. Unlike the soil
            # top node (soil_T[0] = soil-top temperature, below the in-canopy
            # resistance), these are the surface the radiometer sees. Both are
            # reassigned each step, so save() re-points the dict before writing.
            'T_skin': np.asarray(self.sfc_state.temperature),
            'lw_out': np.asarray(self.atm_state.lw_out),
            'seb_res': self.sfc_state.seb_residual,
            'seb_storage': np.asarray(self.atm_state.seb_storage),
            'precip': np.asarray(self.atm_state.precipitation),
            'runoff': self.sfc_state.fluxes.runoff,
            'bypass': getattr(
                self, '_bypass_flux',
                np.zeros_like(np.asarray(self.sfc_state.fluxes.runoff)),
            ),
            'soil_z': self.input.grid.z,
            'soil_type': self.input.soil_type_names,
            'soil_T': self.soil_state.temperature,
            'soil_q': self.soil_state.moisture,
        }
        if getattr(self, 'canopy', None) is not None:
            self.full_output_fields.update({
                'r_s': self.canopy_state.resistance,
                'theta_root': self.canopy_state.theta_root,
                'lhf_soil': self.canopy_state.latent_soil,
                'lhf_veg': self.canopy_state.latent_veg,
                'lhf_wet': self.canopy_state.latent_wet,
                'canopy_water': self.canopy_state.water_storage,
            })
        self.output_fields = self._select_output_fields(self.full_output_fields)
        self.output.set_fields(self.output_fields)
        self.output.save(self.output_fields, 0, 0, initial=True)

    def _select_output_fields(
            self,
            available_fields: dict[str, np.ndarray | list[str]],
    ) -> dict[str, np.ndarray | list[str]]:
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

    def _partition_flux_wq_components(
        self,
        sfc_T: FloatOrArray,
        gnd_q: FloatOrArray,
        atm_q: FloatOrArray,
        atm_p: FloatOrArray,
        ustar: np.ndarray,
        fh: FloatOrArray,
        cols: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Kinematic soil, dry-canopy, and wet-canopy moisture fluxes.

        Bare-soil mode: single-source form
        ``flux_wq = (gnd_q - atm_q) · u* · f_h``.

        Canopy-active mode uses three terms: bare-soil exchange, dry
        stomatal transpiration, and non-stomatal wet-canopy exchange
        against a prognostic water store. The wet term keeps the signed
        vapor gradient so dewfall can enter the SEB as negative LE, but
        evaporation and condensation are capped by the water available
        or storage capacity over the current timestep.
        """
        canopy = getattr(self, 'canopy', None)
        if canopy is None:
            soil = cast(np.ndarray, (gnd_q - atm_q) * ustar * fh)
            return soil, np.zeros_like(soil), np.zeros_like(soil)

        if cols is not None:
            f_veg = canopy.veg_fraction[cols]
            r_s = self.canopy_state.resistance[cols]
            storage = self.canopy_state.water_storage[cols]
            capacity = self.canopy_state.water_capacity[cols]
            wet_cooling_max = canopy.wet_cooling_max[cols]
            lw_in = np.asarray(self.atm_state.lw_in)[cols]
            sw_in = np.asarray(self.atm_state.sw_in)[cols]
        else:
            f_veg = canopy.veg_fraction
            r_s = self.canopy_state.resistance
            storage = self.canopy_state.water_storage
            capacity = self.canopy_state.water_capacity
            wet_cooling_max = canopy.wet_cooling_max
            lw_in = np.asarray(self.atm_state.lw_in)
            sw_in = np.asarray(self.atm_state.sw_in)

        u_fh = ustar * fh
        q_sat = thermo.saturation_specific_humidity(sfc_T, atm_p)
        e_soil = (1.0 - f_veg) * (gnd_q - atm_q) * u_fh

        wet_fraction = np.divide(
            storage,
            capacity,
            out=np.zeros_like(storage),
            where=capacity > 0.0,
        )
        wet_fraction = np.clip(wet_fraction, 0.0, 1.0)

        # 1/(ra + rs) normalized dry transpiration. It is one-way: roots
        # cannot take up negative water, and dewfall belongs to the wet
        # storage path below. Wet leaves also reduce the dry transpiring
        # fraction.
        denom = 1.0 + u_fh * r_s
        t_veg = (
            f_veg
            * (1.0 - wet_fraction)
            * np.maximum(q_sat - atm_q, 0.0)
            * u_fh
            / denom
        )

        SB = c.radiation.STEFAN_BOLTZMANN
        CP = c.thermodynamic.SPECIFIC_HEAT
        emissivity = getattr(self.input.surface, 'emissivity', 1.0)
        lw_cooling_flux = np.maximum(emissivity * SB * sfc_T**4 - lw_in, 0.0)
        night_factor = np.clip(1.0 - sw_in / 50.0, 0.0, 1.0)

        rho = (np.asarray(self.sfc_state.air_density)[cols]
               if cols is not None else np.asarray(self.sfc_state.air_density))
        cooling = np.divide(
            lw_cooling_flux,
            rho * CP * np.maximum(u_fh, 1e-4),
            out=np.zeros_like(sfc_T),
            where=rho > 0.0,
        )
        cooling = np.minimum(cooling * night_factor, wet_cooling_max)
        q_sat_wet = thermo.saturation_specific_humidity(sfc_T - cooling, atm_p)

        wet = f_veg * (q_sat_wet - atm_q) * u_fh
        wet = np.where(wet > 0.0, wet * wet_fraction, wet)

        dt = max(float(getattr(self, 'tstep', 0.0)), 1.0)
        evap_limit = np.divide(
            storage,
            rho * dt,
            out=np.zeros_like(storage),
            where=rho > 0.0,
        )
        dew_limit = np.divide(
            np.maximum(capacity - storage, 0.0),
            rho * dt,
            out=np.zeros_like(storage),
            where=rho > 0.0,
        )
        wet = np.minimum(wet, evap_limit)
        wet = np.maximum(wet, -dew_limit)

        return (
            cast(np.ndarray, e_soil),
            cast(np.ndarray, t_veg),
            cast(np.ndarray, wet),
        )

    def _partition_flux_wq(
        self,
        sfc_T: np.ndarray,
        gnd_q: FloatOrArray,
        atm_q: np.ndarray,
        atm_p: np.ndarray,
        ustar: np.ndarray,
        fh: FloatOrArray,
        cols: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Total kinematic moisture flux with optional canopy partition."""
        e_soil, t_veg, wet = self._partition_flux_wq_components(
            sfc_T, gnd_q, atm_q, atm_p, ustar, fh, cols=cols
        )
        return cast(np.ndarray, e_soil + t_veg + wet)

    def _finalize_canopy_partition(self) -> None:
        """Splits the converged total LH into soil and canopy components.

        Also records per-layer root uptake rate [m^3/m^3/s] used as a
        sink term in the soil moisture RHS. Layer 0's root fraction
        is folded into layer 1 so the surface BC is unaffected by
        transpiration (the SMB already balances the top layer).
        """
        LV = c.thermodynamic.LATENT_HEAT_VAPORIZATION
        RHO_W = c.water.DENSITY

        canopy = getattr(self, 'canopy', None)
        if canopy is None:
            if hasattr(self, 'canopy_state'):
                self.canopy_state.evap_soil[:] = 0.0
                self.canopy_state.transpiration[:] = 0.0
                self.canopy_state.wet_evaporation[:] = 0.0
                self.canopy_state.water_storage[:] = 0.0
                self.canopy_state.water_capacity[:] = 0.0
                self.canopy_state.latent_soil[:] = (
                    self.sfc_state.fluxes.latent_heat)
                self.canopy_state.latent_veg[:] = 0.0
                self.canopy_state.latent_wet[:] = 0.0
                self.canopy_state.root_uptake[:] = 0.0
            return

        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        sfc_T = self.sfc_state.temperature
        sfc_q = self.sfc_state.moisture
        rho_a = self.sfc_state.air_density

        ust = self.sfc_state.turbulence.friction_velocity
        fh = self.sfc.fh(
            self.input.surface.z_s, self.input.surface.z_t,
            self.sfc_state.turbulence.obukhov_length,
        )
        gnd_q = self.soil.surface_specific_humidity(sfc_T, sfc_q, atm_p)

        E_soil_kin, T_veg_kin, W_veg_kin = self._partition_flux_wq_components(
            sfc_T, gnd_q, atm_q, atm_p, ust, fh
        )

        E_soil_mass = rho_a * E_soil_kin          # [kg/m^2/s]
        T_veg_mass = rho_a * T_veg_kin            # [kg/m^2/s]
        W_veg_mass = rho_a * W_veg_kin            # [kg/m^2/s]

        self.canopy_state.evap_soil[:] = E_soil_mass
        self.canopy_state.transpiration[:] = T_veg_mass
        self.canopy_state.wet_evaporation[:] = W_veg_mass
        self.canopy_state.latent_soil[:] = LV * E_soil_mass
        self.canopy_state.latent_veg[:] = LV * T_veg_mass
        self.canopy_state.latent_wet[:] = LV * W_veg_mass

        dt = max(float(getattr(self, 'tstep', 0.0)), 0.0)
        if dt > 0.0:
            self.canopy_state.water_storage[:] = np.clip(
                self.canopy_state.water_storage - W_veg_mass * dt,
                0.0,
                self.canopy_state.water_capacity,
            )

        # Per-layer uptake rate [m^3/m^3/s] folded into layers 1..nz-1.
        nz = self.input.grid.nz
        root_frac = canopy.root_fraction.copy()  # (nz, ncol)
        if nz > 1:
            root_frac[1] += root_frac[0]
            root_frac[0] = 0.0
        # Extraction is availability-weighted (compensatory uptake): each
        # layer's share is root density times the same wilting-to-field-
        # capacity ramp that feeds the f4 stress, renormalized per column.
        # Without the ramp, transpiration sustained by moist deep layers
        # is withdrawn from layers already at or below wilting, drying
        # them indefinitely. When the whole root zone is stressed the
        # plain root profile is used so the column sink still matches the
        # (small, rs_max-limited) surface transpiration flux.
        theta = np.asarray(self.soil_state.moisture, dtype=float)
        if theta.ndim == 1:
            theta = theta[:, None]
        wilt = self.soil.theta_wilt[:, None]
        fc = self.soil.theta_fc[:, None]
        availability = np.clip(
            (theta - wilt) / np.maximum(fc - wilt, 1e-6), 0.0, 1.0
        )
        weights = root_frac * availability
        w_sum = np.sum(weights, axis=0, keepdims=True)
        weights = np.where(
            w_sum > 1e-12,
            np.divide(weights, w_sum, out=np.zeros_like(weights),
                      where=w_sum > 1e-12),
            root_frac,
        )
        dz = self.input.grid.z[0] - self.input.grid.z[1]
        # T_veg_mass broadcast over layers times extraction weight / (rho_w dz).
        self.canopy_state.root_uptake[:] = (
            weights * T_veg_mass[None, :] / (RHO_W * dz)
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

        # Create vectorized brackets around current temperatures. Keep trial
        # temperatures physically valid for thermodynamic helper functions.
        current_T: np.ndarray = np.array(self.sfc_state.temperature, copy=True)
        bracket_width = 1.0
        min_sfc_temp = 100.0
        temp_a: np.ndarray = np.maximum(current_T - bracket_width, min_sfc_temp)
        temp_b: np.ndarray = current_T + bracket_width

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

            temp_a = np.where(
                expand_left, np.maximum(temp_a - step, min_sfc_temp), temp_a
            )
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

        # Single-pass MOST update: refreshes ustar and L so that SMB sees
        # stability consistent with the just-found T_s, without paying for
        # a full iterative MOST solve. The full solve runs once per timestep
        # in _solve_surface_coupling after all Picard iterations converge.
        sfc_q_vec = self._as_column_vector(self.sfc_state.moisture, "sfc_q")
        L_now = np.array(self.sfc_state.turbulence.obukhov_length, copy=True)
        (ustar, _, _, ground_heat, soil_top_T,
         obukhov_l, sensible, latent, _) = self._solve_most(
            self.sfc_state.temperature, sfc_q_vec, L_now, max_iter=1
        )
        self.sfc_state.turbulence.obukhov_length[:] = obukhov_l
        self.sfc_state.turbulence.friction_velocity[:] = ustar
        self.sfc_state.soil_top_temperature = np.asarray(soil_top_T)

        # Refresh outgoing radiation diagnostics from the converged T_s
        # so atm_state.sw_out / lw_out / radiation_net reflect the SEB
        # solution that downstream consumers and output will see.
        self._refresh_radiation_diagnostics(self.sfc_state.temperature)

        # Diagnostic SEB residual. Should be ~ tol_root in magnitude;
        # sustained drift indicates a tolerance or coupling problem.
        self.sfc_state.seb_residual[:] = (
            np.asarray(self.atm_state.radiation_net)
            - np.asarray(self.atm_state.seb_storage)
            - ground_heat
            - sensible
            - latent
        )

    def _solve_most(
        self,
        sfc_T: np.ndarray,
        sfc_q: np.ndarray,
        L_init: np.ndarray,
        max_iter: int,
        cols: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
               np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
        VK = c.physical.VON_KARMAN
        G = c.physical.GRAVITY
        EVT = c.thermodynamic.EPSILON_VIRTUAL_TEMPERATURE
        TOL = self.input.numerics.tolerances.sfc_flux

        canopy = getattr(self, 'canopy', None)
        if cols is not None:
            atm_T: np.ndarray = np.asarray(self.atm_state.temperature)[cols]
            atm_p: np.ndarray = np.asarray(self.atm_state.pressure)[cols]
            atm_q: np.ndarray = np.asarray(self.atm_state.specific_humidity)[cols]
            atm_ws: np.ndarray = np.asarray(self.atm_state.wind_speed)[cols]
            K_mid: np.ndarray = np.asarray(self.solver_state.conductivity_thermal_mid)[cols]
            soil_T1 = self.soil_state.temperature[1, cols]
            if canopy is not None:
                f_veg = canopy.veg_fraction[cols]
                r_g_aero = canopy.r_ground[cols]
            else:
                f_veg = np.zeros_like(soil_T1)
                r_g_aero = np.zeros_like(soil_T1)
        else:
            atm_T = np.asarray(self.atm_state.temperature)
            atm_p = np.asarray(self.atm_state.pressure)
            atm_q = np.asarray(self.atm_state.specific_humidity)
            atm_ws = np.asarray(self.atm_state.wind_speed)
            K_mid = np.asarray(self.solver_state.conductivity_thermal_mid)
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
        rho = (np.asarray(self.sfc_state.air_density)[cols]
               if cols is not None else np.asarray(self.sfc_state.air_density))

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

        obukhov_l = np.array(L_init, copy=True)
        converged = np.zeros_like(obukhov_l, dtype=bool)
        ustar = np.zeros_like(obukhov_l)
        flux_wT = np.zeros_like(obukhov_l)
        flux_wq = np.zeros_like(obukhov_l)

        for _ in range(max_iter):
            fm = self.sfc.fm(z_m, z_o, obukhov_l)
            fh = self.sfc.fh(z_s, z_t, obukhov_l)

            wind_eff = atm_ws
            if gustiness > 0.0:
                gust = np.hypot(atm_ws, gustiness)
                if gustiness_stable_only:
                    wind_eff = np.where(obukhov_l >= 0.0, gust, atm_ws)
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

            diff = np.abs(L_new - obukhov_l)
            converged |= diff <= TOL
            obukhov_l = np.where(converged, obukhov_l, L_new)

            if converged.all():
                break

        sensible = rho * CP * flux_wT
        latent = rho * LV * flux_wq

        return (ustar, flux_wT, flux_wq, ground_heat, soil_top_T,
                obukhov_l, sensible, latent, converged)

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
            sfc_q: np.ndarray = np.asarray(self.sfc_state.moisture)[cols]
            sw_in: np.ndarray = np.asarray(self.atm_state.sw_in)[cols]
            lw_in: np.ndarray = np.asarray(self.atm_state.lw_in)[cols]
            storage_all = np.asarray(self.atm_state.seb_storage)
            storage: np.ndarray = (
                storage_all if storage_all.ndim == 0 else storage_all[cols]
            )
        else:
            sfc_q = np.asarray(self.sfc_state.moisture)
            sw_in = np.asarray(self.atm_state.sw_in)
            lw_in = np.asarray(self.atm_state.lw_in)
            storage = np.asarray(self.atm_state.seb_storage)

        # Outgoing components respond to the trial T_s; this is what
        # supplies the ~4εσT^3 longwave restoring inside Brent.
        sw_out, lw_out = self.rad.compute_outgoing(
            sfc_T, sw_in, lw_in, self.atm_state, self.sfc_state
        )
        rad_net = sw_in - np.asarray(sw_out) + lw_in - np.asarray(lw_out)

        _, _, _, ground_heat, _, _, sensible, latent, _ = self._solve_most(
            sfc_T, sfc_q, initial_L, max_iter=1, cols=cols
        )

        return np.asarray(rad_net - storage - ground_heat - sensible - latent)

    @staticmethod
    def _moisture_face_conductivity(
        K_upper: np.ndarray,
        K_lower: np.ndarray,
        psi_upper: np.ndarray,
        psi_lower: np.ndarray,
        dz: float,
    ) -> np.ndarray:
        """Lagged hydraulic conductivity at a vertical soil interface.

        The high-conductivity wetting-front branch is selected from the
        Darcy driving term, not from volumetric moisture ordering. That
        matters at texture breaks where clay can hold a larger θ than
        silty loam at the same or lower matric potential.
        """
        K_geom = np.sqrt(np.maximum(K_upper * K_lower, 0.0))
        darcy_drive = 1.0 - (psi_lower - psi_upper) / dz
        return np.where(darcy_drive > 0.0, K_upper, K_geom)

    # Green-Ampt event-state constants: floor on cumulative infiltration
    # (caps the early-storm capacity enhancement) and the e-folding time
    # over which the event memory relaxes after rain stops, representing
    # post-storm redistribution of the wetting front.
    _GA_MIN_CUM_INFIL = 1.0e-4   # [m]
    _GA_REDIST_TAU = 6.0 * 3600.0  # [s]

    def _macropore_config(self) -> tuple[float, float, float, float]:
        """Returns (fraction, z_top, z_bottom, e_folding) for bypass flow.

        Reads the soil namelist section with inert defaults so partial
        test rigs (and namelists without macropore keys) run bare-matrix.
        """
        soil_cfg = getattr(self.input, 'soil', None)
        fraction = float(getattr(soil_cfg, 'macropore_fraction', 0.0))
        z_top = float(getattr(soil_cfg, 'macropore_z_top', 0.0))
        z_bottom = float(getattr(soil_cfg, 'macropore_z_bottom', 0.0))
        e_folding = float(getattr(soil_cfg, 'macropore_e_folding', 0.2))
        return min(max(fraction, 0.0), 1.0), z_top, z_bottom, e_folding

    def _wetting_front_suction(self, theta_i: np.ndarray) -> np.ndarray:
        """Green-Ampt effective wetting-front suction [m, >= 0].

        Uses the Neuman (1976) definition psi_f = (1/K_sat) * int K dpsi
        from the antecedent moisture theta_i to saturation, evaluated by
        midpoint quadrature on the surface layer's own retention and
        conductivity curves - so the suction is consistent with whichever
        soil model is configured rather than a tabulated approximation.

        Args:
            theta_i: Antecedent surface moisture [m^3/m^3], shape (ncol,).

        Returns:
            Wetting-front suction head per column [m], shape (ncol,).
        """
        residual = float(self.soil.properties.residual[0])
        porosity = float(self.soil.properties.porosity[0])
        span = max(porosity - residual, 1e-12)
        theta_lo = np.clip(
            np.asarray(theta_i, dtype=float),
            residual + 1e-6 * span,
            porosity - 2e-6 * span,
        )
        theta_hi = np.full_like(theta_lo, porosity - 1e-6 * span)

        n_seg = 24
        psi_prev = np.asarray(
            self.soil.water_potential(theta_lo, level=0), dtype=float
        )
        total = np.zeros_like(theta_lo)
        for j in range(1, n_seg + 1):
            theta_j = theta_lo + (theta_hi - theta_lo) * (j / n_seg)
            theta_mid = theta_lo + (theta_hi - theta_lo) * ((j - 0.5) / n_seg)
            psi_j = np.asarray(
                self.soil.water_potential(theta_j, level=0), dtype=float
            )
            K_mid = np.asarray(
                self.soil.conductivity_moisture(theta_mid, level=0),
                dtype=float,
            )
            total += K_mid * (psi_j - psi_prev)
            psi_prev = psi_j

        K_sat = np.asarray(
            self.soil.conductivity_moisture(theta_hi, level=0), dtype=float
        )
        return np.asarray(
            np.maximum(
                np.divide(total, K_sat, out=np.zeros_like(total),
                          where=K_sat > 0.0),
                0.0,
            ),
            dtype=float,
        )

    def _update_infiltration_history(self) -> None:
        """Advances the Green-Ampt cumulative-infiltration event state.

        While raining, the matrix-infiltrated depth accumulates and the
        suction enhancement of the infiltration capacity decays toward
        the sealed K_sat limit. After rain stops, the wetting front
        redistributes, so the event memory relaxes exponentially with
        ``_GA_REDIST_TAU`` and a later storm again sees dry-soil capacity.
        Called once per timestep (after surface coupling) so the repeated
        SMB calls inside the Picard loop do not multi-count.
        """
        precip = np.asarray(
            self.atm_state.precipitation, dtype=float).reshape(-1)
        infil_depth = (
            np.asarray(self._infiltration_flux, dtype=float).reshape(-1)
            * self.tstep / c.water.DENSITY
        )
        decayed = self._ga_cum_infil * np.exp(
            -self.tstep / self._GA_REDIST_TAU)
        self._ga_cum_infil[:] = np.where(
            precip > 0.0, self._ga_cum_infil + infil_depth, decayed
        )

    def _distribute_bypass(self) -> None:
        """Deposits captured macropore water into the crack-zone layers.

        Bypass water descends the macropore network quickly and is
        absorbed by the crack walls, so deposition is weighted
        ``exp(-(z - z_top)/e_folding)`` within [z_top, z_bottom] and
        capped per layer by the rate that fills the layer to field
        capacity this step (crack-wall uptake stalls once the matrix
        near the wall saturates). Water the zone cannot absorb ponds
        and overflows to runoff. The resulting (nz, ncol) mass-rate
        profile [kg/m^2/s] is consumed by both the moisture solve (as a
        volumetric source) and the rain heat-advection source (deposited
        at rain temperature).
        """
        self._bypass_deposit = None
        M = np.asarray(
            getattr(self, '_bypass_flux', 0.0), dtype=float).reshape(-1)
        if not np.any(M > 0.0):
            return

        _, z_top, z_bottom, e_fold = self._macropore_config()
        nz = self.input.grid.nz
        dz = abs(self.input.grid.z[1] - self.input.grid.z[0])
        depth = np.abs(np.asarray(self.input.grid.z, dtype=float))

        theta = np.asarray(self.soil_state.moisture, dtype=float)
        if theta.ndim == 1:
            theta = theta[:, None]
        ncol = theta.shape[1]
        runoff_extra = np.zeros(ncol)

        zone = (depth >= z_top) & (depth <= z_bottom)
        zone[0] = False  # layer 0 is the Dirichlet surface BC
        if not zone.any():
            runoff_extra = M.copy()
            deposit = np.zeros((nz, ncol))
        else:
            theta_fc = np.asarray(self.soil.theta_fc, dtype=float)[:, None]
            headroom = np.where(
                zone[:, None],
                np.maximum(theta_fc - theta, 0.0)
                * dz * c.water.DENSITY / self.tstep,
                0.0,
            )
            weight = np.where(
                zone, np.exp(-(depth - z_top) / max(e_fold, 1e-6)), 0.0
            )[:, None] * np.ones((1, ncol))

            deposit = np.zeros((nz, ncol))
            remaining = M.astype(float).copy()
            # Proportional fill with per-layer caps; a few passes route
            # any capped layer's share to the rest of the zone.
            for _ in range(4):
                room = headroom - deposit
                active = (room > 1e-15) & (weight > 0.0)
                w_act = np.where(active, weight, 0.0)
                w_sum = w_act.sum(axis=0)
                cols = (remaining > 1e-15) & (w_sum > 0.0)
                if not cols.any():
                    break
                alloc = np.where(
                    cols[None, :],
                    w_act * np.divide(
                        remaining, w_sum,
                        out=np.zeros_like(remaining), where=w_sum > 0.0,
                    )[None, :],
                    0.0,
                )
                alloc = np.minimum(alloc, np.maximum(room, 0.0))
                deposit += alloc
                remaining = np.maximum(remaining - alloc.sum(axis=0), 0.0)
            runoff_extra = remaining

        if np.any(runoff_extra > 0.0):
            self.sfc_state.fluxes.runoff[:] = (
                np.asarray(self.sfc_state.fluxes.runoff) + runoff_extra
            )
            self._bypass_flux[:] = M - runoff_extra
        if np.any(deposit > 0.0):
            self._bypass_deposit = deposit

    def _solve_smb(self) -> None:
        """Solves the Surface Moisture Budget (SMB) for surface moisture.

        Finds θ_sfc that balances evaporative demand and precipitation
        infiltration against the Darcy flux from the subsurface via
        bracketed Brent root-finding on θ_sfc ∈ [θ_residual, θ_porosity].
        The residual is

            R(θ_sfc) = E(θ_sfc, T_sfc)
                       + rho_w K_mid (θ_sfc) [ (ψ(θ_sfc) - ψ_1) / Δz + 1 ]
                       - P_eff
                       + rho_w Δz (θ_sfc - θ_old,0) / Δt

        Sign convention: E > 0 is upward (loss to atmosphere); the Darcy
        term > 0 is downward gravity drainage; P_eff > 0 is downward
        mass input that infiltrates into the soil. Positive storage is
        water retained in the top soil control volume. P_eff enters with
        a minus sign because it opposes E, drainage, and storage in the
        surface mass balance.

        Before any runoff, a configurable macropore fraction of the
        throughfall is captured by the crack network (bypass flow) and
        removed from the matrix input; it is deposited at depth once per
        step by ``_distribute_bypass``.

        Runoff has two components. Infiltration-excess (Hortonian) runoff
        removes any rain rate exceeding the transient Green-Ampt
        infiltration capacity ``rho_w * K_sat * (1 + Δθ ψ_f / F)`` of the
        surface soil, so an intense burst late in a storm cannot be
        buffered as storage in the thin top cell while a dry antecedent
        surface still absorbs it. The remaining rate enters the residual;
        saturation-excess runoff is then diagnosed from the SMB transport
        limit if the wet endpoint still cannot carry it. Both excesses are
        accumulated into ``runoff`` and removed from the effective input
        ``P_eff``.

        In canopy mode, the precipitation entering this residual is
        throughfall after filling the intercepted-water store. In bare-soil
        mode, it is the raw precipitation forcing.

        At θ_sfc = θ_residual: K → 0, ψ → -∞, gnd_q → 0, so
        R ≈ -rho_a atm_q u* f_h ≤ 0 (before the residual-endpoint
        regularization).
        At θ_sfc = θ_porosity, after any runoff adjustment: K → K_sat,
        ψ → 0, gnd_q is saturated, so R ≥ 0. The sign change guarantees
        a bracketed root, which Brent's method finds robustly - including
        the dry-soil regime where the previous ψ-inversion approach was
        ill-conditioned.
        """
        RHO_W = c.water.DENSITY

        z_s = self.input.surface.z_s
        z_t = self.input.surface.z_t
        dz = abs(self.input.grid.z[0] - self.input.grid.z[1])
        dt = float(self.tstep)
        if dt <= 0.0:
            raise ValueError('SMB solve requires a positive timestep.')

        atm_p = self.atm_state.pressure
        atm_q = self.atm_state.specific_humidity
        sfc_T = self.sfc_state.temperature
        obukhov_l = self.sfc_state.turbulence.obukhov_length
        ust = self.sfc_state.turbulence.friction_velocity
        rho_a = self.sfc_state.air_density
        fh = self.sfc.fh(z_s, z_t, obukhov_l)

        # Use canopy throughfall when the timestep has been prepared by
        # run(); direct unit tests and partial callers fall back to raw
        # precipitation.
        precip_source = (
            self._soil_precipitation
            if getattr(self, '_precipitation_prepared', False)
            else self.atm_state.precipitation
        )
        P_total = np.asarray(precip_source, dtype=float)
        runoff = np.zeros_like(P_total, dtype=float)
        P_eff = np.array(P_total, copy=True)

        # Macropore bypass capture: a fixed fraction of throughfall runs
        # into the crack network at the surface and never enters the
        # matrix surface balance. Deposition into the crack zone happens
        # once per step in _distribute_bypass(); here it is only removed
        # from the matrix input. Stored in place so the output reference
        # and the per-step distribution see the final coupling iterate.
        f_mac, _, _, _ = self._macropore_config()
        bypass = np.where(P_eff > 0.0, f_mac * P_eff, 0.0)
        P_eff = P_eff - bypass
        bypass_store = getattr(self, '_bypass_flux', None)
        if bypass_store is None or np.shape(bypass_store) != np.shape(bypass):
            self._bypass_flux = np.array(bypass, dtype=float)
        else:
            bypass_store[:] = bypass

        residual_q = float(self.soil.properties.residual[0])
        porosity = float(self.soil.properties.porosity[0])

        # Subsurface properties (fixed during SMB solve)
        theta_old0: np.ndarray = np.asarray(self.soil_state.moisture)[0]
        psi1: np.ndarray = np.asarray(self.soil.water_potential(self.soil_state.moisture))[1]
        K1: np.ndarray = np.asarray(self.soil.conductivity_moisture(self.soil_state.moisture))[1]

        canopy = getattr(self, 'canopy', None)
        if canopy is not None:
            f_veg = canopy.veg_fraction
        else:
            f_veg = np.zeros_like(sfc_T)

        def smb_capacity(theta: np.ndarray) -> np.ndarray:
            theta_c = np.clip(theta, residual_q, porosity)
            psi0 = np.asarray(
                self.soil.water_potential(theta_c, level=0), dtype=float
            )
            K0 = np.asarray(
                self.soil.conductivity_moisture(theta_c, level=0), dtype=float
            )
            K_mid = np.maximum(
                self._moisture_face_conductivity(K0, K1, psi0, psi1, dz),
                1e-14,
            )
            gnd_q = self.soil.surface_specific_humidity(sfc_T, theta_c, atm_p)
            # SMB is a bare-soil top-layer balance; transpiration is
            # removed from the root-zone moisture budget (diffusion RHS),
            # not the surface flux residual.
            E_soil = (1.0 - f_veg) * rho_a * (gnd_q - atm_q) * ust * fh
            darcy = RHO_W * K_mid * ((psi0 - psi1) / dz + 1.0)
            return np.asarray(E_soil + darcy)

        def smb_residual(theta: np.ndarray) -> np.ndarray:
            storage = RHO_W * dz * (theta - theta_old0) / dt
            return np.asarray(smb_capacity(theta) - P_eff + storage)

        # Shrink the bracket slightly off the physical bounds to keep
        # ψ(θ) and K(θ) finite and well-defined at the endpoints.
        eps = 1e-6
        span = porosity - residual_q
        a = np.full_like(sfc_T, residual_q + eps * span)
        b = np.full_like(sfc_T, porosity - eps * span)

        # Infiltration-excess (Hortonian) runoff. Rain in excess of the
        # surface's transient infiltration capacity ``I_max`` cannot enter
        # the column and is shed as runoff. ``I_max`` is the Green-Ampt
        # capacity ``rho_w * K_sat * (1 + delta_theta * psi_f / F)``: the
        # gravity-limited saturated Darcy flux enhanced by wetting-front
        # suction. Early in a storm on dry soil (cumulative infiltration
        # F small, antecedent deficit delta_theta large) the suction term
        # dominates and the surface swallows intense bursts; as F grows
        # the capacity decays toward the sealed/saturated K_sat limit
        # (raindrop-impact sealing, sub-grid saturation heterogeneity).
        # The suction psi_f is derived from the configured retention
        # model in _wetting_front_suction; F is per-event state advanced
        # in _update_infiltration_history. Evaluating conductivity at the
        # wet endpoint avoids the singular dry-end ψ/K that made an
        # earlier K_sat pre-cap ill-conditioned. The remainder feeds the
        # saturation-excess balance below.
        if np.any(P_eff > 0.0):
            K0_sat = np.asarray(
                self.soil.conductivity_moisture(b, level=0), dtype=float
            )
            F_cum = np.maximum(
                np.asarray(getattr(self, '_ga_cum_infil', 0.0), dtype=float),
                self._GA_MIN_CUM_INFIL,
            )
            theta_ante = np.asarray(theta_old0, dtype=float)
            psi_f = self._wetting_front_suction(theta_ante)
            delta_theta = np.maximum(porosity - theta_ante, 0.0)
            infil_capacity = RHO_W * K0_sat * (
                1.0 + delta_theta * psi_f / F_cum
            )
            infil_excess = np.where(
                P_eff > 0.0, np.maximum(P_eff - infil_capacity, 0.0), 0.0
            )
            runoff = runoff + infil_excess
            P_eff = P_eff - infil_excess

        res_a = smb_residual(a)
        res_b = smb_residual(b)
        too_much_rain = (res_a * res_b > 0.0) & (res_b < 0.0) & (P_eff > 0.0)
        if np.any(too_much_rain):
            extra_runoff = np.minimum(P_eff, -res_b)
            runoff = np.where(too_much_rain, runoff + extra_runoff, runoff)
            P_eff = np.where(too_much_rain, P_eff - extra_runoff, P_eff)

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
        self.sfc_state.fluxes.runoff[:] = runoff
        # Water that actually entered the column this step (rain minus
        # interception minus runoff). Consumed by the rain heat-advection
        # source in the soil heat solve.
        self._infiltration_flux = np.maximum(P_eff, 0.0)

    def _compute_fluxes(self, sfc_T: FloatOrArray, sfc_q: FloatOrArray) -> None:
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
         obukhov_l, sensible, latent, converged) = self._solve_most(
            sfc_T_vec, sfc_q_vec, L_init, max_iter=ITER_MAX
        )

        if not converged.all():
            self.logger.warning(
                'Obukhov length did not converge for %d columns.',
                int(np.sum(~converged)),
            )

        self.sfc_state.turbulence.obukhov_length[:] = obukhov_l
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
        diff_T = np.zeros(1)
        diff_q = np.zeros(1)

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

    def _solve_diffusion(self,state_field: np.ndarray,
                         get_diffusivity: Callable[[np.ndarray], np.ndarray],
                         sfc_boundary: FloatOrArray,field_name: str = 'field',
                         source_term: Optional[np.ndarray] = None,
                         avg_diffusivity: str = 'arithmetic') -> None:
        """Solves a generic 1D diffusion equation using a theta scheme.

        This helper serves the soil heat solve. The moisture equation
        uses a dedicated mixed-form Richards implementation.

        Args:
            state_field: Reference to the field to update (temperature) [NDArray].
            get_diffusivity: Callable that computes diffusivity profile from
                soil moisture [Callable[[NDArray[np.float64]], NDArray[np.float64]]].
            sfc_boundary: Surface boundary value for Dirichlet BC.
            field_name: Name of the field for logging/documentation.
            source_term: Optional per-layer source (or sink, if negative)
                with units of the state field per second, shape (nz, ncol).
                Applied as ``dt · source`` to the RHS for layers 1..nz-1.
            avg_diffusivity: Averaging rule for interface diffusivity. Use
                ``"arithmetic"`` for simple means or ``"geometric"`` for
                multiplicative averaging.

        Physics:
            - Diffusivity always depends on soil moisture
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

        # Compute diffusivity using soil moisture
        D = get_diffusivity(self.soil_state.moisture)
        if avg_diffusivity == 'geometric':
            D_mid = np.sqrt(np.maximum(D[:-1] * D[1:], 0.0))
        else:
            D_mid = 0.5 * (D[:-1] + D[1:])

        # === First soil level below surface (i=0) ===
        cp = dt * D_mid[0] / dz2
        cm = dt * D_mid[1] / dz2

        # Backward (implicit) coefficients
        cb_p = -theta_b * cp
        cb_m = -theta_b * cm
        cb = 1.0 - cb_p - cb_m

        # Forward (explicit) coefficients
        cf_p = theta_f * cp
        cf_m = theta_f * cm
        cf = 1.0 - cf_p - cf_m

        f[0] = cb
        g[0] = cb_m
        r[0] = (cf_p * field[0] + cf * field[1] +
                cf_m * field[2] - cb_p * sfc_boundary_vec)

        # === Interior soil levels (Vectorized) ===
        # Define slices to represent indices i, i+1, and i+2
        # Original loop: for i in range(1, nz - 2)
        # indices: 1, 2, ..., nz-3
        idx     = slice(1, nz - 2)  # corresponds to i
        idx_p1  = slice(2, nz - 1)  # corresponds to i+1
        idx_p2  = slice(3, nz)      # corresponds to i+2

        # Compute diffusion coefficients for all interior points
        # D_mid is size (nz-1), so we slice up to nz-2
        cp = dt * D_mid[idx] / dz2
        cm = dt * D_mid[idx_p1] / dz2

        # Backward (implicit) coefficients
        cb_p = -theta_b * cp
        cb_m = -theta_b * cm
        cb   = 1.0 - cb_p - cb_m

        # Forward (explicit) coefficients
        cf_p = theta_f * cp
        cf_m = theta_f * cm
        cf   = 1.0 - cf_p - cf_m

        # Assign coefficients to tridiagonal matrix arrays
        e[idx] = cb_p
        f[idx] = cb
        g[idx] = cb_m

        # Compute the Right Hand Side (RHS) vector r
        # state_field is size (nz)
        r[idx] = (cf_p * field[idx] +
                  cf   * field[idx_p1] +
                  cf_m * field[idx_p2])

        # === Bottom level (Neumann BC: zero gradient) ===
        j = nz - 2
        cp = dt * D_mid[j] / dz2
        cm = dt * D_mid[j] / dz2

        # Backward (implicit) coefficients
        cb_p = -theta_b * cp
        cb_m = -theta_b * cm
        cb = 1.0 - cb_p - cb_m

        # Forward (explicit) coefficients
        cf_p = theta_f * cp
        cf_m = theta_f * cm
        cf = 1.0 - cf_p - cf_m

        # Assign coefficients to tridiagonal matrix arrays.  The bottom
        # Neumann condition uses the symmetric ghost node
        # field[nz] = field[nz - 2], so the missing lower coefficient folds
        # into the sub-diagonal coupling to the layer above.
        e[j] = cb_p + cb_m
        f[j] = cb

        # Compute the Right Hand Side (RHS) vector r
        # state_field is size (nz)
        r[j] = ((cf_p + cf_m) * field[j] +
                cf * field[j + 1])

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

    def _rain_advection_source(self) -> Optional[np.ndarray]:
        """Sensible-heat source from rain advected through the soil column.

        Liquid precipitation enters the soil near the near-surface air
        temperature, which during convective rain is several K cooler than
        the warm daytime soil. The infiltrated water (post-interception,
        post-runoff) carries that enthalpy deficit into the column, then
        drains downward along the Richards water flux, advecting it deeper.
        Soil heat is otherwise pure conduction, so without this term the
        model under-cools during and after rain - and the conduction-only
        reach previously limited the cooling to the top layer.

        Writing the per-layer enthalpy balance with mass conservation, the
        outgoing advective flux cancels against the heat-capacity-storage
        change, leaving a storage-safe form in which each layer is heated
        only by *incoming* water carrying its upwind source temperature
        relative to the layer's own temperature::

            C_vol_k dT_k/dt = conduction + c_w * sum_in m_in (T_src - T_k)

        The surface-input ("term 1") contribution is layer 1's top-face
        inflow ``c_w * P_infil * (T_rain - T_1)``. The interior
        ("term 2") contributions advect heat between prognostic layers
        along the lagged Darcy mass flux ``rho_w * q_face`` (downward
        positive), with the upwind neighbour supplying ``T_src``. Free
        drainage leaves the bottom layer at its own temperature, so it
        contributes nothing. Because every contribution is relative to the
        layer temperature, the source vanishes identically when the soil
        column is isothermal (no spurious heating from net mass storage).

        The Darcy flux is evaluated from the start-of-step soil state
        (the moisture solve runs after the heat solve), so the interior
        advection is explicit/lagged by one step.

        Macropore bypass water descends the crack network with little
        thermal equilibration, so its per-layer deposition profile also
        deposits enthalpy at the rain temperature:
        ``c_w * m_bypass_k * (T_rain - T_k)``.

        Returns:
            A (nz, ncol) temperature tendency [K/s], or None when there is
            no infiltration this step.
        """
        P_infil = np.asarray(
            getattr(self, '_infiltration_flux', 0.0), dtype=float
        ).reshape(-1)
        bypass_deposit = getattr(self, '_bypass_deposit', None)
        if not np.any(P_infil > 0.0) and bypass_deposit is None:
            return None

        # Specific heat of liquid water [J/kg/K], consistent with the
        # volumetric value used for the soil heat capacity.
        c_w = c.water.VOLUMETRIC_HEAT_CAPACITY / c.water.DENSITY
        rho_w = c.water.DENSITY
        nz = self.input.grid.nz
        dz = abs(self.input.grid.z[1] - self.input.grid.z[0])

        T = np.asarray(self.soil_state.temperature, dtype=float)
        theta = np.asarray(self.soil_state.moisture, dtype=float)
        if T.ndim == 1:
            T = T[:, None]
        if theta.ndim == 1:
            theta = theta[:, None]
        ncol = T.shape[1]
        T_rain = np.broadcast_to(
            np.asarray(self.atm_state.temperature, dtype=float).reshape(-1),
            (ncol,),
        )
        P_infil = np.broadcast_to(np.maximum(P_infil, 0.0), (ncol,))
        c_vol = np.asarray(self.soil.heat_capacity(theta), dtype=float)

        # Lagged downward Darcy mass flux at interior faces [kg/m^2/s].
        # Face j sits between nodes j and j+1; face 0 (surface -> layer 1)
        # is intentionally unused, as layer 1's top-face inflow is the
        # actual infiltration P_infil rather than the internal Darcy flux.
        psi = np.asarray(self.soil.water_potential(theta), dtype=float)
        K_node = np.asarray(
            self.soil.conductivity_moisture(theta), dtype=float)
        K_face = self._moisture_face_conductivity(
            K_node[:-1], K_node[1:], psi[:-1], psi[1:], dz)
        m_face = rho_w * K_face * (1.0 - (psi[1:] - psi[:-1]) / dz)

        # Incoming water flux into each prognostic layer (1..nz-1) and the
        # temperature of its upwind source layer.
        inflow_above = np.zeros((nz, ncol))
        Tsrc_above = np.zeros((nz, ncol))
        inflow_above[1] = P_infil
        Tsrc_above[1] = T_rain

        inflow_below = np.zeros((nz, ncol))
        Tsrc_below = np.zeros((nz, ncol))
        if nz >= 3:
            # Layers 2..nz-1 receive downward flow across faces 1..nz-2.
            inflow_above[2:] = np.maximum(m_face[1:], 0.0)
            Tsrc_above[2:] = T[1:nz - 1]
            # Layers 1..nz-2 receive upward flow across the same faces.
            inflow_below[1:nz - 1] = np.maximum(-m_face[1:], 0.0)
            Tsrc_below[1:nz - 1] = T[2:nz]
        # The bottom layer's lower face is free drainage (downward only),
        # so it contributes no incoming flux.

        heating = c_w * (
            inflow_above * (Tsrc_above - T) + inflow_below * (Tsrc_below - T)
        )                                                    # [W/m^2]
        if bypass_deposit is not None:
            dep = np.asarray(bypass_deposit, dtype=float)
            if dep.ndim == 1:
                dep = dep[:, None]
            heating = np.asarray(
                heating + c_w * dep * (T_rain[None, :] - T), dtype=float
            )
        denom = c_vol * dz
        source = np.asarray(
            np.divide(
                heating, denom, out=np.zeros((nz, ncol)), where=denom > 0.0
            ),
            dtype=float,
        )
        source[0] = 0.0
        return source

    def _solve_diffusion_heat(self) -> None:
        """Solves the soil heat diffusion equation using a theta scheme.

        The Dirichlet BC is the soil-top temperature, not the radiative
        skin: when an in-canopy thermal resistance is configured the
        two differ by ``G * r_canopy_thermal``. With no canopy
        resistance the two are identical, recovering the original
        bare-skin behaviour.

        Infiltrating rain advects sensible heat into the top soil layer
        (see ``_rain_advection_source``); with no precipitation the source
        is omitted and the solve is pure conduction.
        """
        self._ensure_writable_soil_field("temperature")
        self._solve_diffusion(
            state_field=self.soil_state.temperature,
            get_diffusivity=self.soil.diffusivity_thermal,
            sfc_boundary=self.sfc_state.soil_top_temperature,
            field_name='temperature',
            source_term=self._rain_advection_source(),
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

        self._ensure_writable_soil_field("moisture")
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

        _2d = moisture.ndim == 2
        residual = self.soil.properties.residual[:, None] if _2d else self.soil.properties.residual
        porosity = self.soil.properties.porosity[:, None] if _2d else self.soil.properties.porosity
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
            K_face = self._moisture_face_conductivity(
                K_node[:-1], K_node[1:], psi_iter[:-1], psi_iter[1:], dx
            )

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
        balance. Macropore bypass water enters the same way: its
        deposition profile [kg/m^2/s] becomes a volumetric source in the
        crack zone, adding the water the matrix surface never saw.
        """
        source = None
        if getattr(self, 'canopy', None) is not None:
            source = -self.canopy_state.root_uptake
        bypass_deposit = getattr(self, '_bypass_deposit', None)
        if bypass_deposit is not None:
            dz = abs(self.input.grid.z[1] - self.input.grid.z[0])
            bypass_source = (
                np.asarray(bypass_deposit, dtype=float)
                / (c.water.DENSITY * dz)
            )
            if source is None:
                source = bypass_source
            else:
                src_arr = np.asarray(source, dtype=float)
                if src_arr.ndim == 1:
                    src_arr = src_arr[:, None]
                source = src_arr + bypass_source
        self._solve_mixed_moisture(source_term=source)
        tol = max(float(self.input.numerics.tolerances.moisture_bounds), 1e-8)
        self.soil.enforce_moisture_bounds(self.soil_state.moisture, tol)
