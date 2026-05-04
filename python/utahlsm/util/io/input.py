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
"""Handles all input data loading and validation for UtahLSM.

This module defines the `Input` class, which is responsible for reading the
JSON namelist, initial conditions from a NetCDF file, and offline forcing
data. It validates the inputs, populates the configuration and state
dataclasses, and provides a single, clean interface for the main model to
access all setup information.
"""
import json
import logging
from importlib import resources
from typing import Any, Optional, cast

import jsonschema
import netCDF4 as nc
import numpy as np
from numpy.typing import NDArray

from ...data_models import (
    AtmosphericState,
    CanopyConfig,
    ForcingData,
    GeneralConfig,
    GridConfig,
    IterationsConfig,
    NumericsConfig,
    OutputConfig,
    RadiationConfig,
    SoilConfig,
    SoilState,
    SurfaceConfig,
    TimeConfig,
    TolerancesConfig,
)
from ...exceptions import NamelistError
from . import logging_helper
from .soil_properties_loader import SoilPropertiesLoader


class Input:
    """Orchestrates the loading and validation of all model inputs.

    This class reads configuration from a JSON namelist file and initial
    conditions from a NetCDF file. It also handles optional offline forcing
    data. All data is validated and organized into the appropriate dataclasses.

    Attributes:
        logger: A logger for this class.
        general: Dataclass with general simulation settings.
        numerics: Dataclass with numerical scheme parameters.
        time: Dataclass with time-related parameters.
        surface: Dataclass with surface-related parameters.
        soil: Dataclass with soil model configuration.
        radiation: Dataclass with radiation model configuration.
        canopy: Dataclass with canopy / vegetation model configuration.
        output: Dataclass with output file configuration.
        grid: Dataclass with grid and spatial discretization parameters.
        initial: Dataclass holding the initial soil state.
        forcing: Dataclass holding the time-series of offline forcing data,
            or None if not provided.
        soil_properties: Dictionary mapping soil type names to their properties.
        soil_properties_name: Name of the soil property dataset being used.
        soil_type_names: List of soil type names for each layer.
    """

    def __init__(self, namelist_path: str, inputfile: str,
                 offlinefile: Optional[str] = None) -> None:
        """Initializes the Input class and loads all data.

        Args:
            namelist_path: The file path to the JSON namelist.
            inputfile: The file path to the NetCDF initial conditions file.
            offlinefile: The optional file path to the NetCDF offline
                forcing file. Defaults to None.

        Raises:
            FileNotFoundError: If any of the required input files do not exist.
            json.JSONDecodeError: If the namelist JSON file is malformed.
            jsonschema.ValidationError: If the namelist does not conform to the
                expected schema.
        """
        self.logger: logging.Logger = logging_helper.get_logger('Input')
        self.logger.info('Reading %s', namelist_path)
        namelist_data = self._load_and_validate_namelist(namelist_path)
        log_level = namelist_data['general']['log_level']
        logging_helper.finalize_logging(log_level)

        self.logger.info('Reading %s', inputfile)
        nx = namelist_data['grid']['nx']
        ny = namelist_data['grid']['ny']
        nz = namelist_data['grid']['nz']
        init_data = self._load_initial_conditions(inputfile, nx, ny, nz)

        # Load soil properties before creating configs
        self.logger.info('Loading soil properties')
        self._load_soil_properties(namelist_data['soil'], init_data['type'])

        self.general: GeneralConfig = GeneralConfig(**namelist_data['general'])
        iterations_data = namelist_data['numerics']['iterations']
        tolerances_data = namelist_data['numerics']['tolerances']
        numerics_data = namelist_data['numerics']
        self.numerics = NumericsConfig(
            heat_diffusion_back_weight=(
                numerics_data['heat_diffusion_back_weight']),
            warm_start_turbulence=bool(
                numerics_data.get('warm_start_turbulence', False)),
            initialize_surface_temperature_from_seb=bool(
                numerics_data.get('initialize_surface_temperature_from_seb', False)),
            iterations=IterationsConfig(**iterations_data),
            tolerances=TolerancesConfig(**tolerances_data)
        )
        self.time: TimeConfig = TimeConfig(**namelist_data['time'])
        self.surface: SurfaceConfig = SurfaceConfig(**namelist_data['surface'])
        self.soil: SoilConfig = SoilConfig(**namelist_data['soil'])
        rad_data = namelist_data['radiation']
        self.radiation: RadiationConfig = RadiationConfig(**rad_data)
        self.output: OutputConfig = OutputConfig(**namelist_data['output'])
        self.canopy: CanopyConfig = CanopyConfig(
            **namelist_data.get('canopy', {}))
        self.grid: GridConfig = GridConfig(
            nx=nx,
            ny=ny,
            nz=nz,
            z=init_data['z']
        )
        self.initial: SoilState = SoilState(
            temperature=init_data['temperature'],
            moisture=init_data['moisture'],
            type=init_data['type']
        )

        self.forcing: Optional[ForcingData] = None
        if offlinefile:
            self.logger.info('Reading %s', offlinefile)
            self._load_offline_data(offlinefile)

        self._validate_physical_consistency()

    @staticmethod
    def _validate_uniform_soil_z(z: NDArray[np.float64]) -> NDArray[np.float64]:
        """Validates that the soil grid uses uniform vertical spacing."""
        z_arr = np.asarray(z, dtype=float)
        if z_arr.ndim != 1:
            raise ValueError(
                f"soil_z must be one-dimensional after loading, got {z_arr.ndim}D."
            )
        if z_arr.size <= 1:
            return z_arr
        dz = np.diff(z_arr)
        if not np.allclose(dz, dz[0], rtol=1e-10, atol=1e-12):
            raise ValueError(
                "soil_z must have uniform spacing; non-uniform vertical grids "
                "are not supported by the current solver."
            )
        return z_arr

    def _load_and_validate_namelist(self, namelist_path: str) -> dict[str, Any]:
        """Loads and validates the JSON namelist against a schema.

        Args:
            namelist_path: The path to the JSON namelist file.

        Returns:
            A dictionary containing the validated namelist data.

        Raises:
            FileNotFoundError: If the namelist or schema file cannot be
                found.
            json.JSONDecodeError: If the namelist is not valid JSON.
            jsonschema.ValidationError: If the namelist does not match
                the schema.
        """
        try:
            schema_resource = resources.files('utahlsm.util.io').joinpath(
                'schema_namelist.json')
            schema = json.loads(schema_resource.read_text(encoding='utf-8'))
            with open(namelist_path, encoding='utf-8') as f:
                namelist_data = json.load(f)
            jsonschema.validate(instance=namelist_data, schema=schema)
            self.logger.info('--- namelist validation successful')
            return cast(dict[str, Any], namelist_data)
        except (FileNotFoundError, json.JSONDecodeError,
                jsonschema.ValidationError) as e:
            self.logger.error('--- namelist error: %s', e)
            raise

    def _load_initial_conditions(
            self, inputfile: str, nx: int, ny: int,
            nz: int) -> dict[str, NDArray[Any]]:
        """Loads data from the NetCDF initialization file.

        Args:
            inputfile: The path to the NetCDF initial conditions file.
            nx: Number of grid columns in the x direction.
            ny: Number of grid columns in the y direction.
            nz: Number of soil layers expected in the input file.

        Returns:
            A dictionary of NumPy arrays for soil depth, temperature,
            moisture, and type. The 'type' array contains soil type names
            as strings.

        Raises:
            IOError: If the file cannot be read.
            KeyError: If a required variable is missing from the NetCDF
                file.
        """
        try:
            with nc.Dataset(inputfile) as inifile:
                inifile.set_auto_mask(False)
                soil_z_var = inifile.variables['soil_z'][:]
                soil_T_var = inifile.variables['soil_T'][:]
                soil_q_var = inifile.variables['soil_q'][:]
                soil_type_var = inifile.variables['soil_type'][:]

                ncol = nx * ny

                def _decode_soil_type(array: NDArray[Any]) -> NDArray[Any]:
                    if array.dtype.kind == 'S':
                        return np.char.decode(array, 'utf-8').astype(object)
                    if array.dtype.kind in ('U', 'O'):
                        return np.asarray(array, dtype=object).astype(str)
                    raise ValueError(
                        f"soil_type variable must contain strings, "
                        f"got dtype {array.dtype}"
                    )

                def _ensure_z_1d(soil_z: NDArray[Any]) -> NDArray[np.float64]:
                    z = (-1) * soil_z.astype('float')
                    if z.ndim == 1:
                        if z.shape[0] != nz:
                            raise ValueError(
                                f"soil_z length {z.shape[0]} does not match "
                                f"namelist nz={nz}."
                            )
                        return self._validate_uniform_soil_z(z)
                    if z.ndim == 2:
                        if z.shape == (nz, ncol):
                            ref = z[:, 0]
                            if not np.allclose(z, ref[:, None]):
                                raise ValueError(
                                    "soil_z varies across columns; "
                                    "horizontal variation is not supported."
                                )
                            return self._validate_uniform_soil_z(ref)
                        raise ValueError(
                            f"soil_z shape {z.shape} must be (nz, ncol) or "
                            f"(nz,) when using flattened columns."
                        )
                    if z.ndim == 3:
                        if z.shape != (nz, ny, nx):
                            raise ValueError(
                                f"soil_z shape {z.shape} does not match "
                                f"(nz, ny, nx)=({nz}, {ny}, {nx})."
                            )
                        ref = z[:, 0, 0]
                        if not np.allclose(z, ref[:, None, None]):
                            raise ValueError(
                                "soil_z varies across columns; "
                                "horizontal variation is not supported."
                            )
                        return self._validate_uniform_soil_z(ref)
                    raise ValueError(
                        f"soil_z has unsupported dimensions: {z.ndim}."
                    )

                def _reshape_soil_field(
                    field: NDArray[Any], name: str
                ) -> NDArray[np.float64]:
                    data = field.astype('float')
                    if data.ndim == 1:
                        if data.shape[0] != nz:
                            raise ValueError(
                                f"{name} length {data.shape[0]} does not "
                                f"match namelist nz={nz}."
                            )
                        if ncol == 1:
                            return data[:, None]
                        return np.repeat(data[:, None], ncol, axis=1)
                    if data.ndim == 2:
                        if data.shape == (nz, ncol):
                            return data
                        if data.shape == (nz, 1) and ncol > 1:
                            return np.repeat(data, ncol, axis=1)
                        raise ValueError(
                            f"{name} shape {data.shape} must be (nz, ncol) "
                            f"or (nz,) for single-column runs."
                        )
                    if data.ndim == 3:
                        if data.shape != (nz, ny, nx):
                            raise ValueError(
                                f"{name} shape {data.shape} does not match "
                                f"(nz, ny, nx)=({nz}, {ny}, {nx})."
                            )
                        return data.reshape(nz, ncol)
                    raise ValueError(
                        f"{name} has unsupported dimensions: {data.ndim}."
                    )

                def _reshape_soil_type(
                    field: NDArray[Any], name: str
                ) -> NDArray[Any]:
                    data = _decode_soil_type(field)
                    if data.ndim == 1:
                        if data.shape[0] != nz:
                            raise ValueError(
                                f"{name} length {data.shape[0]} does not "
                                f"match namelist nz={nz}."
                            )
                        return data
                    if data.ndim == 2:
                        if data.shape != (nz, ncol):
                            raise ValueError(
                                f"{name} shape {data.shape} must be "
                                f"(nz, ncol) when using flattened columns."
                            )
                        ref = data[:, 0]
                        if not np.all(data == ref[:, None]):
                            raise ValueError(
                                "soil_type varies across columns; "
                                "horizontal variation is not supported."
                            )
                        return ref
                    if data.ndim == 3:
                        if data.shape != (nz, ny, nx):
                            raise ValueError(
                                f"{name} shape {data.shape} does not match "
                                f"(nz, ny, nx)=({nz}, {ny}, {nx})."
                            )
                        ref = data[:, 0, 0]
                        if not np.all(data == ref[:, None, None]):
                            raise ValueError(
                                "soil_type varies across columns; "
                                "horizontal variation is not supported."
                            )
                        return ref
                    raise ValueError(
                        f"{name} has unsupported dimensions: {data.ndim}."
                    )

                init_dict: dict[str, NDArray[Any]] = {
                    'z': _ensure_z_1d(soil_z_var),
                    'temperature': _reshape_soil_field(soil_T_var, "soil_T"),
                    'moisture': _reshape_soil_field(soil_q_var, "soil_q"),
                    'type': _reshape_soil_type(soil_type_var, "soil_type"),
                }
            self.logger.info('--- initial conditions loaded successfully')
            return init_dict
        except (OSError, KeyError) as e:
            self.logger.error('--- initial conditions error: %s', e)
            raise

    def _load_soil_properties(
        self, soil_config: dict[str, Any], soil_type_array: NDArray[np.str_]
    ) -> None:
        """Loads soil properties from JSON files and validates soil types.

        Args:
            soil_config: Dictionary from namelist with 'properties' and 'model'.
            soil_type_array: Array of soil type names (strings) from initial
                conditions.

        Raises:
            NamelistError: If soil property file not found or invalid, or if
                soil types in initial conditions are not found in properties.
        """
        try:
            properties_spec = soil_config['properties']
            self.soil_properties = SoilPropertiesLoader.load(properties_spec)
            self.soil_properties_name: str = properties_spec

            # Validate that all soil types in initial conditions are available
            # in the loaded properties
            self.soil_type_names: list[str] = []
            for soil_type_name in soil_type_array:
                soil_type_lower = soil_type_name.lower()
                if soil_type_lower not in self.soil_properties:
                    available = ', '.join(sorted(self.soil_properties.keys()))
                    raise NamelistError(
                        f"Soil type '{soil_type_name}' from initial conditions "
                        f"not found in properties dataset '{properties_spec}'. "
                        f"Available soil types: {available}"
                    )
                self.soil_type_names.append(soil_type_lower)

            self.logger.info(
                'Soil properties loaded from: %s', properties_spec
            )
        except NamelistError:
            raise
        except Exception as e:
            self.logger.error('Error loading soil properties: %s', e)
            raise NamelistError(
                f'Failed to load soil properties: {e}'
            ) from e

    def _load_offline_data(self, offlinefile: str) -> None:
        """Loads data from the NetCDF offline forcing file.

        Args:
            offlinefile: The path to the NetCDF offline forcing file.

        Raises:
            IOError: If the file cannot be read.
            KeyError: If a required variable is missing from the NetCDF
                file.
        """
        try:
            with nc.Dataset(offlinefile) as metfile:
                metfile.set_auto_mask(False)
                ntime = len(metfile.dimensions['t'])
                tstep = metfile.variables['tstep'][0].astype('float')
                atm_U = metfile.variables['atm_U'][:]
                atm_T = metfile.variables['atm_T'][:]
                atm_q = metfile.variables['atm_q'][:]
                atm_p = metfile.variables['atm_p'][:]
                sw_in = metfile.variables['sw_in'][:]
                sw_out = metfile.variables['sw_out'][:]
                lw_in = metfile.variables['lw_in'][:]
                lw_out = metfile.variables['lw_out'][:]

                ncol = self.grid.nx * self.grid.ny
                ny = self.grid.ny
                nx = self.grid.nx

                def _reshape_forcing(field: NDArray[Any], name: str) -> NDArray[Any]:
                    data = field.astype('float')
                    if data.ndim == 1:
                        if data.shape[0] != ntime:
                            raise ValueError(
                                f"{name} length {data.shape[0]} does not "
                                f"match forcing ntime={ntime}."
                            )
                        if ncol == 1:
                            return data[:, None]
                        return np.repeat(data[:, None], ncol, axis=1)
                    if data.ndim == 2:
                        if data.shape == (ntime, ncol):
                            return data
                        if data.shape == (ntime, 1) and ncol > 1:
                            return np.repeat(data, ncol, axis=1)
                        raise ValueError(
                            f"{name} shape {data.shape} must be (ntime, ncol) "
                            f"or (ntime,) for single-column runs."
                        )
                    if data.ndim == 3:
                        if data.shape != (ntime, ny, nx):
                            raise ValueError(
                                f"{name} shape {data.shape} does not match "
                                f"(ntime, ny, nx)=({ntime}, {ny}, {nx})."
                            )
                        return data.reshape(ntime, ncol)
                    raise ValueError(
                        f"{name} has unsupported dimensions: {data.ndim}."
                    )

                atm_U = _reshape_forcing(atm_U, "atm_U")
                atm_T = _reshape_forcing(atm_T, "atm_T")
                atm_q = _reshape_forcing(atm_q, "atm_q")
                atm_p = _reshape_forcing(atm_p, "atm_p")
                sw_in = _reshape_forcing(sw_in, "sw_in")
                sw_out = _reshape_forcing(sw_out, "sw_out")
                lw_in = _reshape_forcing(lw_in, "lw_in")
                lw_out = _reshape_forcing(lw_out, "lw_out")

                # Net radiation is derived from the four components so the
                # SEB residual and any component-level consumers (e.g. the
                # Jarvis f1 stress factor reading sw_in) stay consistent.
                r_net = sw_in - sw_out + lw_in - lw_out

                # Validate forcing data and clip only minor boundary excursions.
                self._validate_forcing_data(atm_U, atm_T, atm_q, atm_p,
                                            sw_in, sw_out, lw_in, lw_out,
                                            r_net, ntime)

                atm_data = [
                    AtmosphericState(
                        wind_speed=atm_U[i], temperature=atm_T[i],
                        specific_humidity=atm_q[i], pressure=atm_p[i],
                        sw_in=sw_in[i], sw_out=sw_out[i],
                        lw_in=lw_in[i], lw_out=lw_out[i],
                        radiation_net=r_net[i])
                    for i in range(ntime)
                ]

                self.forcing = ForcingData(ntime=ntime, tstep=tstep,
                                           atmos=atm_data)
                self.logger.info(
                    '--- loaded %d timesteps of forcing data', ntime)
        except (OSError, KeyError) as e:
            self.logger.error('--- offline forcing error: %s', e)
            raise

    def _validate_forcing_data(self, atm_U: NDArray[np.float64], atm_T: NDArray[np.float64],
                               atm_q: NDArray[np.float64], atm_p: NDArray[np.float64],
                               sw_in: NDArray[np.float64], sw_out: NDArray[np.float64],
                               lw_in: NDArray[np.float64], lw_out: NDArray[np.float64],
                               r_net: NDArray[np.float64], _ntime: int) -> None:
        """Validates atmospheric forcing data for physical consistency.

        Slight excursions beyond the supported forcing bounds are clipped to
        preserve robustness against boundary-value artifacts. Larger
        violations raise an error so invalid experiments do not continue with
        silently modified forcing.

        Args:
            atm_U: Wind speed array with one value per time step [m/s].
            atm_T: Temperature array with one value per time step [K].
            atm_q: Specific humidity array with one value per time step [kg/kg].
            atm_p: Pressure array with one value per time step [Pa].
            sw_in: Downwelling shortwave radiation [W/m²].
            sw_out: Upwelling (reflected) shortwave radiation [W/m²].
            lw_in: Downwelling longwave radiation [W/m²].
            lw_out: Upwelling (emitted) longwave radiation [W/m²].
            r_net: Derived net radiation [W/m²].
            _ntime: Number of time steps in forcing arrays (unused but kept for
                API compatibility with other validation functions).

        Raises:
            ValueError: If forcing data contains values outside the supported
                bounds by more than the configured clipping tolerance.
        """
        def _clip_or_raise(
                data: NDArray[np.float64], *, name: str, lower: float, upper: float,
                lower_tol: float, upper_tol: float, units: str) -> bool:
            below = data < lower
            above = data > upper
            out_of_range = below | above
            if not np.any(out_of_range):
                return False

            small_below = below & ((lower - data) <= lower_tol)
            small_above = above & ((data - upper) <= upper_tol)
            small_excursions = small_below | small_above
            hard_failures = out_of_range & ~small_excursions

            if np.any(hard_failures):
                bad_values: NDArray[np.float64] = data[hard_failures]
                sample_indices = np.flatnonzero(hard_failures)[:5].tolist()
                raise ValueError(
                    f'Offline forcing {name} contains {bad_values.size} '
                    f'entries outside [{lower:.6g}, {upper:.6g}] {units} by '
                    f'more than the allowed tolerance '
                    f'(-{lower_tol:.6g}/+{upper_tol:.6g} {units}). '
                    f'Observed range {bad_values.min():.6g} to '
                    f'{bad_values.max():.6g} {units}; sample flat indices '
                    f'{sample_indices}.'
                )

            num_clipped = int(np.count_nonzero(small_excursions))
            self.logger.warning(
                'Clipping %d forcing entries for %s to [%f, %f] %s.',
                num_clipped, name, lower, upper, units)
            data[small_excursions] = np.clip(
                data[small_excursions], lower, upper)
            return True

        issues_found = False
        issues_found |= _clip_or_raise(
            atm_T, name='temperature', lower=200.0, upper=350.0,
            lower_tol=0.5, upper_tol=0.5, units='K')
        issues_found |= _clip_or_raise(
            atm_p, name='pressure', lower=50000.0, upper=110000.0,
            lower_tol=100.0, upper_tol=100.0, units='Pa')
        issues_found |= _clip_or_raise(
            atm_q, name='humidity', lower=0.0, upper=0.05,
            lower_tol=1e-4, upper_tol=5e-4, units='kg/kg')
        issues_found |= _clip_or_raise(
            atm_U, name='wind speed', lower=1e-4, upper=50.0,
            lower_tol=1e-4, upper_tol=0.5, units='m/s')
        issues_found |= _clip_or_raise(
            sw_in, name='SW_in', lower=0.0, upper=1400.0,
            lower_tol=5.0, upper_tol=25.0, units='W/m^2')
        issues_found |= _clip_or_raise(
            sw_out, name='SW_out', lower=0.0, upper=1400.0,
            lower_tol=5.0, upper_tol=25.0, units='W/m^2')
        issues_found |= _clip_or_raise(
            lw_in, name='LW_in', lower=100.0, upper=600.0,
            lower_tol=10.0, upper_tol=10.0, units='W/m^2')
        issues_found |= _clip_or_raise(
            lw_out, name='LW_out', lower=100.0, upper=700.0,
            lower_tol=10.0, upper_tol=10.0, units='W/m^2')
        issues_found |= _clip_or_raise(
            r_net, name='net radiation', lower=-200.0, upper=1200.0,
            lower_tol=25.0, upper_tol=25.0, units='W/m^2')

        if issues_found:
            self.logger.info(
                'Forcing data validation: clipped minor boundary excursions. '
                'Please review input data quality.')
        else:
            self.logger.info(
                'Forcing data validation: All variables within '
                'expected ranges')

    def _validate_physical_consistency(self) -> None:
        """Performs validation checks on inter-variable relationships.

        Raises:
            ValueError: If a physical consistency check fails.
        """
        # Grid size validation
        if self.grid.nz < 3:
            raise ValueError(
                f"Grid must have at least 3 soil layers for diffusion "
                f"solvers (3-point stencil), got nz={self.grid.nz}.")

        if self.surface.z_m <= self.surface.z_o:
            raise ValueError(
                f"z_m={self.surface.z_m} must be > "
                f"z_o={self.surface.z_o}.")

        if self.surface.z_s <= self.surface.z_t:
            raise ValueError(
                f"z_s={self.surface.z_s} must be > "
                f"z_t={self.surface.z_t}.")

        self._validate_uniform_soil_z(self.grid.z)

        ncol = self.grid.nx * self.grid.ny
        soil_temp = np.asarray(self.initial.temperature)
        soil_mois = np.asarray(self.initial.moisture)

        if soil_temp.ndim == 1:
            if soil_temp.shape[0] != self.grid.nz:
                raise ValueError(
                    f"Namelist nlevs={self.grid.nz} does not match "
                    f"init file soil_T length of {soil_temp.shape[0]}."
                )
        elif soil_temp.ndim == 2:
            if soil_temp.shape != (self.grid.nz, ncol):
                raise ValueError(
                    f"init file soil_T shape {soil_temp.shape} does not "
                    f"match (nz, ncol)=({self.grid.nz}, {ncol})."
                )
        else:
            raise ValueError(
                f"init file soil_T has unsupported dimensions: "
                f"{soil_temp.ndim}."
            )

        if soil_mois.ndim == 1:
            if soil_mois.shape[0] != self.grid.nz:
                raise ValueError(
                    f"Namelist nlevs={self.grid.nz} does not match "
                    f"init file soil_q length of {soil_mois.shape[0]}."
                )
        elif soil_mois.ndim == 2:
            if soil_mois.shape != (self.grid.nz, ncol):
                raise ValueError(
                    f"init file soil_q shape {soil_mois.shape} does not "
                    f"match (nz, ncol)=({self.grid.nz}, {ncol})."
                )
        else:
            raise ValueError(
                f"init file soil_q has unsupported dimensions: "
                f"{soil_mois.ndim}."
            )

        if len(self.initial.type) != self.grid.nz:
            raise ValueError(
                f"Namelist nlevs={self.grid.nz} does not match "
                f"init file soil_type length of {len(self.initial.type)}.")

        # Validate soil moisture values are within physically possible bounds
        # Moisture must be >= 0 and <= porosity (will be validated
        # against residual later). Issues warnings instead of errors to
        # allow running with imperfect data
        moisture_by_layer = soil_mois
        if moisture_by_layer.ndim == 1:
            moisture_by_layer = moisture_by_layer[:, None]

        for i, soil_type_name in enumerate(self.soil_type_names):
            layer_moisture = moisture_by_layer[i]
            if np.any(layer_moisture < 0):
                self.logger.warning(
                    'Layer %d: soil moisture has negative values. '
                    'Moisture must be >= 0.', i)
            # Get porosity for this soil type to validate upper bound
            soil_type_lower = soil_type_name.lower()
            if soil_type_lower not in self.soil_properties:
                self.logger.warning(
                    'Layer %d: soil type %s not found in properties.',
                    i, soil_type_name)
            else:
                porosity = self.soil_properties[soil_type_lower]['porosity']
                if np.any(layer_moisture > porosity):
                    self.logger.warning(
                        'Layer %d: soil moisture exceeds porosity %f.',
                        i, porosity)

        self.logger.info('Physical consistency checks passed')
