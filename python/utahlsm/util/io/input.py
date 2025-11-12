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
"""Handles all input data loading and validation for UtahLSM.

This module defines the `Input` class, which is responsible for reading the
JSON namelist, initial conditions from a NetCDF file, and offline forcing
data. It validates the inputs, populates the configuration and state
dataclasses, and provides a single, clean interface for the main model to
access all setup information.
"""
from dataclasses import dataclass
import json
import jsonschema
import netCDF4 as nc
import numpy as np
from numpy.typing import NDArray
from typing import Dict, List, Optional
import logging

from ...data_models import (
    GeneralConfig, NumericsConfig, IterationsConfig, TolerancesConfig,
    TimeConfig, GridConfig, SurfaceConfig, SoilConfig, RadiationConfig, 
    OutputConfig, SoilState, ForcingData, AtmosphericState
)
from utahlsm.util.io import logging_helper

class Input(object):
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
        output: Dataclass with output file configuration.
        grid: Dataclass with grid and spatial discretization parameters.
        initial: Dataclass holding the initial soil state.
        forcing: Dataclass holding the time-series of offline forcing data,
            or None if not provided.
    """
    
    def __init__(self, namelist_path: str, inputfile: str, offlinefile: Optional[str] = None) -> None:
        """Initializes the Input class and loads all data.
        
        Args:
            namelist_path: The file path to the JSON namelist.
            inputfile: The file path to the NetCDF initial conditions file.
            offlinefile: The optional file path to the NetCDF offline
                forcing file. Defaults to None.
        """
        self.logger: logging.Logger = logging_helper.get_logger("Input")
        self.logger.info(f"Reading {namelist_path}")
        namelist_data = self._load_and_validate_namelist(namelist_path)
        log_level = namelist_data["general"]["log_level"]
        logging_helper.finalize_logging(log_level)
        
        self.logger.info(f"Reading {inputfile}")
        init_data = self._load_initial_conditions(inputfile)
        
        self.general: GeneralConfig = GeneralConfig(**namelist_data["general"])
        iterations_data = namelist_data["numerics"]["iterations"]
        tolerances_data = namelist_data["numerics"]["tolerances"]
        self.numerics = NumericsConfig(
            diffusion_back_weight=namelist_data["numerics"]["diffusion_back_weight"],
            iterations=IterationsConfig(**iterations_data),
            tolerances=TolerancesConfig(**tolerances_data)
        )
        self.time: TimeConfig = TimeConfig(**namelist_data["time"])
        self.surface: SurfaceConfig = SurfaceConfig(**namelist_data["surface"])
        self.soil: SoilConfig = SoilConfig(**namelist_data["soil"])
        self.radiation: RadiationConfig = RadiationConfig(**namelist_data["radiation"])
        self.output: OutputConfig = OutputConfig(**namelist_data["output"])
        self.grid: GridConfig = GridConfig(nx = namelist_data["grid"]["nx"],
                                           ny = namelist_data["grid"]["ny"],
                                           nz = namelist_data["grid"]["nz"],
                                           z  = init_data["z"]
                                        )
        self.initial: SoilState = SoilState(
            temperature = init_data["temperature"],
            moisture = init_data["moisture"],
            type = init_data["type"]
        )
        
        self.forcing: Optional[ForcingData] = None
        if offlinefile:
            self.logger.info(f"Reading {offlinefile}")
            self._load_offline_data(offlinefile)
            
        self._validate_physical_consistency()

    def _load_and_validate_namelist(self, namelist_path: str) -> Dict:
        """Loads and validates the JSON namelist against a schema.
        
        Args:
            namelist_path: The path to the JSON namelist file.
        
        Returns:
            A dictionary containing the validated namelist data.
        
        Raises:
            FileNotFoundError: If the namelist or schema file cannot be found.
            json.JSONDecodeError: If the namelist is not valid JSON.
            jsonschema.ValidationError: If the namelist does not match the schema.
        """
        schema_path = "utahlsm/util/io/schema_namelist.json"
        
        try:
            with open(schema_path) as f: 
                schema = json.load(f)
            with open(namelist_path) as f: 
                namelist_data = json.load(f)
            jsonschema.validate(instance=namelist_data, schema=schema)
            self.logger.info("--- namelist validation successful")
            return namelist_data
        except (FileNotFoundError, json.JSONDecodeError, jsonschema.ValidationError) as e:
            self.logger.error(f"--- namelist error: {e}")
            raise
    
    def _load_initial_conditions(self, inputfile: str)-> Dict[str, NDArray[np.float64]]:
        """Loads data from the NetCDF initialization file.
        
        Args:
            inputfile: The path to the NetCDF initial conditions file.
        
        Returns:
            A dictionary of NumPy arrays for soil depth, temperature,
            moisture, and type.
        
        Raises:
            IOError: If the file cannot be read.
            KeyError: If a required variable is missing from the NetCDF file.
        """
        try:
            with nc.Dataset(inputfile) as inifile:
                inifile.set_auto_mask(False)
                init_dict = {
                    "z" : (-1)*inifile.variables['soil_z'][:].astype('float'),
                    "temperature" : inifile.variables['soil_T'][:].astype('float'),
                    "moisture" : inifile.variables['soil_q'][:].astype('float'),
                    "type" : inifile.variables['soil_type'][:].astype('int')
                }
            self.logger.info("--- initial conditions loaded successfully")
            return init_dict
        except (IOError, KeyError) as e:
            self.logger.error(f"--- initial conditions error: {e}")
            raise
    
    def _load_offline_data(self, offlinefile: str) -> None:
        """Loads data from the NetCDF offline forcing file.
        
        Args:
            offlinefile: The path to the NetCDF offline forcing file.
        
        Raises:
            IOError: If the file cannot be read.
            KeyError: If a required variable is missing from the NetCDF file.
        """
        try:
            with nc.Dataset(offlinefile) as metfile:
                metfile.set_auto_mask(False)
                ntime    = len(metfile.dimensions['t'])
                tstep    = metfile.variables['tstep'][0].astype('float')
                atm_U    = metfile.variables['atm_U'][:].astype('float')
                atm_T    = metfile.variables['atm_T'][:].astype('float')
                atm_q    = metfile.variables['atm_q'][:].astype('float')
                atm_p    = metfile.variables['atm_p'][:].astype('float')
                r_net    = metfile.variables['R_net'][:].astype('float')

                # Ensure wind speed is never zero or negative (causes numerical issues)
                min_wind_speed = 1.0e-4  # m/s
                zero_wind_mask = atm_U <= 0
                num_zero_wind = np.sum(zero_wind_mask)
                if num_zero_wind > 0:
                    self.logger.warning(
                        f"Found {num_zero_wind} timesteps with zero/negative wind speed. "
                        f"Setting to minimum threshold {min_wind_speed} m/s."
                    )
                    atm_U[zero_wind_mask] = min_wind_speed

                atm_data = [
                    AtmosphericState(wind_speed=atm_U[i], temperature=atm_T[i], specific_humidity=atm_q[i], pressure=atm_p[i], radiation_net=r_net[i])
                    for i in range(ntime)
                ]

                self.forcing = ForcingData(ntime=ntime, tstep=tstep, atmos=atm_data)
                self.logger.info(f"--- loaded {ntime} timesteps of forcing data")
        except (IOError, KeyError) as e:
            self.logger.error(f"--- offline forcing error: {e}")
            raise
    
    def _validate_physical_consistency(self) -> None:
        """Performs validation checks on inter-variable relationships.

        Raises:
            ValueError: If a physical consistency check fails.
        """
        # Grid size validation
        if self.grid.nz < 2:
            raise ValueError(f"Grid must have at least 2 soil layers for diffusion solvers, "
                           f"got nz={self.grid.nz}.")

        if self.surface.z_m <= self.surface.z_o:
            raise ValueError(f"z_m={self.surface.z_m} must be > z_o={self.surface.z_o}.")

        if self.surface.z_s <= self.surface.z_t:
            raise ValueError(f"z_s={self.surface.z_s} must be > z_t={self.surface.z_t}.")

        if len(self.initial.temperature) != self.grid.nz:
            raise ValueError(f"Namelist nlevs={self.grid.nz} does not match "
                             f"init file soil_T length of {len(self.initial.temperature)}.")

        if len(self.initial.moisture) != self.grid.nz:
            raise ValueError(f"Namelist nlevs={self.grid.nz} does not match "
                             f"init file soil_q length of {len(self.initial.moisture)}.")

        if len(self.initial.type) != self.grid.nz:
            raise ValueError(f"Namelist nlevs={self.grid.nz} does not match "
                             f"init file soil_type length of {len(self.initial.type)}.")

        self.logger.info("Physical consistency checks passed")
