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
from dataclasses import dataclass
import json
import jsonschema
import logging
import os
import time
from typing import Dict, List, Optional

import netCDF4 as nc
import numpy as np

from ...data_models import (
    TimeConfig, GridConfig, SurfaceConfig, SoilConfig, RadiationConfig, OutputConfig,
    SoilState, ForcingData, AtmosphericState
)

# local logger
logger = logging.getLogger("IO: Input")

class Input(object):

    def __init__(self, namelist_path: str, inputfile: str, offlinefile: str = None):
        
        # Load and validate the namelist into the configuration dataclasses
        namelist_data = self._load_and_validate_namelist(namelist_path)
            
        # Load initial conditions into the InitialConditions dataclass
        init_data = self._load_initial_conditions(inputfile)
        
        # Assemble the final, structured dataclasses from the raw data
        self.time:      TimeConfig      = TimeConfig(**namelist_data["time"])
        self.surface:   SurfaceConfig   = SurfaceConfig(**namelist_data["surface"])
        self.soil:      SoilConfig      = SoilConfig(**namelist_data["soil"])
        self.radiation: RadiationConfig = RadiationConfig(**namelist_data["radiation"])
        self.output:    OutputConfig    = OutputConfig(**namelist_data["output"])
        self.grid:      GridConfig      = GridConfig(
                                            nx = namelist_data["grid"]["nx"],
                                            ny = namelist_data["grid"]["ny"],
                                            nz = namelist_data["grid"]["nz"],
                                            z  = init_data["z"]
                                        )
        # InitialConditions now only contains dynamic variables
        self.initial: SoilState = SoilState(
            T    = init_data["T"],
            q    = init_data["q"],
            type = init_data["type"]
        )
        
        # Load offline forcing data if provided into the ForcingData dataclass
        self.forcing: Optional[ForcingData] = None
        if offlinefile:
            self._load_offline_data(offlinefile)
            
        # Perform additional semantic and physical validation
        self._validate_physical_consistency()

    def _load_and_validate_namelist(self, namelist_path: str) -> Dict:
        """Loads the JSON namelist, validates it, and populates the config dataclasses."""
        
        schema_path = "utahlsm/util/io/schema_namelist.json"
        
        try:
            with open(schema_path) as f: 
                schema = json.load(f)
            with open(namelist_path) as f: 
                namelist_data = json.load(f)
            
            # validate namelist structure
            jsonschema.validate(instance=namelist_data, schema=schema)
            logger.info("Namelist validation successful")
            return namelist_data
        # raise an error
        except (FileNotFoundError, json.JSONDecodeError, jsonschema.ValidationError) as e:
            logger.error(f"Namelist Error: {e}")
            raise
    
    def _load_initial_conditions(self, inputfile: str)-> Dict[str, np.ndarray]:
        """Loads data from the NetCDF initialization file into a dataclass."""
        try:
            with nc.Dataset(inputfile) as inifile:
                inifile.set_auto_mask(False)
                init_dict = {
                    "z"    : (-1)*inifile.variables['soil_z'][:].astype('float'),
                    "T"    : inifile.variables['soil_T'][:].astype('float'),
                    "q"    : inifile.variables['soil_q'][:].astype('float'),
                    "type" : inifile.variables['soil_type'][:].astype('int')
                }
            logger.info("Initial conditions data loaded successfully")
            return init_dict
        except (IOError, KeyError) as e:
            logger.error(f"Initial Conditions Error: {e}")
            raise
    
    def _load_offline_data(self, offlinefile: str):
        """Loads data from the NetCDF offline forcing file into a dataclass."""
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
                
                atm_data = [
                    AtmosphericState(U=atm_U[i], T=atm_T[i], q=atm_q[i], p=atm_p[i], R_net=r_net[i])
                    for i in range(ntime)
                ]
                
                self.forcing = ForcingData(ntime=ntime, tstep=tstep, atmos=atm_data)
                logger.info(f"Loaded {ntime} timesteps of forcing data")
        except (IOError, KeyError) as e:
            logger.error(f"Offline forcing error: {e}")
            raise
    
    def _validate_physical_consistency(self):
        """Performs validation checks on inter-variable relationships."""
        if self.surface.z_m <= self.surface.z_o:
            raise ValueError(f"z_m={self.surface.z_m} must be > z_o={self.surface.z_o}.")
        
        if self.surface.z_s <= self.surface.z_t:
            raise ValueError(f"z_s={self.surface.z_s} must be > z_t={self.surface.z_t}.")
        
        if len(self.initial.T) != self.grid.nz:
            raise ValueError(f"Namelist nlevs={self.grid.nz} does not match "
                             f"init file soil_T length of {len(self.initial.T)}.")
        logger.info("Physical consistency checks passed")
