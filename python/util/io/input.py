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

import json
import jsonschema
import logging
import os
import time
import netCDF4 as nc
import numpy as np

# local logger
logger = logging.getLogger("IO: Input")

class Input(object):

    def __init__(self, namelist, inputfile, offlinefile=None):
        
        # validation schema file path
        schema = "util/io/schema_namelist.json"
        
        # validation schema file
        try:
            with open(schema) as json_file:
                namelist_schema = json.load(json_file)
        except FileNotFoundError as e:
            logger.error('Error: %s — %s'%(schema,e.strerror))
        except json.decoder.JSONDecodeError as e:
            logger.error('Error parsing %s: %s (see line %s)'%(schema,e.msg,e.lineno))
        
        # namelist json files
        try:
            with open(namelist) as json_file:
                namelist_data = json.load(json_file)
        except FileNotFoundError as e:
            logger.error('Error: %s — %s'%(namelist,e.strerror))
        except json.decoder.JSONDecodeError as e:
            logger.error('Error parsing %s: %s (see line %s)'%(namelist,e.msg,e.lineno))
        else:
            # validate the data against the schema
            try:
                jsonschema.validate(instance=namelist_data, schema=namelist_schema)
            except jsonschema.ValidationError as e:
                logger.error("Namelist validation failed!")
                logger.error(f"Error: {e.message}")
                logger.error(f"Path to error: {list(e.path)}")
            else:
                # time section
                self.step_seb   = namelist_data["time"]["step_seb"]
                self.step_dif   = namelist_data["time"]["step_dif"]
                self.utc_start  = namelist_data["time"]["utc_start"]
                self.julian_day = namelist_data["time"]["julian_day"]
                
                # grid section
                self.nx         = namelist_data["grid"]["nx"]
                self.ny         = namelist_data["grid"]["ny"]
                
                # length section
                self.z_o        = namelist_data["surface"]["z_o"]
                self.z_t        = namelist_data["surface"]["z_t"]
                self.z_m        = namelist_data["surface"]["z_m"]
                self.z_s        = namelist_data["surface"]["z_s"]
                self.albedo     = namelist_data["surface"]["albedo"]
                self.emissivity = namelist_data["surface"]["emissivity"]
                self.sfc_model  = namelist_data["surface"]["model"]
                
                # soil section
                self.nsoil      = namelist_data["soil"]["nsoil"]
                self.soil_param = namelist_data["soil"]["param"]
                self.soil_model = namelist_data["soil"]["model"]
                
                # radiation section
                self.rad_model  = namelist_data["radiation"]["model"]
                self.latitude   = namelist_data["radiation"]["latitude"]
                self.longitude  = namelist_data["radiation"]["longitude"]
                
                # output section
                self.save       = namelist_data["output"]["save"]
                self.fields     = namelist_data["output"]["fields"]
        
        # open and parse the netcdf initialization data
        try:
            inifile = nc.Dataset(inputfile)
        # report file open error to user and exit program
        except (RuntimeError,FileNotFoundError) as e:
            logger.error('There was an issue opening \'%s\'.'%(inputfile))
            logger.error('Error: ',e.strerror)
            raise SystemExit(1)  
        # process the netcdf input file
        else:
            # load initial data from netcdf into local variables
            try:
                self.soil_z    = inifile.variables['soil_z'][:].astype('float')
                self.soil_T    = inifile.variables['soil_T'][:].astype('float')
                self.soil_q    = inifile.variables['soil_q'][:].astype('float')
                self.soil_type = inifile.variables['soil_type'][:].astype('int')
            # report a netcdf dictionary error to user and exit program
            except (KeyError) as e:
                logger.error("There was an issue accessing data from \'%s\'"%inputfile)
                logger.error("Error: The key",e,"does not exist")
                raise SystemExit(1)
        
        # open and parse the netcdf offline data if available
        if (offlinefile):
            try:
                metfile = nc.Dataset(offlinefile)
            # report file open error to user and exit program
            except (RuntimeError,FileNotFoundError) as e:
                logger.error('There was an issue opening \'%s\'.'%(offlinefile))
                logger.error('Error: ',e.strerror)
                raise SystemExit(1)  
            # process the netcdf offline file
            else:
                # load offline data from netcdf into local variables
                try:
                    metfile.set_auto_mask(False)
                    self.ntime = len(metfile.dimensions['t'])
                    self.tstep = metfile.variables['tstep'][0].astype('float')
                    self.atm_U = metfile.variables['atm_U'][:].astype('float')
                    self.atm_T = metfile.variables['atm_T'][:].astype('float')
                    self.atm_q = metfile.variables['atm_q'][:].astype('float')
                    self.atm_p = metfile.variables['atm_p'][:].astype('float')
                    self.r_net = metfile.variables['R_net'][:].astype('float')
                # report a netcdf dictionary error to user and exit program
                except (KeyError) as e:
                    logger.error("There was an issue accessing data from \'%s\'"%inputfile)
                    logger.error("Error: The key",e,"does not exist")
                    raise SystemExit(1)
    