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
import time
import netCDF4 as nc
from . import logging_helper

class Output(object):
            
    def __init__(self,outfile):
        # create output file
        self.logger = logging_helper.get_logger("Output")
        self.logger.info(f"Saving output to {outfile}")
        self.outfile             = nc.Dataset(outfile,'w')
        # self.outfile.description = "UtahLSM output"
        # self.outfile.source      = "Jeremy A. Gibbs"
        # self.outfile.history     = "Created " + time.ctime(time.time())
    
        # dictionary of fields to be saved
        self.fields_time   = {}
        self.fields_static = {}
    
        self.attributes = {
            'time': {
                'dimension':("t",),
                'long_name':'time',
                'units':'s'
            },
            'soil_z': {
                'dimension':("z",),
                'long_name':'z-distance',
                'units':'m'
            },
            'soil_type': {
                'dimension':("z",),
                'long_name':'soil type',
                'units':''
            },
            'soil_T': {
                'dimension':("t","z",),
                'long_name':'soil temperature',
                'units':'K'
            },
            'soil_q': {
                'dimension':("t","z",),
                'long_name':'soil moisture',
                'units':'m3 m-3'
            },
            'ust': {
                'dimension':("t",),
                'long_name':'friction velocity',
                'units':'m s-1'
            },
            'obl': {
                'dimension':("t",),
                'long_name':'Obukhov length',
                'units':'m'
            },
            'shf': {
                'dimension':("t",),
                'long_name':'sensible heat flux',
                'units':'W m-2'
            },
            'lhf': {
                'dimension':("t",),
                'long_name':'latent heat flux',
                'units':'W m-2'
            },
            'ghf': {
                'dimension':("t",),
                'long_name':'ground heat flux',
                'units':'W m-2'
            },
        }
    
    # function to set the dimensions of each variable
    def set_dims(self,dims):
    
        # iterate through keys in dictionary
        for dim in dims:
            size = dims[dim]
            if size==0:
                self.outfile.createDimension(dim)
            else:
                self.outfile.createDimension(dim,size)
    
    # function to create desired output fields
    def set_fields(self,fields):
        
        # add time manually
        dims  = self.attributes['time']['dimension']
        name  = self.attributes['time']['long_name']
        units = self.attributes['time']['units']
        ncvar = self.outfile.createVariable('time', "f8", dims)
        ncvar.units              = units
        ncvar.long_name          = name

        self.fields_time['time'] = ncvar
    
        # iterate through keys in dictionary
        for field in fields:
            dims  = self.attributes[field]['dimension']
            units = self.attributes[field]['units']
            name  = self.attributes[field]['long_name']
            ncvar = self.outfile.createVariable(field, "f8", dims)
            ncvar.units        = units
            ncvar.long_name    = name
            if 't' in dims:
                self.fields_time[field] = ncvar
            else:
                self.fields_static[field] = ncvar
    
    # function to save data to the output file
    def save(self,fields,tidx,time,initial=False):
    
        # save static only for initial time
        if initial:
            for field in self.fields_static:
                self.fields_static[field][:] = fields[field]
        
        # save time
        for field in self.fields_time:
            
            dim = self.attributes[field]['dimension']
    
            if len(dim)==1:
                if field=='time':
                    self.fields_time[field][tidx] = time
                else:
                    self.fields_time[field][tidx] = fields[field]
            else:
                self.fields_time[field][tidx,:] = fields[field]
    
        # sync
        self.outfile.sync()
    
    # function to close the output file
    def close(self):
        self.outfile.close()
    