# 
# UtahLSM
# 
# Copyright (c) 2017–2023 Jeremy A. Gibbs
# Copyright (c) 2017–2023 Rob Stoll
# Copyright (c) 2017–2023 Eric Pardyjak
# Copyright (c) 2017–2023 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

import json
import time
import netCDF4 as nc
import numpy as np
from util import constants

class Input(object):

	def __init__(self, namelist, inputfile, offlinefile=None):
				
		# read the namelist
		try:
			with open(namelist) as json_file:
				data = json.load(json_file)
		except FileNotFoundError as e:
			print('Error: %s — %s'%(namelist,e.strerror))
		except json.decoder.JSONDecodeError as e:
			print('Error parsing %s: %s (see line %s)'%(namelist,e.msg,e.lineno))
		else:
			# time section
			self.step_seb   = data["time"]["step_seb"]
			self.step_dif   = data["time"]["step_dif"]
			
			# grid section
			self.nx         = data["grid"]["nx"]
			self.ny         = data["grid"]["ny"]
			
			# length section
			self.z_o        = data["length"]["z_o"]
			self.z_t        = data["length"]["z_t"]
			self.z_m        = data["length"]["z_m"]
			self.z_s        = data["length"]["z_s"]
			
			# soil section
			self.nsoil      = data["soil"]["nsoil"]
			self.param      = data["soil"]["param"]
			self.model      = data["soil"]["model"]
			
			# radiation section
			self.utc_start  =  data["radiation"]["utc_start"]
			self.comp_rad   =  data["radiation"]["comp_rad"]
			self.albedo     =  data["radiation"]["albedo"]
			self.emissivity =  data["radiation"]["emissivity"]
			self.latitude   =  data["radiation"]["latitude"]
			self.longitude  =  data["radiation"]["longitude"]
			self.julian_day =  data["radiation"]["julian_day"]
			
			# output section
			self.save       = data["output"]["save"]
			self.fields     = data["output"]["fields"]
		
		# open and parse the netcdf initialization data
		try:
			inifile = nc.Dataset(inputfile)
		# report file open error to user and exit program
		except (RuntimeError,FileNotFoundError) as e:
			print('There was an issue opening \'%s\'.'%(inputfile))
			print('Error: ',e.strerror)
			sys.exit(1)  
		# process the netcdf input file
		else:
			# load initial data from netcdf into local variables
			try:
				self.soil_z    = inifile.variables['soil_z'][:]
				self.soil_T    = inifile.variables['soil_T'][:]
				self.soil_q    = inifile.variables['soil_q'][:]
				self.soil_type = inifile.variables['soil_type'][:].astype('int')
			# report a netcdf dictionary error to user and exit program
			except (KeyError) as e:
				print("There was an issue accessing data from \'%s\'"%inputfile)
				print("Error: The key",e,"does not exist")
				sys.exit(1)
		
		# open and parse the netcdf offline data if available
		if (offlinefile):
			try:
				metfile = nc.Dataset(offlinefile)
			# report file open error to user and exit program
			except (RuntimeError,FileNotFoundError) as e:
				print('There was an issue opening \'%s\'.'%(offlinefile))
				print('Error: ',e.strerror)
				sys.exit(1)  
			# process the netcdf offline file
			else:
				# load offline data from netcdf into local variables
				try:
					self.ntime = len(metfile.dimensions['t'])
					self.tstep = metfile.variables['tstep'][:]
					self.atm_U = metfile.variables['atm_U'][:]
					self.atm_T = metfile.variables['atm_T'][:]
					self.atm_q = metfile.variables['atm_q'][:]
					self.atm_p = metfile.variables['atm_p'][:]
					self.r_net = metfile.variables['R_net'][:]
				# report a netcdf dictionary error to user and exit program
				except (KeyError) as e:
					print("There was an issue accessing data from \'%s\'"%inputfile)
					print("Error: The key",e,"does not exist")
					sys.exit(1)

class Output(object):
			
	def __init__(self,outfile):
		
		# create output file
		self.outfile             = nc.Dataset(outfile,'w')
		self.outfile.description = "UtahLSM output"
		self.outfile.source      = "Jeremy A. Gibbs"
		self.outfile.history     = "Created " + time.ctime(time.time())
	
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
		ncvar = self.outfile.createVariable('time', "f4", dims)
		ncvar.long_name          = name
		ncvar.units              = units
		self.fields_time['time'] = ncvar
	
		# iterate through keys in dictionary
		for field in fields:
			dims  = self.attributes[field]['dimension']
			name  = self.attributes[field]['long_name']
			units = self.attributes[field]['units']
			ncvar = self.outfile.createVariable(field, "f4", dims)
			ncvar.long_name    = name
			ncvar.units        = units
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