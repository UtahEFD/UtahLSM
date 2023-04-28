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
			self.utc_start  = data["time"]["utc_start"]
			self.julian_day = data["time"]["julian_day"]
			
			# grid section
			self.nx         = data["grid"]["nx"]
			self.ny         = data["grid"]["ny"]
			
			# length section
			self.z_o        = data["surface"]["z_o"]
			self.z_t        = data["surface"]["z_t"]
			self.z_m        = data["surface"]["z_m"]
			self.z_s        = data["surface"]["z_s"]
			self.albedo     = data["surface"]["albedo"]
			self.emissivity = data["surface"]["emissivity"]
			
			# soil section
			self.nsoil      = data["soil"]["nsoil"]
			self.param      = data["soil"]["param"]
			self.model      = data["soil"]["model"]
			
			# radiation section
			self.comp_rad   =  data["radiation"]["comp_rad"]
			self.latitude   =  data["radiation"]["latitude"]
			self.longitude  =  data["radiation"]["longitude"]
			
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
				self.soil_z    = inifile.variables['soil_z'][:].astype('float')
				self.soil_T    = inifile.variables['soil_T'][:].astype('float')
				self.soil_q    = inifile.variables['soil_q'][:].astype('float')
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
					self.atm_U = metfile.variables['atm_U'][:].astype('float')
					self.atm_T = metfile.variables['atm_T'][:].astype('float')
					self.atm_q = metfile.variables['atm_q'][:].astype('float')
					self.atm_p = metfile.variables['atm_p'][:].astype('float')
					self.r_net = metfile.variables['R_net'][:].astype('float')
				# report a netcdf dictionary error to user and exit program
				except (KeyError) as e:
					print("There was an issue accessing data from \'%s\'"%inputfile)
					print("Error: The key",e,"does not exist")
					sys.exit(1)

class Output(object):
			
	def __init__(self,outfile):
		
		# create output file
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
	
class Logger(object):
	
	def print_number(x,label=False):
		
		if label:
			print()
			print('LOG -> %s: %0.17g'%(label,x))
		else:
			print()
			print('LOG -> %0.17g'%x)