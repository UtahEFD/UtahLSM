import json
import time
import netCDF4 as nc
import numpy as np
from util import constants

class Input(object):

	def __init__(self, namelist, inputfile):
		self.keys = []
		
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
			self.dt        = data["time"]["dt"]
			self.tt        = data["time"]["tt"]
			#self.adaptive  = data["time"]["adaptive"]
			
			# grid section
			self.nz        = data["grid"]["nz"]
			self.lz        = data["grid"]["lz"]
			self.dz        = self.lz/(self.nz-1)
			
			# large-scale forcing section
			self.fc        = data["force"]["fc"]
			self.wls       = data["force"]["wls"]
			# TODO: check if this flag exists
			#       if so, skip next part and extract keys
			#       then check if ug,vg,wls exist
			#self.ls_time   = data["force"]["timedep"]
			# TODO: add option for ug,vg to be time varying
			# TODO: add in w option for large-scale subsidence
			self.ug        = data["force"]["ug"]
			self.vg        = data["force"]["vg"]
			
			# thermodynamic section
			self.Tref      = data["thermo"]["Tref"]
			self.lapse     = data["thermo"]["lapse"]
			self.beta      = constants.grav / self.Tref
			
			# surface section
			self.z0        = data["surface"]["z0"]
			self.z0h       = data["surface"]["z0h"]
			self.sfc_model = data["surface"]["model"]
			self.n_patch   = data["surface"]["n_patch"]
			self.Tdiff     = data["surface"]["Tdiff"]
			#self.sfc_time  = data["surface"]["timedep"]
			
			# landuse section
			self.l1        = data["landuse"]["0"]
			self.l2        = data["landuse"]["1"]
			self.fa        = list(data['landuse'].values())
			for key in data['landuse'].keys(): 
				self.keys.append(int(key)) 
			
			# pbl section
			self.pbl_model = data["pbl"]["model"]
			
			# output section
			self.t_save    = data["output"]["t_save"]
			self.n_save    = np.floor(self.t_save/self.dt)   
			
			# initialization data
			inifile = nc.Dataset(inputfile)
			self.zf = inifile.variables['zf'][:]
			self.zh = inifile.variables['zh'][:]
			self.u  = inifile.variables['u'][:]
			self.v  = inifile.variables['v'][:]
			self.T  = inifile.variables['T'][:]
			# TODO: add ug,vg,wls, and Tsfc if not constant
			
			# TODO: put a check of namelist here
			#if self.tt > 15:
				#self.wls= inifile.variables['wls'][:]
			if self.sfc_model!=1:
				#self.Tse = inifile.variables['Tse'][:]
				self.Ts1 = inifile.variables['Ts1'][:]
				self.Ts2 = inifile.variables['Ts2'][:]
			else:
				self.Ts = inifile.variables['Ts'][:]
			
			# TODO: put a try in here and handle missing tke
			if self.pbl_model==2 or self.pbl_model==3:
				self.tke = inifile.variables['tke'][:]

class Output(object):
			
	def __init__(self,outfile):
		
		# create output file
		self.outfile             = nc.Dataset(outfile,'w')
		self.outfile.description = "NSSL-SCM output"
		self.outfile.source      = "Dominic Candela and Jeremy A. Gibbs"
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
			'zf': {
				'dimension':("zf",),
				'long_name':'z-distance on full levels',
				'units':'m'
			},
			'zh': {
				'dimension':("zh",),
				'long_name':'z-distance on half levels',
				'units':'m'
			},
			'zi': {
				'dimension':("t",),
				'long_name':'PBL depth',
				'units':'m'
			},
			'l': {
				'dimension':("t","zf",),
				'long_name':'PBL length scale',
				'units':'m'
			},
			'L': {
				'dimension':("t",),
				'long_name': 'Obukhov Length',
				'units': 'm'
			},
			'u': {
				'dimension':("t","zh",),
				'long_name':'u-component velocity',
				'units':'m s-1'
			},
			'v': {
				'dimension':("t","zh",),
				'long_name':'v-component velocity',
				'units':'m s-1'
			},
			'T': {
				'dimension':("t","zh",),
				'long_name':'temperature',
				'units':'K'
			},
			'us': {
				'dimension':("t",),
				'long_name':'friction velocity',
				'units':'m s-1'
			},
			'wu': {
				'dimension':("t","zf",),
				'long_name':'kinematic momentum flux',
				'units':'m2 s-2'
			},
			'wT': {
				'dimension':("t","zf",),
				'long_name':'kinematic heat flux',
				'units':'m s-1 K'
			},
			'tke': {
				'dimension':("t","zf",),
				'long_name':'turbulence kinetic energy',
				'units':'m2 s-2'
			},
			'Km': {
				'dimension':("t","zf",),
				'long_name':'eddy diffusivity momentum',
				'units':'m2 s-1'
			},
			'Kh': {
				'dimension':("t","zf",),
				'long_name':'eddy diffusivity heat',
				'units':'m2 s-1'
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