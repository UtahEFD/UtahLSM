import json
import time
import netCDF4 as nc
import numpy as np
from util import constants

class Input(object):

	def __init__(self, namelist, inputfile):
				
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
				self.soil_type = inifile.variables['soil_type'][:]
			# report a netcdf dictionary error to user and exit program
			except (KeyError) as e:
				print("There was an issue accessing data from \'%s\'"%inputfile)
				print("Error: The key",e,"does not exist")
				sys.exit(1)

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