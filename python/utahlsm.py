#!/usr/bin/env python
import argparse
import os
import sys
import time
from util import io, matrix

# custom error message for user case entry
class InvalidCase(Exception):
	pass

# land-surface model class
class UtahLSM:
	"""This is the main UtahLSM class
	
	:param inputLSM: A handle to the :class:`util.io.Input`
	:param outputLSM: A handle to the :class:`util.io.Output`
	"""
	
	# model class initialization
	def __init__(self,inputLSM, outputLSM):
		"""Constructor method
		"""
		# set the input and output fields
		self.input  = inputLSM
		self.output = outputLSM
		
		print("[UtahLSM: Setup] \t Reading input settings")
		# configuration variables
		nz         = self.input.nsoil
		
		print("[UtahLSM: Setup] \t Reading input data")
		# input fields
		self.soil_z    = self.input.soil_z
		self.soil_T    = self.input.soil_T
		self.soil_q    = self.input.soil_q
		self.soil_type = self.input.soil_type
# 		
# 		print("[UtahLSM: Setup] \t Creating surface model")
# 		# choose surface model
# 		self.sfc = SFC.get_model(sfc_model,self.input)
# 	  
# 		print("[UtahLSM: Setup] \t Creating PBL model")
# 		# choose PBL model
# 		self.pbl = PBL.get_model(pbl_model,self.input,self.sfc)
# 
# 		print("[NSSL-SCM: Setup] \t Creating additional fields")
# 		# initialize flux arrays
# 		self.zi    = np.ones(1)
# 		self.obukL = np.zeros(1)
# 		self.ustar = np.zeros(1)
# 		self.tau13 = np.zeros(nz)
# 		self.tau23 = np.zeros(nz)
# 		self.Km    = np.zeros(nz)
# 		self.Kh    = np.zeros(nz)
# 		self.hflux = np.zeros(nz)
# 		self.mflux = np.sqrt(self.tau13**2 + self.tau23**2)
# 		
# 		# Runge-Kutta variables
# 		self.rk_a  = [0,-5/9,-153/128]
# 		self.rk_b  = [1/3, 15/16, 8/15]
# 		self.rk_qu = np.zeros(nz-1)
# 		self.rk_qv = np.zeros(nz-1)
# 		self.rk_qt = np.zeros(nz-1)

		print("[UtahLSM: Setup] \t Creating output file")    
		# set reference to output dimensions
		self.output_dims = {
			't':0,
			'z':nz
		}
		self.output.set_dims(self.output_dims)

		# set reference to output fields
		self.output_fields = {
			'soil_z':self.soil_z,
			'soil_T':self.soil_T,
			'soil_q':self.soil_q,
			'soil_type':self.soil_type
			# 'zh':self.zh[:-1],
			# 'zi':self.zi,
			# 'u':self.U[:-1],
			# 'v':self.V[:-1],
			# 'T':self.T[:-1],
			# 'l':self.pbl.ell,
			# 'L':self.obukL,
			# 'us':self.ustar,
			# 'wu':self.mflux,
			# 'wT':self.hflux,
			# 'Km':self.Km,
			# 'Kh':self.Kh
		}
		self.output.set_fields(self.output_fields)

		# write initial data
		self.output.save(self.output_fields,0,0,initial=True)

# main program to run the LSM
if __name__ == "__main__":
	
	# let's time this thing
	t1 = time.time()

	# a nice welcome message
	print("##############################################################")
	print("#                                                            #")
	print("#                     Welcome to UtahLSM                     #")
	print("#   A land surface model created at the University of Utah   #")
	print("#       and the NOAA National Severe Storms Laboratory       #")
	print("#                                                            #")
	print("##############################################################")
	
	# get case from user
	parser = argparse.ArgumentParser(description="Run a case with UtahLSM")
	parser.add_argument("-c", "--case", dest='case', required=True,
						action='store', type=str, help="Case name")
	parser.add_argument("-o", "--output", dest='outfile', 
						action='store', type=str, help="Output file name")
	args = parser.parse_args()
	case = args.case
	outf = args.outfile
	
	# create Input instance
	try:
		if os.path.exists('../cases/%s/'%case):
			namelist    = '../cases/%s/lsm_namelist.json'%case
			initfile    = '../cases/%s/lsm_init.nc'%case 
			offlinefile = '../cases/%s/lsm_offline.nc'%case  
			inputLSM    = io.Input(namelist,initfile,offlinefile)
		else:
			raise InvalidCase('Error: The folder ../cases/%s does not exist.'%case)
	except InvalidCase as e:
		print(e)
		sys.exit(1)

	# create Output instance
	if not outf:
		outf='lsm_%s_py.nc'%case
	outputLSM = io.Output(outf)
	
	# create UtahLSM instance
	lsm = UtahLSM(inputLSM,outputLSM)
	
	# run the model
	#scm.run()
	
	# time info
	t2 = time.time()
	tt = t2 - t1
	print("\n[UtahLSM: Run] \t Done! Completed in %0.2f seconds"%tt)
	print("##############################################################")