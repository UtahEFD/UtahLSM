#!/usr/bin/env python
import argparse
import os
import sys
import time
from util import io, matrix

# custom error message for user case entry
class InvalidCase(Exception):
	pass 

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
			namelist = '../cases/%s/lsm_namelist.json'%case
			initfile = '../cases/%s/lsm_init.nc'%case  
			inputLSM = io.Input(namelist,initfile)
		else:
			raise InvalidCase('Error: The folder ../cases/%s does not exist.'%case)
	except InvalidCase as e:
		print(e)
		sys.exit(1)
	
	# create Output instance
	if not outf:
		outf='../cases/%s/utahlsm_py.nc'%case
	outputSCM = io.Output(outf)
	
	# run the model
	#scm.run()
	
	# time info
	t2 = time.time()
	tt = t2 - t1
	print("\n[UtahLSM: Run] \t Done! Completed in %0.2f seconds"%tt)
	print("##############################################################")