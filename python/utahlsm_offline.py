#!/usr/bin/env python
import argparse
import sys
import time
from util import io, matrix

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
	parser.add_argument("-c", "--case", dest='case', 
						action='store', type=str, help="Case name")
	parser.add_argument("-o", "--output", dest='outfile', 
						action='store', type=str, help="Output file name")
	args = parser.parse_args()
	case = args.case
	outf = args.outfile
	
	# create Input instance
	namelist = 'cases/%s/namelist.json'%case
	initfile = 'cases/%s/scm_init.nc'%case  
	inputSCM = io.Input(namelist,initfile)
	
	# create Output instance
	if not outf:
		outf='utahlsm.nc'
	outputSCM = io.Output(outf)
	
	# create SCM instance
	#scm = SingleColumnModel(inputSCM,outputSCM)
	
	# run the model
	#scm.run()
	
	# time info
	t2 = time.time()
	tt = t2 - t1
	print("\n[UtahLSM: Run] \t Done! Completed in %0.2f seconds"%tt)
	print("##############################################################")