#!/usr/bin/env python
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

import argparse
import time
import utahlsm
    
# let's time this thing
t1 = time.time()

# get case from user
parser = argparse.ArgumentParser(description="Run a case with UtahLSM")
parser.add_argument("-c", "--case", dest='case', required=True,
                    action='store', type=str, help="Case name")
parser.add_argument("-o", "--output", dest='outfile', 
                    action='store', type=str, help="Output file name")
args = parser.parse_args()
case = args.case
outf = args.outfile

# file paths to input data
namelist    = '../cases/%s/lsm_namelist.json'%case
initfile    = '../cases/%s/lsm_init.nc'%case 
offlinefile = '../cases/%s/lsm_offline.nc'%case

# a nice welcome message
print("##############################################################")
print("#                                                            #")
print("#                     Welcome to UtahLSM                     #")
print("#   A land surface model created at the University of Utah   #")
print("#       and the NOAA National Severe Storms Laboratory       #")
print("#                                                            #")
print("##############################################################")

# create Input instance
try:
    input_lsm = utahlsm.Input(namelist,initfile,offlinefile)
except:
    raise SystemExit(1)

# create Output instance
if not outf:
    outf='lsm_%s_py.nc'%case
try:
    output_lsm = utahlsm.Output(outf)
except:
    raise SystemExit(1)

# Create lsm object from input and object
lsm = utahlsm.UtahLSM(input_lsm,output_lsm)

# Loop through each time
runtime = 0
tstep   = input_lsm.forcing.tstep
for step_count, atm_state in enumerate(input_lsm.forcing.atmos):
    runtime += tstep
    
    # update user-specified fields
    lsm.update(tstep, runtime, atm_state)
    
    # run the model
    lsm.run(step_count, runtime)
    
    # save the data
    lsm.save(step_count,runtime)

# time info
t2 = time.time()
tt = t2 - t1
print("Done! Completed in %0.4f seconds"%tt)
print("##############################################################")