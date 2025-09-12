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
import logging
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

# configure logging
log_format = '{asctime} [{levelname:^8s}] {name:^20s} {message}'
logging.basicConfig(level=logging.INFO,
                    format=log_format,
                    datefmt='%Y-%m-%d %H:%M:%S',
                    style='{',
                    filename='utahlsm.log',
                    filemode='w')    

# create a console handler for printing to the screen
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter(log_format, "%Y-%m-%d %H:%M:%S",style='{',))
logging.getLogger('').addHandler(console_handler)

# local logger
logger = logging.getLogger("UtahLSM")

# a nice welcome message
logger.info("##############################################################")
logger.info("#                                                            #")
logger.info("#                     Welcome to UtahLSM                     #")
logger.info("#   A land surface model created at the University of Utah   #")
logger.info("#       and the NOAA National Severe Storms Laboratory       #")
logger.info("#                                                            #")
logger.info("##############################################################")

# create Input instance
try:
    input_lsm = utahlsm.Input(namelist,initfile,offlinefile)
except:
    raise SystemExit(1)
logger.info("Running offline for the %s case"%case)

   # create Output instance
if not outf:
    outf='lsm_%s_py.nc'%case
output_lsm = utahlsm.Output(outf)

# grid information (not used yet)
nx = input_lsm.grid.nx
ny = input_lsm.grid.ny

# --- Main Time Loop ---
logger.info("Starting simulation time loop...")

# local fluxes to be modified by lsm
lsm = utahlsm.UtahLSM(input_lsm,output_lsm)

# Loop through each time
runtime = 0
tstep = input_lsm.forcing.tstep
for step_count, atm_state in enumerate(input_lsm.forcing.atmos):
    runtime += tstep
    logger.info(f"Running for time: {runtime:8.2f} of {input_lsm.forcing.ntime*tstep:8.2f}")
    
    # update user-specified fields
    lsm.update(tstep, runtime, atm_state)
    lsm.run(runtime)
    lsm.save(step_count,runtime)

# time info
t2 = time.time()
tt = t2 - t1
logger.info("Done! Completed in %0.4f seconds"%tt)
logger.info("##############################################################")