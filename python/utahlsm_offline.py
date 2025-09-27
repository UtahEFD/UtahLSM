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
"""Offline driver for the Utah Land-Surface Model (UtahLSM).

This script serves as the primary entry point for running standalone, offline
simulations with UtahLSM. It handles command-line argument parsing for specifying
a simulation case, sets up the necessary input and output files, initializes the
LSM, and executes the main time-stepping loop.

To run an offline simulation, provide the case name via the command line:
    $ python utahlsm_offline.py -c my_case_name
"""
import argparse
import time
import utahlsm

def main():
    """Parses arguments, runs the simulation, and prints timing information."""
    # Start a timer for the simulation
    t1 = time.time()

    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Run a case with UtahLSM")
    parser.add_argument("-c", "--case", dest='case', required=True,
                        action='store', type=str, help="Case name")
    parser.add_argument("-o", "--output", dest='outfile',
                        action='store', type=str, help="Output file name")
    args = parser.parse_args()
    case = args.case
    outf = args.outfile

    # Define file paths based on the case name
    namelist    = f'../cases/{case}/lsm_namelist.json'
    initfile    = f'../cases/{case}/lsm_init.nc'
    offlinefile = f'../cases/{case}/lsm_offline.nc'

    # Display a welcome message
    print("##############################################################")
    print("#                                                            #")
    print("#                     Welcome to UtahLSM                     #")
    print("#   A land surface model created at the University of Utah   #")
    print("#       and the NOAA National Severe Storms Laboratory       #")
    print("#                                                            #")
    print("##############################################################")

    # Create Input and Output instances
    try:
        input_lsm = utahlsm.Input(namelist, initfile, offlinefile)
        if not outf:
            outf = f'lsm_{case}_py.nc'
        output_lsm = utahlsm.Output(outf)
    except Exception as e:
        print(f"Error during initialization: {e}")
        raise SystemExit(1)

    # Create the main LSM object
    lsm = utahlsm.UtahLSM(input_lsm, output_lsm)

    # --- Main Time-Stepping Loop ---
    runtime = 0
    tstep   = input_lsm.forcing.tstep
    for step_count, atm_state in enumerate(input_lsm.forcing.atmos):
        runtime += tstep

        # Update the model with the latest atmospheric forcing
        lsm.update(tstep, runtime, atm_state)

        # Run the core physics solvers
        lsm.run(step_count, runtime)

        # Save the output for the current time step
        lsm.save(step_count, runtime)

    # Calculate and print the total runtime
    t2 = time.time()
    tt = t2 - t1
    print(f"Done! Completed in {tt:0.4f} seconds")
    print("##############################################################")

if __name__ == "__main__":
    main()