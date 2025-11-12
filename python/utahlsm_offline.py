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
from typing import Optional
import argparse
import time
import utahlsm
from utahlsm.exceptions import UtahLSMError

def main() -> None:
    """Parses arguments, runs the simulation, and prints timing information."""
    
    # Start a timer for the simulation
    t1: float = time.time()

    # Set up command-line argument parsing
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description="Run a case with UtahLSM")
    parser.add_argument("-c", "--case", dest='case', required=True,
                        action='store', type=str, help="Case name")
    parser.add_argument("-o", "--output", dest='outfile',
                        action='store', type=str, help="Output file name")
    args: argparse.Namespace = parser.parse_args()
    case: str = args.case
    outf: Optional[str] = args.outfile

    # Define file paths based on the case name
    namelist: str = f'../cases/{case}/lsm_namelist.json'
    initfile: str = f'../cases/{case}/lsm_init.nc'
    offlinefile: str = f'../cases/{case}/lsm_offline.nc'

    # Display a welcome message
    print("##############################################################")
    print("#                                                            #")
    print("#                     Welcome to UtahLSM                     #")
    print("#   A land surface model created at the University of Utah   #")
    print("#       and the NOAA National Severe Storms Laboratory       #")
    print("#                                                            #")
    print("##############################################################")

    output_lsm: Optional[utahlsm.Output] = None
    try:
        # Create input and output objects
        input_lsm: utahlsm.Input = utahlsm.Input(namelist, initfile, offlinefile)
        if not outf:
            outf = f'lsm_{case}_py.nc'
        output_lsm = utahlsm.Output(outf)

        # Create the main LSM object
        lsm: utahlsm.UtahLSM = utahlsm.UtahLSM(input_lsm, output_lsm)

        # --- Main Time-Stepping Loop ---
        runtime: float = 0
        tstep: float = input_lsm.forcing.tstep
        for step_count, atm_state in enumerate(input_lsm.forcing.atmos):
            runtime += tstep

            # Update the model with the latest atmospheric forcing
            lsm.update(tstep, runtime, atm_state)

            # Run the core physics solvers
            lsm.run(step_count, runtime)

            # Save the output for the current time step
            lsm.save(step_count, runtime)
    except (UtahLSMError) as e:
        print(f"\nAn error occurred: {e} Check namelist settings.")
        print("UtahLSM simulation failed.")
        raise SystemExit(1)
    except (FileNotFoundError, Exception) as e:
        print(f"\nAn error occurred: {e}")
        print("UtahLSM simulation failed.")
        raise SystemExit(1)
    finally:
        # Ensure the output file is properly closed, even if an error occurred
        if output_lsm is not None:
            output_lsm.close()
    
    # Calculate and print the total runtime
    t2: float = time.time()
    tt: float = t2 - t1
    print(f"Done! Completed in {tt:0.4f} seconds")
    print("##############################################################")

if __name__ == "__main__":
    main()