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
simulations with UtahLSM. It handles command-line argument parsing for
specifying
a simulation case, sets up the necessary input and output files, initializes the
LSM, and executes the main time-stepping loop.

To run an offline simulation, provide the case name via the command line:
    $ python utahlsm_offline.py -c my_case_name
"""
import argparse
import time
from pathlib import Path
from typing import Optional

import utahlsm
from utahlsm.exceptions import SolverError, UtahLSMError


def main() -> None:
    """Parses arguments, runs the simulation, and prints timing information.

    Raises:
        ValueError: If the case path contains path traversal attempts.
        UtahLSMError: If an error occurs during model setup or execution.
        FileNotFoundError: If required input files (namelist, init file) cannot
            be found.
    """

    # Start a timer for the simulation
    t1: float = time.time()

    # Set up command-line argument parsing
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Run a case with UtahLSM")
    parser.add_argument("-c", "--case", dest="case", required=True,
                        action="store", type=str, help="Case name")
    parser.add_argument("-o", "--output", dest="outfile", action="store",
                        type=str, help="Output file name")
    args: argparse.Namespace = parser.parse_args()
    case: str = args.case
    outf: Optional[str] = args.outfile

    # Define file paths based on the case name with path traversal protection
    base_path: Path = Path("../cases").resolve()

    try:
        case_path: Path = (base_path / case).resolve()
        case_path.relative_to(base_path)
    except (ValueError, RuntimeError) as e:
        raise ValueError(
            f'Invalid case name {case}: path traversal detected. '
            f'Case directory must be within {base_path}.'
        ) from e

    namelist: str = str(case_path / "lsm_namelist.json")
    initfile: str = str(case_path / "lsm_init.nc")
    offlinefile: str = str(case_path / "lsm_offline.nc")

    # Display a welcome message
    print('##############################################################')
    print('#                                                            #')
    print('#                     Welcome to UtahLSM                     #')
    print('#   A land surface model created at the University of Utah   #')
    print('#       and the NOAA National Severe Storms Laboratory       #')
    print('#                                                            #')
    print('##############################################################')

    output_lsm: Optional[utahlsm.Output] = None
    try:
        # Create input and output objects
        input_lsm: utahlsm.Input = utahlsm.Input(
            namelist, initfile, offlinefile)
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
            lsm.run()

            # Save the output for the current time step
            lsm.save(step_count+1, runtime)
    except SolverError as e:
        print('\n!!! NUMERICAL SOLVER FAILURE !!!')
        print(f'Error details: {e}')
        if lsm is not None:
            crash_file = f'lsm_crash_{case}.nc'
            print(f'Attempting to save debug state to: {crash_file}')
            try:
                # Create a specialized output object for the crash dump
                crash_out = utahlsm.Output(crash_file)
                crash_out.set_dims(lsm.output_dims)
                crash_out.set_fields(lsm.output_fields)
                crash_out.save(lsm.output_fields, step_count, runtime)
                crash_out.close()
                print('>> Crash dump saved successfully.')
            except Exception as dump_e:
                print(f'>> Failed to save crash dump: {dump_e}')
        print('UtahLSM simulation failed.')
        raise SystemExit(1) from e
    except ValueError as e:
        print(f'\nConfiguration error: {e}')
        print('UtahLSM simulation failed.')
        raise SystemExit(1) from e
    except UtahLSMError as e:
        print(f'\nAn error occurred: {e}')
        print('Check namelist settings.')
        print('UtahLSM simulation failed.')
        raise SystemExit(1) from e
    except FileNotFoundError as e:
        print(f'\nFile not found: {e}')
        print('UtahLSM simulation failed.')
        raise SystemExit(1) from e
    finally:
        # Ensure the output file is properly closed, even if an error occurred
        if output_lsm is not None:
            output_lsm.close()

    # Calculate and print the total runtime
    t2: float = time.time()
    tt: float = t2 - t1
    print(f'Done! Completed in {tt:0.4f} seconds')
    print('##############################################################')

if __name__ == "__main__":
    main()
