#!/usr/bin/env python
#
# UtahLSM
#
# Copyright (c) 2017–2026 Jeremy A. Gibbs
# Copyright (c) 2017–2026 Rob Stoll
# Copyright (c) 2017–2026 Eric Pardyjak
# Copyright (c) 2017–2026 Pete Willemsen
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
import sys
import time
from pathlib import Path
from typing import Callable, Optional, TextIO

import numpy as np

import utahlsm
from utahlsm.exceptions import SolverError, UtahLSMError


class _ProgressReporter:
    """Reports model progress without adding an external dependency."""

    # ``[`` + 28 cells + ``]`` matches the 30-character timestamp/level
    # prefix used by the logger.  The following dotted 10-character phase
    # field likewise aligns its colon with logger-name colons.
    _BAR_WIDTH = 28
    _LOG_MILESTONE = 0.25

    def __init__(
        self,
        total: int,
        *,
        enabled: bool = True,
        stream: Optional[TextIO] = None,
        min_interval: float = 0.25,
        clock: Callable[[], float] = time.perf_counter,
    ) -> None:
        """Initializes a terminal bar or sparse redirected-log reporter."""
        self.total = max(int(total), 0)
        self.stream = stream if stream is not None else sys.stdout
        self.enabled = enabled and self.total > 0
        self.min_interval = max(float(min_interval), 0.0)
        self._clock = clock
        self._start = self._clock()
        self._last_render = self._start
        self._next_log_fraction = self._LOG_MILESTONE
        self._last_line_width = 0
        self._started = False
        self._is_tty = bool(
            getattr(self.stream, 'isatty', lambda: False)()
        )

    @staticmethod
    def _format_duration(seconds: float) -> str:
        """Formats a nonnegative duration as MM:SS or H:MM:SS."""
        total_seconds = max(int(seconds), 0)
        minutes, secs = divmod(total_seconds, 60)
        hours, minutes = divmod(minutes, 60)
        if hours > 0:
            return f'{hours:d}:{minutes:02d}:{secs:02d}'
        return f'{minutes:02d}:{secs:02d}'

    def update(self, completed: int, *, phase: str = 'Running') -> None:
        """Updates progress, throttling terminals and sparsifying log output."""
        if not self.enabled:
            return

        completed = min(max(int(completed), 0), self.total)
        fraction = completed / self.total
        now = self._clock()
        is_complete = completed >= self.total

        if self._is_tty:
            if (
                self._started
                and not is_complete
                and now - self._last_render < self.min_interval
            ):
                return
        elif self._started and not is_complete:
            if fraction < self._next_log_fraction:
                return
            while self._next_log_fraction <= fraction:
                self._next_log_fraction += self._LOG_MILESTONE

        elapsed = max(now - self._start, 0.0)
        eta = (
            elapsed * (1.0 - fraction) / fraction
            if fraction > 0.0
            else None
        )
        eta_text = self._format_duration(eta) if eta is not None else '--:--'
        metrics = (
            f'{100.0 * fraction:5.1f}% '
            f'({completed}/{self.total}) '
            f'elapsed {self._format_duration(elapsed)} ETA {eta_text}'
        )

        if self._is_tty:
            filled = int(self._BAR_WIDTH * fraction)
            bar = '*' * filled + '-' * (self._BAR_WIDTH - filled)
            line = f'[{bar}] {phase:.>10s}: {metrics}'
            padded = line.ljust(self._last_line_width)
            self.stream.write(f'\r{padded}')
            self._last_line_width = len(line)
            if is_complete:
                self.stream.write('\n')
        else:
            self.stream.write(f'Progress {phase}: {metrics}\n')
        self.stream.flush()
        self._last_render = now
        self._started = True


def main() -> None:
    """Parses arguments, runs the simulation, and prints timing information.

    Raises:
        ValueError: If the case path contains path traversal attempts.
        UtahLSMError: If an error occurs during model setup or execution.
        FileNotFoundError: If required input files (namelist, init file) cannot
            be found.
    """
    # Start a timer for the simulation
    t1: float = time.perf_counter()

    # Set up command-line argument parsing
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Run a case with UtahLSM")
    parser.add_argument("-c", "--case", dest="case", required=True,
                        action="store", type=str, help="Case name")
    parser.add_argument("-o", "--output", dest="outfile", action="store",
                        type=str, help="Output file name")
    parser.add_argument("-s", "--spinup", dest="spinup", action="store",
                        type=int, default=0,
                        help="Number of diurnal cycles to spin up the soil "
                             "temperature before the scored run (precip off, "
                             "moisture held at the initial condition).")
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable the terminal progress bar and redirected-log milestones.",
    )
    args: argparse.Namespace = parser.parse_args()
    case: str = args.case
    outf: Optional[str] = args.outfile
    spinup: int = max(0, args.spinup)

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
    lsm: Optional[utahlsm.UtahLSM] = None
    step_count: int = 0
    runtime: float = 0.0
    try:
        # Create input and output objects
        input_lsm: utahlsm.Input = utahlsm.Input(
            namelist, initfile, offlinefile)
        if not outf:
            outf = f'lsm_{case}_py.nc'
        output_lsm = utahlsm.Output(outf, enabled=input_lsm.output.save)

        # Create the main LSM object
        lsm = utahlsm.UtahLSM(input_lsm, output_lsm)

        # --- Main Time-Stepping Loop ---
        assert input_lsm.forcing is not None
        tstep: float = input_lsm.forcing.tstep
        atmos = input_lsm.forcing.atmos
        ntime: int = len(atmos)
        step_count = 0
        run_steps = max(ntime - 1, 0)
        total_steps = (spinup + 1) * run_steps
        completed_steps = 0
        progress = _ProgressReporter(
            total_steps, enabled=not args.no_progress
        )
        initial_phase = f'Spin-up 1/{spinup}' if spinup > 0 else 'Running'
        progress.update(0, phase=initial_phase)

        # Optional soil-temperature spin-up. The deep soil equilibrates to the
        # mean surface forcing, but the initial condition is a single-time
        # snapshot, so a cold (or non-equilibrium) start drifts over the scored
        # day. Cycle the diurnal forcing `spinup` times to settle the soil
        # temperature, with precipitation disabled and soil moisture held at its
        # initial condition (the episodic rain cannot be spun up, and we do not
        # want spin-up to alter the moisture profile). Only the soil temperature
        # carries into the scored run; record 0 is then rewritten to reflect the
        # spun-up state so the saved series starts consistently.
        if spinup > 0:
            print(f'Spinning up soil temperature: {spinup} diurnal cycle(s)...')
            moisture_ic = np.array(lsm.soil_state.moisture, copy=True)
            for cycle in range(spinup):
                for k in range(1, ntime):
                    lsm.update(tstep, k * tstep, atmos[k])
                    np.asarray(lsm.atm_state.precipitation)[...] = 0.0
                    lsm.run()
                    np.asarray(lsm.soil_state.moisture)[...] = moisture_ic
                    completed_steps += 1
                    progress.update(
                        completed_steps,
                        phase=f'Spin-up {cycle + 1}/{spinup}',
                    )
            lsm.save(0, 0.0)

        # Record 0 (t=0) was already written during model setup using
        # forcing[0] as the consistent initial state. The forcing series has
        # ntime inclusive samples at times t_k = k*tstep, so there are
        # ntime-1 intervals to integrate. We use the right-endpoint
        # convention: the step landing at t_k is driven by forcing[k], and
        # record k reports forcing[k] at time t_k. This keeps the output
        # record, forcing sample, and time label aligned (no overshoot past
        # the forcing window, no double-counting of forcing diagnostics such
        # as precip). forcing[0] therefore serves only as the t=0 state.
        for step_count in range(1, ntime):
            runtime = step_count * tstep

            # Load the atmospheric forcing valid at this output time.
            lsm.update(tstep, runtime, atmos[step_count])

            # Run the core physics solvers, advancing the state to t = runtime.
            lsm.run()

            # Save the state and forcing under their shared time label.
            lsm.save(step_count, runtime)
            completed_steps += 1
            progress.update(completed_steps)
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
                fields = getattr(lsm, 'full_output_fields', lsm.output_fields)
                crash_out.set_fields(fields)
                crash_out.save(fields, step_count, runtime)
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
    t2: float = time.perf_counter()
    tt: float = t2 - t1
    print(f'Done! Completed in {tt:0.4f} seconds')
    print('##############################################################')

if __name__ == "__main__":
    main()
