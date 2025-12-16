#!/usr/bin/env python3
"""Compare UtahLSM output with GABLS3 observations.

This script creates a multi-panel comparison plot of key surface fluxes
and friction velocity between model output and observations.
"""

import argparse
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def load_observations(obs_file: str) -> dict:
    """Load GABLS3 observation data.

    Args:
        obs_file: Path to observation NetCDF file.

    Returns:
        Dictionary containing time and flux observations.
    """
    with nc.Dataset(obs_file, 'r') as ds:
        data = {
            'time': ds.variables['time'][:],
            'ust': ds.variables['UST'][:],
            'H': ds.variables['H'][:],
            'LE': ds.variables['LE'][:],
            'G': ds.variables['G0'][:],
        }
        # Mask fill values
        for key in ['ust', 'H', 'LE', 'G']:
            data[key] = np.ma.masked_equal(data[key], -9999.0)
    return data


def load_lsm_output(lsm_file: str, time_offset: float = 0.0) -> dict:
    """Load UtahLSM output data.

    Args:
        lsm_file: Path to LSM output NetCDF file.
        time_offset: Offset to add to LSM time to align with observations [hours].

    Returns:
        Dictionary containing time and flux model output.
    """
    with nc.Dataset(lsm_file, 'r') as ds:
        # Convert time to hours if needed
        time = ds.variables['time'][:]
        time_units = ds.variables['time'].units

        # Handle 3D output (squeeze out spatial dimensions for single column)
        def get_var(name):
            var = ds.variables[name][:]
            # If 3D (time, y, x), squeeze to 1D
            if var.ndim == 3:
                return var[:, 0, 0]
            return var

        data = {
            'time': time,
            'ust': get_var('ust'),
            'H': get_var('shf'),   # Sensible heat flux
            'LE': get_var('lhf'),  # Latent heat flux
            'G': get_var('ghf'),   # Ground heat flux
        }

        # Convert time from seconds to hours
        # LSM output is always in seconds since start of simulation
        data['time'] = data['time'] / 3600.0  # Convert to hours

        # Apply time offset to align with observations
        # (obs is in hours since start of month, e.g., 24, 24.1667, etc.)
        data['time'] = data['time'] + time_offset

    return data


def create_comparison_plot(obs: dict, lsm: dict, output_file: str = 'gabls3_comparison.png'):
    """Create multi-panel comparison plot.

    Args:
        obs: Observation data dictionary.
        lsm: LSM output data dictionary.
        output_file: Path to save output figure.
    """
    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    variables = [
        ('ust', r'$u_*$ [m s$^{-1}$]', 'Friction Velocity'),
        ('H', r'$H$ [W m$^{-2}$]', 'Sensible Heat Flux'),
        ('LE', r'$\Lambda E$ [W m$^{-2}$]', 'Latent Heat Flux'),
        ('G', r'$G$ [W m$^{-2}$]', 'Ground Heat Flux'),
    ]

    for ax, (var_name, ylabel, title) in zip(axes, variables):
        # Plot observations
        ax.plot(obs['time'], obs[var_name], 'ko-',
                label='Observations', linewidth=1.5, markersize=4)

        # Plot LSM output
        ax.plot(lsm['time'], lsm[var_name], 'r-',
                label='UtahLSM', linewidth=1.5)

        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold', loc='left')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='best', framealpha=0.9)

        # Add zero line for flux plots
        if var_name != 'ust':
            ax.axhline(0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

    # Set x-axis label only on bottom plot
    axes[-1].set_xlabel('Time [hours since 2006-07-01 00:00 UTC]', fontsize=11)
    axes[-1].set_xlim([obs['time'].min(), obs['time'].max()])

    # Overall title
    fig.suptitle('GABLS3 Comparison: UtahLSM vs. Observations',
                 fontsize=14, fontweight='bold', y=0.995)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved to: {output_file}")
    plt.show()


def compute_statistics(obs: dict, lsm: dict):
    """Compute and print comparison statistics.

    Args:
        obs: Observation data dictionary.
        lsm: LSM output data dictionary.
    """
    print("\n" + "="*60)
    print("COMPARISON STATISTICS")
    print("="*60)

    # Interpolate LSM to observation times (simple nearest neighbor)
    obs_times = obs['time']
    lsm_times = lsm['time']

    variables = ['ust', 'H', 'LE', 'G']
    var_names = ['u*', 'H', 'LE', 'G']

    for var, var_name in zip(variables, var_names):
        # Get valid observation indices
        valid = ~obs[var].mask if hasattr(obs[var], 'mask') else np.ones_like(obs[var], dtype=bool)

        if not np.any(valid):
            print(f"\n{var_name}: No valid observations")
            continue

        # Simple statistics on available data
        obs_mean = np.mean(obs[var][valid])
        lsm_mean = np.mean(lsm[var])

        # If same length, compute bias and RMSE
        if len(obs[var]) == len(lsm[var]):
            bias = np.mean(lsm[var][valid] - obs[var][valid])
            rmse = np.sqrt(np.mean((lsm[var][valid] - obs[var][valid])**2))

            print(f"\n{var_name}:")
            print(f"  Obs mean:  {obs_mean:8.3f}")
            print(f"  LSM mean:  {lsm_mean:8.3f}")
            print(f"  Bias:      {bias:8.3f}")
            print(f"  RMSE:      {rmse:8.3f}")
        else:
            print(f"\n{var_name}:")
            print(f"  Obs mean:  {obs_mean:8.3f}")
            print(f"  LSM mean:  {lsm_mean:8.3f}")
            print(f"  (Different time lengths - no bias/RMSE computed)")

    print("\n" + "="*60 + "\n")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Compare UtahLSM output with GABLS3 observations.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-l', '--lsm',
        type=str,
        default='output.nc',
        help='Path to LSM output NetCDF file'
    )
    parser.add_argument(
        '-o', '--obs',
        type=str,
        default='../cases/gabls3/observations/gabls3_fluxes.nc',
        help='Path to observation NetCDF file'
    )
    parser.add_argument(
        '-p', '--plot',
        type=str,
        default='gabls3_comparison.png',
        help='Output plot filename'
    )
    parser.add_argument(
        '--no-stats',
        action='store_true',
        help='Skip printing statistics'
    )
    parser.add_argument(
        '--time-offset',
        type=float,
        default=None,
        help='Manual time offset to add to LSM time [hours]. If not specified, automatically aligns with obs start time.'
    )

    args = parser.parse_args()

    # Check if files exist
    if not Path(args.obs).exists():
        print(f"Error: Observation file not found: {args.obs}")
        return 1

    if not Path(args.lsm).exists():
        print(f"Error: LSM output file not found: {args.lsm}")
        print("\nRun a simulation first:")
        print("  python utahlsm_offline.py -c GABLS3 -o output.nc")
        return 1

    # Load data
    print(f"Loading observations from: {args.obs}")
    obs = load_observations(args.obs)

    # Determine time offset for LSM to align with observations
    print(f"Observation time range: {obs['time'][0]:.2f} to {obs['time'][-1]:.2f} hours")

    if args.time_offset is not None:
        time_offset = args.time_offset
        print(f"Using manual time offset: {time_offset:.2f} hours")
    else:
        # Auto-align: use observation start time as offset
        time_offset = obs['time'][0]
        print(f"Auto-aligning LSM time with observation start: offset = {time_offset:.2f} hours")

    print(f"Loading LSM output from: {args.lsm}")
    lsm = load_lsm_output(args.lsm, time_offset=time_offset)
    print(f"LSM time range (adjusted): {lsm['time'][0]:.2f} to {lsm['time'][-1]:.2f} hours")

    # Create comparison plot
    print(f"Creating comparison plot...")
    create_comparison_plot(obs, lsm, args.plot)

    # Print statistics
    if not args.no_stats:
        compute_statistics(obs, lsm)

    return 0


if __name__ == '__main__':
    exit(main())
