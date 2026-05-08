#!/usr/bin/env python3
"""
Calculate the total integrated flux (Φ_a) for each resonant mass
by integrating the spectral flux over energy, following the Giannotti paper approach.

Usage:
    python3 integrated_tp_fluxes.py [--results-dir PATH]
"""

import argparse
import numpy as np
from pathlib import Path


def calculate_integrated_tp_fluxes(results_dir=None):
    """
    Calculate the total integrated flux (Φ_a) for each resonant mass by integrating
    the spectral flux over energy.
    
    Φ_a = ∫ (dΦ_a/dω) dω
    
    Args:
        results_dir: Path to directory containing TP_resonant_*.dat files
        
    Returns:
        Dictionary with mass labels and their integrated fluxes (cm^-2 s^-1)
    """
    if results_dir is None:
        results_dir = Path(__file__).resolve().parents[1] / "results"
    
    results_dir = Path(results_dir)
    
    # Dictionary to store integrated fluxes
    fluxes = {}
    
    # Resonant mass files with their labels
    mass_files = {
        "m131 (131 eV)": "TP_resonant_m131.dat",
        "m11 (11 eV)": "TP_resonant_m11.dat",
        "m0 (0 eV)": "TP_resonant_m0.dat",
    }
    
    print("\nIntegrated Transverse Plasmon Flux (Φ_a) for Each Mass")
    print("=" * 75)
    print(f"{'Mass':<20} {'Integrated Flux (Φ_a)':<30} {'Unit':<25}")
    print("-" * 75)
    
    for mass_label, filename in mass_files.items():
        filepath = results_dir / filename
        
        if not filepath.exists():
            print(f"{mass_label:<20} {'FILE NOT FOUND':<30} {'':<25}")
            fluxes[mass_label] = None
            continue
        
        try:
            # Load the spectral flux data
            # Format: energy (keV) in column 0, flux (cm^-2 s^-1 keV^-1) in column 1
            data = np.genfromtxt(str(filepath))
            
            if data.size == 0:
                print(f"{mass_label:<20} {'EMPTY FILE':<30} {'':<25}")
                fluxes[mass_label] = None
                continue
            
            # Handle both 1D and 2D arrays
            if data.ndim == 1:
                print(f"{mass_label:<20} {'INVALID FORMAT':<30} {'':<25}")
                fluxes[mass_label] = None
                continue
            
            # Extract energy (column 0) and flux (column 1)
            energies = data[:, 0]
            flux_spectrum = data[:, 1]
            
            # Integrate flux over energy using trapezoidal rule
            # Φ_a = ∫ (dΦ_a/dω) dω
            integrated_flux = float(np.trapz(flux_spectrum, energies))
            
            fluxes[mass_label] = integrated_flux
            
            # Format output with appropriate precision
            flux_str = f"{integrated_flux:.4e}"
            unit_str = "cm⁻² s⁻¹"
            
            print(f"{mass_label:<20} {flux_str:<30} {unit_str:<25}")
        
        except Exception as e:
            print(f"{mass_label:<20} {'ERROR: ' + str(e)[:25]:<30} {'':<25}")
            fluxes[mass_label] = None
    
    print("=" * 75)
    return fluxes


def main():
    parser = argparse.ArgumentParser(
        description="Calculate integrated TP flux for each resonant mass",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="Path to results directory with TP_resonant_*.dat files",
    )
    args = parser.parse_args()
    
    fluxes = calculate_integrated_tp_fluxes(args.results_dir)
    
    # Print summary
    print("\nSummary of integrated fluxes:")
    for mass_label, flux_value in fluxes.items():
        if flux_value is not None:
            print(f"  {mass_label:<25}: {flux_value:.4e} cm⁻² s⁻¹")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
