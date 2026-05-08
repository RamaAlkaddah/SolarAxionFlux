#!/usr/bin/env python3
"""
Plot spectral transverse plasmon flux for different axion masses
in the style of Giannotti et al., following the Breit-Wigner resonance model.

This recreates a figure similar to the Giannotti paper showing how the 
resonant TP mode flux varies with energy for different axion masses.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def plot_tp_resonances(results_dir=None, output_file=None):
    """
    Create a log-log plot of spectral TP flux for multiple axion masses.
    
    Args:
        results_dir: Path to directory containing TP_resonant_*.dat files
        output_file: Path to save the figure (PDF/PNG)
    """
    if results_dir is None:
        results_dir = Path(__file__).resolve().parents[1] / "results"
    
    results_dir = Path(results_dir)
    
    # Mass configurations to plot
    masses = {
        "m = 131 eV": {
            "file": "TP_resonant_m131.dat",
            "color": "#1f77b4",
            "linestyle": "-",
            "linewidth": 2.5,
            "alpha": 0.9,
        },
        "m = 11 eV": {
            "file": "TP_resonant_m11.dat",
            "color": "#eb327f",
            "linestyle": "--",
            "linewidth": 2.5,
            "alpha": 0.9,
        },
        "m = 0 eV (off-res)": {
            "file": "TP_resonant_m0.dat",
            "color": "#2ca02c",
            "linestyle": ":",
            "linewidth": 2.5,
            "alpha": 0.9,
        },
    }
    
    # Create figure
    # Note: LaTeX rendering may fail on some systems; disable if needed
    try:
        plt.rc('text', usetex=True)
        plt.rc('text.latex', 
               preamble=r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}')
    except:
        plt.rc('text', usetex=False)
        print("Warning: LaTeX rendering disabled, using basic text")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot each mass configuration
    for label, config in masses.items():
        filepath = results_dir / config["file"]
        
        if not filepath.exists():
            print(f"Warning: {config['file']} not found, skipping {label}")
            continue
        
        try:
            # Load data: column 0 = energy (keV), column 1 = flux (cm^-2 s^-1 keV^-1)
            data = np.genfromtxt(str(filepath))
            
            if data.size == 0 or data.ndim != 2:
                print(f"Warning: {config['file']} has invalid format, skipping {label}")
                continue
            
            energy = data[:, 0]
            flux = data[:, 1]
            
            # Plot on log-log scale
            ax.plot(energy, flux,
                   label=label,
                   color=config["color"],
                   linestyle=config["linestyle"],
                   linewidth=config["linewidth"],
                   alpha=config["alpha"])
            
            print(f"✓ Loaded {label}: {len(energy)} energy points")
            
        except Exception as e:
            print(f"Error loading {config['file']}: {e}")
    
    # Set log scales for both axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Labels and title
    ax.set_xlabel(r'Energy $\omega$ [keV]', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'Spectral flux $\Phi_a(\omega)$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]', 
                  fontsize=14, fontweight='bold')
    
    ax.set_title(r'Transverse Plasmon Resonance: Axion Mass Dependence' + '\n' + 
                 r'(AGSS09 solar model, $g_{a\gamma\gamma} = 5 \times 10^{-11}$ GeV$^{-1}$)',
                 fontsize=13, fontweight='bold', pad=15)
    
    # Grid
    ax.grid(True, which='both', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.grid(True, which='minor', alpha=0.15, linestyle=':', linewidth=0.3)
    
    # Legend
    ax.legend(fontsize=12, frameon=True, fancybox=True, shadow=True, loc='best')
    
    # Tick parameters
    ax.tick_params(which='both', direction='in', bottom=True, top=True, 
                   left=True, right=True, labelsize=11)
    ax.tick_params(which='major', length=6, width=1.2)
    ax.tick_params(which='minor', length=3, width=0.8)
    
    # Set reasonable axis limits
    ax.set_xlim([0.05, 20])
    ax.set_ylim([1e4, 1e13])
    
    # Tight layout
    plt.tight_layout()
    
    # Save figure
    if output_file is None:
        output_file = results_dir / "tp_resonance_comparison.pdf"
    
    output_file = Path(output_file)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved to {output_file}")
    
    return fig, ax


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Plot TP resonance spectra for different axion masses"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="Path to results directory with TP_resonant_*.dat files"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output file path for the figure (PDF/PNG)"
    )
    args = parser.parse_args()
    
    plot_tp_resonances(args.results_dir, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
