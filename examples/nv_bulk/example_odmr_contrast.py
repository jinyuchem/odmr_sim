#!/usr/bin/env python
"""
Example: ODMR Contrast Simulation

Demonstrates how to compute ODMR contrast and simulate ODMR spectra.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import get_preset
from odmr_sim.simulations import ODMRSimulation
from odmr_sim.plotting import plot_odmr_spectrum


def main():
    # Create model using NV bulk preset
    model = get_preset('nv_bulk')

    # Create simulation
    sim = ODMRSimulation(model)

    # --- Part 1: Contrast vs excitation rate ---
    print("=== ODMR Contrast vs Excitation Rate ===")
    print("For NV@bulk, the +/-1 levels are degenerate, so MW drives both transitions.")
    print()

    gammas = np.array([0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8])

    # For NV@bulk: driving BOTH +/-1 transitions together (degenerate case)
    _, contrasts_both = sim.sweep_gamma(
        gammas, kmw_minus=1.0, kmw_plus=1.0
    )

    # Also show individual transitions for comparison (non-degenerate case)
    _, contrasts_minus = sim.sweep_gamma(
        gammas, kmw_minus=1.0, kmw_plus=0.0
    )

    print("Gamma (MHz) | Contrast (+/-1 together) | Contrast (-1 only)")
    print("-" * 55)
    for g, cb, cm in zip(gammas, contrasts_both, contrasts_minus):
        print(f"{g:7.1f} | {cb*100:21.2f}% | {cm*100:17.2f}%")
    print()

    # Plot contrast vs gamma
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.semilogx(gammas, contrasts_both * 100, 'o-', label=r'$k_{\mathrm{MW}}^{\pm 1}$ = 1 MHz (NV@bulk)')
    ax.semilogx(gammas, contrasts_minus * 100, 's--', alpha=0.6, label=r'$k_{\mathrm{MW}}^{-1}$ only = 1 MHz')
    ax.set_xlabel(r'Excitation rate $\Gamma$ (MHz)')
    ax.set_ylabel('ODMR Contrast (%)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("odmr_contrast_vs_gamma.pdf", bbox_inches='tight', dpi=300)

    # --- Part 2: ODMR Spectrum ---
    print("=== ODMR Spectrum ===")

    # For NV@bulk: both +/-1 peaks are at the SAME frequency (degenerate at zero field)
    # Peak position from plot_v3.py: mu1 = 3230.78 * 1e-3 GHz
    peak_freq = 3.231  # GHz

    # Compute the ACTUAL contrast (peak height) using the rate model
    # For NV@bulk with degenerate +/-1, drive both transitions
    peak_contrast = sim.compute_contrast(
        gamma=0.1,
        kmw_minus=1.0,  # MW driving both transitions
        kmw_plus=1.0
    )
    print(f"Peak contrast: {peak_contrast*100:.2f}%")
    print()

    # Generate arbitrary Lorentzian lineshape for visualization
    # The linewidth is phenomenological - just for display
    frequencies = np.linspace(0, 6, 1001)  # GHz
    linewidth = 0.05  # Arbitrary linewidth for visualization (GHz)

    # Lorentzian: L(x) = gamma^2 / ((x - x0)^2 + gamma^2)
    lorentzian = linewidth**2 / ((frequencies - peak_freq)**2 + linewidth**2)
    contrast_spectrum = peak_contrast * lorentzian  # Scale to actual peak contrast

    # Plot spectrum
    fig, ax = plt.subplots(figsize=(5, 3))
    plot_odmr_spectrum(
        frequencies,
        contrast_spectrum,
        ax=ax,
        xlim=(0, 6),
        title="NV@Bulk ODMR Spectrum"
    )

    plt.tight_layout()
    plt.savefig("odmr_spectrum_example.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
