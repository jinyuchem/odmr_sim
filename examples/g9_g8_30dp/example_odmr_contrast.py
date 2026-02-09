#!/usr/bin/env python3
"""
ODMR contrast simulation for G9-G8 NV configuration at 30 degrees.

This example demonstrates ODMR simulations using g9_g8_30dp preset.
For this NV configuration, the + and - levels are split (non-degenerate).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import get_preset
from odmr_sim.simulations import ODMRSimulation
from odmr_sim.plotting import plot_odmr_spectrum


def main():
    # Load G9-G8 30dp preset (split +/- levels)
    model = get_preset('g9_g8_30dp')

    sim = ODMRSimulation(model)

    # --- Part 1: Contrast vs excitation rate ---
    print("=== ODMR Contrast vs Excitation Rate ===")
    print("For this NV configuration, the + and - levels are split (non-degenerate).")
    print()

    gammas = np.array([0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8])

    # Contrast when driving ONLY the - transition
    _, contrasts_minus = sim.sweep_gamma(
        gammas, kmw_minus=1.0, kmw_plus=0.0
    )

    # Contrast when driving ONLY the + transition
    _, contrasts_plus = sim.sweep_gamma(
        gammas, kmw_minus=0.0, kmw_plus=1.0
    )

    print("Gamma (MHz) | Contrast (-) | Contrast (+)")
    print("-" * 50)
    for g, cm, cp in zip(gammas, contrasts_minus, contrasts_plus):
        print(f"{g:7.1f}     | {cm*100:11.2f}% | {cp*100:11.2f}%")
    print()

    # Plot contrast vs gamma
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.semilogx(gammas, contrasts_minus * 100, 'o-', label=r'$k_{\mathrm{MW}}^{-}$ = 1 MHz')
    ax.semilogx(gammas, contrasts_plus * 100, 's-', label=r'$k_{\mathrm{MW}}^{+}$ = 1 MHz')
    ax.set_xlabel(r'Excitation rate $\Gamma$ (MHz)')
    ax.set_ylabel('ODMR Contrast (%)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title("G9-G8 30dp ODMR Contrast")
    plt.savefig("g9_g8_30dp_contrast_vs_gamma.pdf", bbox_inches='tight', dpi=300)

    # --- Part 2: ODMR Spectrum with TWO peaks ---
    print("=== ODMR Spectrum ===")

    # For NV configuration: + and - are at DIFFERENT frequencies (non-degenerate)
    # From plot_v3.py: center = 3465.40 MHz, splitting = 166.94 MHz
    peak_freq_minus = (3465.40 - 166.94) * 1e-3  # 3.298 GHz (- transition)
    peak_freq_plus = (3465.40 + 166.94) * 1e-3   # 3.632 GHz (+ transition)

    # Compute contrast for each transition
    contrast_minus = sim.compute_contrast(gamma=0.1, kmw_minus=1.0, kmw_plus=0.0)
    contrast_plus = sim.compute_contrast(gamma=0.1, kmw_minus=0.0, kmw_plus=1.0)

    print(f"Peak contrast (-): {contrast_minus*100:.2f}%")
    print(f"Peak contrast (+): {contrast_plus*100:.2f}%")
    print()

    # Generate spectrum with TWO Lorentzian peaks
    frequencies = np.linspace(0.8, 4.2, 1001)  # GHz
    linewidth = 0.02  # GHz

    # Two separate Lorentzian peaks
    lorentz_minus = linewidth**2 / ((frequencies - peak_freq_minus)**2 + linewidth**2)
    lorentz_plus = linewidth**2 / ((frequencies - peak_freq_plus)**2 + linewidth**2)

    contrast_spectrum = contrast_minus * lorentz_minus + contrast_plus * lorentz_plus

    # Plot spectrum
    fig, ax = plt.subplots(figsize=(5, 3))
    plot_odmr_spectrum(
        frequencies,
        contrast_spectrum,
        ax=ax,
        xlim=(1, 4.2),
        title="G9-G8 30dp ODMR Spectrum (Split Peaks)"
    )

    plt.tight_layout()
    plt.savefig("g9_g8_30dp_odmr_spectrum.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
