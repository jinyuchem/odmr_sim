#!/usr/bin/env python3
"""
ODMR contrast comparison: Steady-state vs Transient methods.

This example demonstrates the difference between:
- Steady-state: True equilibrium (t -> infinity)
- Transient: Final ES population at finite time

For g11_gm9_30dp, these give different results due to the asymmetric ISC rates.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import get_preset
from odmr_sim.simulations import ODMRSimulation


def main():
    # Load G11-Gm9 30dp preset
    model = get_preset('g11_gm9_30dp')
    sim = ODMRSimulation(model)

    print("=== G11-Gm9 30dp: Steady-State vs Transient Contrast ===")
    print()
    print("Rate parameters (key asymmetries):")
    print(f"  k47 = {model.k47} MHz (ES|0> -> SS)")
    print(f"  k57 = {model.k57} MHz (ES|-> -> SS)")
    print(f"  k67 = {model.k67} MHz (ES|+> -> SS)")
    print(f"  k71 = {model.k71} MHz (SS -> GS|0>)")
    print(f"  k72 = {model.k72} MHz (SS -> GS|->)")
    print(f"  k73 = {model.k73} MHz (SS -> GS|+>)")
    print()

    gammas = np.array([0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8])
    t_integration = 0.1  # 0.1 seconds

    # --- Steady-State Method ---
    print("=== Steady-State Method (t -> infinity) ===")
    _, contrasts_ss_minus = sim.sweep_gamma(gammas, kmw_minus=1.0, kmw_plus=0.0, method='steady_state')
    _, contrasts_ss_plus = sim.sweep_gamma(gammas, kmw_minus=0.0, kmw_plus=1.0, method='steady_state')

    print("Gamma (MHz) | Contrast (-) | Contrast (+)")
    print("-" * 50)
    for g, cm, cp in zip(gammas, contrasts_ss_minus, contrasts_ss_plus):
        print(f"{g:7.1f}     | {cm*100:11.2f}% | {cp*100:11.2f}%")
    print()

    # --- Transient Method ---
    print(f"=== Transient Method (t = {t_integration} s) ===")
    contrasts_tr_minus = []
    contrasts_tr_plus = []
    for g in gammas:
        # For transient, we scale kmw with gamma
        kmw = g * 10  # kmw = gamma * 10
        c_minus = sim.compute_contrast(gamma=g, kmw_minus=kmw, kmw_plus=0.0,
                                        method='transient', t_integration=t_integration)
        c_plus = sim.compute_contrast(gamma=g, kmw_minus=0.0, kmw_plus=kmw,
                                       method='transient', t_integration=t_integration)
        contrasts_tr_minus.append(c_minus)
        contrasts_tr_plus.append(c_plus)

    contrasts_tr_minus = np.array(contrasts_tr_minus)
    contrasts_tr_plus = np.array(contrasts_tr_plus)

    print("Gamma (MHz) | Contrast (-) | Contrast (+)")
    print("-" * 50)
    for g, cm, cp in zip(gammas, contrasts_tr_minus, contrasts_tr_plus):
        print(f"{g:7.1f}     | {cm*100:11.2f}% | {cp*100:11.2f}%")
    print()

    print("Note: Negative contrast means PL INCREASES with MW (inverted peak)")
    print()

    # --- Plot comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # |-⟩ transition
    axes[0].semilogx(gammas, contrasts_ss_minus * 100, 'o-', label='Steady-state')
    axes[0].semilogx(gammas, contrasts_tr_minus * 100, 's--', label='Transient')
    axes[0].axhline(y=0, color='k', linestyle=':', alpha=0.5)
    axes[0].set_xlabel(r'Excitation rate $\Gamma$ (MHz)')
    axes[0].set_ylabel('ODMR Contrast (%)')
    axes[0].set_title(r'$|-\rangle$ Transition', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # |+⟩ transition
    axes[1].semilogx(gammas, contrasts_ss_plus * 100, 'o-', label='Steady-state')
    axes[1].semilogx(gammas, contrasts_tr_plus * 100, 's--', label='Transient')
    axes[1].axhline(y=0, color='k', linestyle=':', alpha=0.5)
    axes[1].set_xlabel(r'Excitation rate $\Gamma$ (MHz)')
    axes[1].set_ylabel('ODMR Contrast (%)')
    axes[1].set_title(r'$|+\rangle$ Transition', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.suptitle('G11-Gm9 30dp: Steady-State vs Transient Contrast', fontsize=12)
    plt.tight_layout()
    plt.savefig("g11_gm9_30dp_contrast_comparison.pdf", bbox_inches='tight', dpi=300)

    # --- ODMR Spectrum Comparison ---
    print("=== ODMR Spectrum Comparison ===")

    # Peak frequencies from plot_v3.py
    peak_freq_minus = (3439.97 - 480.98) * 1e-3  # 2.959 GHz
    peak_freq_plus = (3439.97 + 480.98) * 1e-3   # 3.921 GHz

    # Use gamma=0.1 for spectrum
    gamma = 0.1
    kmw = gamma * 10  # For transient method

    # Steady-state contrasts
    c_ss_minus = sim.compute_contrast(gamma=gamma, kmw_minus=1.0, kmw_plus=0.0, method='steady_state')
    c_ss_plus = sim.compute_contrast(gamma=gamma, kmw_minus=0.0, kmw_plus=1.0, method='steady_state')

    # Transient contrasts
    c_tr_minus = sim.compute_contrast(gamma=gamma, kmw_minus=kmw, kmw_plus=0.0,
                                       method='transient', t_integration=t_integration)
    c_tr_plus = sim.compute_contrast(gamma=gamma, kmw_minus=0.0, kmw_plus=kmw,
                                      method='transient', t_integration=t_integration)

    print(f"Steady-state: Contrast(-) = {c_ss_minus*100:.2f}%, Contrast(+) = {c_ss_plus*100:.2f}%")
    print(f"Transient:    Contrast(-) = {c_tr_minus*100:.2f}%, Contrast(+) = {c_tr_plus*100:.2f}%")
    print()

    # Generate Lorentzian lineshapes
    frequencies = np.linspace(2.5, 4.5, 1001)
    linewidth = 0.02  # GHz

    lorentz_minus = linewidth**2 / ((frequencies - peak_freq_minus)**2 + linewidth**2)
    lorentz_plus = linewidth**2 / ((frequencies - peak_freq_plus)**2 + linewidth**2)

    # Create spectra
    spectrum_ss = c_ss_minus * lorentz_minus + c_ss_plus * lorentz_plus
    spectrum_tr = c_tr_minus * lorentz_minus + c_tr_plus * lorentz_plus

    # Plot ODMR spectra comparison
    fig, ax = plt.subplots(figsize=(6, 4))

    # Steady-state spectrum (negative for dip convention)
    ax.plot(frequencies, -spectrum_ss * 100, 'b-', linewidth=2, label='Steady-state')

    # Transient spectrum (negative for dip convention)
    ax.plot(frequencies, -spectrum_tr * 100, 'r--', linewidth=2, label='Transient')

    ax.axhline(y=0, color='k', linestyle=':', alpha=0.5)
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('ODMR Signal (%)')
    ax.set_title(f'G11-Gm9 30dp ODMR Spectrum ($\\Gamma$ = {gamma} MHz)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add annotations for peak positions
    ax.axvline(x=peak_freq_minus, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=peak_freq_plus, color='gray', linestyle=':', alpha=0.5)
    ax.text(peak_freq_minus, -0.25, r'$|-\rangle$', ha='center', fontsize=12)
    ax.text(peak_freq_plus, 0.0, r'$|+\rangle$', ha='center', fontsize=12)

    plt.tight_layout()
    plt.savefig("g11_gm9_30dp_spectrum_comparison.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
