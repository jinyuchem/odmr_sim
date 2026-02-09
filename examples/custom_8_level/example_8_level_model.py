#!/usr/bin/env python3
"""
Custom 8-level model example with two metastable states.

This model extends the 7-level NV model by adding a second metastable state.
The level structure is:
  - GS|0>, GS|->, GS|+>: Ground states (indices 0, 1, 2)
  - ES|0>, ES|->, ES|+>: Excited states (indices 3, 4, 5)
  - MS1 (singlet 1): First metastable state (index 6), connected to ES via ISC
  - MS2 (singlet 2): Second metastable state (index 7), connected to GS via ISC

The decay pathway is: ES -> MS1 -> MS2 -> GS
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import RateModel
from odmr_sim.solvers import RateSolver
from scipy.linalg import expm


def create_8_level_model():
    """Create an 8-level model with two metastable states."""

    # State labels
    labels = [
        r"GS|$0\rangle$",
        r"GS|$-\rangle$",
        r"GS|$+\rangle$",
        r"ES|$0\rangle$",
        r"ES|$-\rangle$",
        r"ES|$+\rangle$",
        "MS1",  # First metastable state (connected to ES)
        "MS2",  # Second metastable state (connected to GS)
    ]

    model = RateModel(n_states=8, state_labels=labels)

    # Rate parameters (same as NV bulk for ES and GS)
    k41 = 62.5   # ES|0> -> GS|0> (radiative)
    k52 = 62.5   # ES|-> -> GS|-> (radiative)
    k63 = 62.5   # ES|+> -> GS|+> (radiative)

    # Upper ISC: ES -> MS1 (first metastable)
    k47 = 10.5   # ES|0> -> MS1
    k57 = 76.9   # ES|-> -> MS1
    k67 = 76.9   # ES|+> -> MS1

    # MS1 -> MS2 decay
    k78 = 1000.0  # MS1 -> MS2

    # Lower ISC: MS2 -> GS (second metastable)
    k81 = 3.0  # MS2 -> GS|0>
    k82 = 2.63 # MS2 -> GS|->
    k83 = 2.63 # MS2 -> GS|+>

    # Set static rates
    model.set_rates({
        # Radiative decay
        (3, 0): k41,  # ES|0> -> GS|0>
        (4, 1): k52,  # ES|-> -> GS|->
        (5, 2): k63,  # ES|+> -> GS|+>
        # Upper ISC: ES -> MS1
        (3, 6): k47,  # ES|0> -> MS1
        (4, 6): k57,  # ES|-> -> MS1
        (5, 6): k67,  # ES|+> -> MS1
        # MS1 -> MS2 decay
        (6, 7): k78,
        # Lower ISC: MS2 -> GS
        (7, 0): k81,
        (7, 1): k82,
        (7, 2): k83,
    })

    # Dynamic rates for excitation and microwave
    # Excitation: GS -> ES
    model.add_dynamic_rate('gamma', from_state=0, to_state=3)  # GS|0> -> ES|0>
    model.add_dynamic_rate('gamma', from_state=1, to_state=4)  # GS|-> -> ES|->
    model.add_dynamic_rate('gamma', from_state=2, to_state=5)  # GS|+> -> ES|+>

    # Microwave: GS|0> <-> GS|+/->
    model.add_dynamic_rate('kmw_minus', from_state=0, to_state=1)
    model.add_dynamic_rate('kmw_minus', from_state=1, to_state=0)
    model.add_dynamic_rate('kmw_plus', from_state=0, to_state=2)
    model.add_dynamic_rate('kmw_plus', from_state=2, to_state=0)

    return model, {
        'k41': k41, 'k52': k52, 'k63': k63,
        'k47': k47, 'k57': k57, 'k67': k67,
        'k78': k78,
        'k81': k81, 'k82': k82, 'k83': k83,
    }


def run_initialization(model, gammas, t_max=1e-1, n_points=1000):
    """Run initialization simulation."""
    solver = RateSolver()
    t_eval = np.logspace(-9, np.log10(t_max), n_points)

    # Initial state: equal population in ground states
    P0 = np.array([1/3, 1/3, 1/3, 0, 0, 0, 0, 0])

    results = []
    for gamma in gammas:
        W = model.build_rate_matrix(gamma=gamma, kmw_minus=0.0, kmw_plus=0.0)
        populations = solver.solve_expm(W, P0, t_eval)
        results.append({
            't': t_eval * 1e9,  # Convert to ns
            'populations': populations,
            'labels': model.state_labels,
        })

    return results


def run_readout(model, gamma, t_max=1e-5, n_points=1000):
    """Run readout simulation comparing different initial states."""
    solver = RateSolver()
    t_eval = np.linspace(0, t_max, n_points)

    # Different initial states
    initial_states = {
        r'GS|$0\rangle$': np.array([1, 0, 0, 0, 0, 0, 0, 0]),
        r'GS|$-\rangle$': np.array([0, 1, 0, 0, 0, 0, 0, 0]),
        r'GS|$+\rangle$': np.array([0, 0, 1, 0, 0, 0, 0, 0]),
    }

    W = model.build_rate_matrix(gamma=gamma, kmw_minus=0.0, kmw_plus=0.0)

    results = {}
    for name, P0 in initial_states.items():
        populations = solver.solve_expm(W, P0, t_eval)
        # ES total = ES|0> + ES|-> + ES|+>
        es_total = populations[:, 3] + populations[:, 4] + populations[:, 5]
        results[name] = {
            't_ns': t_eval * 1e9,
            'es_total': es_total,
        }

    return results


def compute_odmr_contrast(model, gamma, kmw_minus, kmw_plus):
    """Compute ODMR contrast using steady-state populations."""
    solver = RateSolver()

    # Without microwave
    W_nomw = model.build_rate_matrix(gamma=gamma, kmw_minus=0.0, kmw_plus=0.0)
    P_ss_nomw = solver.solve_steady_state(W_nomw)

    # With microwave
    W_mw = model.build_rate_matrix(gamma=gamma, kmw_minus=kmw_minus, kmw_plus=kmw_plus)
    P_ss_mw = solver.solve_steady_state(W_mw)

    # PL intensity proportional to ES populations (assuming equal radiative rates)
    I_nomw = P_ss_nomw[3] + P_ss_nomw[4] + P_ss_nomw[5]
    I_mw = P_ss_mw[3] + P_ss_mw[4] + P_ss_mw[5]

    if I_nomw == 0:
        return 0.0

    return (I_nomw - I_mw) / I_nomw


def main():
    print("=== 8-Level Model with Two Metastable States ===")
    print()

    # Create model
    model, rates = create_8_level_model()

    print("Level structure:")
    print("  GS|0>, GS|->, GS|+>  (indices 0, 1, 2)")
    print("  ES|0>, ES|->, ES|+>  (indices 3, 4, 5)")
    print("  MS1 (singlet 1)      (index 6) - connected to ES")
    print("  MS2 (singlet 2)      (index 7) - connected to GS")
    print()
    print("Decay pathway: ES -> MS1 -> MS2 -> GS")
    print()
    print("Rate parameters:")
    for key, value in rates.items():
        print(f"  {key} = {value} MHz")
    print()

    # --- Initialization ---
    print("=== Initialization Simulation ===")
    gammas = [0.1, 0.3, 1.0]
    results = run_initialization(model, gammas)

    fig, axes = plt.subplots(len(gammas), 1, figsize=(6, 3*len(gammas)), sharex=True)

    colors = ['#4285F4', '#DB4437', '#0F9D58', '#F4B400', '#AB47BC',
              '#FF7043', '#26A69A', '#78909C']

    for i, (gamma, result) in enumerate(zip(gammas, results)):
        for j, label in enumerate(result['labels']):
            axes[i].semilogx(result['t'], result['populations'][:, j],
                            label=label, color=colors[j], linewidth=1.5)

        axes[i].set_ylabel('Population')
        axes[i].set_ylim(-0.05, 1.05)
        axes[i].legend(loc='right', fontsize=8, ncols=2)
        axes[i].grid(True, alpha=0.3)
        axes[i].text(1.02, 0.5, f"$\\Gamma$ = {gamma} MHz",
                    transform=axes[i].transAxes, fontsize=12, va='center')

    axes[-1].set_xlabel('Time (ns)')
    plt.suptitle('8-Level Model Initialization', fontsize=12)
    plt.tight_layout()
    plt.savefig("8_level_initialization.pdf", bbox_inches='tight', dpi=300)

    # Print final populations
    for gamma, result in zip(gammas, results):
        print(f"Gamma = {gamma} MHz - Final populations:")
        for j, label in enumerate(result['labels']):
            print(f"  {label}: {result['populations'][-1, j]:.4f}")
        print()

    # --- Readout ---
    print("=== Readout Simulation ===")
    gamma = 12.8
    readout_results = run_readout(model, gamma=gamma)

    fig, ax = plt.subplots(figsize=(5, 3.5))

    for name, result in readout_results.items():
        ax.plot(result['t_ns'], result['es_total'], label=name, linewidth=1.5)

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('ES Population')
    ax.set_title(f'8-Level Readout ($\\Gamma$ = {gamma} MHz)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-50, 1200)
    plt.tight_layout()
    plt.savefig("8_level_readout.pdf", bbox_inches='tight', dpi=300)

    for name, result in readout_results.items():
        peak_idx = result['es_total'].argmax()
        print(f"{name}: Peak ES = {result['es_total'][peak_idx]:.4f} at {result['t_ns'][peak_idx]:.1f} ns")
    print()

    # --- ODMR Contrast ---
    print("=== ODMR Contrast ===")
    gammas = np.array([0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8])

    contrasts_both = []
    contrasts_minus = []
    for g in gammas:
        c_both = compute_odmr_contrast(model, g, kmw_minus=1.0, kmw_plus=1.0)
        c_minus = compute_odmr_contrast(model, g, kmw_minus=1.0, kmw_plus=0.0)
        contrasts_both.append(c_both)
        contrasts_minus.append(c_minus)

    contrasts_both = np.array(contrasts_both)
    contrasts_minus = np.array(contrasts_minus)

    print("Gamma (MHz) | Both MW | Only minus")
    print("-" * 45)
    for g, cb, cm in zip(gammas, contrasts_both, contrasts_minus):
        print(f"{g:7.1f}     | {cb*100:6.2f}% | {cm*100:6.2f}%")
    print()

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.semilogx(gammas, contrasts_both * 100, 'o-', label=r'Both $k_{\mathrm{MW}}^{\pm}$ = 1 MHz')
    ax.semilogx(gammas, contrasts_minus * 100, 's--', label=r'Only $k_{\mathrm{MW}}^{-}$ = 1 MHz')
    ax.set_xlabel(r'Excitation rate $\Gamma$ (MHz)')
    ax.set_ylabel('ODMR Contrast (%)')
    ax.set_title('8-Level Model ODMR Contrast', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("8_level_odmr_contrast.pdf", bbox_inches='tight', dpi=300)

    # --- ODMR Spectrum ---
    print("=== ODMR Spectrum ===")

    # ZFS frequency (typical for NV-like systems)
    zfs_freq = 2.87  # GHz
    linewidth = 0.02  # GHz

    # Use gamma = 0.1 MHz for spectrum
    gamma = 0.1
    contrast_both = compute_odmr_contrast(model, gamma, kmw_minus=1.0, kmw_plus=1.0)
    contrast_minus = compute_odmr_contrast(model, gamma, kmw_minus=1.0, kmw_plus=0.0)
    contrast_plus = compute_odmr_contrast(model, gamma, kmw_minus=0.0, kmw_plus=1.0)

    print(f"At gamma = {gamma} MHz:")
    print(f"  Contrast (both): {contrast_both*100:.2f}%")
    print(f"  Contrast (minus only): {contrast_minus*100:.2f}%")
    print(f"  Contrast (plus only): {contrast_plus*100:.2f}%")
    print()

    # Generate Lorentzian spectrum (assuming degenerate +/- like NV)
    frequencies = np.linspace(2.4, 3.4, 1001)
    lorentzian = linewidth**2 / ((frequencies - zfs_freq)**2 + linewidth**2)

    # Spectrum with both MW (degenerate case - single peak)
    spectrum = contrast_both * lorentzian

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(frequencies, -spectrum * 100, 'b-', linewidth=2)
    ax.axhline(y=0, color='k', linestyle=':', alpha=0.5)
    ax.axvline(x=zfs_freq, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('ODMR Signal (%)')
    ax.set_title(f'8-Level ODMR Spectrum at ZFS = {zfs_freq} GHz', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.text(zfs_freq, -22, f'ZFS = {zfs_freq} GHz', ha='center', fontsize=10)
    ax.set_ylim(-25, 0)
    plt.tight_layout()
    plt.savefig("8_level_odmr_spectrum.pdf", bbox_inches='tight', dpi=300)

    plt.show()


if __name__ == "__main__":
    main()
