#!/usr/bin/env python3
"""
Initialization simulation for G11-G-9 NV configuration at 30 degrees.

This example demonstrates spin polarization dynamics using g11_gm9_30dp preset.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import get_preset
from odmr_sim.simulations import InitializationSimulation
from odmr_sim.plotting import plot_populations


def main():
    # Load G11-G-9 30dp preset (split ± levels)
    model = get_preset('g11_gm9_30dp')

    print("Rate parameters:")
    print(model.get_rate_summary())
    print()

    sim = InitializationSimulation(model)

    # Run initialization for different excitation rates
    gammas = [0.1, 0.3, 1.0]
    results = sim.run_sweep_gamma(gammas, t_max=1e-1, n_points=1000)

    # Create plot
    fig, axes = plt.subplots(len(gammas), 1, figsize=(6, 3*len(gammas)), sharex=True)

    for i, (gamma, result) in enumerate(zip(gammas, results)):
        plot_populations(
            result['t'], result['populations'],
            labels=result['labels'],
            ax=axes[i],
            xlabel="" if i < len(gammas) - 1 else "Time (ns)",
            title=f"G11-G-9 30dp Initialization",
            xlim=(1, 1e6),
            ylim=(-0.05, 1.05),
        )

        # Add annotation
        axes[i].text(
            1.02, 0.5, f"$\\Gamma$ = {gamma} MHz",
            transform=axes[i].transAxes,
            fontsize=12,
            verticalalignment='center'
        )

        # Final populations
        final_pops = result['populations'][-1]
        print(f"Gamma = {gamma} MHz - Final populations:")
        for j, label in enumerate(result['labels']):
            print(f"  {label}: {final_pops[j]:.4f}")
        print()

    plt.tight_layout()
    plt.savefig("g11_gm9_30dp_initialization.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
