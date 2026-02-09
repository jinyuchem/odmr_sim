#!/usr/bin/env python
"""
Example: Initialization Simulation

Demonstrates how to simulate spin initialization (polarization) dynamics
starting from a mixed ground state.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Arial', 'Helvetica']
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import SevenLevelModel, get_preset
from odmr_sim.simulations import InitializationSimulation
from odmr_sim.plotting import plot_populations


def main():
    # Create model using NV bulk preset
    model = get_preset('nv_bulk')

    print("Rate parameters:")
    print(model.get_rate_summary())
    print()

    # Create simulation
    sim = InitializationSimulation(model)

    # Run for different excitation rates
    gammas = [0.1, 0.3, 1.0]  # MHz

    fig, axes = plt.subplots(len(gammas), 1, figsize=(8, 8), sharex=True)

    for i, gamma in enumerate(gammas):
        result = sim.run(
            gamma=gamma,
            kmw_minus=0.0,  # No microwave
            kmw_plus=0.0,
            t_max=1e-1,     # 100 ms
            t_min=1e-9,     # 1 ns
            n_points=1000,
        )

        # Plot
        plot_populations(
            result['t'], result['populations'],
            labels=result['labels'],
            ax=axes[i],
            xlim=(1, 1e6),
            show_legend=(i == 0),
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

    axes[-1].set_xlabel("Time (ns)")
    fig.supylabel("Population", x=0.02)
    plt.tight_layout()
    plt.savefig("initialization_example.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
