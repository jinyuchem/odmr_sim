#!/usr/bin/env python3
"""
Readout simulation for G9-G8 NV configuration at 30 degrees.

This example demonstrates optical readout dynamics using g9_g8_30dp preset.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import get_preset
from odmr_sim.simulations import ReadoutSimulation
from odmr_sim.plotting import plot_readout_comparison


def main():
    # Load G9-G8 30dp preset (split +/- levels)
    model = get_preset('g9_g8_30dp')

    sim = ReadoutSimulation(model)

    # Run readout with high excitation rate
    gamma = 12.8  # MHz
    results = sim.run_comparison(
        gamma=gamma,
        t_max=1e-5,  # 10 us
        n_points=1000,
    )

    print(f"Readout simulation with Gamma = {gamma} MHz")
    print()

    for state_name, result in results.items():
        if result['es_total'] is not None:
            peak_idx = result['es_total'].argmax()
            print(f"Initial state {state_name}:")
            print(f"  Peak ES population: {result['es_total'][peak_idx]:.4f}")
            print(f"  Time to peak: {result['t_ns'][peak_idx]:.1f} ns")
            print()

    # Plot
    fig, ax = plt.subplots(figsize=(5, 3.5))
    plot_readout_comparison(
        results, ax=ax,
        xlim=(-50, 1200),
        ylim=(0, 0.2),
        title=f"G9-G8 30dp Readout Dynamics ($\\Gamma$ = {gamma} MHz)"
    )

    plt.tight_layout()
    plt.savefig("g9_g8_30dp_readout.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
