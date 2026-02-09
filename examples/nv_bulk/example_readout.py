#!/usr/bin/env python
"""
Example: Readout Simulation

Demonstrates how to simulate optical readout dynamics from different
initial spin states.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Arial', 'Helvetica']
plt.rcParams.update({'font.size': 12})
from odmr_sim.models import get_preset
from odmr_sim.simulations import ReadoutSimulation
from odmr_sim.plotting import plot_readout_comparison


def main():
    # Create model using NV bulk preset
    model = get_preset('nv_bulk')

    # Create simulation
    sim = ReadoutSimulation(model)

    # Run comparison for all three ground states
    gamma = 12.8  # MHz - typical readout excitation rate

    results = sim.run_comparison(
        gamma=gamma,
        t_max=1e-5,     # 10 μs
        n_points=1000,
    )

    print(f"Readout simulation with Gamma = {gamma} MHz")
    print()

    for state_name, result in results.items():
        peak_es = result['es_total'].max()
        t_peak = result['t_ns'][np.argmax(result['es_total'])]
        print(f"Initial state {state_name}:")
        print(f"  Peak ES population: {peak_es:.4f}")
        print(f"  Time to peak: {t_peak:.1f} ns")
        print()

    # Plot comparison
    fig, ax = plt.subplots(figsize=(5, 3.5))
    plot_readout_comparison(
        results, ax=ax,
        xlim=(-50, 1200),
        ylim=(0, 0.2),
        title=f"Readout Dynamics ($\\Gamma$ = {gamma} MHz)"
    )

    plt.tight_layout()
    plt.savefig("readout_example.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
