#!/usr/bin/env python
"""
Example: Custom N-Level Model

Demonstrates how to create a custom rate equation model with an
arbitrary number of states.
"""

import numpy as np
import matplotlib.pyplot as plt
from odmr_sim.models import RateModel
from odmr_sim.solvers import RateSolver
from odmr_sim.plotting import plot_populations


def main():
    # --- Example 1: Simple 3-level system ---
    print("=== 3-Level System ===")

    model_3 = RateModel(
        n_states=3,
        state_labels=["Ground", "Excited", "Metastable"]
    )

    # Set transition rates (MHz)
    model_3.set_rates({
        (1, 0): 50.0,    # Excited -> Ground (radiative)
        (1, 2): 10.0,    # Excited -> Metastable (ISC)
        (2, 0): 5.0,     # Metastable -> Ground
    })

    # Add excitation as dynamic rate
    model_3.add_dynamic_rate('gamma', from_state=0, to_state=1)

    # Build and solve
    W = model_3.build_rate_matrix(gamma=1.0)
    P0 = np.array([1.0, 0.0, 0.0])  # Start in ground
    t_eval = np.logspace(-9, -3, 500)

    solver = RateSolver()
    populations = solver.solve_expm(W, P0, t_eval)

    print("Final populations:")
    for i, label in enumerate(model_3.state_labels):
        print(f"  {label}: {populations[-1, i]:.4f}")
    print()

    # --- Example 2: 5-level system (like simplified NV) ---
    print("=== 5-Level System ===")

    model_5 = RateModel(
        n_states=5,
        state_labels=["GS|0>", "GS|+/->", "ES|0>", "ES|+/->", "Singlet"]
    )

    # Transition rates similar to NV but simplified
    model_5.set_rates({
        # Radiative
        (2, 0): 62.5,   # ES|0> -> GS|0>
        (3, 1): 62.5,   # ES|+/-> -> GS|+/->
        # ISC up
        (2, 4): 4.0,    # ES|0> -> Singlet
        (3, 4): 30.0,   # ES|+/-> -> Singlet (stronger)
        # ISC down
        (4, 0): 1000.0, # Singlet -> GS|0> (strongly polarizing)
        (4, 1): 5.0,    # Singlet -> GS|+/->
    })

    # Dynamic rates
    model_5.add_dynamic_rate('gamma', 0, 2)  # GS|0> -> ES|0>
    model_5.add_dynamic_rate('gamma', 1, 3)  # GS|+/-> -> ES|+/->

    W = model_5.build_rate_matrix(gamma=0.1)
    P0 = np.array([0.5, 0.5, 0.0, 0.0, 0.0])  # Mixed ground state
    t_eval = np.logspace(-9, -1, 1000)

    populations = solver.solve_expm(W, P0, t_eval)

    print("Final populations:")
    for i, label in enumerate(model_5.state_labels):
        print(f"  {label}: {populations[-1, i]:.4f}")
    print()

    # --- Plot both ---
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # 3-level
    plot_populations(
        np.logspace(-9, -3, 500),
        solver.solve_expm(model_3.build_rate_matrix(gamma=1.0),
                          np.array([1.0, 0.0, 0.0]),
                          np.logspace(-9, -3, 500)),
        labels=model_3.state_labels,
        ax=axes[0],
        title="3-Level System",
        colors=['#4285F4', '#DB4437', '#0F9D58'],
        linestyles=['-', '-', '-'],
    )

    # 5-level
    plot_populations(
        np.logspace(-9, -1, 1000),
        solver.solve_expm(model_5.build_rate_matrix(gamma=0.1),
                          np.array([0.5, 0.5, 0.0, 0.0, 0.0]),
                          np.logspace(-9, -1, 1000)),
        labels=model_5.state_labels,
        ax=axes[1],
        title="5-Level System (Simplified NV)",
        colors=['#4285F4', '#DB4437', '#4285F4', '#DB4437', '#0F9D58'],
        linestyles=['-', '-', ':', ':', '--'],
    )

    plt.tight_layout()
    plt.savefig("custom_model_example.pdf", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
