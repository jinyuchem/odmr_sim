"""
Solvers for rate equation systems.

Provides matrix exponential and ODE integration methods for solving
the rate equation dP/dt = W @ P.
"""

import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from typing import Optional, Tuple, Literal


class RateSolver:
    """Solver for rate equation systems.

    Solves the master equation dP/dt = W @ P where:
    - P is the population vector (probability in each state)
    - W is the rate matrix

    Two methods are available:
    1. Matrix exponential: P(t) = exp(W*t) @ P0 (exact, fast for small systems)
    2. ODE integration: scipy.integrate.solve_ivp (flexible, handles stiff systems)

    Examples
    --------
    >>> from odmr_sim.models import SevenLevelModel
    >>> from odmr_sim.solvers import RateSolver
    >>>
    >>> model = SevenLevelModel()
    >>> W = model.build_rate_matrix(gamma=0.1)
    >>> P0 = model.get_ground_state_mixed()
    >>> t_eval = np.logspace(-9, -1, 1000)
    >>>
    >>> solver = RateSolver()
    >>> populations = solver.solve_expm(W, P0, t_eval)
    """

    def solve_expm(
        self,
        W: np.ndarray,
        P0: np.ndarray,
        t_eval: np.ndarray
    ) -> np.ndarray:
        """Solve using matrix exponential method.

        This is the exact solution: P(t) = exp(W*t) @ P0

        Parameters
        ----------
        W : np.ndarray
            Rate matrix with shape (n_states, n_states) in units of 1/s.
        P0 : np.ndarray
            Initial population vector with shape (n_states,).
        t_eval : np.ndarray
            Time points at which to evaluate the solution, in seconds.

        Returns
        -------
        populations : np.ndarray
            Population array with shape (len(t_eval), n_states).
            populations[i, j] is the population of state j at time t_eval[i].
        """
        n_states = len(P0)
        n_times = len(t_eval)

        populations = np.zeros((n_times, n_states))

        for i, t in enumerate(t_eval):
            populations[i] = expm(W * t) @ P0

        return populations

    def solve_ivp(
        self,
        W: np.ndarray,
        P0: np.ndarray,
        t_span: Tuple[float, float],
        t_eval: Optional[np.ndarray] = None,
        method: Literal['RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'] = 'RK45',
        **kwargs
    ) -> np.ndarray:
        """Solve using ODE integration (scipy.integrate.solve_ivp).

        This method is more flexible and can handle stiff systems with
        appropriate method choice (e.g., 'Radau', 'BDF').

        Parameters
        ----------
        W : np.ndarray
            Rate matrix with shape (n_states, n_states) in units of 1/s.
        P0 : np.ndarray
            Initial population vector with shape (n_states,).
        t_span : tuple
            (t0, tf) - start and end times in seconds.
        t_eval : np.ndarray, optional
            Time points at which to store the solution. If None,
            the solver chooses the points.
        method : str
            Integration method. Default is 'RK45'.
            Use 'Radau' or 'BDF' for stiff problems.
        **kwargs
            Additional arguments passed to scipy.integrate.solve_ivp.

        Returns
        -------
        populations : np.ndarray
            Population array with shape (len(t_eval), n_states).
        """
        def rate_equations(t, P):
            return W @ P

        solution = solve_ivp(
            rate_equations,
            t_span,
            P0,
            method=method,
            t_eval=t_eval,
            **kwargs
        )

        if not solution.success:
            raise RuntimeError(f"ODE solver failed: {solution.message}")

        # Transpose to get (n_times, n_states) shape
        return solution.y.T

    def solve_steady_state(self, W: np.ndarray) -> np.ndarray:
        """Find the steady-state population (null space of W).

        The steady state satisfies W @ P_ss = 0.

        Parameters
        ----------
        W : np.ndarray
            Rate matrix with shape (n_states, n_states).

        Returns
        -------
        P_ss : np.ndarray
            Normalized steady-state population vector.
        """
        # Find eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eig(W)

        # Find the eigenvector corresponding to eigenvalue ≈ 0
        zero_idx = np.argmin(np.abs(eigenvalues))

        P_ss = np.real(eigenvectors[:, zero_idx])

        # Normalize to sum to 1
        P_ss = P_ss / np.sum(P_ss)

        # Ensure non-negative (should be, but fix numerical issues)
        P_ss = np.maximum(P_ss, 0)
        P_ss = P_ss / np.sum(P_ss)

        return P_ss

    def compute_photon_emission(
        self,
        populations: np.ndarray,
        radiative_rates: np.ndarray,
        excited_state_indices: np.ndarray
    ) -> np.ndarray:
        """Compute photon emission rate over time.

        The photon emission rate is proportional to the excited state
        populations weighted by their radiative decay rates.

        Parameters
        ----------
        populations : np.ndarray
            Population array with shape (n_times, n_states).
        radiative_rates : np.ndarray
            Radiative decay rates for excited states in MHz.
        excited_state_indices : np.ndarray
            Indices of the excited states in the population vector.

        Returns
        -------
        emission_rate : np.ndarray
            Photon emission rate at each time point (arbitrary units).
        """
        n_es = len(excited_state_indices)
        if len(radiative_rates) != n_es:
            raise ValueError(
                f"Number of radiative rates ({len(radiative_rates)}) must match "
                f"number of excited states ({n_es})"
            )

        emission_rate = np.zeros(len(populations))

        for i, es_idx in enumerate(excited_state_indices):
            emission_rate += populations[:, es_idx] * radiative_rates[i]

        return emission_rate
