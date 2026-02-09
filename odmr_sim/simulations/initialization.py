"""
Initialization (spin polarization) simulations.
"""

import numpy as np
from typing import Optional, List, Union
from odmr_sim.models.base import RateModel
from odmr_sim.models.seven_level import SevenLevelModel
from odmr_sim.solvers.rate_solver import RateSolver


class InitializationSimulation:
    """Simulation for spin initialization (polarization) dynamics.

    Simulates how the spin system evolves under optical pumping,
    typically starting from a thermally mixed state.

    Parameters
    ----------
    model : RateModel or SevenLevelModel
        The rate equation model to use.

    Examples
    --------
    >>> from odmr_sim.models import SevenLevelModel
    >>> from odmr_sim.simulations import InitializationSimulation
    >>>
    >>> model = SevenLevelModel()
    >>> sim = InitializationSimulation(model)
    >>> result = sim.run(gamma=0.1, t_max=1e-1, n_points=1000)
    >>> print(f"Final GS|0> population: {result['populations'][-1, 0]:.3f}")
    """

    def __init__(self, model: Union[RateModel, SevenLevelModel]):
        self.model = model
        self.solver = RateSolver()

    def run(
        self,
        gamma: float = 0.1,
        kmw_minus: float = 0.0,
        kmw_plus: float = 0.0,
        P0: Optional[np.ndarray] = None,
        t_max: float = 1e-1,
        t_min: float = 1e-9,
        n_points: int = 1000,
        use_log_time: bool = True,
        **model_kwargs
    ) -> dict:
        """Run initialization simulation.

        Parameters
        ----------
        gamma : float
            Optical excitation rate in MHz.
        kmw_minus : float
            Microwave rate for |0> <-> |-> in MHz.
        kmw_plus : float
            Microwave rate for |0> <-> |+> in MHz.
        P0 : np.ndarray, optional
            Initial population. If None, uses mixed ground state for
            SevenLevelModel or equal population for generic RateModel.
        t_max : float
            Maximum simulation time in seconds.
        t_min : float
            Minimum simulation time in seconds (for log scale).
        n_points : int
            Number of time points.
        use_log_time : bool
            If True, use logarithmic time spacing.
        **model_kwargs
            Additional keyword arguments passed to model.build_rate_matrix().

        Returns
        -------
        result : dict
            Dictionary containing:
            - 't': time array in seconds
            - 't_ns': time array in nanoseconds
            - 'populations': population array (n_times, n_states)
            - 'labels': state labels
            - 'model': the model used
            - 'params': simulation parameters
        """
        # Build rate matrix
        if isinstance(self.model, SevenLevelModel):
            W = self.model.build_rate_matrix(
                gamma=gamma,
                kmw_minus=kmw_minus,
                kmw_plus=kmw_plus
            )
            if P0 is None:
                P0 = self.model.get_ground_state_mixed()
        else:
            W = self.model.build_rate_matrix(
                gamma=gamma,
                kmw_minus=kmw_minus,
                kmw_plus=kmw_plus,
                **model_kwargs
            )
            if P0 is None:
                # Equal population in first half of states (assume ground states)
                n_gs = self.model.n_states // 2
                P0 = np.zeros(self.model.n_states)
                P0[:n_gs] = 1.0 / n_gs

        # Create time array
        if use_log_time:
            t_eval = np.logspace(np.log10(t_min), np.log10(t_max), n_points)
        else:
            t_eval = np.linspace(t_min, t_max, n_points)

        # Solve
        populations = self.solver.solve_expm(W, P0, t_eval)

        return {
            't': t_eval,
            't_ns': t_eval * 1e9,
            'populations': populations,
            'labels': self.model.state_labels,
            'model': self.model,
            'params': {
                'gamma': gamma,
                'kmw_minus': kmw_minus,
                'kmw_plus': kmw_plus,
                'P0': P0,
            }
        }

    def run_sweep_gamma(
        self,
        gammas: List[float],
        kmw_minus: float = 0.0,
        kmw_plus: float = 0.0,
        P0: Optional[np.ndarray] = None,
        t_max: float = 1e-1,
        t_min: float = 1e-9,
        n_points: int = 1000,
    ) -> List[dict]:
        """Run initialization for multiple excitation rates.

        Parameters
        ----------
        gammas : list of float
            List of excitation rates to sweep.

        Returns
        -------
        results : list of dict
            List of result dictionaries, one for each gamma value.
        """
        results = []
        for gamma in gammas:
            result = self.run(
                gamma=gamma,
                kmw_minus=kmw_minus,
                kmw_plus=kmw_plus,
                P0=P0,
                t_max=t_max,
                t_min=t_min,
                n_points=n_points,
            )
            results.append(result)
        return results
