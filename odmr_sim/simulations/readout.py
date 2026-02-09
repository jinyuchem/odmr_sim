"""
Optical readout simulations.
"""

import numpy as np
from typing import Optional, List, Union, Dict
from odmr_sim.models.base import RateModel
from odmr_sim.models.seven_level import SevenLevelModel
from odmr_sim.solvers.rate_solver import RateSolver


class ReadoutSimulation:
    """Simulation for optical readout dynamics.

    Simulates photoluminescence readout from different initial spin states.
    This is used to determine spin-state-dependent fluorescence contrast.

    Parameters
    ----------
    model : RateModel or SevenLevelModel
        The rate equation model to use.

    Examples
    --------
    >>> from odmr_sim.models import SevenLevelModel
    >>> from odmr_sim.simulations import ReadoutSimulation
    >>>
    >>> model = SevenLevelModel()
    >>> sim = ReadoutSimulation(model)
    >>> result = sim.run(gamma=12.8, initial_state='gs0')
    >>> print(f"Peak ES population: {result['es_total'].max():.3f}")
    """

    def __init__(self, model: Union[RateModel, SevenLevelModel]):
        self.model = model
        self.solver = RateSolver()

    def run(
        self,
        gamma: float = 12.8,
        initial_state: Union[str, int, np.ndarray] = 'gs0',
        t_max: float = 1e-5,
        t_min: float = 0.0,
        n_points: int = 1000,
        use_log_time: bool = False,
        **model_kwargs
    ) -> dict:
        """Run readout simulation from a specific initial state.

        Parameters
        ----------
        gamma : float
            Optical excitation rate in MHz. Typically higher for readout.
        initial_state : str, int, or np.ndarray
            Initial state specification:
            - 'gs0', 'gs_minus', 'gs_plus': Pure ground states (7-level model)
            - int: State index to populate
            - np.ndarray: Custom initial population vector
        t_max : float
            Maximum simulation time in seconds.
        t_min : float
            Minimum simulation time in seconds.
        n_points : int
            Number of time points.
        use_log_time : bool
            If True, use logarithmic time spacing.
        **model_kwargs
            Additional keyword arguments for model.build_rate_matrix().

        Returns
        -------
        result : dict
            Dictionary containing:
            - 't': time array in seconds
            - 't_ns': time array in nanoseconds
            - 'populations': population array
            - 'es_total': total excited state population (if applicable)
            - 'initial_state': the initial state used
        """
        # Build rate matrix (no MW during readout)
        if isinstance(self.model, SevenLevelModel):
            W = self.model.build_rate_matrix(
                gamma=gamma,
                kmw_minus=0.0,
                kmw_plus=0.0
            )
        else:
            W = self.model.build_rate_matrix(gamma=gamma, **model_kwargs)

        # Parse initial state
        P0 = self._parse_initial_state(initial_state)

        # Create time array
        if use_log_time and t_min > 0:
            t_eval = np.logspace(np.log10(t_min), np.log10(t_max), n_points)
        else:
            t_eval = np.linspace(t_min, t_max, n_points)

        # Solve
        populations = self.solver.solve_expm(W, P0, t_eval)

        # Calculate excited state total (for 7-level model)
        es_total = None
        if isinstance(self.model, SevenLevelModel):
            es_total = (
                populations[:, SevenLevelModel.ES_0] +
                populations[:, SevenLevelModel.ES_MINUS] +
                populations[:, SevenLevelModel.ES_PLUS]
            )

        return {
            't': t_eval,
            't_ns': t_eval * 1e9,
            'populations': populations,
            'es_total': es_total,
            'labels': self.model.state_labels,
            'initial_state': initial_state,
            'params': {
                'gamma': gamma,
                'P0': P0,
            }
        }

    def run_comparison(
        self,
        gamma: float = 12.8,
        t_max: float = 1e-5,
        n_points: int = 1000,
    ) -> Dict[str, dict]:
        """Run readout for all three ground states to compare.

        For 7-level model, runs readout starting from GS|0>, GS|->, and GS|+>.

        Returns
        -------
        results : dict
            Dictionary mapping initial state names to result dicts.
        """
        if not isinstance(self.model, SevenLevelModel):
            raise ValueError("run_comparison is only available for SevenLevelModel")

        states = ['gs0', 'gs_minus', 'gs_plus']
        results = {}

        for state in states:
            results[state] = self.run(
                gamma=gamma,
                initial_state=state,
                t_max=t_max,
                n_points=n_points,
            )

        return results

    def _parse_initial_state(self, initial_state) -> np.ndarray:
        """Convert initial state specification to population vector."""
        if isinstance(initial_state, np.ndarray):
            return initial_state

        if isinstance(initial_state, int):
            return self.model.get_initial_state(initial_state)

        if isinstance(initial_state, str):
            if not isinstance(self.model, SevenLevelModel):
                raise ValueError(
                    "String initial states only supported for SevenLevelModel"
                )

            state_map = {
                'gs0': SevenLevelModel.GS_0,
                'gs_0': SevenLevelModel.GS_0,
                'gs_minus': SevenLevelModel.GS_MINUS,
                'gs-': SevenLevelModel.GS_MINUS,
                'gs_plus': SevenLevelModel.GS_PLUS,
                'gs+': SevenLevelModel.GS_PLUS,
                '100': SevenLevelModel.GS_0,
                '010': SevenLevelModel.GS_MINUS,
                '001': SevenLevelModel.GS_PLUS,
            }

            state_lower = initial_state.lower()
            if state_lower not in state_map:
                raise ValueError(f"Unknown initial state: {initial_state}")

            return self.model.get_initial_state(state_map[state_lower])

        raise ValueError(f"Invalid initial_state type: {type(initial_state)}")
