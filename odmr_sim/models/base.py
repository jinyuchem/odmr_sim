"""
Base rate equation model with configurable N-level support.
"""

import numpy as np
from typing import Optional, Dict, Tuple, List


class RateModel:
    """Configurable N-level rate equation model.

    This base class allows you to define arbitrary rate equation models
    with any number of states and custom transition rates.

    Parameters
    ----------
    n_states : int
        Number of states in the model.
    state_labels : list[str], optional
        Human-readable labels for each state. If not provided,
        defaults to ["State 0", "State 1", ...].

    Examples
    --------
    Create a custom 5-level model:

    >>> model = RateModel(n_states=5, state_labels=["GS", "ES1", "ES2", "MS1", "MS2"])
    >>> model.set_rates({
    ...     (0, 1): 10.0,   # GS -> ES1
    ...     (1, 0): 62.5,   # ES1 -> GS
    ...     (1, 3): 5.0,    # ES1 -> MS1
    ...     (3, 0): 100.0,  # MS1 -> GS
    ... })
    >>> W = model.build_rate_matrix()

    Notes
    -----
    - Rates are specified in MHz
    - The rate matrix is built using the convention W[j,i] = rate from state i to state j
    - Diagonal elements are automatically computed to conserve probability
    """

    def __init__(
        self,
        n_states: int,
        state_labels: Optional[List[str]] = None
    ):
        if n_states < 2:
            raise ValueError("Model must have at least 2 states")

        self.n_states = n_states
        self.state_labels = state_labels or [f"State {i}" for i in range(n_states)]

        if len(self.state_labels) != n_states:
            raise ValueError(
                f"Number of labels ({len(self.state_labels)}) must match "
                f"number of states ({n_states})"
            )

        # Dictionary of (from_state, to_state) -> rate (in MHz)
        self._rates: Dict[Tuple[int, int], float] = {}

        # Dictionary for dynamic rate parameters (e.g., excitation rate, MW rate)
        self._dynamic_rate_specs: Dict[str, List[Tuple[int, int, float]]] = {}

    def set_rate(self, from_state: int, to_state: int, rate: float) -> None:
        """Set a fixed transition rate from one state to another.

        Parameters
        ----------
        from_state : int
            Index of the source state (0-indexed).
        to_state : int
            Index of the target state (0-indexed).
        rate : float
            Transition rate in MHz.
        """
        self._validate_state_index(from_state, "from_state")
        self._validate_state_index(to_state, "to_state")

        if from_state == to_state:
            raise ValueError("Cannot set rate from a state to itself")

        self._rates[(from_state, to_state)] = rate

    def set_rates(self, rates: Dict[Tuple[int, int], float]) -> None:
        """Set multiple transition rates at once.

        Parameters
        ----------
        rates : dict
            Dictionary mapping (from_state, to_state) tuples to rate values in MHz.
        """
        for (from_state, to_state), rate in rates.items():
            self.set_rate(from_state, to_state, rate)

    def add_dynamic_rate(
        self,
        param_name: str,
        from_state: int,
        to_state: int,
        coefficient: float = 1.0
    ) -> None:
        """Add a dynamic rate that will be multiplied by a parameter at build time.

        This allows rates like excitation (Gamma) or microwave (k_MW) to be set
        when building the rate matrix rather than at model definition.

        Parameters
        ----------
        param_name : str
            Name of the parameter (e.g., 'gamma', 'kmw_minus').
        from_state : int
            Index of the source state.
        to_state : int
            Index of the target state.
        coefficient : float, optional
            Multiplier for the parameter value. Default is 1.0.

        Examples
        --------
        >>> model.add_dynamic_rate('gamma', from_state=0, to_state=3)  # GS|0> -> ES|0>
        >>> W = model.build_rate_matrix(gamma=0.1)  # Gamma = 0.1 MHz
        """
        self._validate_state_index(from_state, "from_state")
        self._validate_state_index(to_state, "to_state")

        if param_name not in self._dynamic_rate_specs:
            self._dynamic_rate_specs[param_name] = []

        self._dynamic_rate_specs[param_name].append((from_state, to_state, coefficient))

    def build_rate_matrix(self, **kwargs) -> np.ndarray:
        """Build the N×N rate matrix W.

        Parameters
        ----------
        **kwargs : float
            Dynamic rate parameters (e.g., gamma=0.1, kmw_minus=1.0).
            These are multiplied by their coefficients and added to
            the rate matrix.

        Returns
        -------
        W : np.ndarray
            The rate matrix with shape (n_states, n_states) in units of 1/s.
            Convention: W[j, i] = rate from state i to state j.
            Diagonal elements ensure probability conservation (row sums = 0).
        """
        W = np.zeros((self.n_states, self.n_states))

        # Add fixed rates
        for (from_state, to_state), rate in self._rates.items():
            W[to_state, from_state] += rate

        # Add dynamic rates
        for param_name, transitions in self._dynamic_rate_specs.items():
            param_value = kwargs.get(param_name, 0.0)
            for from_state, to_state, coefficient in transitions:
                W[to_state, from_state] += param_value * coefficient

        # Set diagonal elements for probability conservation
        # Diagonal = negative sum of outgoing rates from that state
        for i in range(self.n_states):
            W[i, i] = -np.sum(W[:, i]) + W[i, i]

        # Convert from MHz to Hz (1/s)
        return W * 1e6

    def get_initial_state(self, state_index: int) -> np.ndarray:
        """Get initial probability vector with population in a single state.

        Parameters
        ----------
        state_index : int
            Index of the state to populate.

        Returns
        -------
        P0 : np.ndarray
            Initial probability vector.
        """
        self._validate_state_index(state_index, "state_index")
        P0 = np.zeros(self.n_states)
        P0[state_index] = 1.0
        return P0

    def get_mixed_initial_state(self, populations: Dict[int, float]) -> np.ndarray:
        """Get initial probability vector with specified populations.

        Parameters
        ----------
        populations : dict
            Dictionary mapping state indices to populations.
            Values should sum to 1.0 for a normalized state.

        Returns
        -------
        P0 : np.ndarray
            Initial probability vector.
        """
        P0 = np.zeros(self.n_states)
        for state_index, population in populations.items():
            self._validate_state_index(state_index, "state_index")
            P0[state_index] = population
        return P0

    def _validate_state_index(self, index: int, name: str) -> None:
        """Validate that a state index is within bounds."""
        if not (0 <= index < self.n_states):
            raise ValueError(
                f"{name}={index} is out of bounds for model with "
                f"{self.n_states} states (valid: 0 to {self.n_states - 1})"
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(n_states={self.n_states}, "
            f"labels={self.state_labels})"
        )
