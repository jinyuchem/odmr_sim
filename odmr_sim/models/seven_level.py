"""
7-level rate equation model for NV-like defect systems.

The model includes:
- 3 ground states: GS|0>, GS|->, GS|+>
- 3 excited states: ES|0>, ES|->, ES|+>
- 1 metastable singlet state
"""

import numpy as np
from typing import Optional
from odmr_sim.models.base import RateModel


class SevenLevelModel(RateModel):
    """7-level rate equation model for NV-like spin systems.

    State indexing:
        0: GS|0>  - Ground state, m_s = 0
        1: GS|->  - Ground state, m_s = -1 (or -)
        2: GS|+>  - Ground state, m_s = +1 (or +)
        3: ES|0>  - Excited state, m_s = 0
        4: ES|->  - Excited state, m_s = -1 (or -)
        5: ES|+>  - Excited state, m_s = +1 (or +)
        6: SS     - Singlet (metastable) state

    Rate parameters (in MHz):
        k41, k52, k63 : Radiative decay (ES -> GS)
        k47, k57, k67 : Upper ISC (ES -> Singlet)
        k71, k72, k73 : Lower ISC (Singlet -> GS)

    Dynamic parameters (passed to build_rate_matrix):
        gamma     : Optical excitation rate (GS -> ES)
        kmw_minus : Microwave rate for |0> <-> |-> transition
        kmw_plus  : Microwave rate for |0> <-> |+> transition

    Parameters
    ----------
    k41 : float, optional
        Radiative decay rate ES|0> -> GS|0> (default: 62.5 MHz)
    k52 : float, optional
        Radiative decay rate ES|-> -> GS|-> (default: same as k41)
    k63 : float, optional
        Radiative decay rate ES|+> -> GS|+> (default: same as k41)
    k47 : float, optional
        Upper ISC rate ES|0> -> SS (default: 4.4 MHz)
    k57 : float, optional
        Upper ISC rate ES|-> -> SS (default: 0.005 MHz)
    k67 : float, optional
        Upper ISC rate ES|+> -> SS (default: 44.1 MHz)
    k71 : float, optional
        Lower ISC rate SS -> GS|0> (default: 2336 MHz)
    k72 : float, optional
        Lower ISC rate SS -> GS|-> (default: 3.1 MHz)
    k73 : float, optional
        Lower ISC rate SS -> GS|+> (default: 0.001 MHz)
    use_degenerate_labels : bool, optional
        If True, use |+/-1> labels (degenerate case). If False, use |+/-> labels (split case).
        Default is False.

    References
    ----------
    - L. Robledo et al., New J. Phys. 13, 025013 (2011)
    - S. Ahmadi et al., Phys. Rev. Applied 8, 034001 (2017)
    - Y. Masuyama et al., Sci. Rep. 14, 18135 (2024)

    Examples
    --------
    >>> model = SevenLevelModel()  # Default NV-like rates
    >>> W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)
    >>> print(W.shape)
    (7, 7)
    """

    # State labels
    LABELS_DEGENERATE = [r"GS$|0\rangle$", r"GS$|-1\rangle$", r"GS$|+1\rangle$", r"ES$|0\rangle$", r"ES$|-1\rangle$", r"ES$|+1\rangle$", "SS"]
    LABELS_SPLIT = [r"GS$|0\rangle$", r"GS$|-\rangle$", r"GS$|+\rangle$", r"ES$|0\rangle$", r"ES$|-\rangle$", r"ES$|+\rangle$", "SS"]

    # State indices as class constants for convenience
    GS_0 = 0
    GS_MINUS = 1
    GS_PLUS = 2
    ES_0 = 3
    ES_MINUS = 4
    ES_PLUS = 5
    SINGLET = 6

    def __init__(
        self,
        k41: float = 62.5,
        k52: Optional[float] = None,
        k63: Optional[float] = None,
        k47: float = 4.4,
        k57: float = 0.005,
        k67: float = 44.1,
        k71: float = 2336.0,
        k72: float = 3.1,
        k73: float = 0.001,
        use_degenerate_labels: bool = False
    ):
        labels = self.LABELS_DEGENERATE if use_degenerate_labels else self.LABELS_SPLIT
        super().__init__(n_states=7, state_labels=labels)

        # Store rate parameters
        self.k41 = k41
        self.k52 = k52 if k52 is not None else k41
        self.k63 = k63 if k63 is not None else k41
        self.k47 = k47
        self.k57 = k57
        self.k67 = k67
        self.k71 = k71
        self.k72 = k72
        self.k73 = k73

        # Set up fixed rates
        self._setup_fixed_rates()

        # Set up dynamic rates (gamma, kmw)
        self._setup_dynamic_rates()

    def _setup_fixed_rates(self) -> None:
        """Set up the fixed transition rates."""
        # Radiative decay (ES -> GS)
        self.set_rate(self.ES_0, self.GS_0, self.k41)
        self.set_rate(self.ES_MINUS, self.GS_MINUS, self.k52)
        self.set_rate(self.ES_PLUS, self.GS_PLUS, self.k63)

        # Upper ISC (ES -> Singlet)
        self.set_rate(self.ES_0, self.SINGLET, self.k47)
        self.set_rate(self.ES_MINUS, self.SINGLET, self.k57)
        self.set_rate(self.ES_PLUS, self.SINGLET, self.k67)

        # Lower ISC (Singlet -> GS)
        self.set_rate(self.SINGLET, self.GS_0, self.k71)
        self.set_rate(self.SINGLET, self.GS_MINUS, self.k72)
        self.set_rate(self.SINGLET, self.GS_PLUS, self.k73)

    def _setup_dynamic_rates(self) -> None:
        """Set up dynamic rate specifications for gamma and kmw."""
        # Excitation (GS -> ES), controlled by 'gamma'
        self.add_dynamic_rate('gamma', self.GS_0, self.ES_0)
        self.add_dynamic_rate('gamma', self.GS_MINUS, self.ES_MINUS)
        self.add_dynamic_rate('gamma', self.GS_PLUS, self.ES_PLUS)

        # Microwave driving GS|0> <-> GS|->, controlled by 'kmw_minus'
        self.add_dynamic_rate('kmw_minus', self.GS_0, self.GS_MINUS)
        self.add_dynamic_rate('kmw_minus', self.GS_MINUS, self.GS_0)

        # Microwave driving GS|0> <-> GS|+>, controlled by 'kmw_plus'
        self.add_dynamic_rate('kmw_plus', self.GS_0, self.GS_PLUS)
        self.add_dynamic_rate('kmw_plus', self.GS_PLUS, self.GS_0)

    def build_rate_matrix(
        self,
        gamma: float = 0.0,
        kmw_minus: float = 0.0,
        kmw_plus: float = 0.0
    ) -> np.ndarray:
        """Build the 7×7 rate matrix.

        Parameters
        ----------
        gamma : float
            Optical excitation rate in MHz.
        kmw_minus : float
            Microwave rate for |0> <-> |-> transition in MHz.
        kmw_plus : float
            Microwave rate for |0> <-> |+> transition in MHz.

        Returns
        -------
        W : np.ndarray
            The 7×7 rate matrix in units of 1/s.
        """
        return super().build_rate_matrix(
            gamma=gamma,
            kmw_minus=kmw_minus,
            kmw_plus=kmw_plus
        )

    def get_ground_state_mixed(self) -> np.ndarray:
        """Get initial state with equal population in ground states.

        Returns
        -------
        P0 : np.ndarray
            Initial probability vector with P = [1/3, 1/3, 1/3, 0, 0, 0, 0].
        """
        P0 = np.zeros(7)
        P0[self.GS_0] = 1/3
        P0[self.GS_MINUS] = 1/3
        P0[self.GS_PLUS] = 1/3
        return P0

    def get_rate_summary(self) -> str:
        """Get a formatted summary of all rate parameters."""
        return (
            f"Radiative rates: k41={self.k41}, k52={self.k52}, k63={self.k63} MHz\n"
            f"Upper ISC rates: k47={self.k47}, k57={self.k57}, k67={self.k67} MHz\n"
            f"Lower ISC rates: k71={self.k71}, k72={self.k72}, k73={self.k73} MHz"
        )
