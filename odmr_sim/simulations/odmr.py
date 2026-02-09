"""
ODMR contrast simulations.
"""

import numpy as np
from typing import Optional, Union, Tuple
from odmr_sim.models.base import RateModel
from odmr_sim.models.seven_level import SevenLevelModel
from odmr_sim.solvers.rate_solver import RateSolver


class ODMRSimulation:
    """Simulation for ODMR (Optically Detected Magnetic Resonance) contrast.

    Computes the change in photoluminescence intensity when microwave
    fields are applied at resonance with spin transitions.

    Parameters
    ----------
    model : RateModel or SevenLevelModel
        The rate equation model to use.

    Examples
    --------
    >>> from odmr_sim.models import SevenLevelModel
    >>> from odmr_sim.simulations import ODMRSimulation
    >>>
    >>> model = SevenLevelModel()
    >>> sim = ODMRSimulation(model)
    >>> contrast = sim.compute_contrast(gamma=0.1, kmw_minus=1.0)
    >>> print(f"ODMR contrast: {contrast*100:.2f}%")
    """

    def __init__(self, model: Union[RateModel, SevenLevelModel]):
        self.model = model
        self.solver = RateSolver()

    def compute_contrast(
        self,
        gamma: float = 0.1,
        kmw_minus: float = 0.0,
        kmw_plus: float = 0.0,
        P0: Optional[np.ndarray] = None,
        t_integration: float = 1e-3,
        method: str = 'steady_state'
    ) -> float:
        """Compute ODMR contrast for given parameters.

        Contrast is defined as: (I_noMW - I_MW) / I_noMW
        where I is the integrated photoluminescence intensity.

        Parameters
        ----------
        gamma : float
            Optical excitation rate in MHz.
        kmw_minus : float
            Microwave rate for |0> <-> |-> in MHz.
        kmw_plus : float
            Microwave rate for |0> <-> |+> in MHz.
        P0 : np.ndarray, optional
            Initial population. If None, uses steady state without MW.
        t_integration : float
            Integration time in seconds (for time-integrated method).
        method : str
            'steady_state': Use null-space steady state (t -> infinity)
            'transient': Use final ES population at t_integration
            'time_integrated': Use trapz-integrated fluorescence over [0, t_integration]

        Returns
        -------
        contrast : float
            ODMR contrast (fractional change in PL).
        """
        if not isinstance(self.model, SevenLevelModel):
            raise ValueError("compute_contrast requires SevenLevelModel")

        if method == 'steady_state':
            return self._contrast_steady_state(gamma, kmw_minus, kmw_plus)
        elif method == 'transient':
            return self._contrast_transient(
                gamma, kmw_minus, kmw_plus, P0, t_integration
            )
        elif method == 'time_integrated':
            return self._contrast_time_integrated(
                gamma, kmw_minus, kmw_plus, P0, t_integration
            )
        else:
            raise ValueError(f"Unknown method: {method}")

    def _contrast_steady_state(
        self,
        gamma: float,
        kmw_minus: float,
        kmw_plus: float
    ) -> float:
        """Compute contrast using steady-state populations."""
        # Without microwave
        W_nomw = self.model.build_rate_matrix(
            gamma=gamma, kmw_minus=0.0, kmw_plus=0.0
        )
        P_ss_nomw = self.solver.solve_steady_state(W_nomw)

        # With microwave
        W_mw = self.model.build_rate_matrix(
            gamma=gamma, kmw_minus=kmw_minus, kmw_plus=kmw_plus
        )
        P_ss_mw = self.solver.solve_steady_state(W_mw)

        # PL intensity ∝ sum of excited state populations × radiative rates
        # For simplicity, assume equal radiative rates
        I_nomw = (
            P_ss_nomw[SevenLevelModel.ES_0] * self.model.k41 +
            P_ss_nomw[SevenLevelModel.ES_MINUS] * self.model.k52 +
            P_ss_nomw[SevenLevelModel.ES_PLUS] * self.model.k63
        )
        I_mw = (
            P_ss_mw[SevenLevelModel.ES_0] * self.model.k41 +
            P_ss_mw[SevenLevelModel.ES_MINUS] * self.model.k52 +
            P_ss_mw[SevenLevelModel.ES_PLUS] * self.model.k63
        )

        if I_nomw == 0:
            return 0.0

        return (I_nomw - I_mw) / I_nomw

    def _contrast_transient(
        self,
        gamma: float,
        kmw_minus: float,
        kmw_plus: float,
        P0: Optional[np.ndarray],
        t_integration: float
    ) -> float:
        """Compute contrast using final ES population at t_integration.

        Evolve from initial state P0 to time t_integration and compare final ES populations.
        """
        if P0 is None:
            # Start from mixed ground state (1/3, 1/3, 1/3, 0, 0, 0, 0)
            P0 = self.model.get_ground_state_mixed()

        # Evolve to final time
        t_eval = np.array([0.0, t_integration])

        # Without microwave
        W_nomw = self.model.build_rate_matrix(
            gamma=gamma, kmw_minus=0.0, kmw_plus=0.0
        )
        pop_nomw = self.solver.solve_expm(W_nomw, P0, t_eval)

        # With microwave
        W_mw = self.model.build_rate_matrix(
            gamma=gamma, kmw_minus=kmw_minus, kmw_plus=kmw_plus
        )
        pop_mw = self.solver.solve_expm(W_mw, P0, t_eval)

        # Get final ES populations (at t_integration)
        ES_nomw = (
            pop_nomw[-1, SevenLevelModel.ES_0] +
            pop_nomw[-1, SevenLevelModel.ES_MINUS] +
            pop_nomw[-1, SevenLevelModel.ES_PLUS]
        )
        ES_mw = (
            pop_mw[-1, SevenLevelModel.ES_0] +
            pop_mw[-1, SevenLevelModel.ES_MINUS] +
            pop_mw[-1, SevenLevelModel.ES_PLUS]
        )

        if ES_nomw == 0:
            return 0.0

        # Contrast = 1 - ES_mw / ES_nomw
        return 1.0 - ES_mw / ES_nomw

    def _contrast_time_integrated(
        self,
        gamma: float,
        kmw_minus: float,
        kmw_plus: float,
        P0: Optional[np.ndarray],
        t_integration: float
    ) -> float:
        """Compute contrast using time-integrated fluorescence."""
        if P0 is None:
            # Start from mixed ground state
            P0 = self.model.get_ground_state_mixed()

        # Time array
        t_eval = np.linspace(0, t_integration, 10000)

        # Without microwave
        W_nomw = self.model.build_rate_matrix(
            gamma=gamma, kmw_minus=0.0, kmw_plus=0.0
        )
        pop_nomw = self.solver.solve_expm(W_nomw, P0, t_eval)

        # With microwave
        W_mw = self.model.build_rate_matrix(
            gamma=gamma, kmw_minus=kmw_minus, kmw_plus=kmw_plus
        )
        pop_mw = self.solver.solve_expm(W_mw, P0, t_eval)

        # Compute time-integrated fluorescence
        def compute_pl(populations):
            return (
                populations[:, SevenLevelModel.ES_0] * self.model.k41 +
                populations[:, SevenLevelModel.ES_MINUS] * self.model.k52 +
                populations[:, SevenLevelModel.ES_PLUS] * self.model.k63
            )

        I_nomw = np.trapz(compute_pl(pop_nomw), t_eval)
        I_mw = np.trapz(compute_pl(pop_mw), t_eval)

        if I_nomw == 0:
            return 0.0

        return (I_nomw - I_mw) / I_nomw

    def sweep_gamma(
        self,
        gammas: np.ndarray,
        kmw_minus: float = 0.0,
        kmw_plus: float = 0.0,
        method: str = 'steady_state'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Sweep excitation rate and compute contrast.

        Parameters
        ----------
        gammas : np.ndarray
            Array of excitation rates to sweep (in MHz).
        kmw_minus, kmw_plus : float
            Microwave rates.
        method : str
            Contrast computation method.

        Returns
        -------
        gammas : np.ndarray
            Input gamma values.
        contrasts : np.ndarray
            Computed contrasts.
        """
        contrasts = np.array([
            self.compute_contrast(gamma, kmw_minus, kmw_plus, method=method)
            for gamma in gammas
        ])
        return gammas, contrasts

    def run_spectrum(
        self,
        gamma: float = 0.1,
        freq_center: float = 3.0,  # GHz
        freq_width: float = 0.5,   # GHz
        n_points: int = 101,
        peak_freq_minus: Optional[float] = None,
        peak_freq_plus: Optional[float] = None,
        linewidth: float = 0.05,   # GHz
        kmw_amplitude: float = 1.0,  # MHz
    ) -> dict:
        """Simulate ODMR spectrum.

        This simulates what you would measure when sweeping the microwave
        frequency across the spin resonances.

        Parameters
        ----------
        gamma : float
            Optical excitation rate in MHz.
        freq_center : float
            Center frequency of sweep in GHz.
        freq_width : float
            Width of frequency sweep in GHz.
        n_points : int
            Number of frequency points.
        peak_freq_minus : float, optional
            Resonance frequency for |-> transition in GHz.
        peak_freq_plus : float, optional
            Resonance frequency for |+> transition in GHz.
        linewidth : float
            Lorentzian linewidth in GHz.
        kmw_amplitude : float
            Maximum microwave rate at resonance in MHz.

        Returns
        -------
        result : dict
            Dictionary with 'frequencies' and 'contrast' arrays.
        """
        frequencies = np.linspace(
            freq_center - freq_width,
            freq_center + freq_width,
            n_points
        )

        # Default peak positions (symmetric around center)
        if peak_freq_minus is None:
            peak_freq_minus = freq_center - freq_width / 4
        if peak_freq_plus is None:
            peak_freq_plus = freq_center + freq_width / 4

        contrasts = []

        for freq in frequencies:
            # Lorentzian profiles for MW coupling efficiency
            lorentz_minus = self._lorentzian(freq, peak_freq_minus, linewidth)
            lorentz_plus = self._lorentzian(freq, peak_freq_plus, linewidth)

            kmw_minus = kmw_amplitude * lorentz_minus
            kmw_plus = kmw_amplitude * lorentz_plus

            contrast = self.compute_contrast(
                gamma=gamma,
                kmw_minus=kmw_minus,
                kmw_plus=kmw_plus,
            )
            contrasts.append(contrast)

        return {
            'frequencies': frequencies,
            'contrast': np.array(contrasts),
            'params': {
                'gamma': gamma,
                'peak_freq_minus': peak_freq_minus,
                'peak_freq_plus': peak_freq_plus,
                'linewidth': linewidth,
                'kmw_amplitude': kmw_amplitude,
            }
        }

    @staticmethod
    def _lorentzian(x: float, x0: float, gamma: float) -> float:
        """Lorentzian function normalized to peak = 1."""
        return gamma**2 / ((x - x0)**2 + gamma**2)
