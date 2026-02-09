"""Tests for ODMR simulation solvers."""

import pytest
import numpy as np
from odmr_sim.models import SevenLevelModel
from odmr_sim.solvers import RateSolver


class TestRateSolver:
    """Tests for the RateSolver class."""

    @pytest.fixture
    def model_and_solver(self):
        """Create a model and solver for testing."""
        model = SevenLevelModel()
        solver = RateSolver()
        return model, solver

    def test_solve_expm_shape(self, model_and_solver):
        """Test matrix exponential solver returns correct shape."""
        model, solver = model_and_solver
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)
        P0 = model.get_ground_state_mixed()
        t_eval = np.linspace(0, 1, 100)

        populations = solver.solve_expm(W, P0, t_eval)
        assert populations.shape == (100, 7)

    def test_solve_expm_conservation(self, model_and_solver):
        """Test probability is conserved during evolution."""
        model, solver = model_and_solver
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)
        P0 = model.get_ground_state_mixed()
        t_eval = np.linspace(0, 1, 100)

        populations = solver.solve_expm(W, P0, t_eval)
        total_prob = np.sum(populations, axis=1)
        np.testing.assert_array_almost_equal(total_prob, 1.0)

    def test_solve_expm_initial_condition(self, model_and_solver):
        """Test initial condition is preserved at t=0."""
        model, solver = model_and_solver
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)
        P0 = model.get_ground_state_mixed()
        t_eval = np.linspace(0, 1, 100)

        populations = solver.solve_expm(W, P0, t_eval)
        np.testing.assert_array_almost_equal(populations[0], P0)

    def test_steady_state_normalization(self, model_and_solver):
        """Test steady state is normalized."""
        model, solver = model_and_solver
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)

        P_ss = solver.solve_steady_state(W)
        np.testing.assert_almost_equal(np.sum(P_ss), 1.0)

    def test_steady_state_is_null_space(self, model_and_solver):
        """Test steady state satisfies W @ P_ss = 0."""
        model, solver = model_and_solver
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)

        P_ss = solver.solve_steady_state(W)
        residual = W @ P_ss
        np.testing.assert_array_almost_equal(residual, 0.0, decimal=6)

    def test_steady_state_positive(self, model_and_solver):
        """Test steady state has non-negative populations."""
        model, solver = model_and_solver
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)

        P_ss = solver.solve_steady_state(W)
        assert np.all(P_ss >= -1e-10)  # Allow small numerical errors
