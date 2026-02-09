"""Tests for ODMR simulations."""

import pytest
import numpy as np
from odmr_sim.models import get_preset
from odmr_sim.simulations import ODMRSimulation, InitializationSimulation


class TestODMRSimulation:
    """Tests for ODMRSimulation class."""

    @pytest.fixture
    def simulation(self):
        """Create simulation for testing."""
        model = get_preset('nv_bulk')
        return ODMRSimulation(model)

    def test_compute_contrast_steady_state(self, simulation):
        """Test steady-state contrast computation."""
        contrast = simulation.compute_contrast(
            gamma=0.1, kmw_minus=1.0, kmw_plus=1.0,
            method='steady_state'
        )
        assert isinstance(contrast, float)
        assert -1 <= contrast <= 1  # Contrast should be in valid range

    def test_compute_contrast_transient(self, simulation):
        """Test transient contrast computation."""
        contrast = simulation.compute_contrast(
            gamma=0.1, kmw_minus=1.0, kmw_plus=1.0,
            method='transient', t_integration=0.1
        )
        assert isinstance(contrast, float)

    def test_contrast_positive_for_nv(self, simulation):
        """Test NV bulk gives positive contrast."""
        contrast = simulation.compute_contrast(
            gamma=0.1, kmw_minus=1.0, kmw_plus=1.0,
            method='steady_state'
        )
        assert contrast > 0  # NV should have positive contrast

    def test_sweep_gamma(self, simulation):
        """Test gamma sweep returns correct shapes."""
        gammas = np.array([0.1, 0.2, 0.4])
        result_gammas, contrasts = simulation.sweep_gamma(
            gammas, kmw_minus=1.0, kmw_plus=0.0
        )
        assert len(result_gammas) == 3
        assert len(contrasts) == 3

    def test_invalid_method_raises(self, simulation):
        """Test invalid method raises error."""
        with pytest.raises(ValueError):
            simulation.compute_contrast(
                gamma=0.1, kmw_minus=1.0, kmw_plus=1.0,
                method='invalid_method'
            )


class TestInitializationSimulation:
    """Tests for InitializationSimulation class."""

    @pytest.fixture
    def simulation(self):
        """Create simulation for testing."""
        model = get_preset('nv_bulk')
        return InitializationSimulation(model)

    def test_run_basic(self, simulation):
        """Test basic initialization run."""
        result = simulation.run(gamma=0.1, t_max=1e-3, n_points=100)
        assert 't' in result
        assert 'populations' in result
        assert 'labels' in result

    def test_population_shapes(self, simulation):
        """Test population array has correct shape."""
        result = simulation.run(gamma=0.1, t_max=1e-3, n_points=100)
        assert result['populations'].shape == (100, 7)

    def test_sweep_gamma(self, simulation):
        """Test gamma sweep returns list of results."""
        gammas = [0.1, 0.3, 1.0]
        results = simulation.run_sweep_gamma(gammas, t_max=1e-3, n_points=50)
        assert len(results) == 3
        for result in results:
            assert 'populations' in result
