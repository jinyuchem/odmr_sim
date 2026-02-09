"""Tests for ODMR simulation models."""

import pytest
import numpy as np
from odmr_sim.models import RateModel, SevenLevelModel, get_preset


class TestRateModel:
    """Tests for the base RateModel class."""

    def test_creation(self):
        """Test basic model creation."""
        model = RateModel(n_states=3, state_labels=["A", "B", "C"])
        assert model.n_states == 3
        assert model.state_labels == ["A", "B", "C"]

    def test_set_rates(self):
        """Test setting transition rates."""
        model = RateModel(n_states=3)
        model.set_rates({(0, 1): 10.0, (1, 0): 5.0})
        W = model.build_rate_matrix()
        # Check that rates are set (may be in MHz converted to Hz)
        assert W[1, 0] > 0  # Rate from 0 to 1 is positive
        assert W[0, 1] > 0  # Rate from 1 to 0 is positive

    def test_rate_matrix_conservation(self):
        """Test that columns sum to zero (probability conservation)."""
        model = RateModel(n_states=3)
        model.set_rates({(0, 1): 10.0, (1, 2): 5.0, (2, 0): 3.0})
        W = model.build_rate_matrix()
        column_sums = np.sum(W, axis=0)
        np.testing.assert_array_almost_equal(column_sums, 0.0)


class TestSevenLevelModel:
    """Tests for the SevenLevelModel class."""

    def test_creation(self):
        """Test 7-level model creation."""
        model = SevenLevelModel()
        assert model.n_states == 7

    def test_default_rates(self):
        """Test default NV-like rates are set."""
        model = SevenLevelModel()
        assert model.k41 == 62.5
        assert model.k52 == 62.5

    def test_rate_matrix_shape(self):
        """Test rate matrix has correct shape."""
        model = SevenLevelModel()
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)
        assert W.shape == (7, 7)

    def test_rate_matrix_conservation(self):
        """Test probability conservation."""
        model = SevenLevelModel()
        W = model.build_rate_matrix(gamma=0.1, kmw_minus=1.0, kmw_plus=1.0)
        column_sums = np.sum(W, axis=0)
        np.testing.assert_array_almost_equal(column_sums, 0.0)

    def test_ground_state_mixed(self):
        """Test mixed ground state initialization."""
        model = SevenLevelModel()
        P0 = model.get_ground_state_mixed()
        assert len(P0) == 7
        np.testing.assert_almost_equal(np.sum(P0), 1.0)


class TestPresets:
    """Tests for preset configurations."""

    @pytest.mark.parametrize("preset_name", [
        "nv_bulk",
        "g4_g9_90dp",
        "g9_g8_30dp",
        "g11_gm9_30dp",
        "gp7_gp3_30dp",
    ])
    def test_preset_loading(self, preset_name):
        """Test all presets can be loaded."""
        model = get_preset(preset_name)
        assert model is not None
        assert model.n_states == 7

    def test_invalid_preset(self):
        """Test invalid preset raises error."""
        with pytest.raises(ValueError):
            get_preset("nonexistent_preset")

    def test_nv_bulk_is_seven_level(self):
        """Test NV bulk is a 7-level model."""
        model = get_preset("nv_bulk")
        assert model.n_states == 7

    def test_nv_config_is_seven_level(self):
        """Test NV configuration presets are 7-level models."""
        model = get_preset("g4_g9_90dp")
        assert model.n_states == 7
