"""
odmr_sim - ODMR simulation using configurable N-level rate equation models.

A Python package for simulating Optically Detected Magnetic Resonance (ODMR)
using rate equation models with support for arbitrary numbers of states.
"""

__version__ = "0.1.0"
__author__ = "Yu Jin"

from odmr_sim.models import RateModel, SevenLevelModel
from odmr_sim.solvers import RateSolver

__all__ = ["RateModel", "SevenLevelModel", "RateSolver"]
