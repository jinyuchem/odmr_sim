# ODMR_SIM

A Python package for simulating **Optically Detected Magnetic Resonance (ODMR)** using configurable N-level rate equation models.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

- **Configurable N-level models**: Support for arbitrary number of states
- **Built-in 7-level model**: Pre-configured for NV-like defect systems
- **Multiple solver methods**: Matrix exponential (fast) and ODE integration
- **Preset configurations**: NV bulk, G9-G8@30°DP, G11-Gm9@30°DP, and more
- **Simulation types**: Initialization, readout, and ODMR contrast

## Installation

```bash
# Clone the repository
git clone https://github.com/jinyuchem/odmr_sim.git
cd odmr_sim

# Install in development mode
pip install -e .
```

## Quick Start

### Using the 7-level model

```python
import numpy as np
from odmr_sim.models import SevenLevelModel
from odmr_sim.solvers import RateSolver

# Create model with default NV-like rates
model = SevenLevelModel()

# Build rate matrix (Gamma = 0.1 MHz excitation)
W = model.build_rate_matrix(gamma=0.1, kmw_minus=0.0, kmw_plus=0.0)

# Initial state: equal population in ground states
P0 = np.array([1/3, 1/3, 1/3, 0, 0, 0, 0])

# Time points (in seconds)
t_eval = np.logspace(-9, -1, 1000)

# Solve
solver = RateSolver()
populations = solver.solve_expm(W, P0, t_eval)

print(f"Final GS|0> population: {populations[-1, 0]:.3f}")
```

### Creating a custom N-level model

```python
from odmr_sim.models import RateModel

# Create a 5-level system
model = RateModel(n_states=5, state_labels=["GS", "ES1", "ES2", "MS1", "MS2"])

# Define transition rates (in MHz)
model.set_rates({
    (0, 1): 10.0,   # GS -> ES1
    (1, 0): 62.5,   # ES1 -> GS (radiative)
    (1, 3): 5.0,    # ES1 -> MS1 (ISC)
    (3, 0): 100.0,  # MS1 -> GS
})

W = model.build_rate_matrix()
```

## Rate Parameters

The 7-level model includes:
- **Radiative transitions**: $k_{41}$, $k_{52}$, $k_{63}$ (ES $\to$ GS)
- **Upper ISC rates**: $k_{47}$, $k_{57}$, $k_{67}$ (ES $\to$ Singlet)
- **Lower ISC rates**: $k_{71}$, $k_{72}$, $k_{73}$ (Singlet $\to$ GS)
- **Excitation rate**: $\Gamma$ (GS $\to$ ES)
- **Microwave rates**: $k_\text{MW}^\pm$ (GS $|0\rangle \leftrightarrow$ GS $|\pm\rangle$)

## Documentation

See the [docs/](docs/) folder for:
- [Theory background](docs/theory.md)
- [Initialization tutorial](docs/tutorial_initialization.md)
- [Readout tutorial](docs/tutorial_readout.md)
- [ODMR contrast tutorial](docs/tutorial_odmr.md)
- [Custom model tutorial](docs/tutorial_custom_model.md)

## References

1. L. Robledo et al., New J. Phys. **13**, 025013 (2011)
2. S. Ahmadi et al., Phys. Rev. Applied **8**, 034001 (2017)
3. Y. Masuyama, C. Shinei, S. Ishii et al., Sci. Rep. **14**, 18135 (2024)
4. C. Zhang, V. W.-z. Yu, Y. Jin et al., npj Comput. Mater. **12**, 77 (2026)

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments
This package was developed with assistance from [Claude](https://www.anthropic.com/claude), an AI assistant by Anthropic.
