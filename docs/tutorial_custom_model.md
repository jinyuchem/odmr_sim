# Tutorial: Custom N-Level Models

Create rate equation models with any number of states.

## Basic Custom Model

```python
from odmr_sim.models import RateModel
import numpy as np

# Create 3-level system
model = RateModel(
    n_states=3,
    state_labels=["Ground", "Excited", "Metastable"]
)

# Set transition rates (MHz)
model.set_rates({
    (1, 0): 50.0,  # Excited -> Ground
    (1, 2): 10.0,  # Excited -> Metastable
    (2, 0): 5.0,   # Metastable -> Ground
})

# Add dynamic rate (set at build time)
model.add_dynamic_rate('gamma', from_state=0, to_state=1)

# Build matrix
W = model.build_rate_matrix(gamma=1.0)
```

## Solving

```python
from odmr_sim.solvers import RateSolver

solver = RateSolver()
P0 = np.array([1.0, 0.0, 0.0])
t_eval = np.logspace(-9, -3, 500)

populations = solver.solve_expm(W, P0, t_eval)
```

## Complex Models

```python
# 9-level with multiple singlet states
model = RateModel(
    n_states=9,
    state_labels=[
        "GS0", "GS-", "GS+",
        "ES0", "ES-", "ES+",
        "S1", "S2", "S3"
    ]
)

model.set_rates({
    # Radiative
    (3, 0): 62.5,
    (4, 1): 62.5,
    (5, 2): 62.5,
    # ISC to singlets
    (3, 6): 4.0,
    (4, 7): 0.5,
    (5, 8): 40.0,
    # Singlet decay
    (6, 0): 1000.0,
    # ... add more
})
```

## Initial States

```python
# Pure state
P0 = model.get_initial_state(0)

# Mixed
P0 = model.get_mixed_initial_state({0: 0.5, 1: 0.25, 2: 0.25})
```
