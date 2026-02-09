# Tutorial: Readout Simulation

This tutorial shows how to simulate optical readout dynamics.

## Basic Usage

```python
from odmr_sim.models import SevenLevelModel
from odmr_sim.simulations import ReadoutSimulation

model = SevenLevelModel()
sim = ReadoutSimulation(model)

# Run from GS|0>
result = sim.run(
    gamma=12.8,           # Higher Gamma for readout
    initial_state='gs0',  # Start in GS|0>
    t_max=1e-5,           # 10 μs
)

# Excited state total (proportional to fluorescence)
es_total = result['es_total']
```

## Initial States

```python
# Different ways to specify initial state
result_0 = sim.run(initial_state='gs0')      # GS|0>
result_m = sim.run(initial_state='gs_minus') # GS|->
result_p = sim.run(initial_state='gs_plus')  # GS|+>

# By index
result_0 = sim.run(initial_state=0)

# Custom
import numpy as np
P0 = np.array([0.8, 0.1, 0.1, 0, 0, 0, 0])
result = sim.run(initial_state=P0)
```

## Comparison

```python
# Compare all three ground states
results = sim.run_comparison(gamma=12.8, t_max=1e-5)

for state, res in results.items():
    peak = res['es_total'].max()
    print(f"{state}: peak ES = {peak:.3f}")
```

## Plotting

```python
from odmr_sim.plotting import plot_readout_comparison

plot_readout_comparison(results, xlim=(-50, 600))
```
