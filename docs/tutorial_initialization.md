# Tutorial: Initialization Simulation

This tutorial shows how to simulate spin polarization (initialization) dynamics.

## Basic Usage

```python
from odmr_sim.models import SevenLevelModel
from odmr_sim.simulations import InitializationSimulation

# Create model
model = SevenLevelModel()

# Create and run simulation
sim = InitializationSimulation(model)
result = sim.run(
    gamma=0.1,        # Excitation rate (MHz)
    kmw_minus=0.0,    # No microwave
    kmw_plus=0.0,
    t_max=1e-1,       # 100 ms
    n_points=1000,
)

# Access results
t = result['t']              # Time in seconds
t_ns = result['t_ns']        # Time in nanoseconds
populations = result['populations']  # Shape: (n_times, 7)
labels = result['labels']    # State labels
```

## With Microwave

```python
# With MW driving the |-> transition
result_mw = sim.run(
    gamma=0.1,
    kmw_minus=1.0,  # MW on
    kmw_plus=0.0,
)
```

## Sweeping Excitation Rate

```python
gammas = [0.1, 0.3, 1.0]
results = sim.run_sweep_gamma(gammas)

for i, result in enumerate(results):
    print(f"Gamma = {gammas[i]} MHz")
    print(f"Final GS|0>: {result['populations'][-1, 0]:.3f}")
```

## Plotting

```python
from odmr_sim.plotting import plot_populations

plot_populations(
    result['t'], result['populations'],
    labels=result['labels'],
    xlim=(1, 1e6),
    use_log_x=True,
)
```
