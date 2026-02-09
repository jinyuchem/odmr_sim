# Tutorial: ODMR Simulation

This tutorial shows how to compute ODMR contrast and spectra.

## Basic Contrast

```python
from odmr_sim.models import SevenLevelModel
from odmr_sim.simulations import ODMRSimulation

model = SevenLevelModel()
sim = ODMRSimulation(model)

# Single contrast calculation
contrast = sim.compute_contrast(
    gamma=0.1,
    kmw_minus=1.0,  # Drive |-> transition
    kmw_plus=0.0,
)
print(f"Contrast: {contrast*100:.2f}%")
```

## Sweep Excitation Rate

```python
import numpy as np

gammas = np.logspace(-1, 1, 20)
gammas, contrasts = sim.sweep_gamma(gammas, kmw_minus=1.0)

# Plot
import matplotlib.pyplot as plt
plt.semilogx(gammas, contrasts * 100, 'o-')
plt.xlabel('Gamma (MHz)')
plt.ylabel('Contrast (%)')
```

## ODMR Spectrum

```python
spectrum = sim.run_spectrum(
    gamma=0.1,
    freq_center=3.0,      # GHz
    freq_width=1.0,       # GHz
    peak_freq_minus=2.7,  # Resonance for |->
    peak_freq_plus=3.3,   # Resonance for |+>
    linewidth=0.05,       # GHz
    kmw_amplitude=10.0,   # MHz at resonance
)

# Plot
from odmr_sim.plotting import plot_odmr_spectrum
plot_odmr_spectrum(spectrum['frequencies'], spectrum['contrast'])
```

## Contrast Methods

```python
# Steady-state (default, fast)
c_ss = sim.compute_contrast(gamma=0.1, kmw_minus=1.0, method='steady_state')

# Time-integrated (more accurate for pulsed)
c_ti = sim.compute_contrast(
    gamma=0.1, kmw_minus=1.0,
    method='time_integrated',
    t_integration=1e-3
)
```
