# Documentation

Welcome to the `odmr_sim` documentation!

## Contents

- [Theory](theory.md) - Background on the rate equation approach
- [Tutorial: Initialization](tutorial_initialization.md) - Simulating spin polarization
- [Tutorial: Readout](tutorial_readout.md) - Simulating optical readout
- [Tutorial: ODMR](tutorial_odmr.md) - Computing ODMR contrast
- [Tutorial: Custom Models](tutorial_custom_model.md) - Creating custom N-level models

## Quick Start

```python
from odmr_sim.models import SevenLevelModel
from odmr_sim.simulations import InitializationSimulation

model = SevenLevelModel()
sim = InitializationSimulation(model)
result = sim.run(gamma=0.1)
```

See the [examples/](../examples/) directory for complete working examples.
