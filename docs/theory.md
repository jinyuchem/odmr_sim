# Theory: Rate Equation Model for ODMR

## Overview

This package simulates ODMR phenomena using the **rate equation approach**, which describes the dynamics of spin state populations under optical excitation and microwave driving.

## The Master Equation

The time evolution of state populations P(t) is governed by:

```
dP/dt = W · P
```

where W is the **rate matrix** containing all transition rates, and P is the population vector.

The solution is:

```
P(t) = exp(W·t) · P₀
```

where P₀ is the initial population and exp(W·t) is the matrix exponential.

## 7-Level Model

The standard NV-like defect is modeled with 7 states:

```
         ┌────────────────┐
         │  Excited State │
         │ ES|0> ES|-> ES|+>│
         └───────┬────────┘
           k_rad │ ↓  ↘ k_ISC (upper)
                 │     ↘
         ┌───────┴───┐ ┌───┐
         │Ground State│ │ SS│ (Singlet)
         │GS|0> GS|-> GS|+>│ └─┬─┘
         └───────────┘   │
                         └──-> k_ISC (lower)
```

### States

| Index | State | Description |
|-------|-------|-------------|
| 0 | GS\|0> | Ground state, mₛ = 0 |
| 1 | GS\|-> | Ground state, mₛ = -1 |
| 2 | GS\|+> | Ground state, mₛ = +1 |
| 3 | ES\|0> | Excited state, mₛ = 0 |
| 4 | ES\|-> | Excited state, mₛ = -1 |
| 5 | ES\|+> | Excited state, mₛ = +1 |
| 6 | SS | Metastable singlet |

### Rate Parameters

| Parameter | Description | Typical Value (MHz) |
|-----------|-------------|---------------------|
| k₄₁, k₅₂, k₆₃ | Radiative decay (ES -> GS) | 62.5 |
| k₄₇ | Upper ISC (ES\|0> -> SS) | 4.4 |
| k₅₇ | Upper ISC (ES\|-> -> SS) | 0.005 |
| k₆₇ | Upper ISC (ES\|+> -> SS) | 44.1 |
| k₇₁ | Lower ISC (SS -> GS\|0>) | 2336 |
| k₇₂ | Lower ISC (SS -> GS\|->) | 3.1 |
| k₇₃ | Lower ISC (SS -> GS\|+>) | 0.001 |
| Gamma | Optical excitation (GS -> ES) | 0.1 - 10 |
| k_MW | Microwave driving (GS\|0> <-> GS\|+/->) | 0 - 10 |

## Key Phenomena

### Spin Initialization

Under optical pumping (with Gamma > 0), the system preferentially populates GS|0> due to:
1. Spin-selective ISC from ES to singlet (k₆₇ >> k₅₇)
2. Preferential decay from singlet to GS|0> (k₇₁ >> k₇₂, k₇₃)

### Optical Readout

Fluorescence intensity is proportional to excited state population. Different ground states give different fluorescence due to:
- GS|0> -> bright (slower ISC path)
- GS|+/-> -> dimmer (faster ISC path through singlet)

### ODMR Contrast

Contrast = (I_noMW - I_MW) / I_noMW

When MW drives the |0> <-> |+/-> transition, population redistributes, changing the average fluorescence intensity.

## References

1. L. Robledo et al., New J. Phys. **13**, 025013 (2011)
2. S. Ahmadi et al., Phys. Rev. Applied **8**, 034001 (2017)
3. Y. Masuyama, C. Shinei, S. Ishii et al., Sci. Rep. **14**, 18135 (2024)
