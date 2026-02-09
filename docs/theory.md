# Theory: Rate Equation Model for ODMR

## Overview

This package simulates ODMR phenomena using the **rate equation approach**, which describes the dynamics of spin state populations under optical excitation and microwave driving.

## The Master Equation

The time evolution of state populations $P(t)$ is governed by:

$$\frac{dP}{dt} = W \cdot P$$

where $W$ is the **rate matrix** containing all transition rates, and $P$ is the population vector.

The solution is:

$$P(t) = \exp(W \cdot t) \cdot P_0$$

where $P_0$ is the initial population and $\exp(W \cdot t)$ is the matrix exponential.

## 7-Level Model

The standard NV-like defect is modeled with 7 states:

### States

| Index | State | Description |
|-------|-------|-------------|
| 0 | GS $|0\rangle$ | Ground state, $m_s = 0$ |
| 1 | GS $|-\rangle$ | Ground state, $m_s = -1$ |
| 2 | GS $|+\rangle$ | Ground state, $m_s = +1$ |
| 3 | ES $|0\rangle$ | Excited state, $m_s = 0$ |
| 4 | ES $|-\rangle$ | Excited state, $m_s = -1$ |
| 5 | ES $|+\rangle$ | Excited state, $m_s = +1$ |
| 6 | SS | Metastable singlet |

### Rate Parameters

| Parameter | Description | Typical Value (MHz) |
|-----------|-------------|---------------------|
| $k_{41}$, $k_{52}$, $k_{63}$ | Radiative decay (ES $\to$ GS) | 62.5 |
| $k_{47}$ | Upper ISC (ES $|0\rangle \to$ SS) | 4.4 |
| $k_{57}$ | Upper ISC (ES $|-\rangle \to$ SS) | 0.005 |
| $k_{67}$ | Upper ISC (ES $|+\rangle \to$ SS) | 44.1 |
| $k_{71}$ | Lower ISC (SS $\to$ GS $|0\rangle$) | 2336 |
| $k_{72}$ | Lower ISC (SS $\to$ GS $|-\rangle$) | 3.1 |
| $k_{73}$ | Lower ISC (SS $\to$ GS $|+\rangle$) | 0.001 |
| $\Gamma$ | Optical excitation (GS $\to$ ES) | 0.1 - 10 |
| $k_\text{MW}$ | Microwave driving (GS $|0\rangle \leftrightarrow$ GS $|\pm\rangle$) | 0 - 10 |

## Key Phenomena

### Spin Initialization

Under optical pumping (with $\Gamma > 0$), the system preferentially populates GS $|0\rangle$ due to:
1. Spin-selective ISC from ES to singlet ($k_{67} \gg k_{57}$)
2. Preferential decay from singlet to GS $|0\rangle$ ($k_{71} \gg k_{72}, k_{73}$)

### Optical Readout

Fluorescence intensity is proportional to excited state population. Different ground states give different fluorescence due to:
- GS $|0\rangle \to$ bright (slower ISC path)
- GS $|\pm\rangle \to$ dimmer (faster ISC path through singlet)

### ODMR Contrast

$$\text{Contrast} = \frac{I_{\text{noMW}} - I_{\text{MW}}}{I_{\text{noMW}}}$$

When MW drives the $|0\rangle \leftrightarrow |\pm\rangle$ transition, population redistributes, changing the average fluorescence intensity.

## References

1. L. Robledo et al., New J. Phys. **13**, 025013 (2011)
2. S. Ahmadi et al., Phys. Rev. Applied **8**, 034001 (2017)
3. Y. Masuyama, C. Shinei, S. Ishii et al., Sci. Rep. **14**, 18135 (2024)
