# ODMR Simulation Examples

This directory contains example scripts demonstrating the `odmr_sim` package.

## Directory Structure

```
examples/
├── README.md                    # This file
├── nv_bulk/                     # NV center in bulk diamond (7-level)
├── g4_g9_90dp/                  # NV G4-G9 configuration at 90° (7-level)
├── g9_g8_30dp/                  # NV G9-G8 configuration at 30° (7-level)
├── g11_gm9_30dp/                # NV G11-Gm9 configuration at 30° (7-level)
├── gp7_gp3_30dp/                # NV GP7-GP3 configuration at 30° (7-level)
└── custom_8_level/              # Custom 8-level model with two metastable states
```

## Configuration Examples

### nv_bulk/ - NV Center in Bulk Diamond
- **Model**: 7-level with degenerate ±1 spin states
- **Contrast**: ~17% at γ = 0.1 MHz (both MW)
- Files:
  - `example_initialization.py` - Spin polarization dynamics
  - `example_readout.py` - Optical readout comparison
  - `example_odmr_contrast.py` - ODMR contrast vs excitation rate

### g4_g9_90dp/ - NV G4-G9 Configuration at 90°
- **Model**: 7-level with split (non-degenerate) ±1 states
- **Peak frequencies**: 1.640 GHz (−), 3.257 GHz (+)
- **Contrast**: ~10% (−), ~7% (+) at γ = 0.1 MHz
- Files: Same structure as nv_bulk

### g9_g8_30dp/ - NV G9-G8 Configuration at 30°
- **Model**: 7-level with split ±1 states
- **Peak frequencies**: 3.298 GHz (−), 3.632 GHz (+)
- **Contrast**: ~0.4% (−), ~0.1% (+) at γ = 0.1 MHz
- Files: Same structure as nv_bulk

### g11_gm9_30dp/ - NV G11-Gm9 Configuration at 30°
- **Model**: 7-level with split ±1 states
- **Peak frequencies**: 2.959 GHz (−), 3.921 GHz (+)
- **Contrast**: ~0.15% (−), ~1.7% (+) at γ = 0.1 MHz (steady-state)
- Files:
  - Standard examples (initialization, readout, odmr_contrast)
  - `example_contrast_comparison.py` - Compares steady-state vs transient methods

### gp7_gp3_30dp/ - NV GP7-GP3 Configuration at 30°
- **Model**: 7-level with split ±1 states
- **Peak frequencies**: 2.970 GHz (−), 3.710 GHz (+)
- **Contrast**: ~1% (−), ~15% (+) at γ = 0.1 MHz
- Files: Same structure as nv_bulk

### custom_8_level/ - Custom 8-Level Model
- **Model**: 8-level with two metastable states (MS1, MS2)
- **Level structure**:
  - GS|0⟩, GS|−⟩, GS|+⟩ (ground states)
  - ES|0⟩, ES|−⟩, ES|+⟩ (excited states)
  - MS1 (singlet 1) - connected to ES via ISC
  - MS2 (singlet 2) - connected to GS via ISC
- **Decay pathway**: ES → MS1 → MS2 → GS
- **Contrast**: Inverted (negative, ~5-8%)
- Files:
  - `example_custom_model.py` - Original custom model demo
  - `example_8_level_model.py` - Complete 8-level simulation

## Running Examples

```bash
# Navigate to an example directory
cd examples/nv_bulk

# Run an example
python3 example_odmr_contrast.py
```

## Output

Each example generates:
- Console output with numerical results
- PDF figures saved to the example directory

## Contrast Calculation Methods

The package supports multiple methods for computing ODMR contrast:

| Method | Description |
|--------|-------------|
| `steady_state` | True equilibrium at t → ∞ (null-space of rate matrix) |
| `transient` | Final ES population at finite time |
| `time_integrated` | Integrated fluorescence over [0, t_integration] |

See `g11_gm9_30dp/example_contrast_comparison.py` for a comparison of these methods.
