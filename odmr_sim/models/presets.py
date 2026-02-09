"""
Preset configurations for different defect systems.

Each preset provides rate parameters based on experimental or theoretical values.
"""

from typing import Dict, Any
from odmr_sim.models.seven_level import SevenLevelModel


# Preset configurations dictionary
# Each entry contains rate parameters for SevenLevelModel
PRESETS: Dict[str, Dict[str, Any]] = {
    "nv_bulk": {
        "description": "NV center in bulk diamond",
        "k41": 62.5,
        "k52": 62.5,
        "k63": 62.5,
        "k47": 10.5,
        "k57": 76.9,
        "k67": 76.9,
        "k71": 3.0,
        "k72": 2.63,
        "k73": 2.63,
        "use_degenerate_labels": True,
        "references": [
            "L. Robledo et al., New J. Phys. 13, 025013 (2011)",
            "S. Ahmadi et al., Phys. Rev. Applied 8, 034001 (2017)",
            "Y. Masuyama et al., Sci. Rep. 14, 18135 (2024)",
        ],
    },
    "g9_g8_30dp": {
        "description": "NV center G9-G8 configuration at 30 degrees",
        "k41": 62.5,
        "k52": 62.5,
        "k63": 62.5,
        "k47": 0.007,
        "k57": 0.5,
        "k67": 0.15,
        "k71": 279.1,
        "k72": 0.3,
        "k73": 21.3,
        "use_degenerate_labels": False,
        "references": [
            "L. Robledo et al., New J. Phys. 13, 025013 (2011)",
            "S. Ahmadi et al., Phys. Rev. Applied 8, 034001 (2017)",
            "Y. Masuyama et al., Sci. Rep. 14, 18135 (2024)",
        ],
    },
    "g4_g9_90dp": {
        "description": "NV center G4-G9 configuration at 90 degrees",
        "k41": 62.5,
        "k52": 62.5,
        "k63": 62.5,
        "k47": 31.2,
        "k57": 1.7,
        "k67": 12.3,
        "k71": 126.0,
        "k72": 3.3,
        "k73": 0.001,
        "use_degenerate_labels": False,
        "references": [
            "L. Robledo et al., New J. Phys. 13, 025013 (2011)",
            "S. Ahmadi et al., Phys. Rev. Applied 8, 034001 (2017)",
            "Y. Masuyama et al., Sci. Rep. 14, 18135 (2024)",
        ],
    },
    "gp7_gp3_30dp": {
        "description": "NV center G+7-G+3 configuration at 30 degrees",
        "k41": 62.5,
        "k52": 62.5,
        "k63": 62.5,
        "k47": 7.1,
        "k57": 9.4,
        "k67": 69.1,
        "k71": 717.0,
        "k72": 188.0,
        "k73": 23.3,
        "use_degenerate_labels": False,
        "references": [
            "L. Robledo et al., New J. Phys. 13, 025013 (2011)",
            "S. Ahmadi et al., Phys. Rev. Applied 8, 034001 (2017)",
            "Y. Masuyama et al., Sci. Rep. 14, 18135 (2024)",
        ],
    },
    "g11_gm9_30dp": {
        "description": "NV center G11-G-9 configuration at 30 degrees",
        "k41": 62.5,
        "k52": 62.5,
        "k63": 62.5,
        "k47": 4.4,
        "k57": 0.005,
        "k67": 44.1,
        "k71": 2336.0,
        "k72": 3.1,
        "k73": 0.001,
        "use_degenerate_labels": False,
        "references": [
            "L. Robledo et al., New J. Phys. 13, 025013 (2011)",
            "S. Ahmadi et al., Phys. Rev. Applied 8, 034001 (2017)",
            "Y. Masuyama et al., Sci. Rep. 14, 18135 (2024)",
        ],
    },
}


def list_presets() -> Dict[str, str]:
    """List all available presets with their descriptions.

    Returns
    -------
    presets : dict
        Dictionary mapping preset names to descriptions.

    Examples
    --------
    >>> from odmr_sim.models import list_presets
    >>> for name, desc in list_presets().items():
    ...     print(f"{name}: {desc}")
    """
    return {name: config["description"] for name, config in PRESETS.items()}


def get_preset(name: str) -> SevenLevelModel:
    """Get a SevenLevelModel with preset rate parameters.

    Parameters
    ----------
    name : str
        Name of the preset (e.g., 'nv_bulk', 'g11_gm9_30dp').
        Use list_presets() to see available options.

    Returns
    -------
    model : SevenLevelModel
        Configured model with the preset's rate parameters.

    Raises
    ------
    ValueError
        If the preset name is not recognized.

    Examples
    --------
    >>> from odmr_sim.models import get_preset
    >>> model = get_preset('nv_bulk')
    >>> W = model.build_rate_matrix(gamma=0.1)
    """
    name_lower = name.lower().replace("-", "_").replace("@", "_").replace(" ", "_")

    if name_lower not in PRESETS:
        available = ", ".join(PRESETS.keys())
        raise ValueError(
            f"Unknown preset '{name}'. Available presets: {available}"
        )

    config = PRESETS[name_lower]

    # Extract rate parameters (exclude metadata fields)
    rate_params = {
        k: v for k, v in config.items()
        if k not in ("description", "references")
    }

    return SevenLevelModel(**rate_params)


def get_preset_info(name: str) -> Dict[str, Any]:
    """Get detailed information about a preset.

    Parameters
    ----------
    name : str
        Name of the preset.

    Returns
    -------
    info : dict
        Dictionary with description, rate parameters, and references.
    """
    name_lower = name.lower().replace("-", "_").replace("@", "_").replace(" ", "_")

    if name_lower not in PRESETS:
        available = ", ".join(PRESETS.keys())
        raise ValueError(
            f"Unknown preset '{name}'. Available presets: {available}"
        )

    return PRESETS[name_lower].copy()
