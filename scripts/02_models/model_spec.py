"""Single source of truth for transmission parameters and units."""
from __future__ import annotations

import json
import os
from pathlib import Path

SPEC_PATH = Path(os.environ.get("VPDS_TRANSMISSION_SPEC", str(Path(__file__).with_name("transmission_model_spec.json"))))


def load_spec(path: Path = SPEC_PATH) -> dict:
    with path.open(encoding="utf-8") as handle:
        spec = json.load(handle)
    if spec["time"]["unit"] != "day" or spec["time"]["dt_days"] != 1.0:
        raise ValueError("The active transmission kernel requires a one-day time step")
    required = {"R0", "latent_days", "infectious_days", "acquisition_ratio", "clinical_fraction", "reporting"}
    for disease, values in spec["diseases"].items():
        missing = required.difference(values)
        if missing:
            raise ValueError(f"{disease} missing specification fields: {sorted(missing)}")
        if not (0 < values["R0"] and 0 < values["latent_days"] and 0 < values["infectious_days"]):
            raise ValueError(f"Invalid natural history for {disease}")
        for key in ("acquisition_ratio", "clinical_fraction", "reporting"):
            if not 0 <= values[key] <= 1:
                raise ValueError(f"{disease}.{key} must be in [0, 1]")
    return spec


def protection_factors(gap, acquisition_ratio, source_reduction=0.0):
    """Use already-scaled clinical protection mass; keep source effects separate."""
    import numpy as np
    gap = np.clip(np.asarray(gap, dtype=float), 0.0, 1.0)
    protected = 1.0 - gap
    susceptible = 1.0 - acquisition_ratio * protected
    denominator = np.maximum(susceptible, 1e-12)
    protected_among_infected = protected * (1.0 - acquisition_ratio) / denominator
    source = 1.0 - source_reduction * protected_among_infected
    clinical = gap / denominator
    return susceptible, np.clip(source, 0.0, 1.0), np.clip(clinical, 0.0, 1.0)


SPEC = load_spec()
DISEASES = list(SPEC["diseases"])
AGE_LABELS = list(SPEC["age"]["bands"])
