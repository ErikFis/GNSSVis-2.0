"""Synthetic GNSS waterfall simulations built on GNSSVis.

Depends on the ``gnssvis`` library (never the reverse). The scenario JSON is
the single source of truth for grids, site, pointing, beam and references.
"""

from .config import Scenario, ScenarioError, load_scenario, verify_reference_hashes
from .geometry import (
    GeometryResult,
    Passage,
    SatelliteRecord,
    angular_separation_deg,
    beam_masks,
    compute_geometry,
    find_passages,
    geometry_passages,
    load_reference_satellites,
    offline_sky_analyzer,
)

__all__ = [
    "Scenario", "ScenarioError", "load_scenario", "verify_reference_hashes",
    "GeometryResult", "Passage", "SatelliteRecord", "angular_separation_deg", "beam_masks",
    "compute_geometry", "find_passages", "geometry_passages",
    "load_reference_satellites", "offline_sky_analyzer",
    "simulate_realistic", "save_realistic", "RealisticResult",
]


def __getattr__(name):
    # Keep the optional wideband CLI import-safe and avoid loading it twice
    # when invoked with python -m gnssvis_simulation.realistic.
    if name in {"simulate_realistic", "save_realistic", "RealisticResult"}:
        from . import realistic
        return getattr(realistic, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
