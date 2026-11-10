"""Simplified relative satellite emission on the (frequency, time) grid.

For satellite i and in-band nominal carrier c::

    A_i(t)   = P_ref * (r_ref / r_i(t))**2                 range scaling
    G_i(t)   = exp(-4 ln2 (sep_i(t) / full_width)**2)      beam angular gain
    F_c(f)   = exp(-0.5 ((f - f_c) / sigma)**2), 0 beyond cutoff*sigma
    P_ic     = mask_i(t) * A_i(t) * G_i(t) * F_c(f)

``mask_i`` is the geometric passage mask (inside the beam gate and above the
horizon). Contributions are summed in linear relative power. Which satellites
and carriers contribute follows only from the geometry, the signal catalog
(constellation, applicable blocks) and the receiver band.

This is deterministic (no random numbers). It is not a calibrated received
power: no modulation, Doppler, polarization, multipath or receiver model.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .config import Scenario, ScenarioError
from .geometry import GeometryResult, SatelliteRecord

FOUR_LN2 = 4.0 * np.log(2.0)


@dataclass(frozen=True)
class Carrier:
    """A nominal carrier from the signal catalog."""

    signal_id: str
    constellation: str
    carrier_hz: float
    applicable_blocks: Tuple[str, ...]

    def applies_to(self, record: SatelliteRecord) -> bool:
        return (record.constellation == self.constellation
                and record.block in self.applicable_blocks)


@dataclass(frozen=True)
class GaussianCarrierProfile:
    """Relative Gaussian carrier power profile for one satellite block."""

    peak_power: float  # P_ref, relative power at r_ref and zero separation
    reference_range_m: float  # r_ref
    sigma_hz: float
    cutoff_sigma: float

    def __post_init__(self):
        for name in ("peak_power", "reference_range_m", "sigma_hz", "cutoff_sigma"):
            if not getattr(self, name) > 0:
                raise ScenarioError(f"emission profile {name} must be positive")

    def range_scaling(self, range_m: np.ndarray) -> np.ndarray:
        return (self.reference_range_m / np.asarray(range_m, dtype=np.float64)) ** 2

    def support(self, frequencies_hz: np.ndarray, carrier_hz: float) -> np.ndarray:
        """Boolean spectral support: |f - f_c| <= cutoff * sigma."""
        return np.abs(np.asarray(frequencies_hz) - carrier_hz) <= self.cutoff_sigma * self.sigma_hz

    def spectral_shape(self, frequencies_hz: np.ndarray, carrier_hz: float) -> np.ndarray:
        """Weights at the actual channel centers (peak 1 at f_c), 0 outside support."""
        f = np.asarray(frequencies_hz, dtype=np.float64)
        shape = np.exp(-0.5 * ((f - carrier_hz) / self.sigma_hz) ** 2)
        return np.where(self.support(f, carrier_hz), shape, 0.0)


def angular_gain(separation_deg: np.ndarray, full_width_deg: float) -> np.ndarray:
    """Gaussian beam response with half power at separation = full_width / 2."""
    return np.exp(-FOUR_LN2 * (np.asarray(separation_deg, dtype=np.float64) / full_width_deg) ** 2)


@dataclass
class Contribution:
    """One satellite/carrier term: outer(spectral, envelope)."""

    norad_id: int
    signal_id: str
    carrier_hz: float
    spectral: np.ndarray  # (n_freq,) float64 weights, 0 outside support
    support: np.ndarray  # (n_freq,) bool
    envelope: np.ndarray  # (n_time,) float64 = mask * A * G
    range_scaling: np.ndarray  # (n_time,) (r_ref/r)^2, unmasked
    gain: np.ndarray  # (n_time,) angular gain, unmasked
    mask: np.ndarray  # (n_time,) bool passage mask

    @property
    def active(self) -> bool:
        return bool(self.mask.any() and self.support.any())

    def power(self) -> np.ndarray:
        return np.outer(self.spectral, self.envelope)


# -- configuration --------------------------------------------------------

def load_signal_catalog(scenario: Scenario) -> List[Carrier]:
    path = scenario.reference_path("signal_catalog_file")
    with open(path, encoding="utf-8", newline="") as handle:
        return [Carrier(signal_id=row["signal_id"], constellation=row["constellation"],
                        carrier_hz=float(row["carrier_hz"]),
                        applicable_blocks=tuple(b for b in row["applicable_blocks"].split("|") if b))
                for row in csv.DictReader(handle)]


def carrier_in_band(carrier_hz: float, scenario: Scenario) -> bool:
    """Receiver rule: band_min_hz <= carrier_hz <= band_max_hz (inclusive)."""
    return scenario.band_min_hz <= carrier_hz <= scenario.band_max_hz


def profiles_from_scenario(scenario: Scenario) -> Dict[str, GaussianCarrierProfile]:
    cfg = scenario.section("emission")
    if cfg.get("model") != "relative_gaussian_carrier_power":
        raise ScenarioError(f"Unsupported emission model {cfg.get('model')!r}")
    if cfg.get("profile_key") != "satellite_catalog.block":
        raise ScenarioError("Only emission.profile_key = 'satellite_catalog.block' is supported")
    return {block: GaussianCarrierProfile(
                peak_power=float(p["peak_power_relative_at_reference_range"]),
                reference_range_m=float(p["reference_range_m"]),
                sigma_hz=float(p["spectral_sigma_hz"]),
                cutoff_sigma=float(p["spectral_cutoff_sigma"]))
            for block, p in cfg["profiles"].items()}


def satellite_signal_ids(record: SatelliteRecord, carriers: Sequence[Carrier],
                         scenario: Scenario) -> List[str]:
    """In-band signal IDs applicable to a satellite (catalog order)."""
    return [c.signal_id for c in carriers
            if c.applies_to(record) and carrier_in_band(c.carrier_hz, scenario)]


# -- synthesis ------------------------------------------------------------

def plan_contributions(scenario: Scenario, geometry: GeometryResult,
                       carriers: Optional[Sequence[Carrier]] = None,
                       profiles: Optional[Dict[str, GaussianCarrierProfile]] = None
                       ) -> List[Contribution]:
    """All satellite x applicable in-band carrier terms (including inactive ones)."""
    carriers = list(carriers) if carriers is not None else load_signal_catalog(scenario)
    profiles = profiles if profiles is not None else profiles_from_scenario(scenario)
    full_width = float(scenario.section("beam")["full_width_deg"])
    freqs = scenario.frequencies_hz
    if geometry.selected.shape != (len(geometry.satellites), scenario.n_times):
        raise ValueError("geometry does not match the scenario time grid")

    out: List[Contribution] = []
    for row, record in enumerate(geometry.satellites):
        applicable = [c for c in carriers
                      if c.applies_to(record) and carrier_in_band(c.carrier_hz, scenario)]
        if not applicable:
            continue
        if record.block not in profiles:
            raise ScenarioError(f"No emission profile for block {record.block!r} "
                                f"(NORAD {record.norad_id})")
        profile = profiles[record.block]
        mask = np.asarray(geometry.selected[row], dtype=bool)
        scaling = profile.range_scaling(geometry.range_m[row])
        gain = angular_gain(geometry.separation_deg[row], full_width)
        envelope = np.where(mask, profile.peak_power * scaling * gain, 0.0)
        for carrier in applicable:
            out.append(Contribution(
                norad_id=record.norad_id, signal_id=carrier.signal_id,
                carrier_hz=carrier.carrier_hz,
                spectral=profile.spectral_shape(freqs, carrier.carrier_hz),
                support=profile.support(freqs, carrier.carrier_hz),
                envelope=envelope, range_scaling=scaling, gain=gain, mask=mask))
    return out


def compute_signal_power(scenario: Scenario, geometry: GeometryResult,
                         carriers: Optional[Sequence[Carrier]] = None,
                         profiles: Optional[Dict[str, GaussianCarrierProfile]] = None
                         ) -> Tuple[np.ndarray, List[Contribution]]:
    """Summed signal power, float32 (frequency, time), and its contributions.

    Terms are accumulated in float64 and cast once to float32.
    """
    contributions = plan_contributions(scenario, geometry, carriers, profiles)
    total = np.zeros(scenario.matrix_shape, dtype=np.float64)
    for term in contributions:
        if not term.active:
            continue
        rows = np.flatnonzero(term.support)
        cols = np.flatnonzero(term.mask)
        total[np.ix_(rows, cols)] += np.outer(term.spectral[rows], term.envelope[cols])
    return total.astype(np.float32), contributions
