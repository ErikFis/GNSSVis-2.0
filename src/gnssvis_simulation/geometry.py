"""Offline satellite geometry, beam selection masks and passages.

Reuses ``gnssvis.SkyAnalyzer`` in offline mode: local TLEs, SGP4, the WGS84
observer and Skyfield's bundled timescale. No planetary ephemeris (DE421),
Earth-orientation update or network access is needed.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional

import numpy as np

from gnssvis.sky import SkyAnalyzer

from .config import Scenario, ScenarioError, format_utc, parse_utc


@dataclass(frozen=True)
class SatelliteRecord:
    """A selected satellite: catalog metadata plus its SGP4 model."""

    norad_id: int
    name: str
    constellation: str
    block: str
    prn: int
    tle_epoch_utc: datetime
    satellite: object  # skyfield.api.EarthSatellite


@dataclass(frozen=True)
class Passage:
    """A maximal contiguous run of selected samples for one satellite."""

    norad_id: int
    start_index: int  # first included sample
    stop_index: int  # first following excluded sample, or n_times if clipped
    entry_offset_seconds: int
    exit_offset_seconds: int  # exclusive
    duration_seconds: int
    sample_count: int
    entry_clipped: bool
    exit_clipped: bool


@dataclass
class GeometryResult:
    """Per-satellite, per-sample geometry on the scenario time grid.

    All 2-D arrays have shape (n_satellites, n_times); rows follow the
    configured NORAD ID order and columns the scenario time offsets.
    """

    satellites: List[SatelliteRecord]
    time_offsets_seconds: np.ndarray  # int64, (n_times,)
    start_utc: datetime
    step_seconds: int
    alt_deg: np.ndarray
    az_deg: np.ndarray
    range_m: np.ndarray
    separation_deg: np.ndarray
    above_horizon: np.ndarray
    in_beam: np.ndarray
    selected: np.ndarray  # passage mask: in_beam & above_horizon
    range_rate_m_s: Optional[np.ndarray] = None  # positive for recession

    @property
    def norad_ids(self) -> List[int]:
        return [s.norad_id for s in self.satellites]

    @property
    def timestamps_utc(self) -> List[datetime]:
        return [self.start_utc + timedelta(seconds=int(s)) for s in self.time_offsets_seconds]


# --------------------------------------------------------------------------
def offline_sky_analyzer(scenario: Scenario) -> SkyAnalyzer:
    """SkyAnalyzer at the scenario site in offline mode with the local TLEs."""
    site = scenario.section("observatory")
    return SkyAnalyzer(
        float(site["latitude_deg"]),
        float(site["longitude_deg"]),
        float(site["elevation_m"]),
        offline=True,
        tle_file=scenario.reference_path("tle_file"),
    )


def load_satellite_catalog(scenario: Scenario) -> Dict[int, Dict[str, str]]:
    """Read satellites.csv as a mapping from NORAD ID to row."""
    path = scenario.reference_path("satellite_catalog_file")
    with open(path, encoding="utf-8", newline="") as handle:
        return {int(row["norad_id"]): row for row in csv.DictReader(handle)}


def load_reference_satellites(scenario: Scenario,
                              sky: Optional[SkyAnalyzer] = None) -> List[SatelliteRecord]:
    """Load the configured satellites from the local TLE file.

    Checks that every configured NORAD ID exists in the TLE file and catalog,
    that the catalog epoch matches the epoch encoded in the TLE, and that the
    whole observation window is within the configured epoch offset limit.
    """
    sky = sky or offline_sky_analyzer(scenario)
    ids = scenario.norad_ids
    satellites = sky.load_local_tles(scenario.reference_path("tle_file"), norad_ids=ids)
    catalog = load_satellite_catalog(scenario)
    limit_h = float(scenario.config["references"]["max_absolute_tle_epoch_offset_hours"])

    records: List[SatelliteRecord] = []
    for norad_id, sat in zip(ids, satellites):
        if norad_id not in catalog:
            raise ScenarioError(f"NORAD ID {norad_id} missing from the satellite catalog")
        row = catalog[norad_id]
        epoch = sat.epoch.utc_datetime()
        catalog_epoch = parse_utc(row["tle_epoch_utc"])
        if abs((epoch - catalog_epoch).total_seconds()) > 1e-3:
            raise ScenarioError(
                f"Catalog epoch {row['tle_epoch_utc']} differs from TLE epoch "
                f"{format_utc(epoch)} for NORAD {norad_id}")
        worst_h = max(abs((scenario.start_utc - epoch).total_seconds()),
                      abs((scenario.end_utc - epoch).total_seconds())) / 3600.0
        if worst_h > limit_h:
            raise ScenarioError(
                f"NORAD {norad_id}: observation window is up to {worst_h:.1f} h from the "
                f"TLE epoch, exceeding the {limit_h:g} h limit")
        records.append(SatelliteRecord(
            norad_id=norad_id,
            name=row["name"],
            constellation=row["constellation"],
            block=row["block"],
            prn=int(row["prn"]),
            tle_epoch_utc=epoch,
            satellite=sat,
        ))
    return records


def angular_separation_deg(alt1_deg, az1_deg, alt2_deg, az2_deg) -> np.ndarray:
    """Great-circle separation between alt/az directions (Vincenty formula)."""
    a1, z1, a2, z2 = (np.radians(np.asarray(v, dtype=np.float64))
                      for v in (alt1_deg, az1_deg, alt2_deg, az2_deg))
    dz = z2 - z1
    num = np.hypot(np.cos(a2) * np.sin(dz),
                   np.cos(a1) * np.sin(a2) - np.sin(a1) * np.cos(a2) * np.cos(dz))
    den = np.sin(a1) * np.sin(a2) + np.cos(a1) * np.cos(a2) * np.cos(dz)
    return np.degrees(np.arctan2(num, den))


def beam_masks(separation_deg, alt_deg, full_width_deg: float, minimum_altitude_deg: float):
    """Return (above_horizon, in_beam, selected) boolean masks.

    in_beam uses the inclusive rule ``separation <= full_width / 2``;
    above_horizon uses ``altitude >= minimum_altitude``.
    """
    above = np.asarray(alt_deg) >= minimum_altitude_deg
    in_beam = np.asarray(separation_deg) <= full_width_deg / 2.0
    return above, in_beam, above & in_beam


def skyfield_times(scenario: Scenario, sky: SkyAnalyzer):
    """Skyfield Time array for the scenario sample grid (UTC)."""
    start = scenario.start_utc
    seconds = start.second + start.microsecond / 1e6 + scenario.time_offsets_seconds
    return sky.ts.utc(start.year, start.month, start.day, start.hour, start.minute, seconds)


def compute_geometry(scenario: Scenario,
                     sky: Optional[SkyAnalyzer] = None,
                     satellites: Optional[List[SatelliteRecord]] = None) -> GeometryResult:
    """Propagate the configured satellites over the scenario time grid.

    The passage mask is ``separation <= full_width/2`` and
    ``altitude >= minimum_altitude``, evaluated only at the sample times
    (no subsample interpolation).
    """
    sky = sky or offline_sky_analyzer(scenario)
    satellites = satellites if satellites is not None else load_reference_satellites(scenario, sky)
    times = skyfield_times(scenario, sky)
    n_sat, n_t = len(satellites), scenario.n_times

    alt = np.empty((n_sat, n_t))
    az = np.empty((n_sat, n_t))
    rng = np.empty((n_sat, n_t))
    rates = np.empty((n_sat, n_t))
    for row, record in enumerate(satellites):
        values = sky.satellite_altaz(record.satellite, times)
        alt[row], az[row], rng[row] = values["alt_deg"], values["az_deg"], values["range_m"]
        rates[row] = values["range_rate_m_s"]

    pointing = scenario.section("pointing")
    beam = scenario.section("beam")
    separation = angular_separation_deg(alt, az, float(pointing["altitude_deg"]),
                                        float(pointing["azimuth_deg"]))
    above, in_beam, selected = beam_masks(separation, alt, float(beam["full_width_deg"]),
                                          float(beam["minimum_altitude_deg"]))
    if beam.get("model") == "chromatic_elliptical_gaussian":
        from .beam import beam_offsets_deg, beam_gain
        x, y = beam_offsets_deg(alt, az, float(pointing["altitude_deg"]),
                                float(pointing["azimuth_deg"]))
        in_beam = beam_gain(x, y, float(beam["reference_frequency_hz"]), beam) >= 0.5
        selected = above & in_beam
    return GeometryResult(
        satellites=satellites,
        time_offsets_seconds=scenario.time_offsets_seconds,
        start_utc=scenario.start_utc,
        step_seconds=scenario.step_seconds,
        alt_deg=alt,
        az_deg=az,
        range_m=rng,
        separation_deg=separation,
        above_horizon=above,
        in_beam=in_beam,
        selected=selected,
        range_rate_m_s=rates,
    )


def find_passages(mask: np.ndarray, step_seconds: int, norad_id: int = 0) -> List[Passage]:
    """Split one satellite's boolean sample mask into passages.

    Entry is the first included sample; exit is the first following excluded
    sample, or the observation end when clipped. Duration is the number of
    included samples times the step.
    """
    mask = np.asarray(mask, dtype=bool)
    n = mask.size
    padded = np.concatenate(([False], mask, [False])).astype(np.int8)
    edges = np.diff(padded)
    starts = np.flatnonzero(edges == 1)
    stops = np.flatnonzero(edges == -1)
    passages = []
    for start, stop in zip(starts, stops):
        count = int(stop - start)
        passages.append(Passage(
            norad_id=int(norad_id),
            start_index=int(start),
            stop_index=int(stop),
            entry_offset_seconds=int(start) * step_seconds,
            exit_offset_seconds=int(stop) * step_seconds,
            duration_seconds=count * step_seconds,
            sample_count=count,
            entry_clipped=bool(start == 0),
            exit_clipped=bool(stop == n),
        ))
    return passages


def geometry_passages(geometry: GeometryResult) -> List[Passage]:
    """All passages ordered by entry time, then the configured NORAD order."""
    passages = []
    for row, record in enumerate(geometry.satellites):
        passages.extend(find_passages(geometry.selected[row], geometry.step_seconds,
                                      record.norad_id))
    order = {n: i for i, n in enumerate(geometry.norad_ids)}
    return sorted(passages, key=lambda p: (p.start_index, order[p.norad_id]))
