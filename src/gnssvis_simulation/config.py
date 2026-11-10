"""Scenario loading and grid validation for the offline simulation.

The scenario JSON (for example ``examples/scenarios/offline.json``) is the
single source of truth for the observation parameters. Reference paths inside
it are resolved relative to the repository root, never the working directory.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

PathLike = Union[str, Path]


class ScenarioError(ValueError):
    """Invalid or unsupported scenario configuration."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ScenarioError(message)


def _positive_int(value: Any, name: str) -> int:
    _require(isinstance(value, (int, float)) and not isinstance(value, bool),
             f"{name} must be a number, got {value!r}")
    _require(float(value).is_integer(), f"{name} must be an integer, got {value!r}")
    _require(value > 0, f"{name} must be positive, got {value!r}")
    return int(value)


def parse_utc(text: str) -> datetime:
    """Parse an ISO 8601 UTC timestamp with a ``Z`` suffix."""
    _require(isinstance(text, str) and text.endswith("Z"),
             f"UTC timestamps must end with 'Z', got {text!r}")
    parsed = datetime.fromisoformat(text[:-1])
    _require(parsed.tzinfo is None, f"Unexpected UTC offset in {text!r}")
    return parsed.replace(tzinfo=timezone.utc)


def format_utc(moment: datetime) -> str:
    """Format an aware UTC datetime as ISO 8601 with a ``Z`` suffix."""
    moment = moment.astimezone(timezone.utc)
    fmt = "%Y-%m-%dT%H:%M:%S.%fZ" if moment.microsecond else "%Y-%m-%dT%H:%M:%SZ"
    return moment.strftime(fmt)


def find_repo_root(start: PathLike) -> Path:
    """Repository root for resolving reference paths.

    The first ancestor of ``start`` containing pyproject.toml; for a scenario
    stored outside the repository, the checkout that contains this package.
    """
    start = Path(start).resolve()
    package_dir = Path(__file__).resolve().parent
    for origin in (start, package_dir):
        for candidate in [origin, *origin.parents]:
            if (candidate / "pyproject.toml").is_file() and (candidate / "reference_data").is_dir():
                return candidate
    raise ScenarioError(f"No repository root (pyproject.toml + reference_data/) above {start}")


def sha256_file(path: PathLike) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class Scenario:
    """Validated scenario with derived time and frequency grids."""

    config: Dict[str, Any]
    config_path: Path
    repo_root: Path
    start_utc: datetime
    step_seconds: int
    duration_seconds: int
    channel_count: int
    band_min_hz: float
    band_max_hz: float

    # -- time axis (matrix columns) -------------------------------------
    @property
    def n_times(self) -> int:
        return self.duration_seconds // self.step_seconds

    @property
    def end_utc(self) -> datetime:
        """Exclusive observation end."""
        return self.start_utc + timedelta(seconds=self.duration_seconds)

    @property
    def time_offsets_seconds(self) -> np.ndarray:
        """Integer seconds since start_utc, one per column (int64)."""
        return np.arange(self.n_times, dtype=np.int64) * np.int64(self.step_seconds)

    @property
    def timestamps_utc(self) -> List[datetime]:
        return [self.start_utc + timedelta(seconds=int(s)) for s in self.time_offsets_seconds]

    # -- frequency axis (matrix rows) -----------------------------------
    @property
    def channel_width_hz(self) -> float:
        return (self.band_max_hz - self.band_min_hz) / self.channel_count

    @property
    def frequencies_hz(self) -> np.ndarray:
        """Uniform channel centers in Hz (float64); band limits are edges."""
        i = np.arange(self.channel_count, dtype=np.float64)
        return self.band_min_hz + (i + 0.5) * self.channel_width_hz

    @property
    def matrix_shape(self) -> tuple:
        return (self.channel_count, self.n_times)

    # -- sections and paths ---------------------------------------------
    def section(self, name: str) -> Dict[str, Any]:
        return self.config[name]

    def resolve(self, relative: PathLike) -> Path:
        """Resolve a reference path relative to the repository root."""
        path = Path(relative)
        return path if path.is_absolute() else (self.repo_root / path)

    def reference_path(self, key: str) -> Path:
        return self.resolve(self.config["references"][key])

    @property
    def norad_ids(self) -> List[int]:
        return [int(n) for n in self.config["references"]["norad_ids"]]


def _validate(config: Dict[str, Any]) -> None:
    for key in ("time_grid", "observatory", "pointing", "beam", "receiver",
                "references", "geometry", "frequency_grid"):
        _require(key in config, f"Missing scenario section {key!r}")

    observatory = config["observatory"]
    _require(observatory.get("earth_model") == "WGS84",
             "Only observatory.earth_model = 'WGS84' is supported")
    _require(-90.0 <= float(observatory["latitude_deg"]) <= 90.0, "latitude out of range")

    pointing = config["pointing"]
    _require(pointing.get("mode") == "fixed_altaz",
             "Only pointing.mode = 'fixed_altaz' is supported")
    _require(-90.0 <= float(pointing["altitude_deg"]) <= 90.0, "pointing altitude out of range")

    beam = config["beam"]
    _require(beam.get("model") in ("circular_selection_gate", "chromatic_elliptical_gaussian"),
             "Unsupported beam.model")
    width = float(beam["full_width_deg"])
    _require(0.0 < width <= 360.0, f"beam.full_width_deg must be in (0, 360], got {width}")
    if beam.get("model") == "chromatic_elliptical_gaussian":
        _require(0 < float(beam.get("minor_width_deg", 0)) <= width,
                 "beam.minor_width_deg must be positive and <= full_width_deg")
        _require(float(beam.get("reference_frequency_hz", 0)) > 0,
                 "beam.reference_frequency_hz must be positive")
    _require(beam.get("refraction", "none") == "none", "Only unrefracted geometry is supported")

    geometry = config["geometry"]
    _require(geometry.get("propagator") == "SGP4", "Only the SGP4 propagator is supported")
    _require(geometry.get("position_kind") == "topocentric", "Only topocentric geometry is supported")

    ids = config["references"].get("norad_ids")
    _require(isinstance(ids, list) and len(ids) > 0, "references.norad_ids must be a non-empty list")
    _require(len(set(int(n) for n in ids)) == len(ids), "references.norad_ids contains duplicates")


def load_scenario(config_path: PathLike, repo_root: Optional[PathLike] = None) -> Scenario:
    """Load and validate a scenario JSON file.

    Args:
        config_path: Scenario JSON path.
        repo_root: Repository root for reference paths. Defaults to the first
            ancestor of the config file that contains pyproject.toml.

    Raises:
        ScenarioError: For invalid grids or unsupported options.
    """
    config_path = Path(config_path).resolve()
    with open(config_path, encoding="utf-8") as handle:
        config = json.load(handle)
    _validate(config)

    grid = config["time_grid"]
    _require(grid.get("time_scale", "UTC") == "UTC", "time_grid.time_scale must be UTC")
    start = parse_utc(grid["start_utc"])
    step = _positive_int(grid["step_seconds"], "time_grid.step_seconds")
    duration = _positive_int(grid["duration_seconds"], "time_grid.duration_seconds")
    _require(duration % step == 0,
             f"duration_seconds ({duration}) must be an exact multiple of step_seconds ({step})")

    channels = _positive_int(config["frequency_grid"]["channel_count"],
                             "frequency_grid.channel_count")
    band_min = float(config["receiver"]["band_min_hz"])
    band_max = float(config["receiver"]["band_max_hz"])
    _require(band_max > band_min, "receiver.band_max_hz must exceed band_min_hz")

    root = Path(repo_root).resolve() if repo_root is not None else find_repo_root(config_path.parent)
    return Scenario(
        config=config,
        config_path=config_path,
        repo_root=root,
        start_utc=start,
        step_seconds=step,
        duration_seconds=duration,
        channel_count=channels,
        band_min_hz=band_min,
        band_max_hz=band_max,
    )


def verify_reference_hashes(scenario: Scenario) -> Dict[str, str]:
    """Check reference files against the SHA-256 values in provenance.json.

    Returns:
        Mapping of file name to verified SHA-256.

    Raises:
        ScenarioError: If a listed file is missing or its hash differs.
    """
    provenance_path = scenario.reference_path("provenance_file")
    with open(provenance_path, encoding="utf-8") as handle:
        expected = json.load(handle).get("files_sha256", {})
    verified: Dict[str, str] = {}
    for name, digest in expected.items():
        path = provenance_path.parent / name
        _require(path.is_file(), f"Reference file missing: {path}")
        actual = sha256_file(path)
        _require(actual == digest, f"SHA-256 mismatch for {path}: {actual} != {digest}")
        verified[name] = actual
    return verified
