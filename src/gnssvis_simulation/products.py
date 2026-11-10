"""Array contract, NPZ storage and JSON metadata for simulation products.

Power arrays: float32, shape (frequency, time), relative linear power.
``frequencies_hz``: float64 channel centers (rows). ``time_offsets_seconds``:
int64 seconds since the UTC start (columns). NPZ files contain plain numeric
arrays only and are read back with ``allow_pickle=False``.
"""

from __future__ import annotations

import hashlib
import json
import platform
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Dict, Mapping

import numpy as np

from .config import Scenario, format_utc, sha256_file

POWER_SUFFIX = "_power"
AXIS_DTYPES = {"frequencies_hz": np.float64, "time_offsets_seconds": np.int64}
VERSIONED_PACKAGES = ("numpy", "scipy", "skyfield", "sgp4", "astropy", "matplotlib", "gnssvis")


def validate_arrays(arrays: Mapping[str, np.ndarray], scenario: Scenario) -> None:
    """Check dtypes, ranks and shapes against the matrix contract."""
    n_f, n_t = scenario.matrix_shape
    for name in AXIS_DTYPES:
        if name not in arrays:
            raise ValueError(f"Missing coordinate array {name!r}")
    freqs, times = arrays["frequencies_hz"], arrays["time_offsets_seconds"]
    if freqs.dtype != np.float64 or freqs.shape != (n_f,):
        raise ValueError(f"frequencies_hz must be float64 ({n_f},), got {freqs.dtype} {freqs.shape}")
    if times.dtype != np.int64 or times.shape != (n_t,):
        raise ValueError(f"time_offsets_seconds must be int64 ({n_t},), got {times.dtype} {times.shape}")
    for name, arr in arrays.items():
        if name in AXIS_DTYPES:
            continue
        if not name.endswith(POWER_SUFFIX):
            raise ValueError(f"Unexpected array {name!r}; power arrays end with {POWER_SUFFIX!r}")
        if arr.dtype != np.float32 or arr.shape != (n_f, n_t):
            raise ValueError(f"{name} must be float32 {(n_f, n_t)}, got {arr.dtype} {arr.shape}")
        if not np.all(np.isfinite(arr)) or np.any(arr < 0):
            raise ValueError(f"{name} must be finite, non-negative linear power")


def coordinate_arrays(scenario: Scenario) -> Dict[str, np.ndarray]:
    return {"frequencies_hz": scenario.frequencies_hz,
            "time_offsets_seconds": scenario.time_offsets_seconds}


def save_npz(path: Path, arrays: Mapping[str, np.ndarray], scenario: Scenario) -> Path:
    """Validate and save arrays to an uncompressed NPZ file."""
    validate_arrays(arrays, scenario)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **{k: np.ascontiguousarray(v) for k, v in arrays.items()})
    return path


def load_npz(path: Path) -> Dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def array_sha256(arr: np.ndarray) -> str:
    """SHA-256 of the array bytes (C order) plus dtype and shape."""
    arr = np.ascontiguousarray(arr)
    digest = hashlib.sha256(f"{arr.dtype.str}{arr.shape}".encode())
    digest.update(arr.tobytes())
    return digest.hexdigest()


def dependency_versions() -> Dict[str, str]:
    versions = {"python": platform.python_version()}
    for package in VERSIONED_PACKAGES:
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def grid_metadata(scenario: Scenario) -> Dict[str, Any]:
    return {
        "matrix_shape": list(scenario.matrix_shape),
        "dimension_order": ["frequency", "time"],
        "frequency": {
            "array": "frequencies_hz", "dtype": "float64", "unit": "Hz",
            "placement": "uniform channel centers; band limits are channel edges",
            "band_min_hz": scenario.band_min_hz, "band_max_hz": scenario.band_max_hz,
            "channel_count": scenario.channel_count,
            "channel_width_hz": scenario.channel_width_hz,
            "first_center_hz": float(scenario.frequencies_hz[0]),
            "last_center_hz": float(scenario.frequencies_hz[-1]),
        },
        "time": {
            "array": "time_offsets_seconds", "dtype": "int64",
            "unit": "seconds since start_utc",
            "start_utc": format_utc(scenario.start_utc),
            "end_utc_exclusive": format_utc(scenario.end_utc),
            "step_seconds": scenario.step_seconds,
            "sample_count": scenario.n_times,
            "last_sample_utc": format_utc(scenario.timestamps_utc[-1]),
        },
        "power": {"dtype": "float32", "unit": "relative linear power (not dBm/W)",
                  "combination": "sum in linear power; dB only for display"},
    }


def base_metadata(scenario: Scenario) -> Dict[str, Any]:
    return {
        "scenario_id": scenario.config.get("scenario_id"),
        "config_file": str(scenario.config_path),
        "config_sha256": sha256_file(scenario.config_path),
        "effective_config": scenario.config,
        "grid": grid_metadata(scenario),
        "dependency_versions": dependency_versions(),
    }


def write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")
    return path
