"""End-to-end offline simulation: geometry -> noise + emission -> outputs.

``total_power = noise_power + signal_power`` in linear relative power; only the
figure converts to ``10*log10(total_power / P0)``. All parameters come from the
scenario file; only the output directory may be overridden.
"""

from __future__ import annotations

import csv
import sys
import time
from dataclasses import dataclass, field
from datetime import timedelta
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .config import Scenario, format_utc, load_scenario, verify_reference_hashes
from .emission import (Contribution, carrier_in_band, compute_signal_power, load_signal_catalog,
                       profiles_from_scenario, satellite_signal_ids)
from .figures import plot_final_waterfall
from .geometry import GeometryResult, Passage, compute_geometry, geometry_passages
from .noise import generate_noise, noise_model_from_scenario, scenario_seed
from .products import (array_sha256, base_metadata, coordinate_arrays, save_npz, write_json)

SAMPLE_COLUMNS = ["sample_index", "time_offset_seconds", "timestamp_utc", "norad_id",
                  "altitude_deg", "azimuth_deg", "separation_deg", "range_m",
                  "above_horizon", "in_beam", "in_passage", "signal_ids"]
PASSAGE_COLUMNS = ["passage_index", "norad_id", "name", "prn", "block", "entry_utc",
                   "exit_utc_exclusive", "entry_offset_seconds", "exit_offset_seconds",
                   "duration_seconds", "sample_count", "entry_clipped", "exit_clipped",
                   "min_separation_deg", "min_separation_utc", "peak_signal_power",
                   "peak_signal_utc", "signal_ids"]


@dataclass
class SimulationResult:
    scenario: Scenario
    geometry: GeometryResult
    passages: List[Passage]
    noise_power: np.ndarray
    signal_power: np.ndarray
    total_power: np.ndarray
    contributions: List[Contribution]
    signal_ids: Dict[int, List[str]]
    files: Dict[str, Path] = field(default_factory=dict)
    checks: Dict[str, object] = field(default_factory=dict)
    display: Dict[str, float] = field(default_factory=dict)


def _b(value) -> str:
    """Documented lowercase booleans for CSV."""
    return "true" if bool(value) else "false"


def simulate(scenario: Scenario) -> SimulationResult:
    """Compute all arrays (no files). Noise and signal are independent."""
    if 'realistic_emission' in scenario.config:
        raise ValueError('Use simulate_realistic or python -m gnssvis_simulation.realistic '
                         'for a wideband solar scenario')
    geometry = compute_geometry(scenario)
    passages = geometry_passages(geometry)
    noise = generate_noise(scenario)
    signal, contributions = compute_signal_power(scenario, geometry)
    total = noise + signal  # float32 + float32, elementwise, linear power
    carriers = load_signal_catalog(scenario)
    signal_ids = {s.norad_id: satellite_signal_ids(s, carriers, scenario)
                  for s in geometry.satellites}
    return SimulationResult(scenario, geometry, passages, noise, signal, total, contributions,
                            signal_ids)


# -- tables -------------------------------------------------------------------

def write_samples_csv(result: SimulationResult, path: Path) -> int:
    """One row per (timestamp, satellite), ordered by time then configured ID."""
    g, sc = result.geometry, result.scenario
    stamps = [format_utc(t) for t in g.timestamps_utc]
    rows = 0
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=sc.config["output"].get("csv_delimiter", ","),
                            lineterminator="\n")
        writer.writerow(SAMPLE_COLUMNS)
        ids = [";".join(result.signal_ids[s.norad_id]) for s in g.satellites]
        for j, offset in enumerate(g.time_offsets_seconds):
            for r, sat in enumerate(g.satellites):
                writer.writerow([j, int(offset), stamps[j], sat.norad_id,
                                 f"{g.alt_deg[r, j]:.6f}", f"{g.az_deg[r, j]:.6f}",
                                 f"{g.separation_deg[r, j]:.6f}", f"{g.range_m[r, j]:.3f}",
                                 _b(g.above_horizon[r, j]), _b(g.in_beam[r, j]),
                                 _b(g.selected[r, j]), ids[r]])
                rows += 1
    return rows


def passage_records(result: SimulationResult) -> List[Dict[str, object]]:
    g, sc = result.geometry, result.scenario
    index = {n: i for i, n in enumerate(g.norad_ids)}
    records = []
    for k, p in enumerate(result.passages):
        r = index[p.norad_id]
        sat = g.satellites[r]
        sl = slice(p.start_index, p.stop_index)
        j_sep = p.start_index + int(np.argmin(g.separation_deg[r, sl]))
        sat_signal = np.zeros(sc.n_times)
        for term in result.contributions:
            if term.norad_id == p.norad_id:
                sat_signal += term.envelope * term.spectral.max()
        j_pk = p.start_index + int(np.argmax(sat_signal[sl]))
        peak = float(sat_signal[j_pk])
        records.append({
            "passage_index": k, "norad_id": p.norad_id, "name": sat.name, "prn": sat.prn,
            "block": sat.block,
            "entry_utc": format_utc(sc.start_utc + timedelta(seconds=p.entry_offset_seconds)),
            "exit_utc_exclusive": format_utc(sc.start_utc + timedelta(seconds=p.exit_offset_seconds)),
            "entry_offset_seconds": p.entry_offset_seconds,
            "exit_offset_seconds": p.exit_offset_seconds,
            "duration_seconds": p.duration_seconds, "sample_count": p.sample_count,
            "entry_clipped": p.entry_clipped, "exit_clipped": p.exit_clipped,
            "min_separation_deg": round(float(g.separation_deg[r, j_sep]), 6),
            "min_separation_utc": format_utc(g.timestamps_utc[j_sep]),
            "peak_signal_power": round(peak, 6) if peak > 0 else 0.0,
            "peak_signal_utc": format_utc(g.timestamps_utc[j_pk]) if peak > 0 else "",
            "signal_ids": ";".join(result.signal_ids[p.norad_id]),
        })
    return records


def write_passages_csv(records: List[Dict[str, object]], path: Path, delimiter: str = ",") -> int:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PASSAGE_COLUMNS, delimiter=delimiter,
                                lineterminator="\n")
        writer.writeheader()  # header is kept even with no passages
        for rec in records:
            writer.writerow({k: _b(v) if isinstance(v, (bool, np.bool_)) else v
                             for k, v in rec.items()})
    return len(records)


# -- checks and metadata ----------------------------------------------------------

def consistency_checks(result: SimulationResult) -> Dict[str, object]:
    sc, g = result.scenario, result.geometry
    union = np.zeros(sc.n_times, dtype=bool)
    for term in result.contributions:
        if term.active:
            union |= term.mask
    return {
        "shapes_equal_matrix_contract": all(a.shape == sc.matrix_shape and a.dtype == np.float32
                                            for a in (result.noise_power, result.signal_power,
                                                      result.total_power)),
        "total_equals_noise_plus_signal": bool(np.array_equal(
            result.total_power, result.noise_power + result.signal_power)),
        "signal_zero_outside_passage_masks": bool(not result.signal_power[:, ~union].any()),
        "signal_columns_equal_active_masks": bool(np.array_equal(
            result.signal_power.any(axis=0), union)),
        "active_masks_subset_of_geometry_selection": bool(not (union & ~g.selected.any(axis=0)).any()),
        "noise_outside_signal_equals_total": bool(np.array_equal(
            result.total_power[result.signal_power == 0],
            result.noise_power[result.signal_power == 0])),
    }


def ut1_metadata(scenario: Scenario) -> Dict[str, object]:
    from .geometry import offline_sky_analyzer, skyfield_times

    t = skyfield_times(scenario, offline_sky_analyzer(scenario))
    return {"source": "Skyfield bundled IERS table (builtin=True); no external update",
            "ut1_minus_utc_seconds_at_start": float(t.dut1[0]),
            "ut1_minus_utc_seconds_at_last_sample": float(t.dut1[-1])}


def build_metadata(result: SimulationResult, display: Dict[str, float],
                   files: Dict[str, Path], rows: Dict[str, int], elapsed: float) -> Dict:
    sc, g = result.scenario, result.geometry
    meta = base_metadata(sc)
    emission_cfg = sc.section("emission")
    meta.update({
        "description": sc.config.get("purpose"),
        "command": " ".join([Path(sys.executable).name, "-m", "gnssvis_simulation", *sys.argv[1:]])
                   if sys.argv and sys.argv[0].endswith("__main__.py") else None,
        "runtime_seconds": round(elapsed, 2),
        "reference_hashes_sha256": verify_reference_hashes(sc),
        "satellites": [{"norad_id": s.norad_id, "name": s.name, "prn": s.prn, "block": s.block,
                        "constellation": s.constellation,
                        "tle_epoch_utc": format_utc(s.tle_epoch_utc),
                        "max_abs_epoch_offset_hours": round(max(
                            abs((sc.start_utc - s.tle_epoch_utc).total_seconds()),
                            abs((sc.end_utc - s.tle_epoch_utc).total_seconds())) / 3600, 3),
                        "in_band_signal_ids": result.signal_ids[s.norad_id],
                        "passage_samples": int(g.selected[i].sum())}
                       for i, s in enumerate(g.satellites)],
        "earth_orientation": ut1_metadata(sc),
        "randomness": {"generator": "numpy.random.default_rng", "seed": scenario_seed(sc),
                       "consumer": "noise background only; emission is deterministic"},
        "files": {k: p.name for k, p in files.items()},
        "npz_arrays": {
            "noise_power": "float32 (frequency, time), relative linear power",
            "signal_power": "float32 (frequency, time), relative linear power",
            "total_power": "float32 (frequency, time) = noise_power + signal_power",
            "frequencies_hz": "float64 (frequency,), channel centres in Hz",
            "time_offsets_seconds": "int64 (time,), seconds since time_grid.start_utc",
        },
        "array_sha256": {k: array_sha256(v) for k, v in (
            ("noise_power", result.noise_power), ("signal_power", result.signal_power),
            ("total_power", result.total_power))},
        "counts": {"sample_rows": rows["samples"], "passage_rows": rows["passages"],
                   "satellites": len(g.satellites), "time_samples": sc.n_times,
                   "frequency_channels": sc.channel_count,
                   "active_contributions": sum(t.active for t in result.contributions),
                   "nonzero_signal_cells": int(np.count_nonzero(result.signal_power))},
        "in_band_carriers": [c.signal_id for c in load_signal_catalog(sc)
                             if carrier_in_band(c.carrier_hz, sc)],
        "models": {
            "noise": {**noise_model_from_scenario(sc).describe(),
                      "draw": "Generator.gamma in float64, C order, cast to float32"},
            "emission": {
                "formula": "mask(t) * P_ref * (r_ref/r)**2 * exp(-4 ln2 (sep/full_width)**2) "
                           "* exp(-0.5 ((f-f_c)/sigma)**2), zero beyond cutoff*sigma",
                "profiles": {b: vars(p) for b, p in profiles_from_scenario(sc).items()},
                "calibration": emission_cfg.get("calibration"),
                "not_modelled": emission_cfg.get("measurement_models_excluded"),
            },
            "geometry": {"propagator": "SGP4 (Skyfield EarthSatellite)", "observer": "WGS84",
                         "frame": "unrefracted topocentric alt/az", "sampling": "no subsample interpolation"},
        },
        "csv_conventions": {"booleans": "lowercase true/false", "signal_ids": "semicolon-separated",
                            "timestamps": "ISO 8601 UTC with Z", "samples_order": "time, then configured NORAD order",
                            "sample_columns": SAMPLE_COLUMNS, "passage_columns": PASSAGE_COLUMNS},
        "display": {"quantity": "10*log10(total_power / noise.mean_power_relative)",
                    "unit": "dB relative to configured mean background",
                    "color_limits": "full finite range", **display},
        "checks": result.checks,
    })
    return meta


# -- entry point ----------------------------------------------------------------

def run_simulation(config: Path, output: Optional[Path] = None, verbose: bool = True
                   ) -> SimulationResult:
    """Run the scenario and write PNG, NPZ, CSV and JSON outputs."""
    t0 = time.perf_counter()
    scenario = load_scenario(config)
    out_cfg = scenario.config["output"]
    output = Path(output) if output is not None else scenario.resolve(out_cfg["default_directory"])
    output.mkdir(parents=True, exist_ok=True)

    result = simulate(scenario)
    result.checks = consistency_checks(result)

    files = {
        "waterfall_npz": output / out_cfg["waterfall_file"],
        "samples_csv": output / out_cfg["samples_file"],
        "passages_csv": output / out_cfg["passages_file"],
        "figure_png": output / out_cfg["figure_file"],
        "metadata_json": output / out_cfg["metadata_file"],
    }
    save_npz(files["waterfall_npz"], {"noise_power": result.noise_power,
                                      "signal_power": result.signal_power,
                                      "total_power": result.total_power,
                                      **coordinate_arrays(scenario)}, scenario)
    rows = {"samples": write_samples_csv(result, files["samples_csv"]),
            "passages": write_passages_csv(passage_records(result), files["passages_csv"],
                                           out_cfg.get("csv_delimiter", ","))}
    result.display = plot_final_waterfall(scenario, result.geometry, result.passages,
                                          result.total_power, files["figure_png"],
                                          result.signal_ids)
    write_json(files["metadata_json"], build_metadata(result, result.display, files, rows,
                                                      time.perf_counter() - t0))
    result.files = files
    if verbose:
        print_summary(result, rows, time.perf_counter() - t0)
    return result


def print_summary(result: SimulationResult, rows: Dict[str, int], elapsed: float) -> None:
    sc = result.scenario
    print(f"Scenario {sc.config.get('scenario_id')}: {format_utc(sc.start_utc)} .. "
          f"{format_utc(sc.end_utc)} (end excluded), step {sc.step_seconds} s")
    print(f"Arrays (frequency, time) = {sc.matrix_shape}, float32; seed {scenario_seed(sc)}")
    print(f"Passages: {len(result.passages)}")
    for rec in passage_records(result):
        print(f"  NORAD {rec['norad_id']} PRN{rec['prn']:02d} [{rec['signal_ids'] or 'no in-band carrier'}]"
              f": {rec['entry_utc']} -> {rec['exit_utc_exclusive']} ({rec['sample_count']} samples), "
              f"min sep {rec['min_separation_deg']:.2f} deg, peak signal {rec['peak_signal_power']:.3f}")
    print(f"Display range: {result.display['vmin_db']:.2f} .. {result.display['vmax_db']:.2f} dB re P0")
    failed = [k for k, v in result.checks.items() if v is not True]
    print(f"Consistency checks: {len(result.checks) - len(failed)}/{len(result.checks)} passed"
          + (f"; FAILED: {failed}" if failed else ""))
    print(f"Rows: samples.csv {rows['samples']}, passages.csv {rows['passages']}")
    for key, path in result.files.items():
        print(f"  {key:14s} {path}")
    print(f"Done in {elapsed:.1f} s")
