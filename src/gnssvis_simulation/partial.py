"""Intermediate results: geometry, noise background and satellite emission.

Usage::

    python -m gnssvis_simulation.partial --config examples/scenarios/offline.json

Writes figures, ``noise.npz``/``signal.npz`` and their JSON metadata to ``partial_results/``
(relative to the repository root unless ``--output`` is given).
"""

from __future__ import annotations

import argparse
import time
from datetime import timedelta
from pathlib import Path
from typing import List, Optional

import numpy as np

from .config import Scenario, format_utc, load_scenario, verify_reference_hashes
from .diagnostics import (ascii_timeline, contribution_summary, format_table, noise_statistics,
                          passage_rows)
from .emission import (carrier_in_band, compute_signal_power, load_signal_catalog,
                       profiles_from_scenario)
from .figures import (plot_geometry, plot_noise_statistics, plot_power_waterfall,
                      plot_signal_preview)
from .geometry import compute_geometry, geometry_passages
from .noise import generate_noise, noise_model_from_scenario, scenario_seed
from .products import (array_sha256, base_metadata, coordinate_arrays, load_npz, save_npz,
                       write_json)


def _section(title: str) -> None:
    print(f"\n{'=' * 78}\n {title}\n{'=' * 78}")


def _fmt(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}g}"


def report_scenario(scenario: Scenario) -> None:
    _section("1. SCENARIO AND GRIDS")
    site, point, beam = (scenario.section(k) for k in ("observatory", "pointing", "beam"))
    n_f, n_t = scenario.matrix_shape
    mb = n_f * n_t * 4 / 2 ** 20
    print(format_table(["item", "value"], [
        ["scenario", scenario.config.get("scenario_id")],
        ["site (lat, lon, h)", f"{site['latitude_deg']} deg, {site['longitude_deg']} deg, "
                               f"{site['elevation_m']} m (WGS84)"],
        ["pointing", f"alt {point['altitude_deg']} deg, az {point['azimuth_deg']} deg"],
        ["beam gate", f"separation <= {beam['full_width_deg'] / 2:g} deg, "
                      f"alt >= {beam['minimum_altitude_deg']} deg"],
        ["time", f"{format_utc(scenario.start_utc)} .. {format_utc(scenario.end_utc)} (end excluded)"],
        ["time step / samples", f"{scenario.step_seconds} s / {n_t} (last {format_utc(scenario.timestamps_utc[-1])})"],
        ["band", f"{scenario.band_min_hz / 1e6:g} .. {scenario.band_max_hz / 1e6:g} MHz"],
        ["channels / width", f"{n_f} / {scenario.channel_width_hz:.1f} Hz"],
        ["first / last center", f"{scenario.frequencies_hz[0]:.2f} / {scenario.frequencies_hz[-1]:.2f} Hz"],
        ["matrix (freq, time)", f"{(n_f, n_t)}, float32 = {mb:.1f} MiB per power array"],
    ], numeric=False))


def report_geometry(scenario: Scenario, geometry) -> None:
    _section("2. REFERENCE INPUTS AND GEOMETRY (offline SGP4, local TLEs)")
    hashes = verify_reference_hashes(scenario)
    print("provenance SHA-256 verified:", ", ".join(sorted(hashes)))
    rows = []
    for r, s in enumerate(geometry.satellites):
        offs = [abs((t - s.tle_epoch_utc).total_seconds()) / 3600
                for t in (scenario.start_utc, scenario.end_utc)]
        j = int(np.argmin(geometry.separation_deg[r]))
        rows.append([s.norad_id, f"PRN{s.prn:02d}", s.block, format_utc(s.tle_epoch_utc)[:19],
                     f"{max(offs):.1f}", f"{geometry.alt_deg[r].max():.2f}",
                     f"{geometry.separation_deg[r, j]:.2f}",
                     format_utc(geometry.start_utc + timedelta(seconds=int(geometry.time_offsets_seconds[j])))[11:19],
                     f"{100 * geometry.above_horizon[r].mean():.0f}%",
                     int(geometry.selected[r].sum())])
    print(format_table(["NORAD", "PRN", "block", "TLE epoch (UTC)", "max|dt| h", "max alt",
                        "min sep", "at UTC", "visible", "in-beam"], rows))
    print()
    print(ascii_timeline(geometry))
    passages = geometry_passages(geometry)
    print(f"\nPassages ({len(passages)}):")
    if passages:
        print(format_table(["NORAD", "PRN", "block", "entry (UTC)", "exit excl. (UTC)", "samples",
                            "dur s", "min sep", "at UTC", "range km", "clipped in/out"],
                           passage_rows(geometry, passages)))
        for p in passages:
            r = geometry.norad_ids.index(p.norad_id)
            edge = [p.start_index - 1, p.start_index, p.stop_index - 1, p.stop_index]
            edge = [j for j in edge if 0 <= j < geometry.selected.shape[1]]
            print(f"  NORAD {p.norad_id} boundary samples (sep deg): " + ", ".join(
                f"{format_utc(geometry.timestamps_utc[j])[11:19]}={geometry.separation_deg[r, j]:.5f}"
                f"{'(in)' if geometry.selected[r, j] else '(out)'}" for j in edge))
    else:
        print("  none")


def report_noise(scenario: Scenario, noise: np.ndarray, seed: Optional[int]) -> dict:
    model = noise_model_from_scenario(scenario)
    _section("3. WHITE-NOISE BACKGROUND")
    print(f"model: Gamma(shape=M={model.averages}, scale=P0/M={model.scale_parameter:g}), "
          f"P0={model.mean_power:g}; numpy.random.default_rng(seed={seed})")
    print(f"array: {noise.dtype} {noise.shape} (frequency, time), "
          f"min {noise.min():.4f}, max {noise.max():.4f}")
    t0 = time.perf_counter()
    st = noise_statistics(noise, model)
    print(f"\nMoments ({st['n_cells']} cells; z = (measured - expected) / standard error; "
          f"|z| < 4 expected), computed in {time.perf_counter() - t0:.1f} s:")
    rows = [[name, _fmt(c["measured"]), _fmt(c["expected"]), _fmt(c["se"], 2), f"{c['z']:+.2f}"]
            for name, c in st["checks"].items()]
    for axis, c in st["lag1_correlation"].items():
        rows.append([f"lag-1 corr ({axis})", _fmt(c["measured"], 3), "0", _fmt(c["se"], 2), f"{c['z']:+.2f}"])
    print(format_table(["statistic", "measured", "expected", "std err", "z"], rows))
    print(f"\nKolmogorov-Smirnov vs Gamma: D = {st['ks_test']['statistic']:.2e}, "
          f"p = {st['ks_test']['pvalue']:.3f}  (p > 0.001 expected)")
    cm, sm = st["channel_means"], st["sample_means"]
    print(f"per-channel means: sd {cm['sd_measured']:.5f} (expected {cm['sd_expected']:.5f}), "
          f"max |z| {cm['max_abs_z']:.2f} over {noise.shape[0]} channels")
    print(f"per-sample means:  sd {sm['sd_measured']:.5f} (expected {sm['sd_expected']:.5f}), "
          f"max |z| {sm['max_abs_z']:.2f} over {noise.shape[1]} samples")
    print("\nQuantiles of linear power:")
    print(format_table(["quantile", "measured", "Gamma theory"],
                       [[q, f"{v['measured']:.4f}", f"{v['expected']:.4f}"]
                        for q, v in st["quantiles"].items()]))
    print(f"\nNote: the dB display has mean {st['checks']['mean_db_re_p0']['expected']:.3f} dB, "
          f"not 0 dB, because E[log P] < log E[P]; dB range "
          f"{st['db_min']:.2f} .. {st['db_max']:.2f} dB.")
    return st


def report_emission(scenario: Scenario, geometry, signal: np.ndarray, contributions,
                    seed: Optional[int], noise_hash_before: str) -> dict:
    _section("4. SATELLITE EMISSION (signal only, linear relative power)")
    freqs = scenario.frequencies_hz
    rows = []
    for c in load_signal_catalog(scenario):
        inside = carrier_in_band(c.carrier_hz, scenario)
        if inside:
            i = int(np.argmin(np.abs(freqs - c.carrier_hz)))
            place = [i, f"{freqs[i] / 1e6:.6f}", f"{(freqs[i] - c.carrier_hz) / 1e3:+.2f}"]
        else:
            place = ["-", "-", "-"]
        rows.append([c.signal_id, f"{c.carrier_hz / 1e6:.3f}", "yes" if inside else "no (excluded)",
                     "|".join(c.applicable_blocks), *place])
    print(format_table(["signal", "carrier MHz", "in band", "blocks", "nearest ch",
                        "ch centre MHz", "offset kHz"], rows))
    print(format_table(["block", "P_ref", "r_ref km", "sigma kHz", "cutoff", "support"],
                       [[b, f"{p.peak_power:g}", f"{p.reference_range_m / 1e3:g}",
                         f"{p.sigma_hz / 1e3:g}", f"{p.cutoff_sigma:g} sigma",
                         f"+/-{p.cutoff_sigma * p.sigma_hz / 1e6:g} MHz"]
                        for b, p in profiles_from_scenario(scenario).items()]))
    print(f"angular gain: exp(-4 ln2 (sep/{scenario.section('beam')['full_width_deg']:g})^2), "
          f"= 0.5 at the {scenario.section('beam')['full_width_deg'] / 2:g} deg gate")

    summaries = [contribution_summary(t, scenario, geometry) for t in contributions]
    print(f"\nSatellite x in-band carrier terms: {len(contributions)} "
          f"({sum(s['active'] for s in summaries)} active, i.e. with passage samples)")
    print(format_table(["NORAD", "PRN", "signal", "active", "mask samples", "support ch"],
                       [[s["norad_id"], f"PRN{s['prn']:02d}", s["signal_id"], s["active"],
                         s["mask_samples"], s["support_channels"]] for s in summaries]))
    active = [s for s in summaries if s["active"]]
    for s in active:
        pk = s["peak"]
        print(f"\nNORAD {s['norad_id']} PRN{s['prn']:02d} {s['signal_id']} "
              f"({s['carrier_hz'] / 1e6:.2f} MHz):")
        print(format_table(["interval entry (UTC)", "exit excl. (UTC)", "samples", "dur s",
                            "at entry", "at last"],
                           [[iv["entry_utc"], iv["exit_utc_exclusive"], iv["samples"],
                             iv["duration_seconds"], f"{iv['entry_value']:.3f}",
                             f"{iv['last_value']:.3f}"] for iv in s["intervals"]]))
        print(f"  peak contribution {pk['value']:.4f} (= {10 * np.log10(pk['value']):.2f} dB re P0) "
              f"at {pk['utc']}, channel {pk['channel_index']} ({pk['channel_center_hz'] / 1e6:.6f} MHz)")
        print(f"  = P_ref x (r_ref/r)^2 {pk['range_scaling']:.4f} [r = {pk['range_m'] / 1e3:.1f} km] "
              f"x gain {pk['angular_gain']:.4f} [sep {pk['separation_deg']:.3f} deg] "
              f"x spectral {pk['spectral_weight']:.4f}")

    union = np.zeros(scenario.n_times, dtype=bool)
    for t in contributions:
        if t.active:
            union |= t.mask
    nonzero_cols = signal.any(axis=0)
    summed = np.zeros(signal.shape)
    for t in contributions:
        summed += t.power()
    checks = {
        "signal_dtype_shape": f"{signal.dtype} {signal.shape}",
        "nonzero_cells": int(np.count_nonzero(signal)),
        "nonzero_columns_equal_passage_masks": bool(np.array_equal(nonzero_cols, union)),
        "zero_outside_masks": bool(not signal[:, ~union].any()),
        "nonzero_rows_within_spectral_support": bool(not signal[~np.any(
            [t.support for t in contributions] or [np.zeros(signal.shape[0], bool)], axis=0)].any()),
        "sum_of_terms_max_abs_diff": float(np.max(np.abs(summed - signal))),
        # Regenerated after the emission step: the background must not change.
        "noise_regenerated_unchanged": array_sha256(generate_noise(scenario, seed=seed))
                                       == noise_hash_before,
        "max_signal": float(signal.max()),
    }
    print()
    print(format_table(["check", "result"], [[k, v] for k, v in checks.items()], numeric=False))
    return {"contributions": summaries, "checks": checks}


def report_reproducibility(scenario: Scenario, noise: np.ndarray, seed: Optional[int],
                           geometry, npz_path: Path, signal: np.ndarray) -> dict:
    _section("5. REPRODUCIBILITY")
    again = generate_noise(scenario, seed=seed)
    other_seed = (seed or 0) + 1
    other = generate_noise(scenario, seed=other_seed)
    geometry_other = compute_geometry(scenario)
    loaded = load_npz(npz_path)
    result = {
        "same_seed_identical": bool(np.array_equal(noise, again)),
        "other_seed_differs": bool(not np.array_equal(noise, other)),
        "other_seed_correlation": float(np.corrcoef(noise.ravel(), other.ravel())[0, 1]),
        "geometry_unchanged": bool(np.array_equal(geometry.alt_deg, geometry_other.alt_deg)
                                   and np.array_equal(geometry.selected, geometry_other.selected)),
        "npz_roundtrip_identical": bool(np.array_equal(loaded["noise_power"], noise)
                                        and loaded["noise_power"].dtype == np.float32),
        "signal_deterministic": bool(np.array_equal(
            signal, compute_signal_power(scenario, geometry_other)[0])),
    }
    print(format_table(["check", "result"], [
        [f"seed {seed} twice -> bit-identical", result["same_seed_identical"]],
        [f"seed {other_seed} -> different background", result["other_seed_differs"]],
        [f"correlation seed {seed} vs {other_seed}", f"{result['other_seed_correlation']:+.4f}"],
        ["geometry independent of noise run", result["geometry_unchanged"]],
        ["NPZ reload bit-identical (float32)", result["npz_roundtrip_identical"]],
        ["signal recomputed -> identical (no RNG)", result["signal_deterministic"]],
        [f"sha256(noise_power) seed {seed}", array_sha256(noise)[:16] + "..."],
        [f"sha256(noise_power) seed {other_seed}", array_sha256(other)[:16] + "..."],
    ], numeric=False))
    print("NPZ contents:", ", ".join(f"{k} {v.dtype} {v.shape}" for k, v in loaded.items()))
    return result


def run(config: Path, output: Optional[Path] = None, seed: Optional[int] = None) -> List[Path]:
    scenario = load_scenario(config)
    output = Path(output) if output is not None else scenario.repo_root / "partial_results"
    eff_seed = scenario_seed(scenario, seed)

    report_scenario(scenario)
    t0 = time.perf_counter()
    geometry = compute_geometry(scenario)
    print(f"\n[geometry computed in {time.perf_counter() - t0:.2f} s]")
    report_geometry(scenario, geometry)

    t0 = time.perf_counter()
    noise = generate_noise(scenario, seed=eff_seed)
    print(f"\n[noise generated in {time.perf_counter() - t0:.2f} s]")
    stats = report_noise(scenario, noise, eff_seed)

    noise_hash = array_sha256(noise)
    npz_path = save_npz(output / "noise.npz", {"noise_power": noise, **coordinate_arrays(scenario)},
                        scenario)

    t0 = time.perf_counter()
    signal, contributions = compute_signal_power(scenario, geometry)
    print(f"\n[emission computed in {time.perf_counter() - t0:.2f} s]")
    emission = report_emission(scenario, geometry, signal, contributions, eff_seed, noise_hash)
    signal_npz = save_npz(output / "signal.npz",
                          {"signal_power": signal, **coordinate_arrays(scenario)}, scenario)

    repro = report_reproducibility(scenario, noise, eff_seed, geometry, npz_path, signal)

    figures = [plot_geometry(scenario, geometry, output / "01_geometry_separation.png")]
    limits = plot_power_waterfall(
        scenario, noise, output / "02_noise_waterfall.png",
        f"Background only: Gamma white noise (M={noise_model_from_scenario(scenario).averages}, "
        f"seed {eff_seed})")
    figures.append(output / "02_noise_waterfall.png")
    figures.append(plot_noise_statistics(scenario, noise, noise_model_from_scenario(scenario),
                                         output / "03_noise_statistics.png"))
    figures.append(plot_signal_preview(scenario, geometry, signal, contributions,
                                       output / "04_signal_waterfall.png"))

    metadata = base_metadata(scenario)
    metadata.update({
        "stage": "noise_background_only",
        "npz_file": npz_path.name,
        "npz_arrays": {"noise_power": "float32 (frequency, time), relative linear power",
                       "frequencies_hz": "float64 channel centers, Hz",
                       "time_offsets_seconds": "int64 seconds since start_utc"},
        "randomness": {"generator": "numpy.random.default_rng", "seed": eff_seed,
                       "draw": "Generator.gamma(shape, scale, size=(n_freq, n_time)) in float64, "
                               "C order, cast to float32"},
        "noise_model": noise_model_from_scenario(scenario).describe(),
        "noise_statistics": stats,
        "reproducibility_checks": repro,
        "array_sha256": {"noise_power": array_sha256(noise)},
        "display": {"quantity": "10*log10(noise_power/P0)", **limits},
    })
    meta_path = write_json(output / "noise_metadata.json", metadata)

    signal_meta = base_metadata(scenario)
    signal_meta.update({
        "stage": "satellite_emission_signal_only",
        "npz_file": signal_npz.name,
        "npz_arrays": {"signal_power": "float32 (frequency, time), relative linear power",
                       "frequencies_hz": "float64 channel centers, Hz",
                       "time_offsets_seconds": "int64 seconds since start_utc"},
        "model": {
            "formula": "mask(t) * P_ref * (r_ref/r(t))**2 * exp(-4 ln2 (sep/full_width)**2) "
                       "* exp(-0.5 ((f-f_c)/sigma)**2), zero beyond cutoff*sigma",
            "profiles": {b: vars(p) for b, p in profiles_from_scenario(scenario).items()},
            "mask": "geometry.selected = in beam gate and above horizon, per sample",
            "combination": "sum of satellite/carrier terms in float64, cast to float32",
            "randomness": "none (deterministic)",
        },
        "in_band_carriers": [c.signal_id for c in load_signal_catalog(scenario)
                             if carrier_in_band(c.carrier_hz, scenario)],
        "emission": emission,
        "array_sha256": {"signal_power": array_sha256(signal)},
    })
    signal_meta_path = write_json(output / "signal_metadata.json", signal_meta)

    _section("6. FILES WRITTEN")
    written = [npz_path, meta_path, signal_npz, signal_meta_path, *figures]
    for path in written:
        print(f"  {path}  ({path.stat().st_size / 1024:.0f} KiB)")
    return written


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True, type=Path, help="scenario JSON")
    parser.add_argument("--output", type=Path, default=None,
                        help="output directory (default: <repo>/partial_results)")
    parser.add_argument("--seed", type=int, default=None, help="override the configured seed")
    args = parser.parse_args(argv)
    run(args.config, args.output, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
