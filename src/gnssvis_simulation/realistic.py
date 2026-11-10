"""Optional wideband GNSS + Sun simulation; imports cause no network I/O.

Orbital propagation remains SGP4. This module upgrades *received emissions*:
radial Doppler, channel-averaged spread-spectrum PSD, inverse-square loss and
an elliptical, frequency-dependent beam. It does not invent frequency chirps.
All powers are illustrative linear units; solar amplitudes are not measured.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.special import sici
from skyfield.api import Loader, load_file

from .config import (
    Scenario,
    ScenarioError,
    format_utc,
    load_scenario,
    sha256_file,
    verify_reference_hashes,
)
from .geometry import (
    GeometryResult,
    angular_separation_deg,
    compute_geometry,
    find_passages,
    offline_sky_analyzer,
    skyfield_times,
)
from .beam import beam_offsets_deg, beam_gain
from .noise import generate_noise
from .products import base_metadata, validate_arrays, array_sha256

C_M_S = 299792458.0


def doppler_shift_hz(carrier_hz, range_rate_m_s):
    """First-order one-way Doppler: recession shifts the carrier down."""
    return -np.asarray(carrier_hz) * np.asarray(range_rate_m_s) / C_M_S


def _sinc2_integral(x):
    """Integral from zero to x of numpy.sinc(u)**2 (stable at zero)."""
    x = np.asarray(x, dtype=float)
    si = sici(2 * np.pi * x)[0]
    remainder = np.divide(
        np.cos(2 * np.pi * x) - 1, 2 * np.pi**2 * x, out=np.zeros_like(x), where=x != 0
    )
    return si / np.pi + remainder


def bpsk_channel_psd(frequencies_hz, channel_width_hz, center_hz, chip_rate_hz):
    """Exact boxcar-channel mean of unit-peak BPSK sinc² PSD.

    Finite channel integration preserves power even when Doppler is much
    smaller than a channel. Analytic integration resolves spectral nulls
    without creating aliasing from undersampled channel-center evaluations.
    """
    if channel_width_hz <= 0 or chip_rate_hz <= 0:
        raise ValueError("Channel width and chip rate must be positive")
    x = (np.asarray(frequencies_hz) - np.asarray(center_hz)) / chip_rate_hz
    half = channel_width_hz / (2 * chip_rate_hz)
    return np.maximum(
        (_sinc2_integral(x + half) - _sinc2_integral(x - half)) / (2 * half), 0
    )


def solar_ephemeris(scenario: Scenario, ephemeris_path: Optional[Path] = None):
    """Apparent, unrefracted topocentric Sun, DE421 + WGS84.

    With no explicit file, Loader downloads DE421 into the configured cache
    only when allow_download is true. An explicit missing file always fails.
    Bundled Skyfield timescale data avoid an implicit Earth-orientation fetch.
    """
    cfg = scenario.config["solar"]
    path = (
        Path(ephemeris_path)
        if ephemeris_path is not None
        else scenario.resolve(cfg["ephemeris_file"])
    )
    if not path.is_file():
        if ephemeris_path is not None or not cfg.get("allow_download", False):
            raise FileNotFoundError(f"Ephemeris missing: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        Loader(str(path.parent))("de421.bsp").close()
        if not path.is_file():
            raise FileNotFoundError(f"Download did not create {path}")
    planets = load_file(str(path))
    try:
        sky = offline_sky_analyzer(scenario)
        times = skyfield_times(scenario, sky)
        sun = (
            (planets["earth"] + sky.topos).at(times).observe(planets["sun"]).apparent()
        )
        alt, az, distance = sun.altaz()
        p = scenario.section("pointing")
        x, y = beam_offsets_deg(
            alt.degrees, az.degrees, p["altitude_deg"], p["azimuth_deg"]
        )
        return {
            "alt_deg": alt.degrees,
            "az_deg": az.degrees,
            "distance_m": distance.m,
            "separation_deg": angular_separation_deg(
                alt.degrees, az.degrees, p["altitude_deg"], p["azimuth_deg"]
            ),
            "offset_x_deg": x,
            "offset_y_deg": y,
            "ephemeris_sha256": sha256_file(path),
        }
    finally:
        planets.close()


@dataclass
class RealisticResult:
    scenario: Scenario
    geometry: GeometryResult
    noise_power: np.ndarray
    satellite_power: np.ndarray
    solar_power: np.ndarray
    total_power: np.ndarray
    solar: dict
    passages: list
    doppler_hz: np.ndarray
    half_power_mask: np.ndarray


def validate_realistic(scenario):
    """Reject parameters that would produce a plausible-looking invalid run."""
    if scenario.section("beam")["model"] != "chromatic_elliptical_gaussian":
        raise ScenarioError(
            "Realistic synthesis requires a chromatic_elliptical_gaussian beam"
        )
    cfg = scenario.config.get("realistic_emission", {})
    for name in (
        "reference_range_m",
        "peak_power_relative",
        "broad_chip_rate_hz",
        "narrow_chip_rate_hz",
    ):
        if not np.isfinite(float(cfg.get(name, 0))) or float(cfg.get(name, 0)) <= 0:
            raise ScenarioError(
                f"realistic_emission.{name} must be finite and positive"
            )
    if not 0 <= float(cfg.get("broad_fraction", -1)) <= 1:
        raise ScenarioError("realistic_emission.broad_fraction must be in [0,1]")
    for name in ("quiet_power_relative", "active_power_relative"):
        if (
            not np.isfinite(float(scenario.config["solar"].get(name, -1)))
            or float(scenario.config["solar"][name]) < 0
        ):
            raise ScenarioError(f"solar.{name} must be finite and nonnegative")
    solar = scenario.config["solar"]
    if not np.isfinite(float(solar.get("quiet_received_spectral_index", 0))):
        raise ScenarioError("solar.quiet_received_spectral_index must be finite")
    if (
        not np.isfinite(float(solar.get("electron_energy_index", 3)))
        or float(solar.get("electron_energy_index", 3)) <= 1
    ):
        raise ScenarioError("solar.electron_energy_index must be finite and > 1")


def simulate_realistic(scenario: Scenario, ephemeris_path: Optional[Path] = None):
    """Synthesize GPS L1/L2/L5 powers, returning components separately.

    L1/L2 use illustrative mixtures of narrow and broad BPSK PSDs. L5 is broad.
    No health/status assertion is made for any satellite. This is not an ICD
    composite waveform, EIRP calibration, ionospheric or multipath simulator.
    """
    validate_realistic(scenario)
    verify_reference_hashes(scenario)
    g = compute_geometry(scenario)
    beam = scenario.section("beam")
    p = scenario.section("pointing")
    cfg = scenario.config["realistic_emission"]
    x, y = beam_offsets_deg(g.alt_deg, g.az_deg, p["altitude_deg"], p["azimuth_deg"])
    mask = (
        beam_gain(x, y, beam["reference_frequency_hz"], beam) >= 0.5
    ) & g.above_horizon
    passages = [
        v
        for row, sat in enumerate(g.satellites)
        for v in find_passages(mask[row], scenario.step_seconds, sat.norad_id)
    ]
    f = scenario.frequencies_hz[:, None]
    width = scenario.channel_width_hz
    total_sat = np.zeros(scenario.matrix_shape, dtype=np.float64)
    shifts = doppler_shift_hz(1575.42e6, g.range_rate_m_s)
    carriers = [(1575.42e6, 1.0), (1227.60e6, 0.7), (1176.45e6, 0.85)]
    for row, sat in enumerate(g.satellites):
        gain = (
            beam_gain(x[row][None, :], y[row][None, :], f, beam) * g.above_horizon[row]
        )
        amp = (
            float(cfg["peak_power_relative"])
            * (float(cfg["reference_range_m"]) / g.range_m[row]) ** 2
        )
        for center, weight in carriers:
            # IIR/IIR-M do not provide L5. Restrict L5 to IIF and III blocks.
            if center == 1176.45e6 and sat.block not in ("IIF", "III"):
                continue
            if not scenario.band_min_hz <= center <= scenario.band_max_hz:
                continue
            shifted = center + doppler_shift_hz(center, g.range_rate_m_s[row])[None, :]
            broad = bpsk_channel_psd(
                f, width, shifted, float(cfg["broad_chip_rate_hz"])
            )
            frac = float(cfg["broad_fraction"])
            profile = (
                broad
                if center == 1176.45e6
                else (
                    (1 - frac)
                    * bpsk_channel_psd(
                        f, width, shifted, float(cfg["narrow_chip_rate_hz"])
                    )
                    + frac * broad
                )
            )
            total_sat += weight * profile * gain * amp[None, :]
    solar = solar_ephemeris(scenario, ephemeris_path)
    s = scenario.config["solar"]
    fref = float(beam["reference_frequency_hz"])
    # Relative received spectral power, rather than a claimed solar flux in SFU.
    # Quiet thermal continuum + optional active nonthermal continuum.
    # For optically thin synchrotron S_nu ∝ nu^(-(p-1)/2), p=electron index.
    exponent = -(float(s.get("electron_energy_index", 3)) - 1) / 2
    spectrum = (
        float(s["quiet_power_relative"])
        * (f / fref) ** float(s.get("quiet_received_spectral_index", 0))
        + float(s["active_power_relative"]) * (f / fref) ** exponent
    )
    solar_gain = beam_gain(
        solar["offset_x_deg"][None, :], solar["offset_y_deg"][None, :], f, beam
    )
    solar_power = (spectrum * solar_gain * (solar["alt_deg"] > 0)[None, :]).astype(
        np.float32
    )
    noise = generate_noise(scenario)
    satellite_power = total_sat.astype(np.float32)
    total = noise + satellite_power + solar_power
    return RealisticResult(
        scenario,
        g,
        noise,
        satellite_power,
        solar_power,
        total,
        solar,
        passages,
        shifts,
        mask,
    )


def save_realistic(result: RealisticResult, output: Path):
    """Save compressed arrays, compact geometry, passage CSV and provenance."""
    from .realistic_figures import plot_realistic_figures

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    r = result
    sc = r.scenario
    g = r.geometry
    arrays = dict(
        noise_power=r.noise_power,
        satellite_power=r.satellite_power,
        solar_power=r.solar_power,
        total_power=r.total_power,
        frequencies_hz=sc.frequencies_hz,
        time_offsets_seconds=sc.time_offsets_seconds,
    )
    validate_arrays(arrays, sc)
    np.savez_compressed(output / "waterfall.npz", **arrays)
    np.savez_compressed(
        output / "geometry.npz",
        norad_ids=g.norad_ids,
        alt_deg=g.alt_deg.astype("float32"),
        az_deg=g.az_deg.astype("float32"),
        separation_deg=g.separation_deg.astype("float32"),
        range_m=g.range_m.astype("float32"),
        range_rate_m_s=g.range_rate_m_s.astype("float32"),
        l1_doppler_hz=r.doppler_hz.astype("float32"),
        half_power_mask=r.half_power_mask,
        solar_alt_deg=r.solar["alt_deg"],
        solar_az_deg=r.solar["az_deg"],
        solar_separation_deg=r.solar["separation_deg"],
    )
    with (output / "passages.csv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle)
        w.writerow(
            [
                "norad_id",
                "name",
                "entry_utc",
                "exit_utc_exclusive",
                "duration_seconds",
                "min_separation_deg",
            ]
        )
        names = {s.norad_id: s.name for s in g.satellites}
        for v in sorted(r.passages, key=lambda p: p.start_index):
            row = g.norad_ids.index(v.norad_id)
            from datetime import timedelta

            w.writerow(
                [
                    v.norad_id,
                    names[v.norad_id],
                    format_utc(
                        sc.start_utc + timedelta(seconds=v.entry_offset_seconds)
                    ),
                    format_utc(sc.start_utc + timedelta(seconds=v.exit_offset_seconds)),
                    v.duration_seconds,
                    float(g.separation_deg[row, v.start_index : v.stop_index].min()),
                ]
            )
    sun_peak = int(np.argmax(r.solar_power.sum(axis=0)))
    checks = {
        "finite_positive_total": bool(
            np.all(np.isfinite(r.total_power)) and np.all(r.total_power > 0)
        ),
        "geometry_and_emission_half_power_masks_agree": bool(
            np.array_equal(g.selected, r.half_power_mask)
        ),
        "linear_component_sum": bool(
            np.array_equal(
                r.total_power, r.noise_power + r.satellite_power + r.solar_power
            )
        ),
        "power_float32": all(
            a.dtype == np.float32
            for a in [r.noise_power, r.satellite_power, r.solar_power, r.total_power]
        ),
        "solar_transit_inside_beam": bool(
            beam_gain(
                r.solar["offset_x_deg"][sun_peak],
                r.solar["offset_y_deg"][sun_peak],
                sc.section("beam")["reference_frequency_hz"],
                sc.section("beam"),
            )
            >= 0.5
        ),
        "satellite_power_nonnegative": bool(np.all(r.satellite_power >= 0)),
        "no_satellite_power_when_all_below_horizon": bool(
            not r.satellite_power[:, ~g.above_horizon.any(axis=0)].any()
        ),
    }
    if not all(checks.values()):
        raise RuntimeError(f"Consistency checks failed: {checks}")
    meta = base_metadata(sc)
    meta.update(
        {
            "models": {
                "orbit": "SGP4, WGS84, unrefracted topocentric geometry",
                "emission": "channel-integrated sinc² BPSK proxy with radial Doppler; L1/L2 mixture, L5 broad; no actual composite waveform",
                "beam": sc.section("beam"),
                "solar": sc.config["solar"],
                "not_modelled": [
                    "measured EIRP/solar flux",
                    "polarization",
                    "multipath",
                    "ionosphere",
                    "scattering",
                    "receiver bandpass",
                    "solar bursts",
                    "satellite health/status",
                ],
            },
            "reference_hashes_sha256": verify_reference_hashes(sc),
            "array_sha256": {name: array_sha256(arr) for name, arr in arrays.items()},
            "ephemeris_sha256": r.solar["ephemeris_sha256"],
            "counts": {
                "satellites": len(g.satellites),
                "passages": len(r.passages),
                "satellites_entering_half_power_beam": int(
                    r.half_power_mask.any(axis=1).sum()
                ),
                "maximum_simultaneous_half_power": int(
                    r.half_power_mask.sum(axis=0).max()
                ),
                "maximum_simultaneous_above_horizon": int(
                    g.above_horizon.sum(axis=0).max()
                ),
            },
            "solar_peak_utc": format_utc(sc.timestamps_utc[sun_peak]),
            "solar_min_separation_deg": float(r.solar["separation_deg"].min()),
            "l1_doppler_range_hz_visible": [
                float(r.doppler_hz[g.above_horizon].min()),
                float(r.doppler_hz[g.above_horizon].max()),
            ],
            "sampling": {
                "step_seconds": sc.step_seconds,
                "channel_width_hz": sc.channel_width_hz,
                "shape_frequency_time": list(sc.matrix_shape),
                "note": "Doppler is unresolved in the wideband image; see separate kHz plot. Powers are sample-time channel means, not a voltage/FFT simulation.",
            },
            "satellites": [
                {
                    "norad_id": s.norad_id,
                    "name": s.name,
                    "block": s.block,
                    "tle_epoch_utc": format_utc(s.tle_epoch_utc),
                    "max_epoch_offset_hours": max(
                        abs((sc.start_utc - s.tle_epoch_utc).total_seconds()),
                        abs((sc.end_utc - s.tle_epoch_utc).total_seconds()),
                    )
                    / 3600,
                }
                for s in g.satellites
            ],
            "checks": checks,
            "display": {
                "formula": "10 log10(P / P0)",
                "unit": "relative dB, not dBm",
                "common_limits_db": [float((10*np.log10(r.total_power/sc.config["noise"]["mean_power_relative"])).min()), float((10*np.log10(r.total_power/sc.config["noise"]["mean_power_relative"])).max())],
            },
            "passage_definition": "Elliptical half-power contour at reference frequency, above horizon; emission continues in Gaussian wings.",
        }
    )
    (output / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    (output / "scenario.json").write_text(
        json.dumps(sc.config, indent=2), encoding="utf-8"
    )
    plot_realistic_figures(r, output)
    print(
        json.dumps(
            {
                k: meta[k]
                for k in [
                    "counts",
                    "solar_peak_utc",
                    "l1_doppler_range_hz_visible",
                    "checks",
                ]
            },
            indent=2,
        )
    )
    return meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("outputs/realistic"))
    parser.add_argument("--ephemeris", type=Path)
    args = parser.parse_args()
    save_realistic(
        simulate_realistic(load_scenario(args.config), args.ephemeris), args.output
    )


if __name__ == "__main__":
    main()
