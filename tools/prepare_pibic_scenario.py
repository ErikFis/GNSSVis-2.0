"""Prepare an immutable 32-GPS historical snapshot and choose a solar crossing."""

from pathlib import Path
from skyfield.api import EarthSatellite, load, load_file, wgs84
from gnssvis_simulation import load_scenario, compute_geometry
from gnssvis_simulation.realistic import beam_offsets_deg, beam_gain
from gnssvis_simulation.config import sha256_file, format_utc
import json, csv, re, numpy as np
from urllib.request import urlopen

root = Path(__file__).resolve().parents[1]
cache = root / "data/ephemeris"
cache.mkdir(parents=True, exist_ok=True)
archive = root / "data/gps_snapshot_full.tle"
source_url = "https://raw.githubusercontent.com/ErikFis/gnssvis/4b3f5b75e28dc5346e74544e47e9b1ecaf5c5c84/data/tle_repository/gnss/gnss.txt"
if not archive.is_file():
    with urlopen(source_url, timeout=60) as response:
        archive.write_bytes(response.read())
if not (cache / "de421.bsp").is_file():
    from skyfield.api import Loader

    Loader(str(cache))("de421.bsp").close()
ref = root / "reference_data/pibic_24h"
ref.mkdir(exist_ok=True)
lines = archive.read_text().splitlines()
ts = load.timescale(builtin=True)
sats = []
selected = []
for i in range(0, len(lines), 3):
    if lines[i].startswith("GPS "):
        selected += [lines[i].strip(), lines[i + 1], lines[i + 2]]
        sats.append(EarthSatellite(lines[i + 1], lines[i + 2], lines[i].strip(), ts))
(ref / "satellites.tle").write_text("\n".join(selected) + "\n")
with (ref / "satellites.csv").open("w", newline="") as h:
    w = csv.writer(h)
    w.writerow(["norad_id", "name", "constellation", "block", "prn", "tle_epoch_utc"])
    for sat in sats:
        b = "III" if "BIII-" in sat.name else "IIF" if "BIIF-" in sat.name else "IIR"
        w.writerow(
            [
                sat.model.satnum,
                sat.name,
                "GPS",
                b,
                int(re.search(r"PRN\s+(\d+)", sat.name)[1]),
                format_utc(sat.epoch.utc_datetime()),
            ]
        )
(ref / "signals.csv").write_text(
    (root / "reference_data/signals.csv")
    .read_text()
    .replace("IIF|III", "IIR|IIF|III", 2)
)
prov = {
    "tle": {
        "original_provider": "CelesTrak, via historical upstream snapshot",
        "source_url": "https://raw.githubusercontent.com/ErikFis/gnssvis/4b3f5b75e28dc5346e74544e47e9b1ecaf5c5c84/data/tle_repository/gnss/gnss.txt",
        "source_commit": "4b3f5b75e28dc5346e74544e47e9b1ecaf5c5c84",
        "source_sha256": sha256_file(archive),
        "selection": "All 32 records whose name starts with GPS; health/status unknown; retrospective TLE epochs can follow observation samples.",
        "snapshot_retrieved_utc": __import__("datetime")
        .datetime.now(__import__("datetime").timezone.utc)
        .isoformat(),
        "satellite_count": len(sats),
    },
    "frequency_sources": json.loads(
        (root / "reference_data/provenance.json").read_text()
    )["frequency_sources"],
    "files_sha256": {
        f: sha256_file(ref / f)
        for f in ["satellites.tle", "satellites.csv", "signals.csv"]
    },
}
(ref / "provenance.json").write_text(json.dumps(prov, indent=2))
c = json.loads((root / "examples/scenarios/offline.json").read_text())
c["scenario_id"] = "pibic_gps_sun_24h_20260821"
c["purpose"] = (
    "Illustrative 24-hour simulation for the radio astronomy community; no BINGO affiliation or measured observation."
)
c["time_grid"].update(
    start_utc="2026-08-21T06:00:00Z", duration_seconds=86400, step_seconds=60
)
c["observatory"][
    "name"
] = "Configurable reference site in Paraíba; illustrative coordinates"
c["beam"].update(
    model="chromatic_elliptical_gaussian",
    full_width_deg=50.0,
    minor_width_deg=40.0,
    reference_frequency_hz=1575.42e6,
    position_angle_deg=0,
    selection_note="Elliptical half-power contour at reference frequency; full Gaussian tails synthesized above horizon.",
    inside_rule="gain_at_reference_frequency >= 0.5",
)
c["receiver"].update(
    band_min_hz=1150e6,
    band_max_hz=1600e6,
    frequency_model="channel-integrated BPSK proxy",
    doppler="radial velocity from SGP4 topocentric state",
    occupied_bandwidth="BPSK sinc² narrow/broad mixture; not an exact composite ICD waveform",
)
c["references"].update(
    tle_file="reference_data/pibic_24h/satellites.tle",
    satellite_catalog_file="reference_data/pibic_24h/satellites.csv",
    signal_catalog_file="reference_data/pibic_24h/signals.csv",
    provenance_file="reference_data/pibic_24h/provenance.json",
    norad_ids=[s.model.satnum for s in sats],
)
c["geometry"].update(planetary_ephemeris_required=True, network_access_allowed=True)
c["frequency_grid"]["channel_count"] = 768
c["emission"]["profiles"]["IIR"] = c["emission"]["profiles"]["IIF"].copy()
c["realistic_emission"] = {
    "model": "channel_averaged_bpsk_proxy",
    "peak_power_relative": 24.0,
    "reference_range_m": 20e6,
    "narrow_chip_rate_hz": 1.023e6,
    "broad_chip_rate_hz": 10.23e6,
    "broad_fraction": 0.32,
    "note": "Illustrative equal amplitudes; IIR/IIR-M have no L5 in this model.",
}
c["solar"] = {
    "ephemeris_file": "data/ephemeris/de421.bsp",
    "allow_download": True,
    "quiet_power_relative": 5.0,
    "quiet_received_spectral_index": 0.0,
    "active_power_relative": 0.0,
    "electron_energy_index": 3.0,
    "note": "Default quiet thermal Sun, flat relative received continuum proxy. Optional active continuum uses optically thin synchrotron slope -(p-1)/2; no real solar event or measured flux is claimed. Solar disk treated as unresolved versus 40-50 degree beam.",
}
c["output"]["default_directory"] = "outputs/pibic_24h"
c["passages"][
    "definition"
] = "Runs inside the elliptical half-power contour at reference frequency and above horizon. Emission includes Gaussian wings."
path = root / "examples/scenarios/pibic_24h.json"
path.write_text(json.dumps(c, indent=2))
sc = load_scenario(path)
g = compute_geometry(sc)
eph = load_file(str(root / c["solar"]["ephemeris_file"]))
site = c["observatory"]
obs = eph["earth"] + wgs84.latlon(
    site["latitude_deg"], site["longitude_deg"], elevation_m=site["elevation_m"]
)
times = ts.from_datetimes(sc.timestamps_utc)
alt, az, _ = obs.at(times).observe(eph["sun"]).apparent().altaz()
rows = []
# Choose a reproducible, high-Sun pointing with multiple satellites crossing;
# criterion explicitly balances total passages and simultaneous solar/GPS use.
for j in range(360, 601, 5):  # 12:00..16:00 UTC, in five-minute steps
    if alt.degrees[j] < 40:
        continue
    x, y = beam_offsets_deg(g.alt_deg, g.az_deg, alt.degrees[j], az.degrees[j])
    mask = (beam_gain(x, y, 1575.42e6, c["beam"]) >= 0.5) & g.above_horizon
    counts = mask.sum(axis=0)
    n_pass = int(
        np.count_nonzero(
            np.diff(np.pad(mask.astype(int), ((0, 0), (1, 1))), axis=1) == 1
        )
    )
    score = n_pass + 4 * int(counts[j]) + 2 * int(counts.max())
    rows.append(
        {
            "sample": j,
            "utc": format_utc(sc.timestamps_utc[j]),
            "altitude_deg": float(alt.degrees[j]),
            "azimuth_deg": float(az.degrees[j]),
            "passages": n_pass,
            "simultaneous_at_sun": int(counts[j]),
            "maximum_simultaneous": int(counts.max()),
            "score": score,
        }
    )
best = max(rows, key=lambda v: (v["score"], -abs(v["sample"] - 506)))
c["pointing"].update(
    altitude_deg=best["altitude_deg"],
    azimuth_deg=best["azimuth_deg"],
    azimuth_convention="degrees east of north",
)
c["scenario_selection"] = {
    "criterion": "Maximize N_passages + 4*N_GPS_at_solar_center + 2*max_simultaneous across 12-16 UTC in 5-min steps; ties nearest 14:26 UTC.",
    "selected": best,
    "candidates": rows,
}
path.write_text(json.dumps(c, indent=2))
eph.close()
print(json.dumps(best, indent=2))
