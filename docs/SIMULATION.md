# Offline synthetic GNSS waterfall

`gnssvis_simulation` generates a reproducible synthetic spectro-temporal
waterfall: an averaged white-noise background plus simplified GNSS carrier
emission during satellite beam passages computed with SGP4 from local TLEs.
Everything is driven by one scenario file (`examples/scenarios/offline.json`);
no network access, DE421 ephemeris or observation data is needed.

## Reproduce the example

From the repository root (Python >= 3.9):

```bash
python -m venv .venv                       # optional; .venv/ is git-ignored
.venv/bin/pip install -e ".[dev]"          # numpy, scipy, skyfield, astropy, matplotlib, pytest
.venv/bin/python -m gnssvis_simulation --config examples/scenarios/offline.json
```

The command is headless (Agg backend, no `show()`), prints a summary and exits
with status 0 when all internal consistency checks pass (1 if a check fails,
2 for an invalid scenario). Use `--output DIR` to write somewhere other than
`output.default_directory` (`outputs/offline`, relative to the repository root,
whatever the working directory).

Check everything:

```bash
.venv/bin/python -m pytest -v                # full test suite (about 30 s)
.venv/bin/python -m pytest -s -v             # ... printing measured statistics
.venv/bin/python -m gnssvis_simulation.partial --config examples/scenarios/offline.json
```

The last command writes step-by-step diagnostics to `partial_results/`
(geometry, noise statistics, signal-only preview) with a terminal report.

With the same seed, scenario and dependency versions the arrays are
bit-identical; `metadata.json` records the SHA-256 of each array and the
versions. The seed-42 `noise_power` hash is
`60f04c46fbc48b35309b379c31f83a1c2cd88bbcdda5f71c0d3bc7af1f39a766`
(NumPy 2.3.5).

## Outputs (`outputs/offline/`, git-ignored)

| File | Content |
| --- | --- |
| `waterfall.png` | Upper panel: angular separation of each satellite from the pointing, beam gate line, in-beam samples thick. Lower panel: `10*log10(total_power/P0)`, UTC on X, MHz on Y, dotted labelled in-band carrier lines, horizontal colorbar with the full finite dB range. Right: satellite/colour table with passage times and in-band carriers. Dashed verticals mark sampled entry and exclusive exit in both panels. |
| `waterfall.npz` | `noise_power`, `signal_power`, `total_power`: float32 `(frequency, time)` relative linear power, `total = noise + signal`. `frequencies_hz`: float64 channel centres. `time_offsets_seconds`: int64 seconds since `time_grid.start_utc`. Load with `numpy.load(path, allow_pickle=False)`. |
| `samples.csv` | One row per time sample and satellite (time, then configured NORAD order): altitude, azimuth, separation, range in metres, `above_horizon`/`in_beam`/`in_passage` (lowercase `true`/`false`), semicolon-separated in-band `signal_ids`. |
| `passages.csv` | One row per passage: entry, exclusive exit, offsets, duration, sample count, clipping flags, minimum separation, peak signal power, signal IDs. The header is written even when there are no passages. |
| `metadata.json` | Effective configuration and its hash, reference-file hashes, TLE epochs and epoch offsets, UT1-UTC used, seed, dependency versions, units, array shapes and hashes, row/passage counts, model description, display colour limits and the consistency checks. |

Reference result: only GPS BIIF-8 (PRN 03, NORAD 40294) enters the 15 degree
zenith gate, from 05:01:50 to 05:39:20 UTC (exclusive; 225 samples). Only
GPS L1 lies in the 1500-1600 MHz band. The peak signal is 13.24 times P0
(11.2 dB) at 05:20:50 UTC in channel 386 (1575.488 MHz).

## Model

Grid: `frequencies_hz[i] = band_min + (i + 0.5) * df` (band limits are channel
edges), `time_offsets_seconds[j] = j * step`; the end of the observation is
excluded. All sizes come from the scenario and can be changed there.

Geometry (`geometry.py`, reusing `gnssvis.SkyAnalyzer(offline=True)`): SGP4
propagation of the local TLEs, WGS84 observer, unrefracted topocentric
altitude/azimuth/range, Skyfield's bundled UT1-UTC/leap-second table. The
passage mask is `separation <= full_width/2` and `altitude >= minimum_altitude`
at each sample; no interpolation between samples.

Noise (`noise.py`): each cell is `Gamma(shape=M, scale=P0/M)`, exactly the
mean of M independent complex-Gaussian power samples with mean P0. Drawn with
`numpy.random.default_rng(seed)`; the only random component.

Emission (`emission.py`), for each satellite and each catalog carrier that
applies to its constellation/block and lies inside the receiver band:

```
P(f, t) = mask(t) * P_ref * (r_ref / r(t))**2
          * exp(-4 ln2 (sep(t) / full_width)**2)
          * exp(-0.5 ((f - f_c) / sigma)**2)   for |f - f_c| <= cutoff * sigma, else 0
```

evaluated at the actual channel centres; coefficients per satellite block from
`emission.profiles`. Terms are summed in float64 and stored as float32.

Combination (`pipeline.py`): `total_power = noise_power + signal_power` in
linear power; dB only for display.

## Assumptions and limitations

- Relative, illustrative power. P_ref is the same for IIF and III; no EIRP,
  dBm, antenna gain pattern, system temperature or absolute calibration.
- No Doppler (carrier fixed at nominal frequency), no actual GPS modulation or
  occupied bandwidth, no polarization, multipath, receiver electronics or
  antenna thermal noise. No statement about satellite health or transmission.
- The beam is a hard circular gate with a Gaussian taper inside. The gain is
  0.5 at the gate, so the signal starts and stops abruptly (about 9.6 P0 at
  the peak channel here). This is intentional, not a real beam pattern.
- Sample-level timing: passage entry/exit are known only to the time step
  (10 s). The reference entry sample (05:01:50) is only about 0.2 arcsec
  inside the gate; Skyfield 1.55's bundled UT1-UTC prediction (+0.090 s)
  differs from later IERS values (about +0.007 s) by about 1.3 arcsec of
  Earth rotation, which would move that entry one sample later. Both are
  far below TLE accuracy (about 10 arcsec per km at GPS range).
- TLE accuracy degrades away from the epoch; the loader enforces the scenario
  limit (72 h) over the whole window.
- Only fixed alt/az pointing, WGS84 and the SGP4 propagator are supported;
  other values in the scenario are rejected rather than ignored.
- Noise reproducibility across NumPy versions is not guaranteed by NumPy;
  the statistical tests still validate a changed stream, and the metadata
  records versions and hashes.

## Extension points

- Background: any object with `generate(shape, rng) -> float32 array` satisfies
  `noise.NoiseModel` (for example frequency-dependent system temperature or
  RFI-free bandpass ripple) and can be passed to `generate_noise(model=...)`.
- Spectral profile: `GaussianCarrierProfile` provides `range_scaling`,
  `support` and `spectral_shape`; another profile class with the same methods
  (e.g. a BPSK/BOC sinc-squared shape) can be used per block.
- Catalogs: new constellations or signals are rows in `signals.csv` (and
  `satellites.csv`/TLEs) plus a profile for each new block.
- Time-varying effects such as Doppler need the carrier to vary with time;
  `Contribution` currently stores a separable `outer(spectral, envelope)` term,
  so that would be the place to add a per-sample spectral evaluation.
- Transient sources (pulsars, FRBs) are not implemented; they would be
  additional linear-power contributions added before `total_power`.
