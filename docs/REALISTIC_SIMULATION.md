# GNSSVis 2.0 development — wideband GPS and solar transit

This optional extension adds realistic *structure* to received spectral power.
It reuses GNSSVis `SkyAnalyzer`, SGP4, local TLE verification, WGS84, scenario
loading and the seeded noise generator. It is a simulation for the radio
astronomy community, with no claimed institutional affiliation to BINGO.
The development label “2.0” does not mean a 2.0 release was published on PyPI;
the attached distribution metadata and the public PyPI release remain 0.1.0.

## Import and run

```python
from pathlib import Path
from gnssvis_simulation import load_scenario
from gnssvis_simulation.realistic import simulate_realistic, save_realistic

scenario = load_scenario('examples/scenarios/pibic_24h.json')
result = simulate_realistic(scenario)
save_realistic(result, Path('outputs/pibic_24h'))
```

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -e '.[dev]'
python -m gnssvis_simulation.realistic \
  --config examples/scenarios/pibic_24h.json --output outputs/pibic_24h
python -m pytest -q
```

Imports do not fetch ephemerides. The first run downloads DE421 to
`data/ephemeris/de421.bsp` if absent and `solar.allow_download` is true.
Subsequent runs use the cache. Use `--ephemeris /path/de421.bsp` to choose an
existing kernel; a missing explicitly specified file is an error, never an
implicit download. To forbid downloads set `solar.allow_download: false`.
The DE421 file itself is excluded from the small code archive. The original
fully offline CLI and scenario remain available.

## Frozen demonstration

- Time: 2026-08-21 06:00 UTC to 2026-08-22 06:00 UTC, end excluded.
  In Recife: 21 August 03:00 to 22 August 03:00 (UTC−3).
- 32 GPS satellites; TLE epochs stay within 63.41 h of the observation interval.
  A historical snapshot may contain epochs after an observation sample: this
  is a retrospective demonstration, not a real-time forecast.
- Reference coordinates: latitude −7.211637°, longitude −35.908138°, height
  553 m as WGS84 ellipsoidal height. These are inherited demonstration
  coordinates, not a surveyed instrument or a claimed BINGO configuration.
- Fixed pointing: azimuth 59.819450°, altitude 52.016700°.
- Solar center crosses the boresight at 2026-08-21 12:15 UTC (09:15 Recife).
- Elliptical FWHM: 50° × 40° at 1575.42 MHz. These illustrative widths represent
  a broad receiving beam; no measured antenna pattern or large-dish beam is
  implied. A broad beam intentionally increases the number of passages.
- 17 half-power passages, 16 distinct satellites entering half power,
  up to three simultaneous inside half power at the reference frequency;
  two are present when the Sun crosses the axis. Up to 15 GPS are above the
  horizon, which is a different count from GPS inside the beam.
- The boresight was selected from Sun positions between 12:00 and 16:00 UTC
  at 5 min spacing. Objective: N_passages + 4×N_GPS_at_solar_center +
  2×maximum_simultaneous. The full candidate table is in the scenario JSON.
  The selection favours presentation visibility, not an unbiased survey.
- 1150–1600 MHz, 768 frequency channels, 1440 timestamps at 60 s spacing.
  Channel width 585937.5 Hz. One float32 power plane occupies 4.42 MB;
  all four power planes occupy 17.69 MB before compression. NPZ is compressed.
  At 10 s the same four planes would occupy six times as much memory.

## Physics and numerical choices

The orbital propagator remains **SGP4**: changing the apparent spectral shape
requires an emission/receiving model, not an arbitrary replacement orbit.
`SkyAnalyzer.satellite_altaz` now also returns `range_rate_m_s` from the
observer-relative state vector:

```text
r_dot = (position · velocity) / |position|
Delta_f = - f_c r_dot / c
```

Positive radial velocity means recession and negative Doppler. In the frozen
run the L1 shift for visible GPS ranges from approximately −3.48 to +3.77 kHz.
This is much smaller than a wideband channel. Its physical trajectory is
plotted separately; the waterfall never exaggerates the shift for appearance.

The received signal proxy is:

```text
P_sat(f,t) = sum_s,c A_c P_ref (r_ref/r_s(t))²
                    × G_s(f,t) × channel_mean[PSD_c(f-f_c-Delta_f_s(t))]
```

GPS L1 C/A-like narrow spectra use Rc = 1.023 Mchip/s; broad components use
10.23 Mchip/s. L1/L2 use an illustrative unit-peak mixture with broad weight
0.32 (this is a peak-shape weight, not an integrated-power fraction).
L5 is broad and injected only for IIF and III blocks, never IIR/IIR-M.
Carrier weights L1=1, L2=0.7, L5=0.85 are illustrative, not ICD power ratios.
The profiles are channel averages of sinc², integrated analytically via the
sine integral. Rebinning conserves the spectral power represented by each
channel; no artificial chirp is added. L1C/BOC details, code/data waveforms,
health flags and exact composite satellite spectra are outside this proxy.

Beam offsets use the spherical logarithmic map about boresight, avoiding
azimuth wrapping. Rotated offsets x,y define:

```text
FWHM_x(f) = 50° × f_ref/f; FWHM_y(f) = 40° × f_ref/f
G(f,t) = exp[-4 ln(2) ((x/FWHM_x(f))²+(y/FWHM_y(f))²)]
```

The half-power contour at the reference frequency defines *passage tables*.
Emission continues smoothly outside it wherever the satellite is above the
horizon; the half-power boundary does not truncate the signal. This Gaussian
main-lobe model has no measured sidelobes, ground reflection or aperture
calibration. L5 has a wider beam than L1 in the same model.

The Sun is computed with apparent topocentric DE421 ephemerides and no
atmospheric refraction. The solar disk is unresolved relative to this very
broad beam. Default solar power is an illustrative quiet thermal continuum,
flat in relative received power before beam weighting. It is not a solar
flux measurement, antenna-temperature prediction or a pure synchrotron Sun.
For an explicitly hypothetical active nonthermal contribution, set
`solar.active_power_relative > 0`; the optically thin synchrotron spectral
index is −(p−1)/2, with p=`electron_energy_index`. Solar bursts, gyrosynchrotron
turnover, opacity and real activity on the chosen day are not inferred.

Noise uses Gamma(M, P0/M), with M=16 and seed=42, as in the original pipeline.
Power components are summed linearly, then shown as 10 log10(P/P0).
No units of watts, dBm or measured SFU are claimed. The added solar continuum
represents a mean-power component; source self-noise and temporal correlation
are not modelled. Time samples are instantaneous power expectations plus
synthetic noise, not a time-integrated voltage/FFT receiver calculation.

## Outputs and validation

- `waterfall.npz`: separate noise, satellite, solar and total float32 powers,
  shape (frequency,time), frequency centers in Hz, integer time offsets in s.
  Read using `numpy.load(..., allow_pickle=False)`.
- `geometry.npz`: compact per-satellite geometry, radial velocities, L1 Doppler,
  half-power masks, Sun alt/az and separation.
- `passages.csv`: half-power intervals with exclusive exits at sampled times.
- `metadata.json`: effective parameters, dependency versions, input hashes,
  counts, spectral limits, assumptions and consistency checks.
- `scenario.json`: exact configuration used.
- Five PNGs: overview, components, L1/Doppler detail, spectral/beam models and
  a common-grid comparison against the old fixed Gaussian emission proxy.

Plots use common limits −2 to 20 relative dB for readable comparison.
The default total ranges approximately −6.44 to +16.93 dB: some low-power
background pixels fall below the colour scale, while no high-power peaks are
clipped. Generated arrays themselves are never clipped or smoothed.

Tests cover Doppler sign, radial velocity against a finite range derivative,
beam FWHM and chromatic scaling, angular wrapping, channel rebinning,
integrated BPSK power and nulls, deterministic component synthesis, lack of L5
on IIR, below-horizon suppression and explicit missing-ephemeris behaviour.
The full original offline suite is preserved.

## Sources

- Historical CelesTrak GNSS snapshot in upstream commit
  `4b3f5b75e28dc5346e74544e47e9b1ecaf5c5c84`; exact URL and hashes in
  `reference_data/pibic_24h/provenance.json`.
- GPS frequency plans: IS-GPS-200N §3.3.1.1 (L1/L2), IS-GPS-705J §3.3.1.1 (L5),
  links preserved in the provenance file. PSDs remain documented proxies.
- JPL DE421: https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de421.bsp
- Skyfield satellite and radial-state geometry: https://rhodesmill.org/skyfield/earth-satellites.html
- Quiet/active radio Sun background: Dulk (1985), *Radio emission from the Sun and stars*, https://doi.org/10.1146/annurev.aa.23.090185.001235
- PyPI metadata (independent of this local development): https://pypi.org/pypi/gnssvis/json
