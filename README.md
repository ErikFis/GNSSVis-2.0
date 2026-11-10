# gnssvis

**Atualização GNSSVis 2.0:** [ver gráficos e baixar código/dados](docs/pibic-24h/README.md).

A Python library for analyzing and visualizing radio telescope observations,
with an emphasis on GNSS satellites and celestial objects in the telescope
field of view. Developed for radio astronomy research at the Federal University
of Campina Grande (UFCG), Brazil.

## Installation

From this source checkout:

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e .
```

Dependencies are declared in `pyproject.toml`: NumPy, Matplotlib, Astropy,
SciPy, Skyfield, and Requests. Virtual environments and generated requirement
files are excluded from Git. The published package can also be installed with
`python -m pip install gnssvis`; that installs the published release rather
than the prepared source in this checkout.

## Quick start

```python
from gnssvis import GNSSVis

gv = GNSSVis()
gv.plot_waterfall_with_sky("observation.fits", beamwidth_deg=30.0)
detrended = gv.baseline_subtract(axis=1)
smoothed = gv.gaussian_smooth(sigma_freq=2.0, sigma_time=2.0)
```

The library reads the FITS waterfall, computes planetary and GNSS trajectories,
and plots the waterfall with relative declination trajectories and an object
table. Initial use may download DE421 ephemerides and Celestrak TLEs; these
files are deliberately excluded from this source package.

### Offline GNSS geometry

`SkyAnalyzer(..., offline=True, tle_file=...)` propagates local TLEs with SGP4
using Skyfield's bundled timescale data. It loads no DE421 ephemeris, makes no
network requests and writes no caches; solar system objects are then skipped
with a warning. The default constructor is unchanged.

```python
from gnssvis import SkyAnalyzer

sky = SkyAnalyzer(-7.211637, -35.908138, 553.0, offline=True,
                  tle_file="reference_data/satellites.tle")
sat = sky.load_local_tles("reference_data/satellites.tle", norad_ids=[40294])[0]
geometry = sky.satellite_altaz(sat, times)  # alt_deg, az_deg, range_m (unrefracted)
```

## Offline synthetic waterfall simulation

`gnssvis_simulation` builds a reproducible synthetic waterfall (white noise +
GNSS carriers during computed beam passages) from
`examples/scenarios/offline.json`, fully offline:

```bash
pip install -e ".[dev]"
python -m gnssvis_simulation --config examples/scenarios/offline.json
python -m pytest
```

Outputs go to `outputs/offline/` (git-ignored): `waterfall.png`,
`waterfall.npz`, `samples.csv`, `passages.csv`, `metadata.json`. See
[docs/SIMULATION.md](docs/SIMULATION.md) for reproduction, the model, its
assumptions and limitations, and extension points.

### Intermediate results

The geometry, white-noise and satellite emission stages can also be inspected
separately, without notebooks:

```bash
python -m gnssvis_simulation.partial --config examples/scenarios/offline.json
python -m pytest -s -v        # tests, printing the measured statistics
```

The command prints the grids, satellite geometry, passages, noise statistics and
reproducibility checks, and writes to `partial_results/` (ignored by git):
`01_geometry_separation.png`, `02_noise_waterfall.png`,
`03_noise_statistics.png`, `04_signal_waterfall.png` (signal only, linear
power), `noise.npz` and `signal.npz` (`noise_power`/`signal_power` float32
(frequency, time), `frequencies_hz` float64 Hz, `time_offsets_seconds` int64 s
since start), and `noise_metadata.json`/`signal_metadata.json`. Use `--seed N` to override the seed and `--output DIR`
to change the folder.

## Synthetic demonstration

```python
from gnssvis import GNSSVis

results = GNSSVis().quick_demo()
print(f"Objects found: {len(results['objects'])}")
```

Real observation data is unnecessary for the demonstration. Sky data downloads
may still be needed. The temporary synthetic FITS file is removed afterward.

## Custom observatory

```python
gv = GNSSVis(
    observatory_lat=-22.9,
    observatory_lon=-43.2,
    observatory_alt=10.0,
)
gv.plot_waterfall_with_sky("observation.fits")
```

Latitude and longitude are degrees; altitude is meters above sea level.

## Documentation

- [Offline simulation: reproduction, model and limits](docs/SIMULATION.md)
- [Observation configuration](examples/scenarios/offline.json)
- [Frozen TLE and frequency references](reference_data/README.md)
- [Architecture](docs/ARCHITECTURE.md)
- [Plotting guide](docs/PLOTTING_GUIDE.md)
- [Public API](docs/SIMPLE_API.md)
- [Examples](examples/README.md)
- [Notebooks](notebooks/README.md)

## Contributing

Create a feature branch, make focused changes, and validate them before
opening a pull request. Keep documentation, comments, and messages in English.
Keep generated artifacts and downloaded datasets outside version control.

## License

MIT; see [LICENSE](LICENSE). Original upstream:
https://github.com/ErikFis/gnssvis

## Wideband GPS + Sun extension (2.0 development)

The optional online-capable simulation imports `simulate_realistic` and
`save_realistic` from `gnssvis_simulation.realistic`. It adds topocentric radial
Doppler, channel-integrated spread-spectrum profiles, a smooth chromatic
elliptical beam and apparent DE421 solar ephemerides. Imports do not download
data. See [model, assumptions and reproduction](docs/REALISTIC_SIMULATION.md).

```bash
python -m gnssvis_simulation.realistic --config examples/scenarios/pibic_24h.json --output outputs/pibic_24h
```

The frozen example has 32 GPS, 24 hours, 17 half-power passages and a solar
boresight transit at 2026-08-21 12:15 UTC. All powers are illustrative.
