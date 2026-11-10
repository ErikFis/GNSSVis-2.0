# gnssvis

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

- [Shared offline observation scenario](docs/OBSERVATION_SCENARIO.md)
- [Observation configuration](examples/scenarios/offline.json)
- [Frozen TLE and frequency references](reference_data/README.md)
- [Context for the first Taiga prompt](docs/FIRST_PROMPT.md)
- [Architecture](docs/ARCHITECTURE.md)
- [Plotting guide](docs/PLOTTING_GUIDE.md)
- [Public API](docs/SIMPLE_API.md)
- [Examples](examples/README.md)
- [Notebooks](notebooks/README.md)
- [Task preparation](TASK_PREPARATION.md)

## Contributing

Create a feature branch, make focused changes, and validate them before
opening a pull request. Keep documentation, comments, and messages in English.
Keep generated artifacts and downloaded datasets outside version control.

## License

MIT; see [LICENSE](LICENSE). Original upstream:
https://github.com/ErikFis/gnssvis
