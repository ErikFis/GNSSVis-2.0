# Plotting guide

## Unified interface

```python
from gnssvis import GNSSVis

gv = GNSSVis()
gv.plot_waterfall_with_sky("observation.fits", beamwidth_deg=20.0)
```

`beamwidth_deg` is the full beam width in degrees. The method loads the
observation and pointing, finds solar system and GNSS objects entering the
field of view, and displays a Matplotlib figure. The unified workflow disables
pulsars by default. Catalog downloads are currently disabled in `SkyAnalyzer`;
a persistent pulsar cache can still be read through its direct API.

## Figure layout

The upper left panel shows relative declination against time. The lower left
panel shows the waterfall (frequency in MHz against time, intensity in dB).
The right panel lists object names with their colors. Time limits are shared.
The waterfall uses the `plasma` colormap and the second and 98th percentiles
for display limits. Matching GNSS emission frequencies are drawn as horizontal
lines when they fall within the observed band.

These are declination trajectories, not altitude/elevation curves. Methods
called `analyze_observation`, `plot_waterfall`, or `plot_objects_in_fov` are not
part of the current implementation; use the interfaces shown here.

## Processing the loaded waterfall

```python
detrended = gv.baseline_subtract(axis=1, plot=True)
smoothed = gv.gaussian_smooth(
    data=detrended, sigma_freq=2.0, sigma_time=2.0, plot=True
)
```

For an array shaped `(n_freq, n_time)`, axis 1 subtracts each frequency
channel's median over time. Axis 0 subtracts each time sample's median over
frequency. Gaussian standard deviations are measured in array samples.
Pass `plot=False` to return the processed array without plotting.

## Direct visualization without sky downloads

```python
import numpy as np
import astropy.units as u
from astropy.time import Time
from gnssvis import Visualizer

rng = np.random.default_rng(42)
times = Time("2025-11-12T20:00:00", scale="utc") + np.arange(100) * u.s
Visualizer().plot_waterfall_with_objects(
    data=rng.normal(-87.5, 2.0, (32, 100)),
    times=times,
    frequencies=np.linspace(1200, 1300, 32),
    sky_objects={},
)
```

An empty object dictionary displays the waterfall with a no-objects message.
Use an interactive Matplotlib backend to view plots, or `MPLBACKEND=Agg` for
headless validation. Close figures when processing many observations.

## Coordinates and performance

Telescope pointing is an Astropy `SkyCoord`; time arrays use Astropy `Time`.
GNSS positions use Skyfield SGP4 and topocentric vectors. Planetary positions
use DE421. Long observations and dense time grids increase computation cost.
Keep TLEs current when assessing satellite identification accuracy.
