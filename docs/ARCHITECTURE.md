# Modular architecture

GNSSVis separates FITS input, astronomical calculations, and visualization.
`GNSSVis` in `core.py` coordinates these components through one public interface.

## Modules

| Module | Class | Responsibility |
| --- | --- | --- |
| `core.py` | `GNSSVis` | Observation context and analysis workflow |
| `core.py` | `DataAnalyzer` | Normalization, statistics, baseline subtraction, smoothing |
| `io.py` | `IOManager` | FITS input, pointing metadata, synthetic FITS, NPZ consolidation |
| `sky.py` | `SkyAnalyzer` | Solar system and GNSS trajectories, visibility, catalog caches |
| `plot.py` | `Visualizer` | Synchronized trajectories, waterfall, and object table |
| `tle_utils.py` | Utility functions | TLE archiving, age checks, and archive cleanup |

`SkyAnalyzer` stores downloaded astronomical data under
`data/tle_repository/`. The plotting layer receives data and coordinates;
it does not instantiate an astronomical analyzer.

## Typical workflow

```python
from gnssvis import IOManager, SkyAnalyzer, Visualizer

io = IOManager()
observation = io.load_fits_observation("observation.fits")
pointing = io.get_observer_pointing_from_fits("observation.fits")
sky = SkyAnalyzer()
objects = sky.find_objects_in_fov(
    observation["times"], pointing, beamwidth_deg=15.0,
    include_gnss=True, include_pulsars=False,
)
stats = sky.get_objects_statistics(
    objects, pointing, 15.0, observation["times"]
)
Visualizer().plot_waterfall_with_objects(
    data=observation["data"],
    times=observation["times"],
    frequencies=observation["frequencies"],
    sky_objects=objects,
    telescope_pointing=pointing,
    beamwidth_deg=15.0,
)
```

The observation dictionary contains `data` (frequency x time), `times`
(`astropy.time.Time`), `frequencies` (MHz), and the primary FITS `header`.
Visibility uses spherical angular separation and includes an object if it
enters the beam at any sampled time.

## Extending the package

Add file-format support to `IOManager`, new object calculations to
`SkyAnalyzer`, and plotting methods to `Visualizer`. Keep numerical
transformations in `DataAnalyzer` when they do not need the observation context.
Preserve the public method signatures when making preparation-only changes.

This division allows focused tests, reusable components, and independent
maintenance of scientific calculations and presentation.

## Independent simulation pipeline

The sibling `gnssvis-simulation` repository is a separate Python distribution
with its own Git history and `gnssvis_simulation` namespace. It does not import
or modify the existing modules. Its pipeline implementation is intentionally
left for the programming task. Any future integration should be explicit and
optional rather than introducing a circular dependency.
