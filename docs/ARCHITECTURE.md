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

The application lives in the separate `src/gnssvis_simulation` namespace and
explicitly reuses GNSSVis calculations. Implemented so far:

| Module | Responsibility |
| --- | --- |
| `config.py` | Load/validate the scenario JSON, derive frequency and time grids, resolve reference paths from the repository root, verify provenance hashes |
| `geometry.py` | Offline `SkyAnalyzer` (local TLEs, SGP4, WGS84, bundled timescale), TLE epoch checks, alt/az/range/separation, beam masks, passages |
| `noise.py` | `NoiseModel` interface (`generate(shape, rng) -> float32`) and `AveragedWhiteNoise`, Gamma(M, P0/M) with `default_rng(seed)` |
| `emission.py` | Signal catalog and band selection, per-block `GaussianCarrierProfile`, Gaussian beam gain, masked per-satellite/carrier `Contribution`s summed in linear power |
| `products.py` | Array contract checks, NPZ save/load (no pickle), JSON metadata, dependency versions, array hashes |
| `diagnostics.py` | Noise statistics vs theory (z-scores, KS test), contribution intervals/peaks, text tables, ASCII passage timeline |
| `figures.py` | Headless (Agg) geometry, noise waterfall/statistics and signal-only preview figures |
| `partial.py` | `python -m gnssvis_simulation.partial`: terminal report and `partial_results/` |
| `pipeline.py` | `total_power = noise_power + signal_power`, CSV tables, metadata, final figure, consistency checks |
| `__main__.py` | `python -m gnssvis_simulation --config ... [--output ...]` (headless) |

`figures.plot_final_waterfall` reuses `Visualizer.draw_object_table` and the
GNSSVis colour list; the original `Visualizer.plot_waterfall_with_objects`
output is unchanged. See [SIMULATION.md](SIMULATION.md) for the model and frozen inputs. Keep dependencies directed from the application to the library, without a circular dependency.

## Optional wideband and solar pipeline

`beam.py` supplies spherical offsets and the chromatic elliptical Gaussian
response. `realistic.py` reuses scenario/geometry/noise modules, adds radial
Doppler and channel-integrated spread-spectrum proxies, loads cached or
online DE421 solar ephemerides and returns separate power components.
`realistic_figures.py` renders five common-scale diagnostic/presentation figures.
The public API lazily exposes `simulate_realistic`, `save_realistic` and
`RealisticResult` without an ephemeris download at import. See
[model, assumptions and reproduction](REALISTIC_SIMULATION.md).
