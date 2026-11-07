# Package source

The `gnssvis` namespace exposes `GNSSVis`, `DataAnalyzer`, `IOManager`,
`SkyAnalyzer`, and `Visualizer` through `gnssvis/__init__.py`.

| File | Purpose |
| --- | --- |
| `core.py` | Unified workflow and numerical waterfall operations |
| `io.py` | FITS loading, pointing, synthetic data, and NPZ consolidation |
| `sky.py` | Ephemerides, GNSS trajectories, visibility, and catalog caches |
| `plot.py` | Waterfall and trajectory visualization |
| `tle_utils.py` | TLE age checks, archiving, and cleanup |

Use imports such as `from gnssvis import GNSSVis` after installing the source
with `python -m pip install -e .`. The `src` directory is a packaging layout,
not part of the public import path.

See the [API guide](../../docs/SIMPLE_API.md),
[architecture](../../docs/ARCHITECTURE.md), and
[plotting guide](../../docs/PLOTTING_GUIDE.md). Keep generated data, bytecode,
and installation metadata out of version control.
