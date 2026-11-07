# Examples

Install gnssvis from the repository root:

```bash
python -m pip install -e .
python examples/basic_usage.py
```

| Script | Purpose |
| --- | --- |
| `basic_usage.py` | Synthetic observation demo with UIRAPURU coordinates |
| `simple_usage_example.py` | Demos, custom observatory, processing, direct modules |
| `consolidate_observation.py` | Convert local NPZ observations into FITS |
| `concatenate_func.py` | Standalone NPZ consolidation helper |

The synthetic demos do not require real observations, but initialization and
sky analysis may download DE421 and GNSS TLEs. To use real data, edit the FITS
path in the commented examples. To consolidate NPZ data, provide the local
dataset expected by the script; the dataset is excluded from this package.

```python
from gnssvis import GNSSVis

gv = GNSSVis()
gv.plot_waterfall_with_sky("observation.fits", beamwidth_deg=30.0)
detrended = gv.baseline_subtract(axis=1, plot=False)
smoothed = gv.gaussian_smooth(detrended, sigma_freq=2.0, sigma_time=2.0)
```

Generated FITS files, plots, and caches are ignored. Keep code, messages,
comments, and documentation in English. The future simulation pipeline lives
in the separate sibling `gnssvis-simulation` repository.
