# Local data directory

Observation datasets, FITS files, NumPy arrays, downloaded catalogs, caches,
and ephemerides are intentionally excluded from this prepared repository.
Only Markdown documentation in the data directory is versioned.

Place private observation inputs here locally or use an external directory.
`SkyAnalyzer` downloads astronomical resources into `data/tle_repository/`.
The DE421 ephemeris alone is approximately 17 MB and must remain untracked.

To generate a small local observation without downloading sky data:

```python
from gnssvis import IOManager

IOManager().create_dummy_fits("data/synthetic.fits", n_freq=32, n_time=100)
```

The generated file is ignored by Git. `GNSSVis().quick_demo()` additionally
performs sky analysis and can require downloads.
