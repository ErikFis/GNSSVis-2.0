# Local 2.0 simulation upgrade

- Added `range_rate_m_s` to SkyAnalyzer satellite alt/az results and geometry.
- Added `beam.py`: spherical offsets, rotated elliptical Gaussian main lobe,
  frequency-dependent FWHM and continuous power response.
- Added `realistic.py`: importable simulation API and separate CLI, exact
  channel-averaged BPSK proxy, physical radial Doppler, DE421 apparent Sun,
  quiet thermal continuum and optional active nonthermal continuum.
- Added compact NPZ geometry, compressed separate power components, hashes,
  metadata, half-power passages and five presentation PNGs.
- Added an immutable 32-GPS historical snapshot and selected 24 h scenario.
- Added example, model/reproduction documentation and physical-invariant tests.
- Preserved the original offline scenario/CLI; the legacy CLI rejects a new
  wideband scenario to avoid silently ignoring its Sun and spectral options.
- No PyPI publication or package version change was performed.
