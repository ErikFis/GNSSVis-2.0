# Public API

## GNSSVis

```python
from gnssvis import GNSSVis

gv = GNSSVis(
    observatory_lat=-7.211637,
    observatory_lon=-35.908138,
    observatory_alt=553.0,
)
gv.plot_waterfall_with_sky("observation.fits", beamwidth_deg=30.0)
```

The constructor initializes the UIRAPURU observatory coordinates by default
and loads planetary ephemerides. The visualizer is created lazily.

| Method | Behavior | Return |
| --- | --- | --- |
| `plot_waterfall_with_sky(fits_file, beamwidth_deg=30.0)` | Load, analyze, plot, and store context | `None` |
| `baseline_subtract(axis=1, plot=True)` | Subtract median baseline from loaded data | 2D NumPy array |
| `gaussian_smooth(data=None, sigma_freq=1.0, sigma_time=1.0, plot=True)` | Smooth supplied or loaded data | 2D NumPy array |
| `quick_demo()` | Analyze temporary synthetic FITS data | Results dictionary |

`baseline_subtract` requires a loaded context. `gaussian_smooth` requires a
loaded context when using stored data or plotting trajectories. To smooth an
independent array, supply `data` and use `plot=False`.

`quick_demo()` takes no filename argument. It returns `objects`, `pointing`,
`times`, `frequencies`, and `data`, then removes its temporary FITS file.
It can still download astronomical data.

## IOManager

```python
from gnssvis import IOManager

io = IOManager()
io.create_dummy_fits("synthetic.fits", n_time=100, n_freq=32)
observation = io.load_fits_observation("synthetic.fits")
pointing = io.get_observer_pointing_from_fits("synthetic.fits")
```

`load_fits_observation(fits_path, date_override=None, manual_time=None)`
returns a dictionary containing `data`, `times`, `frequencies`, and `header`.
Frequency is in MHz. `manual_time` supplies ISO UTC timestamps.
The loader handles TIME-extension timestamps, midnight rollovers, and a
SAMPRATE fallback. The pointing loader retains a legacy correction from
declination near -60 degrees to -30 degrees.

`consolidate_observation_from_npz(...)` converts observation NPZ files to
FITS. Its current implementation reads CH1 only, retains the CH2 argument,
interpolates timestamp subdivisions, and converts linear power to dB.

`IOManager.gaussian_smooth(...)` is a numerical utility without plotting.
Supplying an explicit array is recommended; loading a FITS file returns data
and does not automatically populate `IOManager._waterfall`.

## SkyAnalyzer

```python
from gnssvis import SkyAnalyzer

sky = SkyAnalyzer()
objects = sky.find_objects_in_fov(
    observation["times"], pointing, beamwidth_deg=20.0,
    include_gnss=True, include_pulsars=False,
)
sky.print_objects_summary(objects, pointing, 20.0, observation["times"])
```

Direct methods include `get_solar_system_objects`, `get_gnss_trajectories`,
`get_pulsars_catalog`, `categorize_objects`, `get_object_display_info`, and
`get_objects_statistics`. `clear_cache` and `clear_pulsar_cache` clear
in-memory caches; they do not delete persistent files.

## DataAnalyzer and Visualizer

`DataAnalyzer` provides static `normalize`, `compute_statistics`,
`smooth_waterfall`, and `baseline_subtract` methods.

`Visualizer.plot_waterfall_with_objects` accepts the waterfall, timestamps,
frequencies, sky-object coordinates, optional telescope pointing, and beam
width. See the [plotting guide](PLOTTING_GUIDE.md) for a complete example.
