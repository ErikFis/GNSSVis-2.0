#!/usr/bin/env python3
"""GNSSVis usage examples.

Analyze radio telescope observations and detect GNSS satellites.

Author: Erik H L Guimaraes"""

from gnssvis import GNSSVis, IOManager, SkyAnalyzer, Visualizer


def exemplo_demo_rapido():
    """Run a synthetic FITS demonstration without real observation data."""
    print("Example 1: Quick demo")

    gv = GNSSVis()
    results = gv.quick_demo()

    print(f"  Detected objects: {len(results['objects'])}")
    print(f"  Frequency channels: {len(results['frequencies'])}")
    return results


def exemplo_observatorio_customizado():
    """Analyze data using custom observatory coordinates."""
    print("\nExample 2: Custom observatory")

    gv = GNSSVis(
        observatory_lat=-22.9068,
        observatory_lon=-43.1729,
        observatory_alt=10.0,
    )

    # Synthetic demonstration
    results = gv.quick_demo()
    print(f"  Observatory configured.")
    print(f"  Detected objects: {len(results['objects'])}")
    return results


def exemplo_pipeline_processamento():
    """Demonstrate waterfall, baseline subtraction, and smoothing."""
    print("\nExample 3: Processing pipeline (requires real FITS data)")

    gv = GNSSVis()

    # Uncomment to use a real FITS file:
    # gv.plot_waterfall_with_sky("observacao.fits", beamwidth_deg=30.0)
    # detrended = gv.baseline_subtract(axis=1, plot=True)
    # smoothed = gv.gaussian_smooth(sigma_freq=2.0, sigma_time=2.0, plot=True)

    print("  (example commented out - supply a real FITS file)")


def exemplo_modulos_separados():
    """Use the internal modules directly for finer control."""
    print("\nExample 4: Separate modules")

    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import numpy as np

    # Set the observation period
    t0 = Time("2025-10-27T20:00:00", format="isot", scale="utc")
    times = t0 + np.linspace(0, 4, 100) * u.hour

    # Telescope pointing
    pointing = SkyCoord(ra=0.0 * u.deg, dec=-7.0 * u.deg, frame="icrs")

    # Find objects in the field of view
    sky = SkyAnalyzer()
    objects = sky.find_objects_in_fov(times, pointing, beamwidth_deg=30.0)
    print(f"  Objects found: {len(objects)}")
    sky.print_objects_summary(objects, pointing, 30.0, times)


if __name__ == "__main__":
    print("GNSSVis - Usage Examples")
    print("=" * 50)

    exemplo_demo_rapido()
    exemplo_observatorio_customizado()
    exemplo_pipeline_processamento()
    exemplo_modulos_separados()

    print("\nSee notebooks/ for complete interactive examples.")