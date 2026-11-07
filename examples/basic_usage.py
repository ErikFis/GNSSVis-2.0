#!/usr/bin/env python3
"""GNSSVis basic usage example for radio telescope observation analysis."""

from gnssvis import GNSSVis


def main():
    """Demonstrate GNSSVis usage."""

    print("GNSSVis - Usage Examples")
    print("=" * 50)

    # Create an instance at the UIRAPURU observatory (UFCG)
    gv = GNSSVis(
        observatory_lat=-7.211637,
        observatory_lon=-35.908138,
        observatory_alt=553.0,
    )

    # 1. Quick synthetic FITS demo (sky data downloads may be required)
    print("\nRunning a quick demo with synthetic data...")
    resultado = gv.quick_demo()
    print(f"Demo complete!")
    print(f"  Detected objects: {len(resultado['objects'])}")
    print(f"  Timesteps: {len(resultado['times'])}")
    print(f"  Frequency channels: {len(resultado['frequencies'])}")

    # 2. Analyze a real FITS file (uncomment and adjust the path)
    # gv.plot_waterfall_with_sky("observacao.fits", beamwidth_deg=30.0)
    # detrended = gv.baseline_subtract(axis=1)
    # smoothed = gv.gaussian_smooth(sigma_freq=2.0, sigma_time=2.0)

    print("\nSee notebooks/ for interactive examples.")


if __name__ == "__main__":
    main()