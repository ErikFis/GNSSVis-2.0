
# core.py: unified GNSSVis interface.

from typing import Dict, Optional
import numpy as np
from scipy.ndimage import gaussian_filter
from .io import IOManager
from .sky import SkyAnalyzer


class DataAnalyzer:
    """Statistical and signal processing utilities for waterfalls."""

    @staticmethod
    def normalize(data: np.ndarray) -> np.ndarray:
        """Normalize data to the interval [0, 1]."""
        dmin = float(np.nanmin(data))
        dmax = float(np.nanmax(data))
        if dmax == dmin:
            return np.zeros_like(data)
        return (data - dmin) / (dmax - dmin)

    @staticmethod
    def compute_statistics(waterfall_data: np.ndarray) -> Dict[str, float]:
        """Compute the mean, standard deviation, minimum, and maximum."""
        return {
            "mean_db": float(np.mean(waterfall_data)),
            "std_db": float(np.std(waterfall_data)),
            "min_db": float(np.min(waterfall_data)),
            "max_db": float(np.max(waterfall_data)),
        }

    @staticmethod
    def smooth_waterfall(
        waterfall_data: np.ndarray,
        sigma_freq: float = 1.0,
        sigma_time: float = 1.0,
    ) -> np.ndarray:
        """Apply 2D Gaussian smoothing along the frequency and time axes."""
        return gaussian_filter(waterfall_data, sigma=(sigma_freq, sigma_time))

    @staticmethod
    def baseline_subtract(waterfall_data: np.ndarray, axis: int = 0) -> np.ndarray:
        """Subtract the median along an axis from a waterfall.
        
        Args:
            waterfall_data: Two-dimensional waterfall (frequency x time).
            axis: Median reduction axis: 0 reduces frequency, 1 reduces time.
        
        Returns:
            Waterfall with the baseline removed."""
        baseline = np.median(waterfall_data, axis=axis, keepdims=True)
        return waterfall_data - baseline


class GNSSVis:
    """Unified gnssvis interface.
    
    Coordinate the three package layers:
    - IOManager: FITS input and synthetic data generation.
    - SkyAnalyzer: Ephemerides and field-of-view object detection.
    - Visualizer: Synchronized plots.
    
    Example::
    
        from gnssvis import GNSSVis
        gv = GNSSVis()
        gv.plot_waterfall_with_sky("observation.fits", beamwidth_deg=30.0)
        gv.baseline_subtract(axis=1)"""

    def __init__(
        self,
        observatory_lat: float = -7.211637,
        observatory_lon: float = -35.908138,
        observatory_alt: float = 553.0,
    ) -> None:
        """Initialize the observatory and load planetary ephemerides.
        
        Args:
            observatory_lat: Geodetic latitude in degrees (negative south).
            observatory_lon: Geodetic longitude in degrees (negative west).
            observatory_alt: Altitude above sea level in meters."""
        self.io = IOManager()
        self.sky = SkyAnalyzer(observatory_lat, observatory_lon, observatory_alt)
        self._viz = None  # Lazy loading to avoid a circular import

    @property
    def viz(self):
        """Return the Visualizer instance, created on first access."""
        if self._viz is None:
            from .plot import Visualizer
            self._viz = Visualizer()
        return self._viz

    # ------------------------------------------------------------------
    def plot_waterfall_with_sky(
        self, fits_file: str, beamwidth_deg: float = 30.0
    ) -> None:
        """Load a FITS observation, detect field-of-view objects, and plot
        the waterfall synchronized with celestial trajectories.
        
        Store the context for baseline_subtract and gaussian_smooth.
        
        Args:
            fits_file: Path to the observation FITS file.
            beamwidth_deg: Full telescope beam width in degrees."""
        self._fits_file = fits_file
        self._beamwidth_deg = beamwidth_deg
        self._obs = self.io.load_fits_observation(fits_file)
        self._pointing = self.io.get_observer_pointing_from_fits(fits_file)
        self._objects = self.sky.find_objects_in_fov(
            self._obs["times"], self._pointing, beamwidth_deg,
            include_gnss=True, include_pulsars=False,
        )
        self._waterfall = self._obs["data"]
        self._frequencies = self._obs["frequencies"]
        self._times = self._obs["times"]
        self.viz.plot_waterfall_with_objects(
            data=self._waterfall,
            times=self._times,
            frequencies=self._frequencies,
            sky_objects=self._objects,
            telescope_pointing=self._pointing,
            beamwidth_deg=self._beamwidth_deg,
        )
        self._print_detected_objects(self._objects, self._pointing, beamwidth_deg)

    # ------------------------------------------------------------------
    def baseline_subtract(self, axis: int = 1, plot: bool = True) -> np.ndarray:
        """Subtract the median baseline from the last loaded waterfall.
        
        Args:
            axis: Median reduction axis: 0 reduces frequency, 1 reduces time.
            plot: Whether to plot the result with updated trajectories.
        
        Returns:
            Two-dimensional waterfall with the baseline removed.
        
        Raises:
            RuntimeError: If plot_waterfall_with_sky has not loaded a context."""
        if not hasattr(self, "_waterfall"):
            raise RuntimeError(
                "No context loaded. Call plot_waterfall_with_sky first."
            )
        detrended = DataAnalyzer.baseline_subtract(self._waterfall, axis=axis)
        objects = self.sky.find_objects_in_fov(
            self._times, self._pointing, self._beamwidth_deg,
            include_gnss=True, include_pulsars=False,
        )
        if plot:
            self.viz.plot_waterfall_with_objects(
                data=detrended,
                times=self._times,
                frequencies=self._frequencies,
                sky_objects=objects,
                telescope_pointing=self._pointing,
                beamwidth_deg=self._beamwidth_deg,
            )
        return detrended

    # ------------------------------------------------------------------
    def gaussian_smooth(
        self,
        data: Optional[np.ndarray] = None,
        sigma_freq: float = 1.0,
        sigma_time: float = 1.0,
        plot: bool = True,
    ) -> np.ndarray:
        """Apply Gaussian smoothing to a waterfall.
        
        Args:
            data: Two-dimensional array (frequency x time). If None, use the
                last waterfall loaded by plot_waterfall_with_sky.
            sigma_freq: Filter standard deviation along frequency.
            sigma_time: Filter standard deviation along time.
            plot: Whether to plot the result with updated trajectories.
        
        Returns:
            Smoothed two-dimensional waterfall.
        
        Raises:
            RuntimeError: If data is None and no context has been loaded."""
        if data is None:
            if not hasattr(self, "_waterfall"):
                raise RuntimeError(
                    "No data supplied. Pass `data` or call "
                    "plot_waterfall_with_sky first."
                )
            data = self._waterfall

        smoothed = DataAnalyzer.smooth_waterfall(
            np.asarray(data), sigma_freq=sigma_freq, sigma_time=sigma_time
        )

        if plot:
            objects = self.sky.find_objects_in_fov(
                self._times, self._pointing, self._beamwidth_deg,
                include_gnss=True, include_pulsars=False,
            )
            self.viz.plot_waterfall_with_objects(
                data=smoothed,
                times=self._times,
                frequencies=self._frequencies,
                sky_objects=objects,
                telescope_pointing=self._pointing,
                beamwidth_deg=self._beamwidth_deg,
            )

        return smoothed

    # ------------------------------------------------------------------
    def quick_demo(self) -> Dict:
        """Run a quick demonstration with a synthetic FITS observation.
        
        Create a waterfall with Gaussian noise, run the full analysis, and return
        the results. Remove the temporary FITS file afterward. Initialization
        and sky analysis may download ephemerides and GNSS TLE data.
        
        Returns:
            Dictionary with objects, pointing, times, frequencies, and data."""
        import os
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as tmp:
            fits_path = tmp.name

        try:
            self.io.create_dummy_fits(fits_path)
            self.plot_waterfall_with_sky(fits_path, beamwidth_deg=30.0)
            return {
                "objects": self._objects,
                "pointing": self._pointing,
                "times": self._times,
                "frequencies": self._frequencies,
                "data": self._waterfall,
            }
        finally:
            if os.path.exists(fits_path):
                os.remove(fits_path)

    # ------------------------------------------------------------------
    def _print_detected_objects(self, objects, pointing, beamwidth_deg: float) -> None:
        """Print a summary of objects detected within the beam."""
        pointing_dec = pointing.dec.deg
        dec_min = pointing_dec - beamwidth_deg / 2.0
        dec_max = pointing_dec + beamwidth_deg / 2.0
        detected = []
        for name, coord in objects.items():
            try:
                deg = coord.dec.deg
                dec_val = float(np.mean(deg)) if np.ndim(deg) > 0 else float(deg)
            except Exception:
                dec_val = float("nan")
            if dec_min <= dec_val <= dec_max:
                detected.append((name, dec_val))

        print(f"\n{'='*50}")
        print("OBJECTS DETECTED IN THE BEAM")
        print(f"{'='*50}")
        print(f"Pointing declination: {pointing_dec:.2f} deg")
        print(f"Beam: {beamwidth_deg:.1f} deg ({dec_min:.2f} to {dec_max:.2f} deg)")
        print(f"Total: {len(detected)} of {len(objects)} objects\n")
        for i, (name, dec_val) in enumerate(detected, 1):
            print(f"  {i:2d}. {name:30s} Dec = {dec_val:7.2f} deg")
        if not detected:
            print("  No objects within the beam.")
