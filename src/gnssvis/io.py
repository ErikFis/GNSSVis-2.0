"""FITS input and output for gnssvis.

Read and standardize observations, extract RA/DEC pointing, and generate
synthetic FITS observations for testing."""

from __future__ import annotations
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import Dict, Any, Optional, Tuple, Union
from scipy.ndimage import gaussian_filter


class IOManager:
    """Manage FITS file input and output."""

    # ------------------------------------------------------------------
    def load_fits_observation(
        self,
        fits_path: str,
        date_override: Optional[str] = None,
        manual_time: Optional[Union[str, list, np.ndarray]] = None
    ) -> Dict[str, Any]:
        """Load FITS data as a dictionary with data, times, frequencies, and header.
        
        Supported layouts include Horn.fits (PRIMARY, TOD, SPECTRA) and
        observacao_*.fits (PRIMARY, DATA, TIME). Frequencies are returned in MHz.
        Use date_override to replace the observation date and manual_time to
        supply ISO UTC timestamps. Without valid TIME data, derive times from
        SAMPRATE, starting at midnight on the observation date."""
        with fits.open(fits_path) as hdul:
            primary_header = hdul[0].header
            date_obs_raw = primary_header.get("DATE-OBS", "2025-01-01")
            
            # Extract the date (YYYY-MM-DD) if a time component is present
            if date_override:
                date_obs = date_override
            elif "T" in str(date_obs_raw):
                date_obs = str(date_obs_raw).split("T")[0]
            else:
                date_obs = str(date_obs_raw)
            
            # Identify the data extension (TOD or DATA)
            data_ext = None
            for hdu in hdul[1:]:
                if hasattr(hdu, 'name') and hdu.name in ['TOD', 'DATA', 'SPECTRA']:
                    if hdu.name != 'SPECTRA':  # SPECTRA is a table, not an image
                        data_ext = hdu
                        break
            
            if data_ext is None:
                data_ext = hdul[1]  # Fallback
            
            img_header = data_ext.header
            data = data_ext.data.astype(np.float64)
            n_freq, n_time = data.shape if data.ndim == 2 else (1, len(data))

            # Frequencies (Hz to MHz)
            # For Horn.fits: FREQCENT and BANDWID in the primary header
            if "FREQCENT" in primary_header:
                freq_center = primary_header.get("FREQCENT", 1120.0)  # MHz
                bandwidth = primary_header.get("BANDWID", 187.5)     # MHz
                delta_freq = (bandwidth * 1e6) / n_freq
                ref_freq = (freq_center - bandwidth/2) * 1e6
            else:
                # For observacao_*.fits: CRVAL2 and CDELT2
                ref_freq = img_header.get("CRVAL2", 1.2e9)
                delta_freq = img_header.get("CDELT2", 1.0e6)
            
            ref_pix = img_header.get("CRPIX2", 1) - 1
            freqs = (ref_freq + (np.arange(n_freq) - ref_pix) * delta_freq) / 1e6

            # Time axis
            if manual_time is not None:
                times = Time(manual_time, format="isot", scale="utc")
            else:
                try:
                    # Look for the TIME extension
                    time_ext = None
                    for hdu in hdul[1:]:
                        if hasattr(hdu, 'name') and hdu.name in ['TIME']:  # Only TIME, not SPECTRA
                            time_ext = hdu
                            break
                    
                    if time_ext is not None:
                        timestamps = time_ext.data["TIMESTAMP"].astype(str)
                        # Filter empty timestamps
                        timestamps = np.array([ts for ts in timestamps if ts.strip() and ts.strip() != ''])
                    else:
                        raise Exception("Missing TIME extension")
                    
                    # Use valid timestamps when available
                    if len(timestamps) > 0:
                        # Detect a day rollover when timestamps move backward
                        time_list = []
                        current_date = date_obs
                        prev_hour = None
                        
                        for ts in timestamps:
                            try:
                                hour = int(ts.split(':')[0])
                            except Exception:
                                continue
                            
                            # Increment the day when the hour drops sharply (for example, 23 to 00)
                            if prev_hour is not None and hour < prev_hour - 12:
                                # Day rollover detected
                                from datetime import datetime, timedelta
                                dt = datetime.strptime(current_date, '%Y-%m-%d')
                                dt += timedelta(days=1)
                                current_date = dt.strftime('%Y-%m-%d')
                            
                            time_list.append(f"{current_date}T{ts}")
                            prev_hour = hour
                        
                        times = Time(time_list, format="isot", scale="utc")
                        print("[INFO] Timestamps read from the TIME extension")
                        print(f"       First: {time_list[:2]}")
                        print(f"       Last: {time_list[-2:]}")
                    else:
                        raise Exception("Empty extension")
                except Exception as e:
                    print(f"[WARN] FITS has no valid TIME extension ({e}), generating the time axis automatically.")
                    # Use primary header SAMPRATE to calculate the sampling interval
                    samprate = primary_header.get("SAMPRATE", 1.0)
                    dt_seconds = 1.0 / samprate  # Sampling interval in seconds
                    print(f"        SAMPRATE = {samprate} Hz → dt = {dt_seconds} s/sample")
                    print(f"        Total time = {n_time * dt_seconds / 3600:.2f} h")
                    t0 = Time(date_obs + "T00:00:00", format="isot", scale="utc")
                    times = t0 + np.arange(n_time) * dt_seconds * u.s

            if data.shape[0] != len(freqs):
                data = data.T

            return {"data": data, "times": times, "frequencies": freqs, "header": primary_header}

    # ------------------------------------------------------------------
    def get_observer_pointing_from_fits(self, fits_path: str) -> SkyCoord:
        """Extract RA/DEC pointing from the FITS header.
        
        Correct DEC near -60 degrees to -30 degrees for affected observation
        files, following the existing Horn.fits convention."""
        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            ra = hdr.get("RA", hdr.get("OBJCTRA", 0.0))
            dec = hdr.get("DEC", hdr.get("OBJCTDEC", 0.0))
            
            # Correct declination near -60 degrees to -30 degrees (legacy convention)
            if abs(float(dec) - (-60.0)) < 0.1:
                dec = -30.0
        
        # FITS coordinates are already in degrees
        return SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg)

    # ------------------------------------------------------------------
    def create_dummy_fits(
        self,
        filename: str = "dummy_waterfall.fits",
        n_time: int = 3600,
        n_freq: int = 200,
        t0: str = "2025-10-27T00:00:00",
        dt: float = 1.0,
        f_start: float = 1200.0 * 1e6,
        df: float = 1.0 * 1e6,
        noise_db: Tuple[float, float] = (-87.5, 2.0)
    ) -> None:
        """Generate a synthetic FITS observation with timestamps and Gaussian noise.
        
        Args:
            filename: Output file path.
            n_time: Number of time samples.
            n_freq: Number of frequency channels.
            t0: Initial ISO timestamp.
            dt: Sampling interval in seconds.
            f_start: Initial frequency in Hz.
            df: Channel spacing in Hz.
            noise_db: Noise mean and standard deviation in dB."""
        mean, sigma = noise_db
        data = np.random.normal(loc=mean, scale=sigma, size=(n_freq, n_time))

        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header["DATE-OBS"] = t0.split("T")[0]
        primary_hdu.header["TELESCOP"] = "UIRAPURU"
        primary_hdu.header["OBSERVER"] = "Simulated"
        # Astronomical coordinates selected for object detection
        # Point toward the region with detected planets (Saturn/Neptune)
        primary_hdu.header["RA"] = 0.0       # 0h right ascension (meridian)
        primary_hdu.header["DEC"] = -2.0     # -2 degrees declination (planet region center)
        primary_hdu.header["COMMENT"] = "Synthetic FITS for gnssvis testing"
        primary_hdu.header["COMMENT"] = "RA/DEC: realistic astronomical coordinates"

        image_hdu = fits.ImageHDU(data=data, name="DATA")
        image_hdu.header["NAXIS1"] = n_time
        image_hdu.header["NAXIS2"] = n_freq
        image_hdu.header["CRVAL2"] = f_start
        image_hdu.header["CDELT2"] = df
        image_hdu.header["CRPIX2"] = 1
        image_hdu.header["CTYPE2"] = "FREQ"
        image_hdu.header["CUNIT2"] = "Hz"

        # Create the synthetic timestamp column
        seconds = np.arange(n_time) * dt
        hh = (seconds // 3600).astype(int)
        mm = ((seconds % 3600) // 60).astype(int)
        ss = (seconds % 60).astype(int)
        timestamps = np.array([f"{h:02d}:{m:02d}:{s:02d}" for h, m, s in zip(hh, mm, ss)])
        col = fits.Column(name="TIMESTAMP", format="8A", array=timestamps)
        table_hdu = fits.BinTableHDU.from_columns([col], name="TIME")

        fits.HDUList([primary_hdu, image_hdu, table_hdu]).writeto(filename, overwrite=True)
        print(f"[OK] Synthetic FITS created: {filename} ({n_freq}x{n_time})")

    # ------------------------------------------------------------------
    def consolidate_observation_from_npz(
        self,
        data_dir: str = "",
        output_filename: str = "observacao_2025-09-10.fits",
        date_obs: str = "2025-09-10",
        ra: float = 0.0,
        dec: float = -2.0,
        f_start_ch1: float = 1200.0e6,
        f_start_ch2: float = 1400.0e6,
        df: float = 122.07031e3
    ) -> str:
        """Consolidate NPZ observation files into a single FITS file.
        
        Read Averages_CH*.npz files from channel subdirectories or the data
        directory, concatenate spectra, and write timestamps and pointing.
        The current implementation processes CH1 only; the CH2 frequency
        argument is retained for API compatibility.
        
        Args:
            data_dir: Root directory containing observation NPZ files.
            output_filename: Output FITS filename.
            date_obs: Observation date (YYYY-MM-DD).
            ra: Pointing right ascension in degrees.
            dec: Pointing declination in degrees.
            f_start_ch1: Initial CH1 frequency in Hz.
            f_start_ch2: Initial CH2 frequency in Hz (currently unused).
            df: Frequency channel width in Hz.
        
        Returns:
            Path to the generated FITS file."""
        from pathlib import Path
        
        data_dir_path = Path(data_dir)
        
        # Detect channel subdirectories or files directly in the data directory
        ch1_dir = data_dir_path / "ch1"
        has_subdirs = ch1_dir.exists()
        
        if has_subdirs:
            # Legacy layout: ch1/ and ch2/
            search_pattern = "Averages_CH1_*.npz"
            npz_dir = ch1_dir
        else:
            # New layout: files directly in the directory
            search_pattern = "Averages_CH1_*.npz"
            npz_dir = data_dir_path
        
        # Helper to load and process one channel
        def load_channel(channel_dir, channel_num):
            npz_files = sorted(channel_dir.glob(search_pattern))
            if not npz_files:
                print(f"[WARN] No NPZ files found in {channel_dir}")
                return None, None
            
            all_data = []
            all_keys = []
            
            for npz_file in npz_files:
                data = np.load(npz_file)
                # Each file contains multiple timestamp keys
                for key in sorted(data.keys()):
                    all_keys.append(key)
                    # data[key] has shape (n_avg, n_freq)
                    # Apply FFT shift and transpose as in the original notebook
                    spectra = np.fft.fftshift(data[key], axes=1)
                    spectra = np.asarray(spectra).T  # shape: (n_freq, n_avg)
                    all_data.append(spectra)
            
            if not all_data:
                return None, None
            
            # Concatenate all spectra along the time axis
            # Each element has shape (n_freq, 10); concatenate along axis 1
            channel_data = np.concatenate(all_data, axis=1)  # shape: (n_freq, n_time_total)
            
            # Generate interpolated timestamps as in the original notebook
            from datetime import datetime
            timestamps = [datetime.strptime(k, '%H%M%S').strftime('%H:%M:%S') for k in all_keys]
            
            # Interpolate timestamps for each subdivision using NumPy
            all_times = []
            
            # Convert the first timestamp to reference seconds within the day
            h0, m0, s0 = map(int, timestamps[0].split(':'))
            base_day_seconds = h0 * 3600 + m0 * 60 + s0
            
            for i in range(len(timestamps)):
                # Each key has 10 subdivisions
                # Convert the base timestamp to seconds
                h, m, s = map(int, timestamps[i].split(':'))
                current_seconds = h * 3600 + m * 60 + s
                
                # Detect a day rollover when the time precedes the initial time
                if current_seconds < base_day_seconds:
                    current_seconds += 86400  # Add 24 hours
                
                # Interpolate to the next timestamp unless this is the last one
                if i < len(timestamps) - 1:
                    h_next, m_next, s_next = map(int, timestamps[i+1].split(':'))
                    next_seconds = h_next * 3600 + m_next * 60 + s_next
                    
                    # Detect a day rollover in the next timestamp
                    if next_seconds < base_day_seconds:
                        next_seconds += 86400
                    
                    interval = (next_seconds - current_seconds) / 10.0
                    
                    for j in range(10):
                        t_sec = current_seconds + j * interval
                        hh = int(t_sec // 3600) % 24
                        mm = int((t_sec % 3600) // 60)
                        ss = t_sec % 60  # Keep fractional seconds
                        all_times.append(f"{hh:02d}:{mm:02d}:{ss:06.3f}")
                else:
                    # Last timestamp: add the 10 subdivisions
                    for j in range(10):
                        t_sec = current_seconds + j * 0.1  # 100ms increments for last block
                        hh = int(t_sec // 3600) % 24
                        mm = int((t_sec % 3600) // 60)
                        ss = t_sec % 60  # Keep fractional seconds
                        all_times.append(f"{hh:02d}:{mm:02d}:{ss:06.3f}")
            
            if not all_data:
                return None, None
            
            # channel_data has shape (n_freq, n_time_total) after concatenation
            # Transpose to (n_time, n_freq)
            channel_data = channel_data.T
            
            return channel_data, all_times
        
        # Load data
        print(f"[INFO] Loading data from {npz_dir}")
        data_ch1, times_ch1 = load_channel(npz_dir, 1)
        
        # Use CH1 only
        if data_ch1 is None:
            raise ValueError("No data found in channel 1")
        
        data = data_ch1
        times = times_ch1
        n_freq = data.shape[1]
        freqs = f_start_ch1 + np.arange(n_freq) * df
        print("[INFO] Using CH1 only")
        
        # Transpose to FITS layout: (n_freq, n_time)
        data = data.T
        n_freq, n_time = data.shape
        
        print(f"[INFO] Consolidated data: {n_freq} channels x {n_time} time samples")
        
        # Convert to dB, assuming linear power data
        # Comment out this line if the data is already in dB
        data_db = 10 * np.log10(data + 1e-20)  # Avoid log(0)
        
        # Create FITS HDUs
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header["DATE-OBS"] = date_obs
        primary_hdu.header["TELESCOP"] = "UIRAPURU"
        primary_hdu.header["OBSERVER"] = "Consolidated"
        primary_hdu.header["RA"] = ra
        primary_hdu.header["DEC"] = dec
        primary_hdu.header["COMMENT"] = f"Observation consolidated from {data_dir}"
        primary_hdu.header["COMMENT"] = "Data from multiple NPZ files"
        
        # Data HDU
        image_hdu = fits.ImageHDU(data=data_db.astype(np.float32), name="DATA")
        image_hdu.header["NAXIS1"] = n_time
        image_hdu.header["NAXIS2"] = n_freq
        image_hdu.header["CRVAL2"] = float(freqs[0])
        image_hdu.header["CDELT2"] = float(df)
        image_hdu.header["CRPIX2"] = 1
        image_hdu.header["CTYPE2"] = "FREQ"
        image_hdu.header["CUNIT2"] = "Hz"
        
        # Timestamp HDU
        timestamps_array = np.array(times, dtype='U12')  # Preserve fractional seconds (HH:MM:SS.sss)
        col = fits.Column(name="TIMESTAMP", format="12A", array=timestamps_array)
        table_hdu = fits.BinTableHDU.from_columns([col], name="TIME")
        
        # Save FITS
        output_path = str(Path(data_dir).parent / output_filename)
        fits.HDUList([primary_hdu, image_hdu, table_hdu]).writeto(
            output_path, overwrite=True
        )
        
        print(f"[OK] Consolidated FITS created: {output_path}")
        print(f"     Shape: {n_freq} freq × {n_time} time")
        print(f"     Frequency: {freqs[0]/1e6:.2f} - {freqs[-1]/1e6:.2f} MHz")
        print(f"     Time: {times[0]} - {times[-1]}")
        
        return output_path

    def gaussian_smooth(self, data: Optional[np.ndarray] = None,
                        sigma_freq: float = 1.0,
                        sigma_time: float = 1.0) -> np.ndarray:
        """Apply Gaussian smoothing to data (frequency x time), or use
        self._waterfall when data is None.
        
        This utility does not plot or access self.sky or self.viz. For integrated
        plots with synchronized trajectories, use GNSSVis.gaussian_smooth."""
        if data is None:
            if not hasattr(self, '_waterfall'):
                raise RuntimeError('No data supplied. Pass `data` or load a waterfall with `load_fits_observation`.')
            data = self._waterfall

        data = np.asarray(data)
        smoothed = gaussian_filter(data, sigma=(sigma_freq, sigma_time))
        return smoothed
