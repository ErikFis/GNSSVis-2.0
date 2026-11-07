 #Author: Erik Limeira

def consolidate_observation_from_npz(
        data_dir: str = "/home/erik-limeira/gnssvis/data/2025-09-10",
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
        import os
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
                        ss = int(t_sec % 60)
                        all_times.append(f"{hh:02d}:{mm:02d}:{ss:02d}")
                else:
                    # Last timestamp: add the 10 subdivisions
                    for j in range(10):
                        t_sec = current_seconds + j
                        hh = int(t_sec // 3600) % 24
                        mm = int((t_sec % 3600) // 60)
                        ss = int(t_sec % 60)
                        all_times.append(f"{hh:02d}:{mm:02d}:{ss:02d}")
            
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
        print(f"[INFO] Using CH1 only")
        
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
        primary_hdu.header["Author"] = "Erik Limeira"
        primary_hdu.header["RA"] = ra
        primary_hdu.header["DEC"] = dec
        primary_hdu.header["COMMENT"] = f"Observation concatenated from {data_dir}"
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
        timestamps_array = np.array([t.split('.')[0] for t in times], dtype='U8')  # Remove fractional seconds for FITS
        col = fits.Column(name="TIMESTAMP", format="8A", array=timestamps_array)
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
