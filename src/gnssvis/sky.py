"""Ephemerides and astronomical objects with topocentric projection.

Support DE421 solar system ephemerides, GNSS satellites from Celestrak
TLEs, and an ATNF pulsar catalog cache. Coordinates vary with time.

Author: Erik H L Guimaraes"""

from __future__ import annotations
from typing import Dict
import warnings
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, CartesianRepresentation
from skyfield.api import EarthSatellite, Topos
# from psrqpy import QueryATNF  # Temporarily disabled

warnings.filterwarnings("ignore", message="you can only subtract two vectors")


class SkyAnalyzer:
    """Calculate ephemerides and find objects in the field of view."""

    def __init__(self, lat_deg: float = -7.211637, lon_deg: float = -35.908138, elev_m: float = 553.0) -> None:
        """Initialize the observatory and load planetary ephemerides."""
        from pathlib import Path
        from skyfield.api import Loader
        
        self.lat_deg: float = lat_deg
        self.lon_deg: float = lon_deg
        self.elev_m: float = elev_m

        # Set a central directory for astronomical data
        project_root = Path(__file__).parent.parent.parent
        data_dir = project_root / "data" / "tle_repository"
        
        # Create data subdirectories
        ephemeris_dir = data_dir / "ephemeris"
        ephemeris_dir.mkdir(parents=True, exist_ok=True)
        
        # Skyfield: use a custom loader for the central cache
        # Download planetary ephemerides (de421.bsp) here automatically
        self.loader = Loader(str(ephemeris_dir))
        self.ts = self.loader.timescale()
        
        # Automatically download planetary ephemerides into a long-lived cache
        # Skyfield downloads the file if it is missing
        self.planets = self.loader('de421.bsp')
        self.earth = self.planets["earth"]
        self.observatory = self.earth + Topos(latitude_degrees=self.lat_deg,
                                              longitude_degrees=self.lon_deg,
                                              elevation_m=self.elev_m)

        # Astropy: location for AltAz transformations
        self.location = EarthLocation(lat=self.lat_deg * u.deg,
                                      lon=self.lon_deg * u.deg,
                                      height=self.elev_m * u.m)
        
        # Cache to avoid repeated downloads within a session
        self._pulsar_cache = None
        self._gnss_satellites_cache = None

    def clear_cache(self) -> None:
        """Clear in-memory catalog caches to allow new downloads."""
        self._pulsar_cache = None
        self._gnss_satellites_cache = None
        print("Cache cleared - subsequent queries may download data")
        
    def clear_pulsar_cache(self) -> None:
        """Clear only the in-memory pulsar cache."""
        self._pulsar_cache = None
        print("Pulsar cache cleared")

    # ----------------------------------------------------------------------
    def _ensure_time(self, times: Time | list[str]) -> Time:
        """Convert input timestamps to astropy.Time."""
        if not isinstance(times, Time):
            times = Time(times, format="isot", scale="utc")
        return times

    def _to_skyfield_time(self, times: Time):
        """Convert astropy.Time to Skyfield time."""
        return self.ts.from_astropy(times)

    # ----------------------------------------------------------------------
    def get_solar_system_objects(self, times: Time) -> Dict[str, SkyCoord]:
        """Return RA/Dec trajectories of the main planets and the Sun."""
        times = self._ensure_time(times)
        sky_times = self._to_skyfield_time(times)

        names = [
            "sun", "mercury", "venus", "mars",
            "jupiter barycenter", "saturn barycenter",
            "uranus barycenter", "neptune barycenter"
        ]
        traj: Dict[str, SkyCoord] = {}
        for name in names:
            try:
                astrometric = self.observatory.at(sky_times).observe(self.planets[name])
                ra, dec, _ = astrometric.radec()
                traj[name.capitalize()] = SkyCoord(ra=ra.degrees * u.deg,
                                                   dec=dec.degrees * u.deg,
                                                   obstime=times)
            except Exception:
                continue
        return traj


    # ----------------------------------------------------------------------
    def get_gnss_trajectories(
        self,
        times: Time,
        tle_url: str = "https://celestrak.org/NORAD/elements/gp.php?GROUP=gnss&FORMAT=tle",
    ) -> Dict[str, SkyCoord]:
        """Calculate GNSS trajectories with Skyfield SGP4 propagation.
        
        Return observer-relative topocentric RA/Dec positions. Download TLEs
        from the Celestrak gp.php API into data/tle_repository/gnss/. The session
        cache avoids repeated downloads during one run. Local files older than
        seven days require an update.
        
        Args:
            times: Observation timestamps as astropy.Time.
            tle_url: Celestrak API URL for GNSS TLEs.
        
        Returns:
            Dictionary mapping satellite names to SkyCoord positions at each time."""
        import requests
        from datetime import datetime
        from pathlib import Path
        import shutil

        times = self._ensure_time(times)
        sky_times = self._to_skyfield_time(times)
        traj: Dict[str, SkyCoord] = {}

        # TLE repository paths
        project_root = Path(__file__).parent.parent.parent
        tle_repo = project_root / "data" / "tle_repository" / "gnss"
        tle_repo.mkdir(parents=True, exist_ok=True)
        tle_file = tle_repo / "gnss.txt"

        # Check whether the local file is still fresh
        needs_update = True
        if tle_file.exists():
            file_age_days = (
                datetime.now() - datetime.fromtimestamp(tle_file.stat().st_mtime)
            ).days
            needs_update = file_age_days > 7
            if needs_update:
                print(
                    f"[TLE UPDATE] GNSS TLEs are {file_age_days} days old — "
                    "downloading updated data..."
                )
            else:
                print(f"[TLE CACHE] GNSS TLEs are {file_age_days} days old (valid).")

        satellites = None
        try:
            if self._gnss_satellites_cache is not None and not needs_update:
                print("[CACHE] Using GNSS TLEs from the current session cache.")
                satellites = self._gnss_satellites_cache
            else:
                # Download with requests from the Celestrak API (gp.php)
                print(f"[DOWNLOAD] Downloading GNSS TLEs from: {tle_url}")
                response = requests.get(tle_url, timeout=30)
                response.raise_for_status()
                tle_text = response.text

                # Save to the local repository
                tle_file.write_text(tle_text, encoding="utf-8")

                # Parse three-line blocks: name / line 1 / line 2
                lines = [ln.strip() for ln in tle_text.splitlines() if ln.strip()]
                satellites = []
                i = 0
                while i + 2 < len(lines):
                    if lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
                        sat = EarthSatellite(
                            lines[i + 1], lines[i + 2], lines[i], self.ts
                        )
                        satellites.append(sat)
                        i += 3
                    else:
                        i += 1

                self._gnss_satellites_cache = satellites
                print(f"[REPOSITORY] {len(satellites)} GNSS satellites loaded.")

                # Dated copy for the archive
                today = datetime.now().strftime("%Y%m%d")
                tle_dated = tle_repo / f"gnss_{today}.txt"
                if not tle_dated.exists():
                    shutil.copy(tle_file, tle_dated)
                    print(f"[ARCHIVE] Dated backup saved: {tle_dated.name}")

        except Exception as exc:
            print(f"[WARN] Failed to load GNSS TLEs: {exc}")
            return {}

        if not satellites:
            print("[WARN] No GNSS satellites available.")
            return {}

        # Geocentric observer
        observer = self.earth + Topos(
            latitude_degrees=self.lat_deg,
            longitude_degrees=self.lon_deg,
            elevation_m=self.elev_m,
        )

        for sat in satellites:
            try:
                sat_geocentric = sat.at(sky_times).position.au
                earth_pos = self.earth.at(sky_times).position.au
                obs_pos_geocentric = observer.at(sky_times).position.au - earth_pos

                # Topocentric vector (satellite minus geocentric observer)
                topocentric_km = (sat_geocentric - obs_pos_geocentric) * 149_597_870.7

                cart_rep = CartesianRepresentation(
                    x=topocentric_km[0] * u.km,
                    y=topocentric_km[1] * u.km,
                    z=topocentric_km[2] * u.km,
                )
                coord = SkyCoord(
                    cart_rep,
                    frame="icrs",
                    obstime=times,
                    location=self.location,
                )
                traj[sat.name] = coord
            except Exception:
                continue

        return traj

    # ----------------------------------------------------------------------
    def get_pulsars_catalog(self) -> Dict[str, SkyCoord]:
        """Load the ATNF pulsar catalog from a persistent cache.
        
        The long-lived cache is data/tle_repository/pulsars/atnf_catalog.json.
        Reuse the session cache when present. Catalog downloads are currently
        disabled; return an empty dictionary when no usable cache is available."""
        from pathlib import Path
        import json
        
        # Reuse the session cache when available
        if self._pulsar_cache is not None:
            print("[CACHE] Using pulsars catalog from current session")
            return self._pulsar_cache
        
        # Set the persistent cache path
        project_root = Path(__file__).parent.parent.parent
        pulsars_dir = project_root / "data" / "tle_repository" / "pulsars"
        pulsars_dir.mkdir(parents=True, exist_ok=True)
        cache_file = pulsars_dir / "atnf_catalog.json"
        
        # Try loading from the persistent cache
        if cache_file.exists():
            try:
                print(f"[PERSISTENT CACHE] Loading pulsars catalog from {cache_file.name}")
                with open(cache_file, 'r') as f:
                    data = json.load(f)
                
                # Reconstruct SkyCoord objects
                out: Dict[str, SkyCoord] = {}
                for name, coords in data.items():
                    out[name] = SkyCoord(ra=coords['ra'], dec=coords['dec'], 
                                        unit=(u.deg, u.deg))
                
                self._pulsar_cache = out
                print(f"[CACHE] {len(out)} pulsars loaded from persistent cache")
                return out
            except Exception as e:
                print(f"[WARN] Failed to load persistent cache: {e}")
        
        # Download if there is no cache (TEMPORARILY DISABLED)
        print("[PULSARS] Download temporarily disabled - returning empty catalog")
        print("[INFO] To enable: uncomment psrqpy import and download code")
        return {}
        
        # Code for future use (uncomment when psrqpy is available):
        # try:
        #     print("[DOWNLOAD] Downloading ATNF pulsars catalog (large file)...")
        #     from psrqpy import QueryATNF
        #     q = QueryATNF(params=["PSRJ", "RAJ", "DECJ"])
        #     table = q.table
        #     out: Dict[str, SkyCoord] = {}
        #     data_to_save = {}
        #     
        #     for row in table:
        #         name = row["PSRJ"]
        #         ra = row["RAJ"]
        #         dec = row["DECJ"]
        #         coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
        #         out[name] = coord
        #         data_to_save[name] = {'ra': coord.ra.deg, 'dec': coord.dec.deg}
        #     
        #     # Save the persistent cache
        #     with open(cache_file, 'w') as f:
        #         json.dump(data_to_save, f, indent=2)
        #     
        #     self._pulsar_cache = out
        #     print(f"[PERSISTENT CACHE] {len(out)} pulsars saved to {cache_file.name}")
        #     return out
        # except Exception as e:
        #     print(f"[ERROR] Failed to download pulsars: {e}")
        #     return {}

    # ----------------------------------------------------------------------
    def find_objects_in_fov(self,
                            times: Time,
                            pointing: SkyCoord,
                            beamwidth_deg: float = 12.0,
                            include_gnss: bool = True,
                            include_pulsars: bool = True) -> Dict[str, SkyCoord]:
        """Find objects visible within the telescope field of view over time.
        
        Compare object coordinates with the ICRS pointing series, accounting
        for the observation times. Include objects entering the beam at any time."""
        times = self._ensure_time(times)
        beam_radius = (beamwidth_deg / 2.0) * u.deg

        # Collect the basic objects
        objects: Dict[str, SkyCoord] = {}
        objects.update(self.get_solar_system_objects(times))
        if include_gnss:
            objects.update(self.get_gnss_trajectories(times))

        # Pulsars are static: replicate their positions
        if include_pulsars:
            pulsars = self.get_pulsars_catalog()
            for name, coord in pulsars.items():
                ra_arr = np.full(len(times), coord.ra.deg) * u.deg
                dec_arr = np.full(len(times), coord.dec.deg) * u.deg
                objects[name] = SkyCoord(ra=ra_arr, dec=dec_arr, obstime=times)

        # Generate the telescope pointing time series in ICRS.
        # Reuse pointing if it is already a time series with obstime and arrays;
        # otherwise, replicate the fixed pointing over time.
        try:
            if hasattr(pointing, 'obstime') and np.ndim(pointing.ra.deg) > 0 and len(pointing.ra.deg) == len(times):
                pointing_series_icrs = pointing.transform_to('icrs')
            else:
                pointing_series_icrs = SkyCoord(
                    ra=np.full(len(times), float(pointing.ra.deg)) * u.deg,
                    dec=np.full(len(times), float(pointing.dec.deg)) * u.deg,
                    obstime=times,
                    frame='icrs'
                )
        except Exception:
            # Simple fallback: replicate the fixed pointing
            pointing_series_icrs = SkyCoord(
                ra=np.full(len(times), float(getattr(pointing, 'ra', 0.0))) * u.deg,
                dec=np.full(len(times), float(getattr(pointing, 'dec', 0.0))) * u.deg,
                obstime=times,
                frame='icrs'
            )

        visible: Dict[str, SkyCoord] = {}
        for name, coord in objects.items():
            try:
                # Keep objects that enter the beam at any observation time
                sep_series = coord.separation(pointing_series_icrs)
                if np.any(sep_series <= beam_radius):
                    visible[name] = coord
            except Exception:
                continue

        return visible

    # ----------------------------------------------------------------------
    def categorize_objects(self, objects: Dict[str, SkyCoord]) -> Dict[str, Dict[str, SkyCoord]]:
        """Group astronomical objects by type for plotting.
        
        Args:
            objects: Dictionary mapping names to SkyCoord positions.
        
        Returns:
            Dictionary with planets, satellites, and pulsars categories."""
        categories = {
            'planets': {},
            'satellites': {},
            'pulsars': {}
        }
        
        for name, coord in objects.items():
            name_lower = name.lower()
            name_upper = name.upper()
            
            # Planets, Sun, and Moon
            if any(planet in name_lower for planet in [
                'sun', 'sol', 'moon', 'lua', 'mercury', 'venus', 
                'mars', 'jupiter', 'saturn', 'uranus', 'neptune'
            ]):
                categories['planets'][name] = coord
            
            # GNSS satellites
            elif (any(sat in name_upper for sat in ['GPS', 'GLONASS', 'GALILEO', 'BEIDOU']) 
                  or 'gnss' in name_lower):
                categories['satellites'][name] = coord
            
            # Pulsars and other objects
            else:
                categories['pulsars'][name] = coord
        
        return categories

    def get_object_display_info(self, obj_name: str) -> Dict[str, any]:
        """Return plotting color, marker, and size for an object.
        
        Args:
            obj_name: Object name.
        
        Returns:
            Dictionary with color, marker, and size."""
        name_lower = obj_name.lower()
        name_upper = obj_name.upper()
        
        # Sun
        if any(sun in name_lower for sun in ['sun', 'sol']):
            return {'color': 'gold', 'marker': 'o', 'size': 8}
        
        # Moon
        elif any(moon in name_lower for moon in ['moon', 'lua']):
            return {'color': 'silver', 'marker': 'o', 'size': 6}
        
        # Specific planets
        elif 'mercury' in name_lower:
            return {'color': 'orange', 'marker': 'o', 'size': 5}
        elif 'venus' in name_lower:
            return {'color': 'yellow', 'marker': 'o', 'size': 5}
        elif 'mars' in name_lower:
            return {'color': 'red', 'marker': 'o', 'size': 5}
        elif 'jupiter' in name_lower:
            return {'color': 'brown', 'marker': 'o', 'size': 5}
        elif 'saturn' in name_lower:
            return {'color': 'tan', 'marker': 'o', 'size': 5}
        elif 'uranus' in name_lower:
            return {'color': 'cyan', 'marker': 'o', 'size': 5}
        elif 'neptune' in name_lower:
            return {'color': 'blue', 'marker': 'o', 'size': 5}
        
        # GNSS satellites
        elif (any(sat in name_upper for sat in ['GPS', 'GLONASS', 'GALILEO', 'BEIDOU']) 
              or 'gnss' in name_lower):
            return {'color': 'green', 'marker': '^', 'size': 3}
        
        # Pulsars and others
        else:
            return {'color': 'purple', 'marker': '*', 'size': 4}

    def get_objects_statistics(self, objects: Dict[str, SkyCoord], 
                              pointing: SkyCoord, beamwidth_deg: float, 
                              times: Time) -> Dict[str, any]:
        """Compute statistics for detected objects.
        
        Args:
            objects: Detected objects.
            pointing: Telescope pointing direction.
            beamwidth_deg: Full beam width in degrees.
            times: Observation timestamps.
        
        Returns:
            Detailed pointing, observation, object count, and category statistics."""
        categories = self.categorize_objects(objects)
        
        stats = {
            'pointing': {
                'ra_deg': pointing.ra.deg,
                'dec_deg': pointing.dec.deg
            },
            'observation': {
                'beamwidth_deg': beamwidth_deg,
                'start_time': times[0].iso,
                'end_time': times[-1].iso,
                'duration_hours': (times[-1] - times[0]).to(u.hour).value
            },
            'objects': {
                'total': len(objects),
                'planets': len(categories['planets']),
                'satellites': len(categories['satellites']),
                'pulsars': len(categories['pulsars'])
            },
            'categories': categories
        }
        
        return stats

    def print_objects_summary(self, objects: Dict[str, SkyCoord], 
                             pointing: SkyCoord, beamwidth_deg: float, 
                             times: Time) -> None:
        """Print a detailed summary of detected objects."""
        stats = self.get_objects_statistics(objects, pointing, beamwidth_deg, times)
        
        print(f"\n{'='*50}")
        print("  OBJECTS IN THE FIELD OF VIEW - SUMMARY")
        print(f"{'='*50}")
        
        print(f"Pointing: RA={stats['pointing']['ra_deg']:.2f}°, "
              f"Dec={stats['pointing']['dec_deg']:.2f}°")
        print(f"Beam width: {stats['observation']['beamwidth_deg']}°")
        print(f"Period: {stats['observation']['start_time']} → {stats['observation']['end_time']}")
        print(f"Duration: {stats['observation']['duration_hours']:.1f} hours")
        print(f"Total objects: {stats['objects']['total']}")
        
        print("\nBy category:")
        if stats['objects']['planets'] > 0:
            planet_names = list(stats['categories']['planets'].keys())
            print(f"   Planets/Sun/Moon: {stats['objects']['planets']} ({', '.join(planet_names)})")
        
        if stats['objects']['satellites'] > 0:
            print(f"   GNSS satellites: {stats['objects']['satellites']} objects")
            
        if stats['objects']['pulsars'] > 0:
            print(f"   Pulsars: {stats['objects']['pulsars']} objects")
        
        print(f"{'='*50}\n")
