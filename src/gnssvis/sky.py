"""Ephemerides and astronomical objects with topocentric projection.

Support DE421 solar system ephemerides, GNSS satellites from Celestrak
TLEs, and an ATNF pulsar catalog cache. Coordinates vary with time.

Author: Erik H L Guimaraes"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Union
import warnings
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, CartesianRepresentation
from skyfield.api import EarthSatellite, Topos, wgs84
# from psrqpy import QueryATNF  # Temporarily disabled

warnings.filterwarnings("ignore", message="you can only subtract two vectors")


class SkyAnalyzer:
    """Calculate ephemerides and find objects in the field of view."""

    def __init__(self, lat_deg: float = -7.211637, lon_deg: float = -35.908138, elev_m: float = 553.0,
                 *, offline: bool = False,
                 tle_file: Optional[Union[str, Path]] = None) -> None:
        """Initialize the observatory and load planetary ephemerides.

        Args:
            lat_deg: Geodetic latitude in degrees (negative south).
            lon_deg: Geodetic longitude in degrees (east positive).
            elev_m: Height above the WGS84 ellipsoid in meters.
            offline: If True, never access the network or write caches: use
                Skyfield's bundled timescale data and do not load DE421.
                Solar system objects are then unavailable; GNSS geometry
                only needs SGP4 and the observer position.
            tle_file: Optional local three-line TLE file. When given, GNSS
                trajectories are propagated from it instead of Celestrak.
        """
        from skyfield.api import Loader

        self.lat_deg: float = lat_deg
        self.lon_deg: float = lon_deg
        self.elev_m: float = elev_m
        self.offline: bool = bool(offline)
        self.tle_file: Optional[Path] = Path(tle_file) if tle_file is not None else None

        # Geocentric WGS84 observer used for topocentric satellite vectors.
        # It does not require a planetary ephemeris.
        self.topos = wgs84.latlon(self.lat_deg, self.lon_deg, elevation_m=self.elev_m)

        # Astropy: location for AltAz transformations
        self.location = EarthLocation(lat=self.lat_deg * u.deg,
                                      lon=self.lon_deg * u.deg,
                                      height=self.elev_m * u.m)

        # Cache to avoid repeated downloads within a session
        self._pulsar_cache = None
        self._gnss_satellites_cache = None

        if self.offline:
            # Bundled Delta T / UT1 / leap-second tables only: no file I/O.
            from skyfield.api import load as _skyfield_load
            self.loader = None
            self.ts = _skyfield_load.timescale(builtin=True)
            self.planets = None
            self.earth = None
            self.observatory = None
            return

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
        if self.planets is None:
            warnings.warn("Solar system objects need the DE421 ephemeris, which is "
                          "not loaded in offline mode; returning no objects.")
            return {}
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
        tle_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, SkyCoord]:
        """Calculate GNSS trajectories with Skyfield SGP4 propagation.
        
        Return observer-relative topocentric RA/Dec positions. Download TLEs
        from the Celestrak gp.php API into data/tle_repository/gnss/. The session
        cache avoids repeated downloads during one run. Local files older than
        seven days require an update.
        
        If ``tle_file`` (or the ``tle_file`` given to the constructor) is set,
        propagate those local elements instead, without network access or
        cache writes. In offline mode without a local file, only an existing
        repository file is used and nothing is downloaded.

        Args:
            times: Observation timestamps as astropy.Time.
            tle_url: Celestrak API URL for GNSS TLEs.
            tle_file: Optional local three-line TLE file.
        
        Returns:
            Dictionary mapping satellite names to SkyCoord positions at each time."""
        from datetime import datetime
        import shutil

        times = self._ensure_time(times)
        sky_times = self._to_skyfield_time(times)

        local_file = tle_file if tle_file is not None else self.tle_file
        project_root = Path(__file__).parent.parent.parent
        repo_file = project_root / "data" / "tle_repository" / "gnss" / "gnss.txt"
        if local_file is None and self.offline:
            if not repo_file.exists():
                print("[WARN] Offline mode: no local GNSS TLE file available.")
                return {}
            local_file = repo_file

        if local_file is not None:
            try:
                satellites = self.load_local_tles(local_file)
            except Exception as exc:
                print(f"[WARN] Failed to load local GNSS TLEs: {exc}")
                return {}
            print(f"[LOCAL TLE] {len(satellites)} satellites loaded from {Path(local_file).name}.")
            return self._topocentric_skycoords(satellites, sky_times, times)

        # TLE repository paths
        tle_repo = repo_file.parent
        tle_repo.mkdir(parents=True, exist_ok=True)
        tle_file = repo_file

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
            elif not needs_update:
                # A fresh repository file exists: reuse it instead of downloading.
                satellites = self.load_local_tles(tle_file, validate_checksum=False)
                self._gnss_satellites_cache = satellites
                print(f"[REPOSITORY] {len(satellites)} GNSS satellites loaded from {tle_file.name}.")
            else:
                import requests

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

        return self._topocentric_skycoords(satellites, sky_times, times)

    # ----------------------------------------------------------------------
    def load_local_tles(self, tle_file: Union[str, Path],
                        norad_ids: Optional[Iterable[int]] = None,
                        validate_checksum: bool = True) -> List[EarthSatellite]:
        """Load EarthSatellite objects from a local TLE file (no network).

        Args:
            tle_file: Three-line TLE file (name, line 1, line 2).
            norad_ids: Optional NORAD IDs to select; the result follows this
                order and a missing or duplicated ID raises ValueError.
            validate_checksum: Verify line checksums and catalog numbers.

        Returns:
            List of Skyfield EarthSatellite objects bound to this timescale.
        """
        from .tle_utils import read_tle_file

        records = read_tle_file(tle_file, validate_checksum=validate_checksum)
        satellites = [EarthSatellite(l1, l2, name, self.ts) for name, l1, l2 in records]
        if norad_ids is None:
            return satellites
        by_id: Dict[int, EarthSatellite] = {}
        for sat in satellites:
            if sat.model.satnum in by_id:
                raise ValueError(f"Duplicate NORAD ID {sat.model.satnum} in {tle_file}")
            by_id[sat.model.satnum] = sat
        missing = [int(n) for n in norad_ids if int(n) not in by_id]
        if missing:
            raise ValueError(f"NORAD IDs not found in {tle_file}: {missing}")
        return [by_id[int(n)] for n in norad_ids]

    def satellite_altaz(self, satellite: EarthSatellite, times) -> Dict[str, np.ndarray]:
        """Unrefracted topocentric altitude, azimuth and range of a satellite.

        Use SGP4 propagation and the WGS84 observer, without a planetary
        ephemeris or atmospheric refraction.

        Args:
            satellite: Skyfield EarthSatellite.
            times: astropy.Time, ISO strings, or a Skyfield Time.

        Returns:
            Dictionary with ``alt_deg``, ``az_deg`` (east of north) and
            ``range_m`` arrays, one value per timestamp.

        Raises:
            ValueError: If SGP4 reports an error or a non-finite position.
        """
        sky_times = times if hasattr(times, "tt") and hasattr(times, "ts") else \
            self._to_skyfield_time(self._ensure_time(times))
        topocentric = (satellite - self.topos).at(sky_times)
        position = topocentric.position.m
        if not np.all(np.isfinite(position)):
            message = getattr(topocentric, "message", "non-finite position")
            raise ValueError(f"SGP4 failed for {satellite.name}: {message}")
        alt, az, distance = topocentric.altaz()
        # Positive radial velocity means recession (negative received Doppler).
        range_rate = np.sum(position * topocentric.velocity.m_per_s, axis=0) / distance.m
        return {
            "alt_deg": np.atleast_1d(alt.degrees),
            "az_deg": np.atleast_1d(az.degrees),
            "range_m": np.atleast_1d(distance.m),
            "range_rate_m_s": np.atleast_1d(range_rate),
        }

    def _topocentric_skycoords(self, satellites, sky_times, times: Time) -> Dict[str, SkyCoord]:
        """Topocentric RA/Dec SkyCoords (legacy GNSSVis trajectory format)."""
        traj: Dict[str, SkyCoord] = {}
        for sat in satellites:
            try:
                # Satellite minus geocentric observer, in GCRS. Equivalent to
                # the former DE421-based difference, in which Earth cancels.
                topocentric_km = (sat - self.topos).at(sky_times).position.km

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
