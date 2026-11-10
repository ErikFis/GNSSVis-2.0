"""The existing SkyAnalyzer interface, used offline with the local TLEs."""

import inspect

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time

from conftest import REPO_ROOT
from gnssvis import SkyAnalyzer

TLE = REPO_ROOT / "reference_data" / "satellites.tle"


@pytest.fixture
def times():
    # 04:55 to 05:45 UTC, 60 s step: spans the PRN 03 zenith passage.
    return Time("2026-08-22T04:55:00", scale="utc") + np.arange(51) * 60 * u.s


def test_public_signatures_are_backward_compatible():
    init = inspect.signature(SkyAnalyzer.__init__).parameters
    assert list(init)[:4] == ["self", "lat_deg", "lon_deg", "elev_m"]
    assert init["offline"].default is False and init["offline"].kind is inspect.Parameter.KEYWORD_ONLY
    traj = inspect.signature(SkyAnalyzer.get_gnss_trajectories).parameters
    assert list(traj)[:3] == ["self", "times", "tle_url"]
    assert traj["tle_file"].default is None


def test_offline_constructor_loads_no_ephemeris():
    sky = SkyAnalyzer(-7.211637, -35.908138, 553.0, offline=True)
    assert sky.planets is None and sky.loader is None
    with pytest.warns(UserWarning, match="DE421"):
        assert sky.get_solar_system_objects(Time("2026-08-22T00:00:00")) == {}


def test_gnss_trajectories_from_local_tle(times):
    sky = SkyAnalyzer(offline=True)
    traj = sky.get_gnss_trajectories(times, tle_file=TLE)
    assert len(traj) == 6
    coord = traj["GPS BIIF-8  (PRN 03)"]
    assert coord.shape == (51,)
    # Topocentric distance matches the alt/az path.
    sat = sky.load_local_tles(TLE, norad_ids=[40294])[0]
    altaz = sky.satellite_altaz(sat, times)
    np.testing.assert_allclose(coord.distance.to_value(u.m), altaz["range_m"], rtol=1e-12)


def test_find_objects_in_fov_offline_matches_altaz_path(times):
    """Legacy RA/Dec selection agrees with the new alt/az separation."""
    sky = SkyAnalyzer(offline=True, tle_file=TLE)
    ts_times = sky.ts.from_astropy(times)
    zenith = sky.topos.at(ts_times).from_altaz(alt_degrees=90.0, az_degrees=0.0)
    ra, dec, _ = zenith.radec()
    pointing = SkyCoord(ra=ra.degrees * u.deg, dec=dec.degrees * u.deg, obstime=times, frame="icrs")

    with pytest.warns(UserWarning, match="DE421"):
        visible = sky.find_objects_in_fov(times, pointing, beamwidth_deg=30.0,
                                          include_gnss=True, include_pulsars=False)
    assert list(visible) == ["GPS BIIF-8  (PRN 03)"]

    sat = sky.load_local_tles(TLE, norad_ids=[40294])[0]
    alt = sky.satellite_altaz(sat, times)["alt_deg"]
    legacy_sep = visible["GPS BIIF-8  (PRN 03)"].separation(pointing).deg
    np.testing.assert_allclose(legacy_sep, 90.0 - alt, atol=1e-6)
