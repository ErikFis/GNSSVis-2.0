"""Offline geometry stage: scenario grids, local TLEs, SGP4, beam masks, passages."""

from datetime import datetime, timedelta, timezone

import numpy as np
import pytest

from conftest import REPO_ROOT, SCENARIO_PATH
from gnssvis.tle_utils import parse_tle_text, read_tle_file
from gnssvis_simulation import (
    ScenarioError,
    angular_separation_deg,
    beam_masks,
    compute_geometry,
    find_passages,
    geometry_passages,
    load_scenario,
    offline_sky_analyzer,
    verify_reference_hashes,
)

UTC = timezone.utc
PRN03 = 40294


def _snapshot(directory):
    return sorted(str(p.relative_to(directory)) for p in directory.rglob("*"))


# -- scenario grids -------------------------------------------------------

def test_reference_grids(scenario):
    assert scenario.matrix_shape == (512, 2160)
    assert scenario.channel_width_hz == 195312.5
    f = scenario.frequencies_hz
    assert f.dtype == np.float64
    assert f[0] == 1500097656.25 and f[-1] == 1599902343.75
    t = scenario.time_offsets_seconds
    assert t.dtype == np.int64 and t[0] == 0 and t[-1] == 21590
    assert np.all(np.diff(t) == 10)
    assert scenario.timestamps_utc[-1] == datetime(2026, 8, 22, 5, 59, 50, tzinfo=UTC)
    assert scenario.end_utc == datetime(2026, 8, 22, 6, 0, 0, tzinfo=UTC)
    assert scenario.repo_root == REPO_ROOT


def test_scenario_outside_repository_uses_package_checkout(tmp_path):
    copy = tmp_path / "my_scenario.json"
    copy.write_text(SCENARIO_PATH.read_text(encoding="utf-8"), encoding="utf-8")
    assert load_scenario(copy).repo_root == REPO_ROOT


def test_reference_paths_do_not_depend_on_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    scenario = load_scenario(SCENARIO_PATH)
    assert scenario.reference_path("tle_file") == REPO_ROOT / "reference_data" / "satellites.tle"


@pytest.mark.parametrize("mutate, message", [
    (lambda c: c["time_grid"].update(duration_seconds=21605), "exact multiple"),
    (lambda c: c["time_grid"].update(step_seconds=0), "positive"),
    (lambda c: c["time_grid"].update(step_seconds=2.5), "integer"),
    (lambda c: c["frequency_grid"].update(channel_count=0), "positive"),
    (lambda c: c["receiver"].update(band_max_hz=1.4e9), "band_max_hz"),
    (lambda c: c["time_grid"].update(start_utc="2026-08-22T00:00:00"), "end with 'Z'"),
    (lambda c: c["observatory"].update(earth_model="GRS80"), "WGS84"),
    (lambda c: c["references"].update(norad_ids=[40294, 40294]), "duplicates"),
])
def test_invalid_scenarios_are_rejected(write_scenario, mutate, message):
    with pytest.raises(ScenarioError, match=message):
        write_scenario(mutate)


def test_resolution_is_configurable(write_scenario, geometry):
    coarse = write_scenario(lambda c: c["time_grid"].update(step_seconds=60)
                            or c["frequency_grid"].update(channel_count=128))
    assert coarse.matrix_shape == (128, 360)
    assert coarse.channel_width_hz == 781250.0
    g60 = compute_geometry(coarse)
    # Coarse samples coincide with every 6th reference sample.
    np.testing.assert_allclose(g60.alt_deg, geometry.alt_deg[:, ::6], atol=1e-9)
    np.testing.assert_array_equal(g60.selected, geometry.selected[:, ::6])


# -- reference inputs ------------------------------------------------------

def test_reference_hashes_match_provenance(scenario):
    assert set(verify_reference_hashes(scenario)) == {"satellites.tle", "satellites.csv",
                                                      "signals.csv"}


def test_tle_checksums_are_validated(scenario):
    text = scenario.reference_path("tle_file").read_text(encoding="utf-8")
    assert len(parse_tle_text(text)) == 6
    lines = text.splitlines()
    # Change an inclination digit in the first line 2 without fixing its checksum.
    lines[2] = lines[2][:9] + ("5" if lines[2][9] != "5" else "6") + lines[2][10:]
    with pytest.raises(ValueError, match="checksum"):
        parse_tle_text("\n".join(lines))


def test_local_tles_follow_configured_order(scenario):
    sky = offline_sky_analyzer(scenario)
    ids = list(reversed(scenario.norad_ids))
    sats = sky.load_local_tles(scenario.reference_path("tle_file"), norad_ids=ids)
    assert [s.model.satnum for s in sats] == ids
    with pytest.raises(ValueError, match="not found"):
        sky.load_local_tles(scenario.reference_path("tle_file"), norad_ids=[12345])


def test_tle_epoch_limit_is_enforced(write_scenario):
    late = write_scenario(lambda c: c["time_grid"].update(start_utc="2026-09-01T00:00:00Z"))
    with pytest.raises(ScenarioError, match="TLE epoch"):
        compute_geometry(late)


# -- geometry --------------------------------------------------------------

def test_geometry_is_offline_and_writes_nothing(scenario):
    data_dir = REPO_ROOT / "data"
    before = _snapshot(data_dir)
    g = compute_geometry(scenario)  # the autouse fixture blocks sockets
    assert _snapshot(data_dir) == before
    assert g.alt_deg.shape == (6, 2160)
    assert g.norad_ids == scenario.norad_ids
    for arr in (g.alt_deg, g.az_deg, g.range_m, g.separation_deg):
        assert np.all(np.isfinite(arr))
    # GPS MEO: topocentric ranges between roughly 20,000 and 26,000+ km.
    assert 19.5e6 < g.range_m.min() and g.range_m.max() < 30.0e6


def test_zenith_separation_is_90_minus_altitude(geometry):
    np.testing.assert_allclose(geometry.separation_deg, 90.0 - geometry.alt_deg, atol=1e-9)


def test_reference_passages(scenario, geometry):
    """Only PRN 03 crosses the 15 degree zenith gate in this window."""
    passages = geometry_passages(geometry)
    assert [p.norad_id for p in passages] == [PRN03]
    p = passages[0]
    start = scenario.start_utc
    # The 05:01:50 entry sample lies only ~0.2 arcsec inside the gate, so a
    # different UT1-UTC table can legitimately move the entry by one sample.
    assert start + timedelta(seconds=p.entry_offset_seconds) in (
        datetime(2026, 8, 22, 5, 1, 50, tzinfo=UTC), datetime(2026, 8, 22, 5, 2, 0, tzinfo=UTC))
    assert start + timedelta(seconds=p.exit_offset_seconds) == datetime(2026, 8, 22, 5, 39, 20,
                                                                       tzinfo=UTC)
    assert p.duration_seconds == p.sample_count * 10
    assert not (p.entry_clipped or p.exit_clipped)
    row = geometry.norad_ids.index(PRN03)
    sel = geometry.selected[row]
    assert np.all(geometry.separation_deg[row, sel] <= 15.0)
    assert geometry.separation_deg[row, p.start_index - 1] > 15.0
    assert geometry.separation_deg[row, p.stop_index] > 15.0
    assert geometry.selected.sum() == p.sample_count


def test_independent_sgp4_topocentric_crosscheck(scenario, geometry):
    """Raw sgp4 + TEME->ECEF (IAU-82 GMST) + WGS84 ENU, with Skyfield's UT1."""
    from sgp4.api import Satrec

    sky = offline_sky_analyzer(scenario)
    from gnssvis_simulation.geometry import skyfield_times

    t = skyfield_times(scenario, sky)
    jd_utc = t.whole + (t.tai_fraction - t._leap_seconds() / 86400.0)
    ut1 = t.ut1 - 2451545.0
    tcen = ut1 / 36525.0
    gmst_s = (67310.54841 + (876600 * 3600 + 8640184.812866) * tcen
              + 0.093104 * tcen ** 2 - 6.2e-6 * tcen ** 3)
    gmst = np.radians((gmst_s % 86400.0) / 240.0)

    site = scenario.section("observatory")
    lat, lon, h = np.radians(site["latitude_deg"]), np.radians(site["longitude_deg"]), site["elevation_m"]
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    obs = np.array([(n + h) * np.cos(lat) * np.cos(lon), (n + h) * np.cos(lat) * np.sin(lon),
                    (n * (1 - e2) + h) * np.sin(lat)])
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])

    tles = {int(l1[2:7]): (l1, l2) for _, l1, l2 in read_tle_file(scenario.reference_path("tle_file"))}
    for row, norad in enumerate(geometry.norad_ids):
        err, r, _ = Satrec.twoline2rv(*tles[norad]).sgp4_array(np.floor(jd_utc), jd_utc - np.floor(jd_utc))
        assert not err.any()
        r = r.T * 1000.0
        c, s = np.cos(gmst), np.sin(gmst)
        d = np.vstack([c * r[0] + s * r[1], -s * r[0] + c * r[1], r[2]]) - obs[:, None]
        rng = np.linalg.norm(d, axis=0)
        alt = np.degrees(np.arcsin(up @ d / rng))
        assert np.max(np.abs(alt - geometry.alt_deg[row])) * 3600 < 0.05
        assert np.max(np.abs(rng - geometry.range_m[row])) < 1.0


# -- separation, masks and passages -----------------------------------------

def test_angular_separation():
    az = np.linspace(0, 360, 13)
    np.testing.assert_allclose(angular_separation_deg(75.0, az, 90.0, 0.0), 15.0, atol=1e-12)
    assert angular_separation_deg(0.0, 0.0, 0.0, 90.0) == pytest.approx(90.0)
    assert angular_separation_deg(10.0, 350.0, 10.0, 10.0) == pytest.approx(
        angular_separation_deg(10.0, 10.0, 10.0, 350.0))
    assert angular_separation_deg(45.0, 123.0, 45.0, 123.0) == pytest.approx(0.0, abs=1e-12)


def test_beam_boundary_is_inclusive():
    sep = np.array([14.999, 15.0, 15.001, 5.0])
    alt = np.array([75.001, 75.0, 74.999, -1.0])
    above, in_beam, selected = beam_masks(sep, alt, 30.0, 0.0)
    np.testing.assert_array_equal(in_beam, [True, True, False, True])
    np.testing.assert_array_equal(above, [True, True, True, False])
    np.testing.assert_array_equal(selected, [True, True, False, False])
    assert beam_masks(0.0, 0.0, 30.0, 0.0)[0]  # altitude 0 counts as visible


@pytest.mark.parametrize("mask, expected", [
    ([], []),
    ([0, 0, 0, 0], []),
    ([1, 1, 1, 1], [(0, 4, True, True)]),
    ([0, 1, 1, 0, 0, 1, 0], [(1, 3, False, False), (5, 6, False, False)]),
    ([1, 0, 0, 1], [(0, 1, True, False), (3, 4, False, True)]),
    ([0, 0, 1, 0], [(2, 3, False, False)]),
])
def test_find_passages(mask, expected):
    passages = find_passages(np.array(mask, dtype=bool), 10, norad_id=1)
    assert [(p.start_index, p.stop_index, p.entry_clipped, p.exit_clipped)
            for p in passages] == expected
    for p in passages:
        assert p.entry_offset_seconds == 10 * p.start_index
        assert p.exit_offset_seconds == 10 * p.stop_index
        assert p.duration_seconds == 10 * p.sample_count == p.exit_offset_seconds - p.entry_offset_seconds


def test_no_passage_case(write_scenario):
    # Narrow beam aimed at the south horizon: nothing passes.
    s = write_scenario(lambda c: c["pointing"].update(altitude_deg=5.0, azimuth_deg=180.0)
                       or c["beam"].update(full_width_deg=1.0))
    g = compute_geometry(s)
    assert not g.selected.any()
    assert geometry_passages(g) == []
