"""Satellite emission: scaling, frequency placement, passage mask, overlaps.

Run with ``pytest -s`` to see the printed reference contribution.
"""

import numpy as np
import pytest

from conftest import REPO_ROOT
from gnssvis_simulation import ScenarioError, beam_masks, compute_geometry
from gnssvis_simulation.emission import (Carrier, GaussianCarrierProfile, angular_gain,
                                         carrier_in_band, compute_signal_power,
                                         load_signal_catalog, plan_contributions,
                                         profiles_from_scenario)
from gnssvis_simulation.geometry import GeometryResult, SatelliteRecord
from gnssvis_simulation.noise import generate_noise
from gnssvis_simulation.products import array_sha256

L1, L2, L5 = 1575.42e6, 1227.6e6, 1176.45e6
NOISE_SEED42_SHA = "60f04c46fbc48b35309b379c31f83a1c2cd88bbcdda5f71c0d3bc7af1f39a766"


# -- helpers ----------------------------------------------------------------

@pytest.fixture
def coarse(write_scenario):
    """Reference band/channels, 36 samples of 600 s."""
    return write_scenario(lambda c: c["time_grid"].update(step_seconds=600))


def synthetic_geometry(scenario, range_m, separation_deg, alt_deg=None, blocks=("IIF",)):
    """GeometryResult with prescribed arrays; masks use the production rule."""
    n_sat, n_t = len(blocks), scenario.n_times
    shape = (n_sat, n_t)
    rng = np.broadcast_to(np.asarray(range_m, dtype=float), shape).copy()
    sep = np.broadcast_to(np.asarray(separation_deg, dtype=float), shape).copy()
    alt = (90.0 - sep) if alt_deg is None else np.broadcast_to(np.asarray(alt_deg, float), shape).copy()
    beam = scenario.section("beam")
    above, in_beam, selected = beam_masks(sep, alt, float(beam["full_width_deg"]),
                                          float(beam["minimum_altitude_deg"]))
    sats = [SatelliteRecord(norad_id=90000 + k, name=f"TEST {k}", constellation="GPS", block=b,
                            prn=k, tle_epoch_utc=scenario.start_utc, satellite=None)
            for k, b in enumerate(blocks)]
    return GeometryResult(satellites=sats, time_offsets_seconds=scenario.time_offsets_seconds,
                          start_utc=scenario.start_utc, step_seconds=scenario.step_seconds,
                          alt_deg=alt, az_deg=np.zeros(shape), range_m=rng, separation_deg=sep,
                          above_horizon=above, in_beam=in_beam, selected=selected)


def on_bin_carrier(scenario, index, signal_id="TEST"):
    """A carrier exactly on a channel centre, so the spectral weight there is 1."""
    return Carrier(signal_id, "GPS", float(scenario.frequencies_hz[index]), ("IIF", "III"))


# -- reference scenario -------------------------------------------------------

@pytest.fixture(scope="module")
def reference_signal(scenario, geometry):
    return compute_signal_power(scenario, geometry)


def test_format_matches_noise_grid(scenario, reference_signal):
    signal, _ = reference_signal
    assert signal.dtype == np.float32 and signal.shape == scenario.matrix_shape == (512, 2160)
    assert np.all(np.isfinite(signal)) and signal.min() == 0.0


def test_reference_contribution_follows_geometry(scenario, geometry, reference_signal):
    """Only PRN 03 at L1 -- derived from masks, catalog and band, not hardcoded."""
    signal, terms = reference_signal
    active = [(t.norad_id, t.signal_id) for t in terms if t.active]
    expected_active = [(n, "GPS_L1") for n, row in zip(geometry.norad_ids, geometry.selected)
                       if row.any()]
    assert active == expected_active == [(40294, "GPS_L1")]
    row = geometry.norad_ids.index(40294)
    np.testing.assert_array_equal(signal.any(axis=0), geometry.selected[row])
    assert not signal[:, ~geometry.selected.any(axis=0)].any()
    support = np.abs(scenario.frequencies_hz - L1) <= 1.5e6
    np.testing.assert_array_equal(signal.any(axis=1), support)
    i, j = np.unravel_index(np.argmax(signal), signal.shape)
    print(f"\nreference peak {signal[i, j]:.4f} at channel {i} "
          f"({scenario.frequencies_hz[i] / 1e6:.6f} MHz), {scenario.timestamps_utc[j]:%H:%M:%S} UTC; "
          f"{np.count_nonzero(signal)} nonzero cells")


def test_reference_peak_equals_formula(scenario, geometry, reference_signal):
    signal, _ = reference_signal
    row = geometry.norad_ids.index(40294)
    j = np.flatnonzero(geometry.selected[row])
    f = scenario.frequencies_hz[:, None]
    expected = (20.0 * (2e7 / geometry.range_m[row, j]) ** 2
                * np.exp(-4 * np.log(2) * (geometry.separation_deg[row, j] / 30.0) ** 2)
                * np.where(np.abs(f - L1) <= 1.5e6, np.exp(-0.5 * ((f - L1) / 5e5) ** 2), 0.0))
    np.testing.assert_allclose(signal[:, j], expected, rtol=1e-6, atol=1e-6)


def test_selection_is_not_hardcoded(write_scenario, geometry):
    """Point the beam at another satellite: that satellite now contributes."""
    row, j = geometry.norad_ids.index(39166), 500
    s = write_scenario(lambda c: c["pointing"].update(
        altitude_deg=float(geometry.alt_deg[row, j]), azimuth_deg=float(geometry.az_deg[row, j])))
    g = compute_geometry(s)
    signal, terms = compute_signal_power(s, g)
    active = {t.norad_id for t in terms if t.active}
    assert 39166 in active
    np.testing.assert_array_equal(signal.any(axis=0), g.selected.any(axis=0))
    assert signal[:, j].max() > 0


def test_noise_unchanged_and_emission_has_no_randomness(scenario, geometry):
    np.random.seed(3)
    expected = np.random.random(2)
    np.random.seed(3)
    compute_signal_power(scenario, geometry)
    assert np.array_equal(np.random.random(2), expected)
    assert array_sha256(generate_noise(scenario)) == NOISE_SEED42_SHA


# -- scaling ------------------------------------------------------------------

def test_inverse_square_range_scaling(coarse):
    k = 200
    ranges = np.full(coarse.n_times, 2e7)
    ranges[:4] = [2e7, 4e7, 1e7, 3e7]
    g = synthetic_geometry(coarse, ranges, 0.0)
    signal, _ = compute_signal_power(coarse, g, carriers=[on_bin_carrier(coarse, k)])
    np.testing.assert_allclose(signal[k, :4], [20.0, 5.0, 80.0, 20.0 / 2.25], rtol=1e-6)


def test_angular_gain_profile(coarse):
    np.testing.assert_allclose(angular_gain([0.0, 7.5, 15.0, 30.0], 30.0),
                               [1.0, 2 ** -0.25, 0.5, 2 ** -4], rtol=1e-12)
    k = 200
    sep = np.zeros(coarse.n_times)
    sep[:4] = [0.0, 7.5, 15.0, 15.000001]  # last one is just outside the gate
    g = synthetic_geometry(coarse, 2e7, sep)
    signal, _ = compute_signal_power(coarse, g, carriers=[on_bin_carrier(coarse, k)])
    np.testing.assert_allclose(signal[k, :4], [20.0, 20.0 * 2 ** -0.25, 10.0, 0.0], rtol=1e-6)


def test_coefficients_come_from_configuration(write_scenario, coarse):
    def mutate(c):
        c["time_grid"].update(step_seconds=600)
        c["emission"]["profiles"]["IIF"].update(peak_power_relative_at_reference_range=40.0,
                                                spectral_sigma_hz=250000.0)
    doubled = write_scenario(mutate)
    carriers = [on_bin_carrier(coarse, 200)]
    base, _ = compute_signal_power(coarse, synthetic_geometry(coarse, 2e7, 0.0), carriers)
    new, _ = compute_signal_power(doubled, synthetic_geometry(doubled, 2e7, 0.0), carriers)
    assert new[200].max() == pytest.approx(2 * base[200].max())
    assert np.count_nonzero(new[:, 0]) < np.count_nonzero(base[:, 0])  # narrower sigma


# -- frequency placement --------------------------------------------------------

def test_profile_evaluated_at_actual_channel_centres(scenario):
    prof = profiles_from_scenario(scenario)["IIF"]
    f = scenario.frequencies_hz
    w = prof.spectral_shape(f, L1)
    i = int(np.argmax(w))
    assert i == int(np.argmin(np.abs(f - L1))) == 386
    # Not snapped to the nearest bin: channel 386 is 68.28 kHz from L1.
    assert w[i] == pytest.approx(np.exp(-0.5 * ((f[i] - L1) / 5e5) ** 2)) and w[i] < 1.0
    support = np.flatnonzero(w > 0)
    assert support.tolist() == list(range(378, 394))
    assert np.all(np.abs(f[support] - L1) <= 1.5e6)
    centroid = np.sum(w * f) / np.sum(w)
    assert abs(centroid - L1) < scenario.channel_width_hz / 2


def test_spectral_cutoff_is_inclusive(scenario):
    prof = GaussianCarrierProfile(20.0, 2e7, 5e5, 3.0)
    f = scenario.frequencies_hz
    carrier = f[100] + 1.5e6  # channel 100 sits exactly at 3 sigma
    w = prof.spectral_shape(f, carrier)
    assert w[100] == pytest.approx(np.exp(-4.5))
    assert w[99] == 0.0
    assert prof.support(f, carrier).sum() == np.count_nonzero(w)


def test_carriers_placed_independently_in_wide_band(write_scenario):
    """A band containing L1, L2 and L5: each appears at its own frequency."""
    s = write_scenario(lambda c: c["receiver"].update(band_min_hz=1.1e9, band_max_hz=1.7e9)
                       or c["frequency_grid"].update(channel_count=3072)
                       or c["time_grid"].update(step_seconds=600))
    g = synthetic_geometry(s, 2e7, 0.0)
    signal, terms = compute_signal_power(s, g)
    assert sorted(t.signal_id for t in terms) == ["GPS_L1", "GPS_L2", "GPS_L5"]
    f = s.frequencies_hz
    for carrier in (L1, L2, L5):
        near = np.abs(f - carrier) <= 1.5e6
        assert signal[near, 0].max() > 19.0
        assert int(np.argmax(np.where(near, signal[:, 0], 0))) == int(np.argmin(np.abs(f - carrier)))
    far = np.all([np.abs(f - c) > 1.5e6 for c in (L1, L2, L5)], axis=0)
    assert not signal[far].any()


# -- out of band ----------------------------------------------------------------

def test_out_of_band_carriers_are_excluded(scenario, reference_signal):
    _, terms = reference_signal
    catalog = {c.signal_id: c for c in load_signal_catalog(scenario)}
    assert not carrier_in_band(catalog["GPS_L2"].carrier_hz, scenario)
    assert not carrier_in_band(catalog["GPS_L5"].carrier_hz, scenario)
    assert {t.signal_id for t in terms} == {"GPS_L1"}


def test_band_edges_are_inclusive(scenario):
    assert carrier_in_band(scenario.band_min_hz, scenario)
    assert carrier_in_band(scenario.band_max_hz, scenario)
    assert not carrier_in_band(scenario.band_max_hz + 1.0, scenario)


def test_band_without_l1_gives_zero_signal(write_scenario):
    s = write_scenario(lambda c: c["receiver"].update(band_min_hz=1.3e9, band_max_hz=1.4e9))
    signal, terms = compute_signal_power(s, compute_geometry(s))
    assert terms == [] and not signal.any()


# -- passage mask -----------------------------------------------------------------

def test_zero_outside_beam_and_below_horizon(coarse):
    n = coarse.n_times
    sep = np.full(n, 5.0)
    alt = np.full(n, 85.0)
    sep[10:15] = 20.0  # outside the 15 deg gate
    alt[20:25] = -0.5  # inside the gate by separation but below the horizon
    g = synthetic_geometry(coarse, 2e7, sep, alt)
    signal, _ = compute_signal_power(coarse, g, carriers=[on_bin_carrier(coarse, 200)])
    assert not signal[:, 10:15].any() and not signal[:, 20:25].any()
    np.testing.assert_array_equal(signal.any(axis=0), g.selected[0])
    assert g.selected[0].sum() == n - 10


def test_each_contribution_is_masked(coarse, geometry, scenario):
    """The per-term envelope/power (used by reports) is zero outside the mask too."""
    sep = np.full(coarse.n_times, 5.0)
    sep[::3] = 20.0
    g = synthetic_geometry(coarse, 2e7, sep)
    for scen, geo in ((coarse, g), (scenario, geometry)):
        for term in plan_contributions(scen, geo):
            assert not term.envelope[~term.mask].any()
            assert not term.power()[:, ~term.mask].any()
            assert np.all(term.envelope[term.mask] > 0)


def test_no_passage_no_signal(coarse):
    g = synthetic_geometry(coarse, 2e7, 40.0)
    signal, terms = compute_signal_power(coarse, g)
    assert not signal.any() and terms and not any(t.active for t in terms)


# -- overlapping contributions -------------------------------------------------------

def test_two_satellites_on_same_carrier_add(coarse):
    ranges = np.stack([np.full(coarse.n_times, 2e7), np.full(coarse.n_times, 2.5e7)])
    seps = np.stack([np.full(coarse.n_times, 3.0), np.full(coarse.n_times, 9.0)])
    g = synthetic_geometry(coarse, ranges, seps, blocks=("IIF", "III"))
    both, terms = compute_signal_power(coarse, g)
    singles = []
    for k in range(2):
        gk = synthetic_geometry(coarse, ranges[k], seps[k], blocks=(g.satellites[k].block,))
        singles.append(compute_signal_power(coarse, gk)[0].astype(np.float64))
    assert len([t for t in terms if t.active]) == 2
    np.testing.assert_allclose(both, singles[0] + singles[1], rtol=1e-6)
    assert both.max() > max(s.max() for s in singles)  # summed, not overwritten


def test_spectrally_overlapping_carriers_add(coarse):
    c1, c2 = on_bin_carrier(coarse, 200, "A"), on_bin_carrier(coarse, 203, "B")  # 0.59 MHz apart
    g = synthetic_geometry(coarse, 2e7, 0.0)
    both, _ = compute_signal_power(coarse, g, carriers=[c1, c2])
    a, _ = compute_signal_power(coarse, g, carriers=[c1])
    b, _ = compute_signal_power(coarse, g, carriers=[c2])
    overlap = (a[:, 0] > 0) & (b[:, 0] > 0)
    assert overlap.sum() >= 10
    np.testing.assert_allclose(both, a.astype(np.float64) + b, rtol=1e-6)
    assert np.all(both[overlap, 0] > np.maximum(a[overlap, 0], b[overlap, 0]))


def test_applicability_by_constellation_and_block(coarse):
    g = synthetic_geometry(coarse, 2e7, 0.0, blocks=("IIF",))
    other = [Carrier("GAL_E1", "Galileo", L1, ("FOC",)), Carrier("GPS_III_ONLY", "GPS", L1, ("III",))]
    signal, terms = compute_signal_power(coarse, g, carriers=other)
    assert terms == [] and not signal.any()


# -- configuration errors ---------------------------------------------------------

def test_missing_block_profile_is_an_error(write_scenario):
    s = write_scenario(lambda c: c["emission"]["profiles"].pop("IIF"))
    with pytest.raises(ScenarioError, match="No emission profile for block 'IIF'"):
        plan_contributions(s, compute_geometry(s))


@pytest.mark.parametrize("mutate, message", [
    (lambda c: c["emission"].update(model="calibrated_eirp"), "Unsupported emission model"),
    (lambda c: c["emission"].update(profile_key="satellite_catalog.prn"), "profile_key"),
    (lambda c: c["emission"]["profiles"]["III"].update(spectral_sigma_hz=0), "sigma_hz"),
])
def test_invalid_emission_configuration(write_scenario, mutate, message):
    with pytest.raises(ScenarioError, match=message):
        profiles_from_scenario(write_scenario(mutate))
