"""Physical invariants for the optional wideband model (no network)."""

import numpy as np
import pytest
from gnssvis_simulation.realistic import (
    doppler_shift_hz,
    bpsk_channel_psd,
    beam_offsets_deg,
    beam_gain,
    solar_ephemeris,
    validate_realistic,
)
from gnssvis_simulation import load_scenario, ScenarioError

BEAM = {
    "full_width_deg": 50.0,
    "minor_width_deg": 40.0,
    "reference_frequency_hz": 1575.42e6,
}


def test_recession_and_approach_have_correct_doppler_sign():
    assert doppler_shift_hz(1575.42e6, 1000) == pytest.approx(-5255.03546857)
    assert doppler_shift_hz(1575.42e6, -1000) > 0
    assert doppler_shift_hz(1575.42e6, 0) == 0


def test_beam_half_power_widths_and_chromaticity():
    assert beam_gain(0, 0, 1575.42e6, BEAM) == 1
    assert beam_gain(25, 0, 1575.42e6, BEAM) == pytest.approx(0.5)
    assert beam_gain(0, 20, 1575.42e6, BEAM) == pytest.approx(0.5)
    assert beam_gain(25 * 1575.42 / 1176.45, 0, 1176.45e6, BEAM) == pytest.approx(0.5)
    assert 0 < beam_gain(35, 0, 1575.42e6, BEAM) < 0.5  # no hard truncation


def test_spherical_offsets_wrap_and_have_correct_norm():
    x, y = beam_offsets_deg(np.array([50.0, 50.0]), np.array([359.0, 1.0]), 50.0, 0.0)
    assert x[0] == pytest.approx(-x[1])
    assert y[0] == pytest.approx(y[1])
    from gnssvis_simulation.geometry import angular_separation_deg

    assert np.hypot(x, y) == pytest.approx(
        angular_separation_deg(50.0, [359.0, 1.0], 50.0, 0.0)
    )
    assert np.hypot(*beam_offsets_deg(90.0, 0.0, 90.0, 100.0)) < 1e-10


def test_channel_integration_is_conservative_under_rebinning():
    fine = np.arange(-20e6, 20e6, 10000) + 5000
    coarse = fine.reshape(-1, 10).mean(axis=1)
    for shift in [-4000.0, 0.0, 4000.0]:
        psd = bpsk_channel_psd(fine, 10000, shift, 1.023e6)
        rebinned = bpsk_channel_psd(coarse, 100000, shift, 1.023e6)
        np.testing.assert_allclose(
            psd.reshape(-1, 10).mean(axis=1), rebinned, rtol=1e-10, atol=1e-12
        )
        assert np.all(psd >= 0)


def test_bpsk_power_and_lobe_nulls():
    f = np.arange(-200e6, 200e6, 10000) + 5000
    p = bpsk_channel_psd(f, 10000, 0, 1.023e6)
    # Integral of unit-peak sinc² is Rc; finite window loses only its tails.
    assert p.sum() * 10000 == pytest.approx(1.023e6, rel=0.001)
    null = bpsk_channel_psd(np.array([0, 1.023e6]), 1000, 0, 1.023e6)
    assert null[0] > 0.999
    assert null[1] < 1e-6


def test_doppler_channel_psd_moves_peak_without_changing_total():
    f = np.arange(-10e6, 10e6, 1000) + 500
    p = bpsk_channel_psd(f, 1000, 4000, 1.023e6)
    assert abs(f[p.argmax()] - 4000) <= 500
    assert p.sum() == pytest.approx(
        bpsk_channel_psd(f, 1000, -4000, 1.023e6).sum(), rel=1e-12
    )


def test_radial_state_velocity_matches_range_derivative(scenario):
    from gnssvis_simulation.geometry import offline_sky_analyzer, skyfield_times

    sky = offline_sky_analyzer(scenario)
    sat = sky.load_local_tles(scenario.reference_path("tle_file"))[0]
    t = skyfield_times(scenario, sky)[100]
    state = sky.satellite_altaz(sat, t)
    before = sky.satellite_altaz(sat, t - 1 / 86400)["range_m"]
    after = sky.satellite_altaz(sat, t + 1 / 86400)["range_m"]
    np.testing.assert_allclose(state["range_rate_m_s"], (after - before) / 2, atol=0.04)


def test_explicit_missing_ephemeris_never_downloads(tmp_path, scenario):
    from dataclasses import replace

    c = dict(scenario.config)
    c["solar"] = {"allow_download": True, "ephemeris_file": "missing.bsp"}
    sc = replace(scenario, config=c)
    with pytest.raises(FileNotFoundError):
        solar_ephemeris(sc, tmp_path / "missing.bsp")


def test_invalid_bandwidth_is_rejected():
    with pytest.raises(ValueError):
        bpsk_channel_psd([0], 0, 0, 1)
    with pytest.raises(ValueError):
        bpsk_channel_psd([0], 1, 0, -1)


def test_component_synthesis_is_reproducible_and_iir_has_no_l5(
    scenario, geometry, monkeypatch
):
    from dataclasses import replace
    from copy import deepcopy
    import gnssvis_simulation.realistic as model

    c = deepcopy(scenario.config)
    c["beam"].update(model="chromatic_elliptical_gaussian", **BEAM)
    c["solar"] = {"quiet_power_relative": 5.0, "active_power_relative": 0.0}
    c["realistic_emission"] = {
        "reference_range_m": 2e7,
        "peak_power_relative": 20.0,
        "broad_chip_rate_hz": 10.23e6,
        "narrow_chip_rate_hz": 1.023e6,
        "broad_fraction": 0.3,
    }
    p = c["pointing"]
    p.update(altitude_deg=60.0, azimuth_deg=0.0)
    sc = replace(
        scenario,
        config=c,
        duration_seconds=20,
        step_seconds=10,
        channel_count=50,
        band_min_hz=1150e6,
        band_max_hz=1200e6,
    )
    sat = replace(geometry.satellites[0], block="IIR")
    g = replace(
        geometry,
        satellites=[sat],
        alt_deg=np.array([[60.0, -30.0]]),
        az_deg=np.zeros((1, 2)),
        range_m=np.full((1, 2), 2e7),
        range_rate_m_s=np.zeros((1, 2)),
        above_horizon=np.array([[True, False]]),
        selected=np.array([[True, False]]),
        in_beam=np.array([[True, False]]),
        separation_deg=np.zeros((1, 2)),
    )
    monkeypatch.setattr(model, "compute_geometry", lambda sc: g)
    sun = {
        "alt_deg": np.array([60.0, -30.0]),
        "offset_x_deg": np.zeros(2),
        "offset_y_deg": np.zeros(2),
    }
    monkeypatch.setattr(model, "solar_ephemeris", lambda *args: sun)
    a = model.simulate_realistic(sc)
    b = model.simulate_realistic(sc)
    assert np.array_equal(a.total_power, b.total_power)
    assert not a.satellite_power.any()  # No L5 on IIR, even while in the beam.
    assert not a.solar_power[:, 1].any()  # Sun below horizon.
    assert np.all(a.solar_power[:, 0] > 0)
    g.satellites = [replace(sat, block="IIF")]
    c = model.simulate_realistic(sc)
    assert np.any(c.satellite_power[:, 0] > 0)
    assert not c.satellite_power[:, 1].any()
    np.testing.assert_array_equal(
        c.total_power, c.noise_power + c.satellite_power + c.solar_power
    )


def test_new_scenario_cannot_silently_run_through_legacy_pipeline(scenario):
    from dataclasses import replace
    from gnssvis_simulation.pipeline import simulate

    c = dict(scenario.config)
    c["realistic_emission"] = {}
    with pytest.raises(ValueError, match="simulate_realistic"):
        simulate(replace(scenario, config=c))


def test_importable_public_api():
    from gnssvis_simulation import simulate_realistic, save_realistic, RealisticResult

    assert callable(simulate_realistic) and callable(save_realistic)
    assert RealisticResult.__name__ == "RealisticResult"
