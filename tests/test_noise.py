"""White-noise background: format, statistics and reproducibility.

Run with ``pytest -s`` to see the measured statistics printed.
"""

import numpy as np
import pytest
from scipy import stats

from gnssvis_simulation import ScenarioError
from gnssvis_simulation.diagnostics import format_table, noise_statistics
from gnssvis_simulation.noise import (AveragedWhiteNoise, generate_noise,
                                      noise_model_from_scenario, scenario_seed)


@pytest.fixture(scope="module")
def noise(scenario):
    return generate_noise(scenario)


def test_format_follows_matrix_contract(scenario, noise):
    assert noise.dtype == np.float32
    assert noise.shape == (512, 2160) == scenario.matrix_shape
    assert noise.flags.c_contiguous
    assert np.all(np.isfinite(noise)) and noise.min() > 0


def test_model_comes_from_configuration(scenario):
    model = noise_model_from_scenario(scenario)
    assert (model.mean_power, model.averages) == (1.0, 16)
    assert model.shape_parameter == 16.0 and model.scale_parameter == 1.0 / 16
    assert scenario_seed(scenario) == 42


def test_statistics_match_gamma_theory(scenario, noise):
    model = noise_model_from_scenario(scenario)
    st = noise_statistics(noise, model)
    print("\n" + format_table(
        ["statistic", "measured", "expected", "z"],
        [[k, f"{c['measured']:.6g}", f"{c['expected']:.6g}", f"{c['z']:+.2f}"]
         for k, c in st["checks"].items()]))
    print(f"KS D={st['ks_test']['statistic']:.2e} p={st['ks_test']['pvalue']:.3f}; "
          f"lag-1 corr time={st['lag1_correlation']['time']['measured']:+.5f} "
          f"freq={st['lag1_correlation']['frequency']['measured']:+.5f}")
    for name, check in st["checks"].items():
        assert abs(check["z"]) < 5, name
    for axis, check in st["lag1_correlation"].items():
        assert abs(check["z"]) < 5, axis
    assert st["ks_test"]["pvalue"] > 1e-3
    # Flat spectrum / stationary background: no channel or sample stands out.
    assert st["channel_means"]["max_abs_z"] < 5
    assert st["sample_means"]["max_abs_z"] < 5.5
    assert st["channel_means"]["sd_measured"] == pytest.approx(
        st["channel_means"]["sd_expected"], rel=0.15)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_statistics_hold_for_other_seeds(write_scenario, seed):
    """The check is not a property of one lucky seed (smaller grid for speed)."""
    s = write_scenario(lambda c: c["frequency_grid"].update(channel_count=128))
    x = generate_noise(s, seed=seed).astype(np.float64)
    n = x.size
    assert abs(x.mean() - 1.0) < 5 * 0.25 / np.sqrt(n)
    assert abs(x.var() - 1 / 16) < 5 * (1 / 16) * np.sqrt(2.375 / n)


def test_gamma_equals_averaged_complex_gaussian_power():
    """Gamma(M, P0/M) == mean of M |complex Gaussian|^2 with mean power P0."""
    rng = np.random.default_rng(123)
    m, p0, cells = 16, 2.5, 100_000
    z = (rng.standard_normal((cells, m)) + 1j * rng.standard_normal((cells, m))) * np.sqrt(p0 / 2)
    averaged = np.mean(np.abs(z) ** 2, axis=1)
    gamma = AveragedWhiteNoise(p0, m).generate((cells,), np.random.default_rng(7))
    ks = stats.ks_2samp(averaged, gamma)
    print(f"\ncomplex-Gaussian average: mean {averaged.mean():.4f} var {averaged.var():.5f}; "
          f"Gamma draw: mean {gamma.mean():.4f} var {gamma.var():.5f}; "
          f"theory mean {p0} var {p0**2/m:.5f}; KS p={ks.pvalue:.3f}")
    assert ks.pvalue > 1e-3
    assert averaged.mean() == pytest.approx(p0, rel=0.01)
    assert gamma.var() == pytest.approx(p0 ** 2 / m, rel=0.03)


def test_mean_power_and_averages_are_not_interchangeable():
    """P0 sets the mean, M sets the spread: swapping them changes the output."""
    rng = np.random.default_rng
    a = AveragedWhiteNoise(1.0, 16).generate((200_000,), rng(1)).astype(np.float64)
    b = AveragedWhiteNoise(16.0, 1).generate((200_000,), rng(1)).astype(np.float64)
    assert a.mean() == pytest.approx(1.0, rel=0.01) and a.var() == pytest.approx(1 / 16, rel=0.03)
    assert b.mean() == pytest.approx(16.0, rel=0.01) and b.var() == pytest.approx(256, rel=0.03)


def test_same_seed_reproduces_bitwise(scenario, noise):
    assert np.array_equal(noise, generate_noise(scenario))
    assert np.array_equal(noise, generate_noise(scenario, seed=42))


def test_different_seed_changes_only_the_background(scenario, noise, geometry):
    from gnssvis_simulation import compute_geometry

    other = generate_noise(scenario, seed=43)
    assert not np.array_equal(noise, other)
    assert abs(np.corrcoef(noise.ravel(), other.ravel())[0, 1]) < 5 / np.sqrt(noise.size)
    again = compute_geometry(scenario)
    assert np.array_equal(again.alt_deg, geometry.alt_deg)
    assert np.array_equal(again.selected, geometry.selected)


def test_global_numpy_random_state_is_untouched(scenario):
    np.random.seed(5)
    expected = np.random.random(3)
    np.random.seed(5)
    generate_noise(scenario)
    assert np.array_equal(np.random.random(3), expected)


def test_seed_42_reference_values(noise):
    """Pins the seed-42 stream (NumPy 2.3 Generator.gamma, float64 -> float32, C order).

    Catches accidental changes to draw order, dtype or shape. If it fails only
    after a NumPy upgrade, the statistical tests above still decide
    correctness; the versions are recorded in the metadata.
    """
    from gnssvis_simulation.products import array_sha256

    first = np.array([1.0564992427825928, 1.1767948865890503, 0.5714749693870544,
                      1.0111337900161743], dtype=np.float32)
    last = np.array([1.0699169635772705, 1.1166319847106934], dtype=np.float32)
    print(f"\nseed 42: first cells {noise[0, :4]}, sha256 {array_sha256(noise)[:16]}...")
    np.testing.assert_array_equal(noise[0, :4], first)
    np.testing.assert_array_equal(noise[-1, -2:], last)
    assert array_sha256(noise) == "60f04c46fbc48b35309b379c31f83a1c2cd88bbcdda5f71c0d3bc7af1f39a766"


def test_configurable_resolution(write_scenario):
    s = write_scenario(lambda c: c["time_grid"].update(step_seconds=60, duration_seconds=3600)
                       or c["frequency_grid"].update(channel_count=64))
    x = generate_noise(s)
    assert x.shape == (64, 60) and x.dtype == np.float32


@pytest.mark.parametrize("mutate, message", [
    (lambda c: c["noise"].update(averages=0), "averages"),
    (lambda c: c["noise"].update(averages=2.5), "averages"),
    (lambda c: c["noise"].update(mean_power_relative=-1.0), "mean power"),
    (lambda c: c["noise"].update(model="pink"), "Unsupported noise model"),
    (lambda c: c["randomness"].update(generator="numpy.random.RandomState"), "default_rng"),
])
def test_invalid_noise_configuration(write_scenario, mutate, message):
    s = write_scenario(mutate)
    with pytest.raises(ScenarioError, match=message):
        generate_noise(s)


def test_randomness_disabled_is_unseeded(write_scenario):
    s = write_scenario(lambda c: c["randomness"].update(enabled=False))
    assert scenario_seed(s) is None
    assert scenario_seed(s, seed=7) == 7
