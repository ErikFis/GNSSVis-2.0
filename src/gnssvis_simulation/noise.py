"""White-noise background for the synthetic waterfall.

Each cell is the power of M averaged independent zero-mean complex Gaussian
white-noise samples with mean power P0. The average of M exponential variables
with mean P0 is exactly Gamma(shape=M, scale=P0/M), so the cells are drawn
directly; no voltage waveform or FFT is needed. Power pixels are therefore
Gamma distributed, not Gaussian.

A background model only needs ``generate(shape, rng) -> float32 array``;
other models (e.g. a frequency-dependent system temperature) can be added
later without changing the pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Protocol

import numpy as np

from .config import Scenario, ScenarioError

POWER_DTYPE = np.float32


class NoiseModel(Protocol):
    """Interface for background power models (relative linear power)."""

    def generate(self, shape: tuple, rng: np.random.Generator) -> np.ndarray:
        """Return a float32 (frequency, time) array of linear relative power."""


@dataclass(frozen=True)
class AveragedWhiteNoise:
    """Averaged white-noise power: Gamma(shape=M, scale=P0/M) per cell."""

    mean_power: float  # P0, relative linear power
    averages: int  # M, number of averaged power samples per cell

    def __post_init__(self):
        if not self.mean_power > 0:
            raise ScenarioError(f"noise mean power must be positive, got {self.mean_power}")
        if int(self.averages) != self.averages or self.averages < 1:
            raise ScenarioError(f"noise averages must be a positive integer, got {self.averages}")

    @property
    def shape_parameter(self) -> float:
        return float(self.averages)

    @property
    def scale_parameter(self) -> float:
        return self.mean_power / self.averages

    def generate(self, shape: tuple, rng: np.random.Generator) -> np.ndarray:
        # Drawn in float64 by Generator.gamma, then stored as float32.
        draws = rng.gamma(self.shape_parameter, self.scale_parameter, size=shape)
        return draws.astype(POWER_DTYPE)

    def expected_moments(self) -> Dict[str, float]:
        """Theoretical moments of the linear power and of 10*log10(P/P0)."""
        from scipy.special import digamma, polygamma

        m, p0 = float(self.averages), self.mean_power
        db = 10.0 / np.log(10.0)
        return {
            "mean": p0,
            "variance": p0 ** 2 / m,
            "std": p0 / np.sqrt(m),
            "skewness": 2.0 / np.sqrt(m),
            "excess_kurtosis": 6.0 / m,
            # E[ln X] = digamma(M) + ln(P0/M); Var[ln X] = trigamma(M)
            "mean_db_re_p0": float(db * (digamma(m) - np.log(m))),
            "std_db": float(db * np.sqrt(polygamma(1, m))),
        }

    def describe(self) -> Dict[str, object]:
        return {
            "model": "independent_averaged_white_noise_power",
            "distribution": "Gamma",
            "shape": self.shape_parameter,
            "scale": self.scale_parameter,
            "mean_power_relative": self.mean_power,
            "averages": self.averages,
        }


def noise_model_from_scenario(scenario: Scenario) -> AveragedWhiteNoise:
    cfg = scenario.section("noise")
    if cfg.get("model") != "independent_averaged_white_noise_power":
        raise ScenarioError(f"Unsupported noise model {cfg.get('model')!r}")
    return AveragedWhiteNoise(mean_power=float(cfg["mean_power_relative"]),
                              averages=cfg["averages"])


def scenario_seed(scenario: Scenario, seed: Optional[int] = None) -> Optional[int]:
    """Effective seed: explicit override, else the configured seed.

    Returns None only when randomness.enabled is false (non-reproducible).
    """
    if seed is not None:
        return int(seed)
    randomness = scenario.section("randomness")
    if randomness.get("generator", "numpy.random.default_rng") != "numpy.random.default_rng":
        raise ScenarioError("Only numpy.random.default_rng is supported")
    return int(randomness["seed"]) if randomness.get("enabled", True) else None


def generate_noise(scenario: Scenario, seed: Optional[int] = None,
                   model: Optional[NoiseModel] = None) -> np.ndarray:
    """Background power on the scenario (frequency, time) grid, float32.

    A fresh ``numpy.random.default_rng(seed)`` is used, so the background is
    independent of geometry and of any other pipeline stage.
    """
    model = model or noise_model_from_scenario(scenario)
    rng = np.random.default_rng(scenario_seed(scenario, seed))
    power = model.generate(scenario.matrix_shape, rng)
    if power.shape != scenario.matrix_shape or power.dtype != POWER_DTYPE:
        raise ValueError(f"Noise model returned {power.dtype} {power.shape}, "
                         f"expected float32 {scenario.matrix_shape}")
    return power
