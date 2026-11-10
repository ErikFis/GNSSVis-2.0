"""NPZ/JSON storage contract and the partial-results command."""

import json

import numpy as np
import pytest

from gnssvis_simulation.noise import generate_noise
from gnssvis_simulation.products import (coordinate_arrays, load_npz, save_npz, validate_arrays)


@pytest.fixture
def small(write_scenario):
    return write_scenario(lambda c: c["time_grid"].update(step_seconds=60)
                          or c["frequency_grid"].update(channel_count=128))


def test_npz_roundtrip_preserves_dtypes_and_axes(small, tmp_path):
    noise = generate_noise(small)
    path = save_npz(tmp_path / "noise.npz", {"noise_power": noise, **coordinate_arrays(small)}, small)
    data = load_npz(path)  # allow_pickle=False
    assert data["noise_power"].dtype == np.float32 and data["noise_power"].shape == (128, 360)
    assert data["frequencies_hz"].dtype == np.float64 and data["frequencies_hz"].shape == (128,)
    assert data["time_offsets_seconds"].dtype == np.int64
    assert data["time_offsets_seconds"][:3].tolist() == [0, 60, 120]
    assert np.array_equal(data["noise_power"], noise)
    # Rows are frequency: row i is centred on frequencies_hz[i].
    assert data["frequencies_hz"][0] == small.band_min_hz + small.channel_width_hz / 2


@pytest.mark.parametrize("bad, message", [
    (lambda a: a.update(noise_power=a["noise_power"].astype(np.float64)), "float32"),
    (lambda a: a.update(noise_power=np.ascontiguousarray(a["noise_power"].T)), "float32"),
    (lambda a: a.update(noise_power=-a["noise_power"]), "non-negative"),
    (lambda a: a.update(frequencies_hz=a["frequencies_hz"] / 1e6 + 0j), "float64"),
    (lambda a: a.update(time_offsets_seconds=a["time_offsets_seconds"].astype(float)), "int64"),
    (lambda a: a.update(noise_db=a["noise_power"]), "Unexpected array"),
    (lambda a: a.pop("frequencies_hz"), "Missing coordinate"),
])
def test_contract_violations_are_rejected(small, bad, message):
    arrays = {"noise_power": generate_noise(small), **coordinate_arrays(small)}
    bad(arrays)
    with pytest.raises(ValueError, match=message):
        validate_arrays(arrays, small)


def test_partial_results_command(small, tmp_path, capsys):
    from gnssvis_simulation.partial import run

    files = run(small.config_path, output=tmp_path)
    out = capsys.readouterr().out
    names = {p.name for p in files}
    assert names == {"noise.npz", "noise_metadata.json", "signal.npz", "signal_metadata.json",
                     "01_geometry_separation.png", "02_noise_waterfall.png",
                     "03_noise_statistics.png", "04_signal_waterfall.png"}
    for p in files:
        assert p.stat().st_size > 0
    for p in files:
        if p.suffix == ".png":
            assert p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    meta = json.loads((tmp_path / "noise_metadata.json").read_text(encoding="utf-8"))
    assert meta["grid"]["matrix_shape"] == [128, 360]
    assert meta["randomness"]["seed"] == 42
    assert meta["noise_model"]["averages"] == 16
    assert meta["dependency_versions"]["numpy"] == np.__version__
    assert meta["effective_config"]["time_grid"]["step_seconds"] == 60
    assert all(meta["reproducibility_checks"][k] for k in
               ("same_seed_identical", "other_seed_differs", "geometry_unchanged",
                "npz_roundtrip_identical"))
    assert meta["display"]["vmin_db"] < 0 < meta["display"]["vmax_db"]
    sig = json.loads((tmp_path / "signal_metadata.json").read_text(encoding="utf-8"))
    assert sig["in_band_carriers"] == ["GPS_L1"]
    assert all(v is True for k, v in sig["emission"]["checks"].items()
               if k in ("nonzero_columns_equal_passage_masks", "zero_outside_masks",
                        "noise_regenerated_unchanged"))
    signal = load_npz(tmp_path / "signal.npz")["signal_power"]
    assert signal.dtype == np.float32 and signal.shape == (128, 360)
    for heading in ("SCENARIO AND GRIDS", "GEOMETRY", "WHITE-NOISE BACKGROUND",
                    "SATELLITE EMISSION", "REPRODUCIBILITY", "FILES WRITTEN"):
        assert heading in out
