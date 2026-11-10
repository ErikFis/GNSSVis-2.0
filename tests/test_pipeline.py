"""End-to-end pipeline: combined arrays, tables, metadata, figure and CLI."""

import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pytest

from conftest import REPO_ROOT, SCENARIO_PATH
from gnssvis_simulation import compute_geometry, load_scenario
from gnssvis_simulation.emission import compute_signal_power
from gnssvis_simulation.noise import generate_noise
from gnssvis_simulation.pipeline import PASSAGE_COLUMNS, SAMPLE_COLUMNS, run_simulation
from gnssvis_simulation.products import load_npz

NOISE_SEED42_SHA = "60f04c46fbc48b35309b379c31f83a1c2cd88bbcdda5f71c0d3bc7af1f39a766"


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


@pytest.fixture(scope="module")
def reference_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("offline")
    result = run_simulation(SCENARIO_PATH, output=out, verbose=False)
    return result, out


# -- reference outputs ----------------------------------------------------------

def test_output_files_follow_configuration(reference_run):
    result, out = reference_run
    cfg = result.scenario.config["output"]
    expected = {cfg[k] for k in ("waterfall_file", "figure_file", "samples_file",
                                 "passages_file", "metadata_file")}
    assert {p.name for p in out.iterdir()} == expected
    assert (out / cfg["figure_file"]).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_npz_arrays_and_linear_combination(reference_run):
    result, out = reference_run
    data = load_npz(out / "waterfall.npz")
    assert sorted(data) == sorted(result.scenario.config["output"]["npz_arrays"])
    for name in ("noise_power", "signal_power", "total_power"):
        assert data[name].dtype == np.float32 and data[name].shape == (512, 2160)
    assert data["frequencies_hz"].dtype == np.float64
    assert data["time_offsets_seconds"].dtype == np.int64
    assert data["time_offsets_seconds"][-1] == 21590
    np.testing.assert_array_equal(data["total_power"], data["noise_power"] + data["signal_power"])
    # Components are exactly the stand-alone stage outputs (unchanged by combining).
    sc = load_scenario(SCENARIO_PATH)
    np.testing.assert_array_equal(data["noise_power"], generate_noise(sc))
    np.testing.assert_array_equal(data["signal_power"], compute_signal_power(sc, compute_geometry(sc))[0])


def test_display_uses_full_db_range(reference_run):
    result, out = reference_run
    meta = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
    db = 10 * np.log10(result.total_power.astype(np.float64) / 1.0).astype(np.float32)
    assert meta["display"]["vmin_db"] == pytest.approx(float(db.min()))
    assert meta["display"]["vmax_db"] == pytest.approx(float(db.max()))
    # The transient L1 signal sets the top of the scale (not clipped at a percentile).
    i, j = np.unravel_index(np.argmax(db), db.shape)
    assert result.signal_power[i, j] > 0 and meta["display"]["vmax_db"] > 10.0


def test_samples_csv(reference_run):
    result, out = reference_run
    header, rows = read_csv(out / "samples.csv")
    assert header == SAMPLE_COLUMNS
    assert len(rows) == 2160 * 6 == 12960
    ids = result.scenario.norad_ids
    assert [int(r["norad_id"]) for r in rows[:12]] == ids * 2
    assert [int(r["sample_index"]) for r in rows[::6]] == list(range(2160))
    for col in ("above_horizon", "in_beam", "in_passage"):
        assert {r[col] for r in rows} <= {"true", "false"}
    assert {r["signal_ids"] for r in rows} == {"GPS_L1"}
    assert rows[0]["timestamp_utc"] == "2026-08-22T00:00:00Z"
    # CSV passage flags match the signal columns in the NPZ.
    in_passage = np.zeros(2160, dtype=bool)
    for r in rows:
        if r["in_passage"] == "true":
            in_passage[int(r["sample_index"])] = True
    np.testing.assert_array_equal(in_passage, result.signal_power.any(axis=0))


def test_passages_csv(reference_run):
    result, out = reference_run
    header, rows = read_csv(out / "passages.csv")
    assert header == PASSAGE_COLUMNS
    assert len(rows) == 1
    p = rows[0]
    assert p["norad_id"] == "40294" and p["signal_ids"] == "GPS_L1"
    assert p["entry_utc"] in ("2026-08-22T05:01:50Z", "2026-08-22T05:02:00Z")
    assert p["exit_utc_exclusive"] == "2026-08-22T05:39:20Z"
    assert int(p["duration_seconds"]) == 10 * int(p["sample_count"])
    assert (p["entry_clipped"], p["exit_clipped"]) == ("false", "false")
    assert float(p["peak_signal_power"]) == pytest.approx(float(result.signal_power.max()), rel=1e-6)
    # Signal columns in the NPZ span exactly entry .. exclusive exit.
    cols = np.flatnonzero(result.signal_power.any(axis=0))
    assert cols[0] * 10 == int(p["entry_offset_seconds"])
    assert (cols[-1] + 1) * 10 == int(p["exit_offset_seconds"])


def test_signal_at_l1_during_passage_only(reference_run):
    result, _ = reference_run
    sc = result.scenario
    rows = np.flatnonzero(result.signal_power.any(axis=1))
    assert np.all(np.abs(sc.frequencies_hz[rows] - 1575.42e6) <= 1.5e6)
    i, j = np.unravel_index(np.argmax(result.signal_power), result.signal_power.shape)
    assert i == int(np.argmin(np.abs(sc.frequencies_hz - 1575.42e6)))
    peak_time = sc.timestamps_utc[j]
    assert datetime(2026, 8, 22, 5, 1, 50, tzinfo=timezone.utc) <= peak_time \
        < datetime(2026, 8, 22, 5, 39, 20, tzinfo=timezone.utc)


def test_metadata(reference_run):
    result, out = reference_run
    meta = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
    assert meta["effective_config"] == json.loads(SCENARIO_PATH.read_text(encoding="utf-8"))
    assert meta["randomness"]["seed"] == 42
    assert meta["array_sha256"]["noise_power"] == NOISE_SEED42_SHA
    assert all(v is True for v in meta["checks"].values())
    assert meta["counts"]["sample_rows"] == 12960 and meta["counts"]["passage_rows"] == 1
    assert meta["grid"]["matrix_shape"] == [512, 2160]
    assert set(meta["reference_hashes_sha256"]) == {"satellites.tle", "satellites.csv", "signals.csv"}
    assert len(meta["satellites"]) == 6 and all(s["tle_epoch_utc"].endswith("Z") for s in meta["satellites"])
    for pkg in ("numpy", "skyfield", "sgp4", "matplotlib"):
        assert meta["dependency_versions"][pkg] != "not installed"
    assert 0.0 < meta["earth_orientation"]["ut1_minus_utc_seconds_at_start"] < 0.9
    assert "Doppler" in meta["models"]["emission"]["not_modelled"]


# -- variations -------------------------------------------------------------------

def test_seed_change_only_changes_background(write_scenario, tmp_path):
    def coarse(seed):
        return lambda c: (c["time_grid"].update(step_seconds=60)
                          or c["randomness"].update(seed=seed))
    a_s = write_scenario(coarse(42))
    a = run_simulation(a_s.config_path, output=tmp_path / "a", verbose=False)
    b_s = write_scenario(coarse(7))
    b = run_simulation(b_s.config_path, output=tmp_path / "b", verbose=False)
    assert not np.array_equal(a.noise_power, b.noise_power)
    assert not np.array_equal(a.total_power, b.total_power)
    np.testing.assert_array_equal(a.signal_power, b.signal_power)
    assert a.signal_power.any()
    assert (tmp_path / "a" / "samples.csv").read_bytes() == (tmp_path / "b" / "samples.csv").read_bytes()
    assert (tmp_path / "a" / "passages.csv").read_bytes() == (tmp_path / "b" / "passages.csv").read_bytes()


def test_no_passage_scenario(write_scenario, tmp_path):
    s = write_scenario(lambda c: c["pointing"].update(altitude_deg=5.0, azimuth_deg=180.0)
                       or c["beam"].update(full_width_deg=1.0)
                       or c["time_grid"].update(step_seconds=60))
    result = run_simulation(s.config_path, output=tmp_path, verbose=False)
    header, rows = read_csv(tmp_path / "passages.csv")
    assert header == PASSAGE_COLUMNS and rows == []
    assert not result.signal_power.any()
    np.testing.assert_array_equal(result.total_power, result.noise_power)
    assert (tmp_path / "waterfall.png").stat().st_size > 0


def test_command_line_is_headless_and_exits(write_scenario, tmp_path):
    s = write_scenario(lambda c: c["time_grid"].update(step_seconds=120)
                       or c["frequency_grid"].update(channel_count=128))
    env = {k: v for k, v in os.environ.items() if k not in ("DISPLAY", "MPLBACKEND")}
    # Any socket use inside the subprocess raises: proves the run is fully offline.
    guard = tmp_path / "guard"
    guard.mkdir()
    (guard / "sitecustomize.py").write_text(
        "import socket\n"
        "def _blocked(*a, **k):\n    raise RuntimeError('NETWORK ACCESS BLOCKED')\n"
        "socket.socket.connect = _blocked\nsocket.create_connection = _blocked\n"
        "socket.getaddrinfo = _blocked\n", encoding="utf-8")
    env["PYTHONPATH"] = os.pathsep.join([str(guard), str(REPO_ROOT / "src")])
    proc = subprocess.run([sys.executable, "-m", "gnssvis_simulation", "--config", str(s.config_path),
                           "--output", str(tmp_path / "out")], capture_output=True, text=True,
                          timeout=600, env=env, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert "NETWORK ACCESS BLOCKED" not in proc.stdout + proc.stderr
    assert "Consistency checks: 6/6 passed" in proc.stdout
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == [
        "metadata.json", "passages.csv", "samples.csv", "waterfall.npz", "waterfall.png"]
    bad = subprocess.run([sys.executable, "-m", "gnssvis_simulation", "--config",
                          str(tmp_path / "missing.json")], capture_output=True, text=True,
                         timeout=120, env=env)
    assert bad.returncode == 2 and "error" in bad.stderr


def test_generated_outputs_are_git_ignored():
    if shutil.which("git") is None:
        pytest.skip("git not available")
    paths = ["outputs/offline/waterfall.png", "outputs/offline/waterfall.npz",
             "outputs/offline/samples.csv", "outputs/offline/metadata.json",
             "partial_results/04_signal_waterfall.png", "partial_results/noise.npz"]
    proc = subprocess.run(["git", "check-ignore", *paths], cwd=REPO_ROOT, capture_output=True,
                          text=True)
    if proc.returncode == 128:
        pytest.skip("not inside a git work tree")
    assert proc.stdout.split() == paths


def test_legacy_visualizer_still_works(monkeypatch):
    """The refactored table helper keeps the original GNSSVis plot working."""
    import astropy.units as u
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from gnssvis import SkyAnalyzer, Visualizer

    times = Time("2026-08-22T05:00:00", scale="utc") + np.arange(60) * 60 * u.s
    sky = SkyAnalyzer(offline=True, tle_file=REPO_ROOT / "reference_data" / "satellites.tle")
    objects = sky.get_gnss_trajectories(times)
    shown = []
    monkeypatch.setattr(plt, "show", lambda: shown.append(plt.gcf()))
    Visualizer().plot_waterfall_with_objects(np.zeros((16, 60)), times, np.linspace(1500, 1600, 16),
                                             objects, SkyCoord(ra=0 * u.deg, dec=-7 * u.deg), 30.0)
    assert len(shown) == 1
    plt.close("all")
