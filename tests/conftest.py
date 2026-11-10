import copy
import json
import socket
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_PATH = REPO_ROOT / "examples" / "scenarios" / "offline.json"

# Allow running the tests from a checkout without an editable install.
sys.path.insert(0, str(REPO_ROOT / "src"))


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Fail any test that tries to open a network connection."""

    def guard(*args, **kwargs):
        raise AssertionError("network access attempted during an offline test")

    monkeypatch.setattr(socket.socket, "connect", guard)
    monkeypatch.setattr(socket, "create_connection", guard)
    monkeypatch.setattr(socket, "getaddrinfo", guard)


@pytest.fixture(scope="session")
def scenario():
    from gnssvis_simulation import load_scenario

    return load_scenario(SCENARIO_PATH)


@pytest.fixture(scope="session")
def geometry(scenario):
    from gnssvis_simulation import compute_geometry

    return compute_geometry(scenario)


@pytest.fixture
def write_scenario(tmp_path):
    """Write a modified copy of the reference scenario; return a loader."""

    def _write(mutate):
        config = copy.deepcopy(json.loads(SCENARIO_PATH.read_text(encoding="utf-8")))
        mutate(config)
        path = tmp_path / "scenario.json"
        path.write_text(json.dumps(config), encoding="utf-8")
        from gnssvis_simulation import load_scenario

        return load_scenario(path, repo_root=REPO_ROOT)

    return _write
