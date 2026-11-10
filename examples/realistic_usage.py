"""Run the wideband demonstration from the checkout, with optional download."""

from pathlib import Path

# These imports do not download data or require an offline-only machine.
from gnssvis_simulation import load_scenario
from gnssvis_simulation.realistic import simulate_realistic, save_realistic


def main():
    root = Path(__file__).resolve().parents[1]
    scenario = load_scenario(root / "examples/scenarios/pibic_24h.json")
    result = simulate_realistic(scenario)  # DE421 is cached/downloaded if needed.
    save_realistic(result, root / "outputs/pibic_24h")


if __name__ == "__main__":
    main()
