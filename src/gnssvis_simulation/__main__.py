"""Command line: ``python -m gnssvis_simulation --config SCENARIO [--output DIR]``.

Runs the offline synthetic GNSS waterfall and writes waterfall.png,
waterfall.npz, samples.csv, passages.csv and metadata.json. Never opens a
plot window. All parameters come from the scenario file; only the output
directory can be overridden.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m gnssvis_simulation",
                                     description="Offline synthetic GNSS waterfall simulation.")
    parser.add_argument("--config", required=True, type=Path, help="scenario JSON file")
    parser.add_argument("--output", type=Path, default=None,
                        help="output directory (default: output.default_directory in the "
                             "scenario, relative to the repository root)")
    args = parser.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")  # headless; figures are saved, never shown

    from .config import ScenarioError
    from .pipeline import run_simulation

    try:
        result = run_simulation(args.config, args.output)
    except (ScenarioError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0 if all(v is True for v in result.checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
