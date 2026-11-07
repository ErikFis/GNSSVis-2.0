# Task preparation

This is a clean local baseline derived from:

- Upstream: https://github.com/ErikFis/gnssvis
- Source commit: `4b3f5b75e28dc5346e74544e47e9b1ecaf5c5c84`
- Local initial branch: `main`
- Task comparison reference: `task-base`

The upstream repository and its history have not been modified. This checkout
has a new local Git history so that downloaded datasets and caches are absent
from both the files and the local history. The source license is preserved.
The upstream URL is recorded here rather than configured as a push remote.

## Preparation changes

- Exclude observation data, NumPy arrays, ephemerides, downloaded catalogs,
  bytecode, installation metadata, and temporary artifacts.
- Keep documentation in the data directory and clear all notebook outputs.
- Use English documentation, comments, and runtime messages, without emojis.
- Refresh outdated documentation examples to match the actual public API.
- Preserve public API names, existing filenames, and astronomical aliases.
- Preserve `pyproject.toml` as the authoritative dependency declaration.
- Ignore local `requirements*.txt` files and requirements export directories.
- Provide the separate sibling `gnssvis-simulation` repository for future work.

## Start the task

```bash
git switch -c task/simulation-work
```

## Files for task review

Keep source code, relevant examples, documentation, packaging metadata, and
intentional test fixtures. Ignore local environments, bytecode, installation
metadata, downloaded catalogs, ephemerides, observation datasets, generated
simulation outputs, logs, temporary files, editor state, and build artifacts.

The ignore rules are based on file purpose. There is no file-size restriction,
size-check script, or commit hook. Small fixtures under `tests/fixtures/` remain
eligible for version control so reviewers can reproduce task behavior.

The provided repositories have local Git histories and `task-base` tags.
Keep the repositories as siblings rather than nesting one inside the other.

## Dependencies and downloads

Install this library with `python -m pip install -e .`. No requirements export
is needed. The simulation scaffold uses its own `pyproject.toml` and has no
runtime dependency on gnssvis. First use of sky analysis can download large
ephemerides and orbital elements, which remain ignored.
