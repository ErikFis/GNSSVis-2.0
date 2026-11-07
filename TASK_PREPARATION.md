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
git config core.hooksPath .githooks
python tools/check_task_size.py
```

Git hooks are configured in the distributed local repository. If copying only
source files or cloning elsewhere, run the configuration command again.
Git does not propagate local hook configuration through a clone.

Commit the preparation before starting task-specific work if you reconstruct
the repository manually. The provided repositories already have their initial
commits and `task-base` tags. Do not move the baseline tag during the task.

## Size policy

The tool checks the entire delta from `task-base`, including changes in earlier
task commits. It sums the Git binary patch, changed file contents, and
non-ignored untracked file sizes, and rejects totals at or above 10,000,000
bytes. The pre-commit hook checks the staged version of that same aggregate
delta. Run the default command before preparing a submission to include
unstaged and untracked work too.

This is a conservative local check. Alignerr may measure a different artifact;
verify the actual submission format and its final size. A `.gitignore` cannot
enforce a size limit, remove files already tracked, or prevent `git add -f`.
The size hook still checks large files added by force. Hooks can be bypassed,
so always run the check explicitly before submission.

Ignore rules keep regenerable artifacts local. Do not submit the virtual
environment, downloaded data, caches, notebook outputs, or original upstream
Git history. Keep the two repositories as siblings, not as nested repositories.

## Dependencies and downloads

Install this library with `python -m pip install -e .`. No requirements export
is needed. The simulation scaffold uses its own `pyproject.toml` and has no
runtime dependency on gnssvis. First use of sky analysis can download large
ephemerides and orbital elements, which remain ignored.
