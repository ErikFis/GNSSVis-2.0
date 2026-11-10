# Task preparation

## Starting point

This is a clean source baseline derived from:

- Original upstream: https://github.com/ErikFis/gnssvis
- Source commit: `4b3f5b75e28dc5346e74544e47e9b1ecaf5c5c84`
- Prepared public repository: https://github.com/ErikFis/GNSSVis-2.0

The original repository has not been rewritten. The prepared repository has a
fresh source history without the large original data objects. The source
license is preserved. Preparation removed observation data, caches, bytecode,
and installation artifacts, cleared notebook outputs, and used English prose
and messages without changing the original numerical operations.

## Fixed task inputs

The scenario in `examples/scenarios/offline.json` is the canonical source for
observation parameters. `reference_data/` contains six real frozen GPS TLEs,
a satellite catalog, nominal L1/L2/L5 carriers, and provenance with hashes.
These small reference files are intentionally versioned rather than ignored.

Read `docs/OBSERVATION_SCENARIO.md` for the agreed conventions and
`docs/FIRST_PROMPT.md` for context to supply at the start of the Taiga session.
The scenario loader, offline GNSS-only integration, and simulation application
are not implemented yet. They are the work to perform during the task.

## Upload to Taiga

Download the source ZIP of the prepared repository's main branch. GitHub source
ZIPs exclude `.git`. Verify that the upload also excludes dependency folders,
local environments, generated outputs, and caches. Keep `.gitignore`.
Use Mount an archive with extraction destination `/workspace` before sending
the first message. Do not run `git init` or `git clone` in that workspace.

Do not upload the older `gnssvis-task-ready.zip`: it was a local-checkout bundle
containing `.git` directories. It is not the Taiga input archive.

## Review contents

Keep source, meaningful tests, documentation, packaging metadata, and intentional
reference fixtures. Ignore downloaded catalogs, planetary ephemerides, raw
observations, generated output tables, logs, build artifacts, and environments.
There is no size-check script, file-size restriction, or commit hook.

Dependencies remain declared in `pyproject.toml`; local requirements exports
are ignored. Check package availability in the task environment's internal
mirror. The existing default sky-analysis path can still download ephemerides
and TLEs. Adding frozen inputs alone does not make those methods offline.
