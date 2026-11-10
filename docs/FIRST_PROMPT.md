# First Taiga message

Rewrite the following in your own words if desired. It is technical task
context, not a model-behavior rating, flag comment, or comparison justification.
The model does not see the reviewer's task-brief fields, so include this context
in the conversation. Upload the GitHub source ZIP first, without any `.git`
directory, and keep its extraction destination at `/workspace`.

## English starting prompt

I am extending my existing GNSSVis project with a separate application for
analyzing potential GNSS interference during radio telescope observations.
This is an additional use of GNSSVis as a calculation library, not a full radio
waveform simulator or a scientific validation suite for the entire library.

Reuse the existing trajectory and angular-geometry calculations instead of
building a disconnected implementation. Keep the application logic separate
from the existing modules and preserve their default public behavior. Small
backward-compatible changes to enable a GNSS-only offline path are in scope.

The input files are already in this upload:

- `examples/scenarios/offline.json`: the single source of truth for observation
  parameters and processing conventions.
- `reference_data/satellites.tle`: six frozen real GPS TLEs, with their epochs
  and identities in `reference_data/satellites.csv`.
- `reference_data/signals.csv`: nominal GPS L1/L2/L5 carriers in Hz.
- `reference_data/provenance.json`: sources and integrity hashes.
- `docs/OBSERVATION_SCENARIO.md`: decisions, output schema, and limitations.

Read these files first. The scenario is already chosen: a six-hour historical
window starting at 2026-08-22 00:00 UTC, 60-second samples, the UIRAPURU reference
site, a fixed zenith pointing, a circular 30-degree full-width beam, and a
1500-1600 MHz receiver band. Use the JSON values throughout the workflow. Do not
ask me to redefine these parameters or silently replace them with different
library defaults. Raise genuine contradictions or missing essentials before
implementing the affected behavior.

The simulation application and JSON loader are not implemented yet. Start by
inspecting the current code, particularly `SkyAnalyzer.__init__`,
`get_gnss_trajectories`, and `find_objects_in_fov`, and check dependency
availability. These paths currently involve downloads. Implement the smallest
maintainable GNSS-only path that accepts the bundled TLEs and uses no network
access, no mandatory DE421 file, and no external Earth-orientation updates.
Do not use file modification times or today's date in place of TLE epochs.

Then develop the configuration-driven application incrementally: calculate
topocentric geometry in the agreed local frame, identify sampled beam passages,
match nominal signal carriers with the receiver band, and generate
`samples.csv`, `passages.csv`, and `metadata.json`. Follow the units, timestamp
boundaries, clipping rules, and schemas in the scenario guide. Nominal carrier
overlap indicates potential interference only; do not add received-power,
Doppler, occupied-bandwidth, or full waveform models.

The proposed execution interface to implement is:

```bash
python -m gnssvis_simulation --config examples/scenarios/offline.json --output outputs/offline
```

Provide meaningful automated tests and a documented offline example. Check
controlled geometry expectations, beam boundaries, no-passage and out-of-band
cases, invalid inputs, configuration propagation, reproducibility with stated
tolerances, and at least one real local-TLE integration path. Do not use mocks
as the only evidence that the application works offline.

This workspace has no general internet access. Keep generated outputs ignored,
do not run `git init` or `git clone` inside `/workspace`, and avoid unrelated
refactoring. Keep code, comments, documentation, and runtime messages in English
without emojis. Report the checks you actually ran and their results, and stop
after the agreed behavior and checks are complete.

## Reviewer brief reminder

The task description should state the existing offline coupling and the new
application goal. Completion criteria should reference the actual implemented
command, the scenario guide, the 360-timestamp / six-satellite grid, the signal
band decisions, offline integration, and passing meaningful tests. Keep proposed
commands labeled as deliverables until they have actually been implemented.
Write your own behavior evaluations and final assessment from observed evidence.
