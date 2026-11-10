# Frozen offline reference inputs

These small files are intentional task inputs, not runtime caches or observation
results. The first scenario uses GPS as a controlled subset of GNSS.

| File | Purpose |
| --- | --- |
| `satellites.tle` | Six real GPS element sets, in name/line-1/line-2 format |
| `satellites.csv` | NORAD identifiers, names, blocks, PRNs, and actual TLE epochs |
| `signals.csv` | Nominal L1, L2, and L5 carrier frequencies in integer Hz |
| `provenance.json` | Immutable upstream source, standards references, and hashes |

The observation specification is
[`examples/scenarios/offline.json`](../examples/scenarios/offline.json).
See the [scenario guide](../docs/OBSERVATION_SCENARIO.md) for the fixed decisions
and [first prompt](../docs/FIRST_PROMPT.md) for the development context.

## TLE provenance and epoch

The elements are a subset of `data/tle_repository/gnss/gnss.txt` from upstream
GNSSVis commit `4b3f5b75e28dc5346e74544e47e9b1ecaf5c5c84`. Upstream identifies
CelesTrak as its GNSS provider. The original download timestamp is not known;
this pack does not claim a fresh download directly from CelesTrak.

The selected elements are dated August 20-22, 2026, and the scenario starts on
August 22, 2026. Read the epoch embedded in each element set. Do not treat an
archive filename, filesystem modification time, or today's date as the epoch.
Do not update these frozen files automatically during an offline task run.

## Frequency provenance and interpretation

| Carrier | Frequency in Hz | Reference |
| --- | ---: | --- |
| GPS L1 | 1575420000 | IS-GPS-200N, section 3.3.1.1 |
| GPS L2 | 1227600000 | IS-GPS-200N, section 3.3.1.1 |
| GPS L5 | 1176450000 | IS-GPS-705J, section 3.3.1.1 |

Official source URLs and document hashes are recorded in `provenance.json`.
The standards PDFs are not bundled; no network access is required to read the
CSV. These are nominal carrier centers, not measured mean frequencies.

The chosen satellites belong to GPS IIF/III blocks. The scenario assumes their
nominal L1/L2/L5 carrier capability; it makes no assertion about health or actual
transmission at the observation time. The frequency catalog does not specify measured occupied bandwidth,
Doppler, calibrated power, polarization, or interference received by an actual telescope. The application will add illustrative Gaussian carrier profiles and relative-power coefficients from the scenario configuration.

## What is not implemented

The JSON and CSV files specify the future application's input contract. The
current GNSSVis methods do not automatically read this scenario. Implementing
the shared loader, offline geometry path, and application is the programming
task. TLE presence alone does not disable the current automatic downloads.
