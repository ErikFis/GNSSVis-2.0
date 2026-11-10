# Shared observation scenario

The single source of truth is
[`examples/scenarios/offline.json`](../examples/scenarios/offline.json).
It is a complete reference specification for the application to be developed,
not evidence that the existing GNSSVis runtime already consumes JSON.

This historical synthetic scenario is deliberately selected to contain both a
beam passage and satellites that do not enter the beam. It is not a real FITS
observation or a measurement of interference.

## Fixed parameters

| Parameter | Reference value |
| --- | --- |
| Scenario ID | `uirapuru_gps_zenith_20260822` |
| Start | `2026-08-22T00:00:00Z` |
| Duration | 21600 seconds (6 hours) |
| End, excluded | `2026-08-22T06:00:00Z` |
| Time step | 60 seconds |
| Sample grid | 360 timestamps, from 00:00 to 05:59 UTC |
| Recife display equivalent | August 21 at 21:00 through August 22 at 03:00, UTC-03:00 |
| Observatory | UIRAPURU reference site, from GNSSVis defaults |
| Latitude / longitude | -7.211637 / -35.908138 degrees; east-positive longitude |
| Reference height | 553 meters, treated as WGS84 ellipsoidal height in this simulation |
| Pointing | Fixed local altitude 90 degrees, azimuth 0 degrees (zenith) |
| Beam | Circular top-hat; full width 30 degrees, radius 15 degrees |
| Horizon | Altitude at least 0 degrees; no atmospheric refraction |
| Receiver band | 1500000000 to 1600000000 Hz, inclusive |
| Satellites | NORAD 39166, 40294, 40534, 40730, 41328, 55268 |
| Randomness | None; no seed is required |
| Reference input formats | JSON configuration, three-line TLE records, CSV catalogs |
| Generated output formats | CSV time series, CSV passage table, JSON metadata |
| Real FITS input / output | Neither is required for this task |
| Internet / DE421 | Neither may be required by the new GNSS-only path |

The coordinates and height are inherited reference defaults, not independently
surveyed site parameters. Zenith is a fixed ground-based direction; do not
silently reinterpret it as fixed RA/Dec. Azimuth is immaterial exactly at the
zenith. Compute both satellite and pointing geometry in a consistent frame.

## Geometry and passage rules

Use GNSSVis as the reusable calculation backend. Adapt its GNSS path minimally
so it can propagate local TLEs without mandatory planetary ephemerides or HTTP.
Use SGP4, WGS84, and local unrefracted topocentric geometry. Use only bundled
timescale resources; do not fetch Earth-orientation tables or leap-second files.

For the chosen zenith pointing, angular separation equals 90 degrees minus
unrefracted altitude. A sample is in the beam when angular separation is at most
15 degrees and the satellite is above the configured horizon. Equality at the
beam and frequency boundaries is included.

A passage is a maximal contiguous run of in-beam samples for one satellite.
Entry is the first included sample. Exit is the first following excluded
sample, or the observation end when the run reaches the final sample. Duration
is the number of included samples multiplied by 60 seconds. Do not interpolate
subsample crossing times. Record entry/exit clipping at observation boundaries.
Do not discard the out-of-beam samples from the time-series output.

The selected elements remain within 72 hours of their epochs across this grid.
Validate that against the actual TLE epochs, not file ages. This bound is an
input rule for this example, not a guarantee of orbit accuracy.

## Frequency and interference rules

Match the nominal carrier center against the receiver band inclusively. GPS L1
is in band; GPS L2 and L5 are out of band. This gives a positive and two negative
frequency cases without ambiguous channel assignments or signal bandwidths.

Potential interference requires the satellite to be in the visible beam and
have at least one nominal carrier in band. This does not predict received
power, verify operational transmission, or establish actual interference.

## Output contract to implement during the task

Proposed command, not currently implemented:

```bash
python -m gnssvis_simulation --config examples/scenarios/offline.json --output outputs/offline
```

- `samples.csv`: one row per timestamp and selected satellite, including
  `sample_index`, `timestamp_utc`, `norad_id`, `altitude_deg`, `azimuth_deg`,
  `angular_separation_deg`, `above_horizon`, `in_beam`, `in_band_signal_ids`,
  and `potential_interference`. The reference grid gives 2160 rows. Order rows
  by timestamp, then by the scenario's NORAD list. Use lowercase `true`/`false`
  for booleans and semicolons for multiple signal IDs in a single CSV cell.
- `passages.csv`: one row per sampled geometric passage, with `norad_id`,
  `entry_utc`, `exit_utc_exclusive`, `duration_seconds`, `sample_count`,
  `entry_clipped`, `exit_clipped`, and `in_band_signal_ids`. Retain a header when
  there are no passages. Sort by entry time, then NORAD ID.
- `metadata.json`: effective configuration, units, reference hashes,
  dependency versions, input epochs, row/passage counts, and limitations.

The default output directory is ignored. Do not commit generated output tables
as a substitute for a working offline reproduction command. Numerical tests
should use documented tolerances, not bit-for-bit comparisons across library
versions. With these inputs, a geometric spot check using Skyfield 1.55 found a
zenith-beam passage of GPS BIIF-8 / NORAD 40294 around 05:02-05:40 UTC. That is a
scenario-selection sanity check, not an independent scientific ground truth or
a replacement for testing the implemented GNSSVis integration.

## One configuration for the whole future workflow

Implement a shared loader and validation during the task. The application,
library adapter, tests, and run metadata must read this file or receive its
validated values. Do not silently fall back to conflicting defaults in
`GNSSVis` or `SkyAnalyzer`. Derive end time, beam radius, sample counts, and
timestamps from the canonical configuration instead of adding separate defaults.

Resolve reference paths relative to the repository root, not whichever working
directory happens to launch the command. An explicit CLI output-directory
override is allowed. Do not alter the scenario automatically to make a test pass.

## Starting limitations in the existing code

`SkyAnalyzer.__init__` currently loads DE421. `get_gnss_trajectories` still
downloads TLEs when its session cache is empty, even if a local TLE file exists.
`find_objects_in_fov` also asks for solar system objects. The task must address
these couplings for its GNSS-only execution path while preserving default APIs.
No application code has been implemented as part of this reference preparation.
