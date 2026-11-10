# Synthetic waterfall scenario

The single source of truth is `examples/scenarios/offline.json`. It specifies
an application to be implemented during the programming task. Adding this JSON
alone does not make the current GNSSVis library consume it or run offline.

## Final deliverable

Generate a reproducible synthetic spectro-temporal waterfall with a white-noise
background and simplified GNSS emission contributions during satellite beam
passages. Reuse GNSSVis geometry and visual conventions. The figure must have:

- An upper panel showing satellite passage geometry against UTC time.
- A lower waterfall with time on X, frequency in MHz on Y, and relative power
  represented by color.
- Dotted horizontal nominal-carrier markers with frequency labels.
- A right-hand satellite/color table, consistent with the upper panel.
- A labeled colorbar and a time axis shared between the panels.

The upper panel uses angular separation from the actual configured pointing,
with the 15-degree selection boundary. Do not copy a relative-declination axis
that is inconsistent with the fixed-ground zenith pointing. Preserve the layout
and reuse/adapt plotting code without breaking the existing public behavior.
The current visualizer calls `show()` and uses fixed RA/Dec-related assumptions;
consider a small backward-compatible adapter or extension for the new workflow.

## Fixed observation and grids

| Parameter | Value |
| --- | --- |
| Scenario | `uirapuru_gps_zenith_20260822` |
| UTC start | `2026-08-22T00:00:00Z` |
| Duration / end | 21600 seconds / `2026-08-22T06:00:00Z`, end excluded |
| Time step | 10 seconds; 2160 timestamps, last sample 05:59:50 UTC |
| Recife equivalent | August 21, 21:00 to August 22, 03:00, UTC-03:00 |
| Latitude / longitude | -7.211637 / -35.908138 degrees, east-positive longitude |
| Height | 553 m, treated as WGS84 ellipsoidal height for this reference |
| Pointing | Fixed local zenith, altitude 90 degrees; azimuth 0 is arbitrary |
| Beam selection | Circular gate, full width 30 degrees; separation <= 15 degrees |
| Visibility | Unrefracted altitude >= 0 degrees |
| Receiver band | 1500000000 to 1600000000 Hz, inclusive |
| Frequency grid | 512 uniform bin centers; 195312.5 Hz channel spacing |
| Selected NORAD IDs | 39166, 40294, 40534, 40730, 41328, 55268 |
| Random seed | 42, using NumPy default_rng |
| Power-array shape | `(512, 2160)`, frequency then time |
| Real FITS input | Not required; all simulation output is generated locally |

This historical scenario is near the actual epochs of the bundled real TLEs.
Its site coordinates are inherited GNSSVis defaults, not newly surveyed values.
The original TLE acquisition timestamp is unknown. Source URLs and hashes are
in `reference_data/provenance.json`. No HTTP, DE421 or external Earth-orientation
update may be required for the new GNSS-only path. Check TLE epochs rather than
filenames or modification dates; the selected inputs stay within 72 hours of
those epochs across the observation window.

## Geometry and time mask

Use SGP4 with consistent unrefracted topocentric geometry and WGS84. Reuse the
existing GNSSVis calculation interface with the minimum necessary offline
extension. The current initializer loads DE421 and the trajectory path can
still download despite a local file; those couplings are part of the task.

For zenith pointing, angular separation is 90 degrees minus satellite altitude.
The time mask is the in-beam, above-horizon decision at each sample. Injection
must be zero outside that mask. Do not claim continuous exact crossing times:
the reference resolution is 10 seconds, with no subsample interpolation.

A passage is a maximal contiguous run of selected samples. Entry is the first
included sample. Exit is the first following excluded sample, or the observation
end if clipped. Duration is included-sample count times the time step. Record
entry/exit clipping at observation boundaries.

## White-noise background

Use independent averaged white-noise power cells, with mean relative power
`P0 = 1` and `M = 16` averages. A direct Gamma(shape=M, scale=P0/M) draw is an
efficient equivalent of averaging squared magnitudes of independent zero-mean
complex Gaussian white-noise samples. A raw voltage waveform or FFT is not
required. Power pixels follow Gamma statistics, not a Gaussian distribution.

This is a synthetic white background with no antenna-temperature, receiver-noise,
or calibrated radiometer model. The background is independent of satellite
passages and is reproducible for a fixed seed and dependency versions.

## Simplified emission and signal mask

The catalog provides nominal GPS carrier centers L1/L2/L5. With this receiver
band, only L1 at 1575420000 Hz is injected. L2 and L5 provide out-of-band cases.
These are carrier frequencies, not measured per-satellite frequency averages.

For satellite i and nominal carrier c:

```
A_i(t) = A_ref * (r_ref / r_i(t))**2
G_i(t) = exp(-4*ln(2) * (separation_i(t) / full_width_deg)**2)
F_c(f) = exp(-0.5 * ((f - carrier_hz) / sigma_hz)**2)
P_i,c(f,t) = passage_mask_i(t) * A_i(t) * G_i(t) * F_c(f)
```

Use `A_ref=20` relative power units, `r_ref=20000000` meters, `sigma_hz=500000`,
and truncate the spectral support outside three sigma. The geometric gate is
hard; the angular taper is Gaussian inside it, with half gain at 15 degrees.
This intentional truncation can cause an edge discontinuity. It is a simple
illustrative model, not a measured GPS spectral shape or a calibrated antenna.

Profiles are selected by satellite block (`IIF` or `III`). Both have identical
illustrative coefficients for now. This does not imply that real block powers
are identical, that any signal was healthy/active at that time, or that the model
provides absolute received power. Do not invent satellite-specific EIRP values.

Calculate profiles at the actual frequency-bin centers. Do not shift the
frequency markers to a nearest bin. A boolean spectral support mask and a
weighted emission profile are distinct. Sum overlapping satellite contributions
instead of overwriting them. Extend the profile catalog later for other GNSS
constellations; the shipped example is GPS only.

Combine everything in linear power:

```
total_power = noise_power + sum_of_signal_contributions
relative_db = 10*log10(total_power / P0)
```

Convert to dB only for display. Do not add dB arrays, inject a signal outside its
passage, or hardcode a rectangle unrelated to calculated geometry. Do not add
Doppler, multipath, actual modulation, polarization, antenna thermal noise,
receiver electronics, or an absolute calibration in this first task.

## Outputs and reproduction contract

Proposed command to implement:

```bash
python -m gnssvis_simulation --config examples/scenarios/offline.json --output outputs/offline
```

- `waterfall.npz`: `noise_power`, `signal_power`, `total_power` as float32 arrays,
  `frequencies_hz` as float64, and integer `time_offsets_seconds`. Preserve units
  and array order; avoid pickle/object arrays.
- `waterfall.png`: the complete synchronized figure described above, saved in
  a headless environment as well as available for interactive use.
- `samples.csv`: one row per timestamp/satellite (12960 rows), with sample index,
  UTC timestamp, NORAD ID, altitude, azimuth, angular separation, range in meters,
  horizon/beam flags, and in-band signal IDs. Use documented lowercase booleans
  and semicolon-separated signal IDs. Order by time, then the configured ID list.
- `passages.csv`: sampled entry, exclusive exit, duration, sample count, clipping
  flags, and signal IDs for each satellite passage; retain a header if empty.
- `metadata.json`: effective configuration, reference hashes, TLE epochs,
  versions, seed, units, model assumptions, array shape, and row/passage counts.

Generated outputs remain ignored, even if their total exceeds the platform's
capture limit. Version source, tests, and small reference inputs rather than
large generated arrays. The actual final output archive must follow the task
platform's requirements. No local size-enforcement hook is needed.

## Tests and extension boundaries

Verify real local-TLE integration and no required network access. Test the beam
and sampling boundaries, a no-passage case, out-of-band signals, independent
carrier placement, additive overlapping contributions, inverse-square scaling,
angular/spectral profile behavior, and deterministic seed behavior. Noise-mean
and variance checks need statistical tolerances, not a demand that every draw
matches an exact mean. Check that changing the seed changes the background but
not geometry or deterministic signal contributions.

Inspect the generated figure as well as numerical arrays. Use shared geometry
and time axes for the mask, upper panel and waterfall. Leave narrow interfaces
for noise and emission profiles so later Doppler, modulation, or receiver models
can be added without a plugin framework or unrelated architecture rewrite.

All parameters must come from the scenario file or an explicit output-path
override. Resolve reference paths relative to the repository root. The JSON,
this document, and `docs/FIRST_PROMPT.md` define a task specification; no noise,
emission, or simulation application code has been implemented in preparation.
