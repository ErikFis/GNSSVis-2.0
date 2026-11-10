# Context for the first Taiga message

This is technical task context to rewrite or adapt before starting. It is not
an evaluation, flag comment, or preference justification. The reviewer brief is
not visible to the model, so repeat the context in your first conversation turn.
Upload the GitHub main-branch source ZIP without `.git`, mounted at `/workspace`.

## English starting message

I wrote GNSSVis to analyze radio telescope observations. I now want a separate
simulation application that uses its geometry and keeps the same visual style.
The result I want is a synthetic waterfall with a white-noise background and
GNSS signal profiles appearing at their nominal frequencies while each satellite
is passing through the telescope beam. Above it, I want the satellite passage
plot; on the right, the satellite/color list, like the existing GNSSVis figure.
Please include dotted carrier-frequency lines and a properly labeled colorbar.

The inputs and observation parameters are already chosen. Read
`examples/scenarios/offline.json`, `docs/OBSERVATION_SCENARIO.md`, and the files in
`reference_data/`. They specify a six-hour historical UIRAPURU zenith observation,
six real GPS TLEs, 10-second samples, 512 frequency channels, and the nominal
L1/L2/L5 frequencies. Only L1 is inside the selected 1500-1600 MHz band. Use the
configuration throughout rather than asking me to choose those values again.

I would like a simple emission model rather than just painting rectangles on
random data: a finite-width Gaussian profile at each in-band carrier, a relative
amplitude that changes with topocentric distance and position in the beam, and a
time mask from the calculated passage. Add contributions in linear power, then
convert to relative dB for display. The specification contains illustrative
coefficients; do not present them as measured GPS transmission powers.

The background should be reproducible synthetic white noise, with no calibrated
antenna or receiver-noise model. Keep the background generator, emission model,
and plotting separate enough to extend later. We are not modelling Doppler,
actual GPS modulation, or absolute received power in this first version.

Please start by reading the existing geometry and plotting code and checking
which parts can be reused. The current SkyAnalyzer initialization and trajectory
loading can trigger downloads, so the new application needs a real local-TLE,
GNSS-only path with no mandatory DE421 file or network access. Preserve the
existing default APIs and avoid a separate duplicate geometry implementation.

Work through this incrementally: offline geometry, seeded background, emission
profiles and passage masks, the combined waterfall, and then the complete figure
and exports. Keep the upper panel and signal mask on exactly the same sampled
time grid. The time resolution is 10 seconds; do not claim exact continuous
crossing times. Use an angular-separation axis consistent with the fixed zenith
pointing rather than blindly copying a fixed-RA/Dec plot.

The proposed command to implement is:

```bash
python -m gnssvis_simulation --config examples/scenarios/offline.json --output outputs/offline
```

It should generate `waterfall.png`, `waterfall.npz`, `samples.csv`, `passages.csv`,
and `metadata.json` as described in the scenario guide. These are deliverables,
not commands or outputs that already exist. Include focused meaningful tests,
a real offline end-to-end example, and documentation of the assumptions and
future extension points. The task environment has no general internet access;
do not run `git init` or `git clone` there. Keep code, comments and documentation
in English without emojis, leave generated outputs ignored, and report the
checks you actually ran. Stop after the agreed result has been verified.

## Matrix convention

The `matrix_contract` section in the shared JSON defines rows, columns and cell
values explicitly. Use `(frequency, time)` with shape `(512, 2160)` for this
scenario. Store frequency centers in Hz, time offsets as integer seconds since
the UTC start, and noise/signal/total cells in linear relative power. In the
figure, display time on X, frequency in MHz on Y, and total power as a relative
dB color scale. Do not interpret the color quantity as a third array dimension.
