# Starting context and follow-up prompts

These are technical task messages, not behavior evaluations, flag comments or
comparison justifications. The reviewer brief is not visible to the model.
The shared specification remains `examples/scenarios/offline.json`, explained
in `docs/OBSERVATION_SCENARIO.md`. The application is to be implemented during
the task, using the existing GNSSVis library.

## First message

I built GNSSVis to help identify possible satellite interference in radio telescope observations. I now want to add a separate simulation pipeline, because real observations are scarce and can be very large. The result I have in mind is a synthetic waterfall with background noise and satellite signals appearing during their beam passages, using the same plot layout as GNSSVis.

I have included the TLEs, signal frequencies and a reference scenario in examples/scenarios/offline.json. The time and frequency resolution should come from that file so I can change them later. For this example, the power arrays will be float32 with 512 frequency rows in Hz and 2160 time columns, 10 seconds apart. Time offsets are integer seconds from the UTC start, and each cell contains linear relative power.

The planned signal model is a Gaussian frequency profile with inverse-square distance scaling and a Gaussian angular response inside the beam. I want to build this in stages, starting with the geometry. Could you look at the existing SGP4 code and make it work with the local TLEs, without needing DE421 or network downloads? I would like to keep the current GNSSVis interfaces working. Let's check this part before adding the noise and emission. Please do not run git init or git clone in this workspace.

## Background noise

Now I would like to generate the background noise on the configured frequency and time grid. Since the matrix represents averaged power, let's use Gamma(shape=M, scale=P0/M), where P0 is the configured mean power and M is the number of averages. I want to use NumPy default_rng with the saved seed so I can reproduce it. The output should be float32 in linear relative power, with frequency rows in Hz and time columns stored as integer seconds from the UTC start.

## Satellite emission

Next I want to add the satellite emission. For each carrier inside the band, let's use a Gaussian profile centered on its nominal frequency. Its amplitude should follow P_ref * (r_ref/r)**2, with the Gaussian angular gain from the scenario. The passage mask should turn the contribution off outside the beam. I would like the widths and coefficients to come from the configuration, and signals at the same frequency to add together rather than replace each other.

## Waterfall and application outputs

Now I would like to see the combined waterfall. Let's add the noise and signal contributions in linear power, then plot 10*log10(total_power/P0). I want UTC time on X, frequency in MHz on Y, and relative dB as the color scale, with dotted carrier-frequency lines.

The layout should look like GNSSVis, with the passage plot above and the satellite/color list on the right. For our zenith pointing, the upper panel should show angular separation. Could you connect this to python -m gnssvis_simulation --config examples/scenarios/offline.json and save the PNG, NPZ, tables and metadata in outputs/offline? The command should finish without waiting for me to close a plot window.

## End-to-end verification

I would like to check the full example now, including the saved figure and arrays. Can you verify that the signal appears at the right frequency and passage times, stays zero outside the mask, and follows the distance and beam response? I also want to check overlapping signals and out-of-band carriers. With a different seed, only the background should change. Please run the offline example and tests, show me the results, and fix anything that still fails.

## Documentation and completion

To finish, I would like clear instructions for reproducing the example and a short explanation of the model's assumptions and limits. The metadata should include the configuration, seed and dependency versions. I would also like the noise and emission functions to be easy to extend later, without adding those future models now. Could you check that the agreed outputs and tests are complete and summarize what has actually been verified?

## Units and scope

The reference arrays have `(frequency, time)` shape `(512, 2160)`, but resolution
must remain configurable through the shared scenario. Frequency centers are in
Hz; time offsets are integer seconds from the UTC start; power cells are float32
in relative linear units. The plot displays UTC time on X, MHz on Y, and relative
dB as color. Future pulsar/FRB models are outside the present implementation.
No simulation implementation is included as part of preparing these messages.
