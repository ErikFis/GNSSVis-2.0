"""Headless figures for intermediate simulation results.

Figures are built with ``matplotlib.figure.Figure`` and saved through the Agg
canvas, so nothing opens a window or calls ``show()``.
"""

from __future__ import annotations

from datetime import timedelta, timezone
from pathlib import Path
from typing import Dict

import matplotlib.dates as mdates
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from gnssvis.plot import COLORS_DISTINCT

from .config import Scenario
from .geometry import GeometryResult
from .noise import AveragedWhiteNoise

DB_LABEL = "Relative power [dB re configured mean background]"


def satellite_colors(geometry: GeometryResult) -> Dict[int, str]:
    return {s.norad_id: COLORS_DISTINCT[i % len(COLORS_DISTINCT)]
            for i, s in enumerate(geometry.satellites)}


def _utc_axis(ax, with_date: bool = True) -> None:
    fmt = "%Y-%m-%d\n%H:%M" if with_date else "%H:%M"
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt, tz=timezone.utc))


def _times_num(scenario: Scenario) -> np.ndarray:
    return mdates.date2num(scenario.timestamps_utc)


def _save(fig: Figure, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    FigureCanvasAgg(fig)
    fig.savefig(path, dpi=130)
    return path


def waterfall_extent(scenario: Scenario):
    """imshow extent: column j spans [t_j, t_j + step); rows span channel edges."""
    return [mdates.date2num(scenario.start_utc), mdates.date2num(scenario.end_utc),
            scenario.band_min_hz / 1e6, scenario.band_max_hz / 1e6]


def _hour_axis(ax, with_date: bool = False) -> None:
    ax.xaxis.set_major_locator(mdates.HourLocator(tz=timezone.utc))
    _utc_axis(ax, with_date)


def plot_geometry(scenario: Scenario, geometry: GeometryResult, path: Path) -> Path:
    """Angular separation from the pointing vs UTC, with passage edge detail."""
    from .geometry import geometry_passages

    colors = satellite_colors(geometry)
    half = float(scenario.section("beam")["full_width_deg"]) / 2.0
    t = _times_num(scenario)
    step = scenario.step_seconds

    fig = Figure(figsize=(14, 9.5), layout="constrained")
    top, bottom = fig.subfigures(2, 1, height_ratios=[1.25, 1])
    ax = top.subplots()
    for row, sat in enumerate(geometry.satellites):
        sep = np.where(geometry.above_horizon[row], geometry.separation_deg[row], np.nan)
        c = colors[sat.norad_id]
        ax.plot(t, sep, color=c, lw=1.0, label=f"{sat.name} [{sat.norad_id}]")
        sel = geometry.selected[row]
        if sel.any():
            ax.plot(t[sel], sep[sel], ".", color=c, ms=4)
    ax.axhspan(0, half, color="green", alpha=0.08)
    ax.axhline(half, color="green", ls="--", lw=1, label=f"beam gate {half:g} deg")
    ax.set_ylim(90, 0)
    ax.set_xlim(t[0], mdates.date2num(scenario.end_utc))
    ax.set_ylabel("Angular separation from pointing [deg]")
    ax.set_xlabel(f"Time [UTC], {scenario.start_utc:%Y-%m-%d}")
    ax.set_title("Satellite geometry: above-horizon samples only; dots = in-beam samples")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    _hour_axis(ax)

    passages = geometry_passages(geometry)
    if not passages:
        ax0 = bottom.subplots()
        ax0.text(0.5, 0.5, "No passages in this observation", ha="center",
                 transform=ax0.transAxes)
        ax0.set_axis_off()
        return _save(fig, path)

    p = passages[0]
    row = geometry.norad_ids.index(p.norad_id)
    sat = geometry.satellites[row]
    c = colors[sat.norad_id]
    sep, sel = geometry.separation_deg[row], geometry.selected[row]
    a_in, a_mid, a_out = bottom.subplots(1, 3, width_ratios=[1, 2, 1])
    bottom.suptitle(f"Passage detail: {sat.name}, {step} s sampling, no subsample interpolation"
                    + (f" (first of {len(passages)} passages)" if len(passages) > 1 else ""),
                    fontsize=11)

    def edge(axis, centre, title):
        lo, hi = max(centre - 6, 0), min(centre + 6, t.size)
        sl = slice(lo, hi)
        axis.plot(t[sl], sep[sl], "-", color=c, lw=0.6)
        axis.plot(t[sl][~sel[sl]], sep[sl][~sel[sl]], "o", mfc="white", color=c, ms=6,
                  label="excluded sample")
        axis.plot(t[sl][sel[sl]], sep[sl][sel[sl]], "o", color=c, ms=6, label="included sample")
        axis.axhline(half, color="green", ls="--", lw=1)
        axis.set_title(title, fontsize=10)
        axis.set_ylabel("Separation [deg]")
        axis.grid(alpha=0.3)
        axis.legend(fontsize=8)
        axis.xaxis.set_major_locator(mdates.SecondLocator(bysecond=[0, 30], tz=timezone.utc))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S", tz=timezone.utc))
        axis.tick_params(axis="x", labelrotation=45)

    edge(a_in, p.start_index, f"Entry: first included sample "
         f"{geometry.timestamps_utc[p.start_index]:%H:%M:%S}")
    exit_j = min(p.stop_index, t.size - 1)
    edge(a_out, exit_j, f"Exit: first excluded sample {geometry.timestamps_utc[exit_j]:%H:%M:%S}")

    lo = max(p.start_index - int(300 / step), 0)
    hi = min(p.stop_index + int(300 / step), t.size)
    sl = slice(lo, hi)
    a_mid.axvspan(t[p.start_index], mdates.date2num(geometry.start_utc + timedelta(
        seconds=p.exit_offset_seconds)), color="green", alpha=0.08, label="passage support")
    a_mid.plot(t[sl], sep[sl], color=c, lw=1.5, label="separation")
    a_mid.axhline(half, color="green", ls="--", lw=1)
    a_mid.set_ylabel("Separation [deg]")
    a_mid.grid(alpha=0.3)
    a_r = a_mid.twinx()
    a_r.plot(t[sl], geometry.range_m[row, sl] / 1e3, ":", color="gray", lw=1.5, label="range")
    a_r.set_ylabel("Topocentric range [km]", color="gray")
    h1, l1 = a_mid.get_legend_handles_labels()
    h2, l2 = a_r.get_legend_handles_labels()
    a_mid.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper center")
    a_mid.set_title(f"{p.sample_count} samples = {p.duration_seconds} s in beam", fontsize=10)
    a_mid.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))
    return _save(fig, path)


def plot_power_waterfall(scenario: Scenario, power: np.ndarray, path: Path,
                         title: str) -> Dict[str, float]:
    """Waterfall of 10*log10(power/P0): UTC on X, MHz on Y, full finite dB range."""
    p0 = float(scenario.section("noise")["mean_power_relative"])
    db = (10.0 * np.log10(power.astype(np.float64) / p0)).astype(np.float32)
    finite = db[np.isfinite(db)]
    vmin, vmax = float(finite.min()), float(finite.max())

    fig = Figure(figsize=(13, 7), layout="constrained")
    ax = fig.subplots()
    im = ax.imshow(db, origin="lower", aspect="auto", interpolation="none", cmap="plasma",
                   extent=waterfall_extent(scenario), vmin=vmin, vmax=vmax)
    ax.set_xlabel(f"Time [UTC], {scenario.start_utc:%Y-%m-%d}")
    ax.set_ylabel("Frequency [MHz]")
    ax.set_title(title)
    _hour_axis(ax)
    fig.colorbar(im, ax=ax, orientation="horizontal", label=DB_LABEL, fraction=0.06, pad=0.02)
    _save(fig, path)
    return {"vmin_db": vmin, "vmax_db": vmax}


def plot_noise_statistics(scenario: Scenario, power: np.ndarray, model: AveragedWhiteNoise,
                          path: Path) -> Path:
    """Histograms vs Gamma theory and channel/sample mean flatness."""
    from scipy import stats

    x = power.astype(np.float64)
    k, theta, p0 = model.shape_parameter, model.scale_parameter, model.mean_power
    exp = model.expected_moments()
    fig = Figure(figsize=(13, 8.5), layout="constrained")
    (a1, a2), (a3, a4) = fig.subplots(2, 2)

    bins = np.linspace(0, np.quantile(x, 0.99999), 200)
    a1.hist(x.ravel(), bins=bins, density=True, color="C0", alpha=0.6, label="simulated cells")
    a1.plot(bins, stats.gamma.pdf(bins, k, scale=theta), "k-", lw=1.2,
            label=f"Gamma(M={model.averages}, P0/M={theta:g})")
    a1.plot(bins, stats.norm.pdf(bins, p0, exp["std"]), "r--", lw=1,
            label="Gaussian, same mean/var (not the model)")
    a1.set_yscale("log")
    a1.set_ylim(1e-5, None)
    a1.set_xlabel("Linear relative power")
    a1.set_ylabel("Probability density")
    a1.set_title("Cell power distribution")
    a1.legend(fontsize=8)

    db = 10 * np.log10(x.ravel() / p0)
    dbins = np.linspace(db.min(), db.max(), 200)
    lin = p0 * 10 ** (dbins / 10)
    pdf_db = stats.gamma.pdf(lin, k, scale=theta) * lin * np.log(10) / 10  # change of variables
    a2.hist(db, bins=dbins, density=True, color="C1", alpha=0.6, label="simulated cells")
    a2.plot(dbins, pdf_db, "k-", lw=1.2, label="Gamma theory in dB")
    a2.axvline(exp["mean_db_re_p0"], color="k", ls=":", label=f"E[dB] = {exp['mean_db_re_p0']:.3f}")
    a2.set_xlabel(DB_LABEL)
    a2.set_title("Displayed quantity (dB)")
    a2.legend(fontsize=8)

    f_mhz = scenario.frequencies_hz / 1e6
    n_f, n_t = x.shape
    sd_row = exp["std"] / np.sqrt(n_t)
    a3.plot(f_mhz, x.mean(axis=1), lw=0.8, color="C2", label=f"mean over {n_t} samples")
    a3.fill_between(f_mhz, p0 - 3 * sd_row, p0 + 3 * sd_row, color="gray", alpha=0.3, label="P0 +/- 3 sigma")
    a3.set_xlabel("Frequency [MHz]")
    a3.set_ylabel("Mean power")
    a3.set_title("Per-channel mean (flat spectrum expected)")
    a3.legend(fontsize=8, loc="upper right", framealpha=0.9)

    t = _times_num(scenario)
    sd_col = exp["std"] / np.sqrt(n_f)
    a4.plot(t, x.mean(axis=0), lw=0.6, color="C4", label=f"mean over {n_f} channels")
    a4.fill_between(t, p0 - 3 * sd_col, p0 + 3 * sd_col, color="gray", alpha=0.3, label="P0 +/- 3 sigma")
    a4.set_xlabel(f"Time [UTC], {scenario.start_utc:%Y-%m-%d}")
    a4.set_ylabel("Mean power")
    a4.set_title("Per-sample mean (stationary background expected)")
    a4.legend(fontsize=8)
    _hour_axis(a4)
    fig.suptitle(f"White-noise background statistics: {n_f} x {n_t} cells, "
                 f"seed {scenario.section('randomness')['seed']}")
    return _save(fig, path)


SIGNAL_LABEL = "Signal power [relative linear units; background P0 = 1]"


def plot_signal_preview(scenario: Scenario, geometry: GeometryResult, signal: np.ndarray,
                        contributions, path: Path) -> Path:
    """Signal-only waterfall in linear power, with zoom and profile cuts."""
    from matplotlib import colormaps
    from .emission import carrier_in_band, load_signal_catalog

    f_mhz = scenario.frequencies_hz / 1e6
    t = _times_num(scenario)
    step_days = scenario.step_seconds / 86400.0
    vmax = float(signal.max()) if signal.max() > 0 else 1.0
    cmap = colormaps["viridis"].with_extremes(bad="white")
    shown = np.ma.masked_equal(signal, 0.0)  # exact zeros drawn white, not as a colour bin

    fig = Figure(figsize=(14, 10), layout="constrained")
    top, bottom = fig.subfigures(2, 1, height_ratios=[1, 1.05])
    ax = top.subplots()
    im = ax.imshow(shown, origin="lower", aspect="auto", interpolation="none", cmap=cmap,
                   extent=waterfall_extent(scenario), vmin=0.0, vmax=vmax)
    in_band = [c for c in load_signal_catalog(scenario) if carrier_in_band(c.carrier_hz, scenario)]
    for mhz in sorted({c.carrier_hz / 1e6 for c in in_band}):
        ids = ", ".join(c.signal_id for c in in_band if c.carrier_hz / 1e6 == mhz)
        ax.axhline(mhz, color="gray", ls=":", lw=0.8)
        ax.text(1.002, mhz, f"{ids}\n{mhz:.2f} MHz", transform=ax.get_yaxis_transform(),
                va="center", fontsize=8)
    ax.set_ylabel("Frequency [MHz]")
    ax.set_xlabel(f"Time [UTC], {scenario.start_utc:%Y-%m-%d}")
    ax.set_title("Signal only (no noise), linear power; white = exactly zero")
    _hour_axis(ax)
    top.colorbar(im, ax=ax, orientation="vertical", label=SIGNAL_LABEL, pad=0.08)

    active = [c for c in contributions if c.active]
    if not active:
        a0 = bottom.subplots()
        a0.text(0.5, 0.5, "No active satellite/carrier contributions", ha="center",
                transform=a0.transAxes)
        a0.set_axis_off()
        return _save(fig, path)

    # Zoom around the strongest contribution.
    term = max(active, key=lambda c: c.envelope.max() * c.spectral.max())
    row = geometry.norad_ids.index(term.norad_id)
    sat = geometry.satellites[row]
    cols = np.flatnonzero(term.mask)
    rows = np.flatnonzero(term.support)
    pad_t, pad_f = max(int(300 / scenario.step_seconds), 3), 6
    c0, c1 = max(cols[0] - pad_t, 0), min(cols[-1] + pad_t + 1, t.size)
    r0, r1 = max(rows[0] - pad_f, 0), min(rows[-1] + pad_f + 1, f_mhz.size)
    df = scenario.channel_width_hz / 1e6
    a_zoom, a_spec, a_time = bottom.subplots(1, 3, width_ratios=[1.3, 1, 1.3])
    a_zoom.imshow(shown[r0:r1, c0:c1], origin="lower", aspect="auto", interpolation="none",
                  cmap=cmap, vmin=0.0, vmax=vmax,
                  extent=[t[c0], t[c1 - 1] + step_days, f_mhz[r0] - df / 2, f_mhz[r1 - 1] + df / 2])
    a_zoom.axhline(term.carrier_hz / 1e6, color="gray", ls=":", lw=1)
    a_zoom.set_title(f"Zoom: PRN{sat.prn:02d} {term.signal_id}", fontsize=10)
    a_zoom.set_ylabel("Frequency [MHz]")
    a_zoom.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))

    # Frequency cut at the peak time: channel values vs continuous profile.
    j = int(np.argmax(term.envelope))
    prof = scenario.section("emission")["profiles"][sat.block]
    sigma, cut = prof["spectral_sigma_hz"], prof["spectral_cutoff_sigma"]
    fine = np.linspace(f_mhz[r0], f_mhz[r1 - 1], 800)
    dense = np.exp(-0.5 * ((fine * 1e6 - term.carrier_hz) / sigma) ** 2)
    dense[np.abs(fine * 1e6 - term.carrier_hz) > cut * sigma] = 0.0
    a_spec.plot(fine, term.envelope[j] * dense, "k-", lw=0.8, label="continuous Gaussian (truncated)")
    a_spec.plot(f_mhz[r0:r1], signal[r0:r1, j], "o", color="C3", ms=4, label="channel centres (array)")
    a_spec.axvline(term.carrier_hz / 1e6, color="gray", ls=":", label=f"carrier {term.carrier_hz / 1e6:.2f} MHz")
    for s in (-1, 1):
        a_spec.axvline((term.carrier_hz + s * cut * sigma) / 1e6, color="gray", ls="--", lw=0.7)
    a_spec.set_xlabel("Frequency [MHz]")
    a_spec.set_ylabel("Linear power")
    a_spec.set_title(f"Spectrum at peak sample {geometry.timestamps_utc[j]:%H:%M:%S}", fontsize=10)
    a_spec.legend(fontsize=7)

    # Time cut at the peak channel: masked array vs unmasked envelope.
    i = int(np.argmax(term.spectral))
    unmasked = (float(prof["peak_power_relative_at_reference_range"]) * term.range_scaling
                * term.gain * term.spectral[i])
    sl = slice(c0, c1)
    a_time.axvspan(t[cols[0]], t[cols[-1]] + step_days, color="green", alpha=0.08, label="passage mask")
    a_time.plot(t[sl], unmasked[sl], "--", color="gray", lw=1, label="A*G*F without mask")
    a_time.plot(t[sl], signal[i, sl], "-", color="C3", lw=1.2, drawstyle="steps-post",
                label=f"array, channel {f_mhz[i]:.3f} MHz")
    a_time.set_ylabel("Linear power")
    a_time.set_title("Time profile at peak channel (hard gate, Gaussian taper)", fontsize=10)
    a_time.legend(fontsize=7, loc="lower center")
    a_time.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))
    return _save(fig, path)



def plot_final_waterfall(scenario: Scenario, geometry: GeometryResult, passages,
                         total_power: np.ndarray, path: Path,
                         signal_ids_by_norad: Dict[int, list]) -> Dict[str, float]:
    """GNSSVis-style figure: passages above, combined waterfall below, table right.

    The upper panel shows the angular separation from the configured pointing
    on the same UTC sample grid used for synthesis; the lower panel shows
    10*log10(total_power/P0) with the full finite dB range as colour limits.
    """
    from gnssvis.plot import Visualizer
    from .config import format_utc
    from .emission import carrier_in_band, load_signal_catalog

    p0 = float(scenario.section("noise")["mean_power_relative"])
    db = (10.0 * np.log10(total_power.astype(np.float64) / p0)).astype(np.float32)
    finite = db[np.isfinite(db)]
    vmin, vmax = float(finite.min()), float(finite.max())
    colors = satellite_colors(geometry)
    half = float(scenario.section("beam")["full_width_deg"]) / 2.0
    t = _times_num(scenario)
    t_end = mdates.date2num(scenario.end_utc)

    fig = Figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[3, 1.2], height_ratios=[1, 3],
                          hspace=0.12, wspace=0.28, left=0.06, right=0.98, top=0.92, bottom=0.07)
    ax_wf = fig.add_subplot(gs[1, 0])
    ax_traj = fig.add_subplot(gs[0, 0], sharex=ax_wf)
    ax_table = fig.add_subplot(gs[:, 1])

    # Upper panel: separation from the pointing, above-horizon samples only.
    for row, sat in enumerate(geometry.satellites):
        c = colors[sat.norad_id]
        sep = np.where(geometry.above_horizon[row], geometry.separation_deg[row], np.nan)
        ax_traj.plot(t, sep, color=c, lw=0.8, alpha=0.45)
        inside = np.where(geometry.selected[row], sep, np.nan)
        ax_traj.plot(t, inside, color=c, lw=2.6)
    for p in passages:
        c = colors[p.norad_id]
        x0 = t[p.start_index]
        x1 = mdates.date2num(scenario.start_utc + timedelta(seconds=p.exit_offset_seconds))
        for ax in (ax_traj, ax_wf):
            ax.axvline(x0, color=c, ls="--", lw=0.9, alpha=0.8)
            ax.axvline(x1, color=c, ls="--", lw=0.9, alpha=0.8)
    ax_traj.axhspan(0, half, color="green", alpha=0.07)
    ax_traj.axhline(half, color="green", ls="--", lw=1.0)
    ax_traj.text(0.002, half, f" beam gate {half:g} deg", transform=ax_traj.get_yaxis_transform(),
                 va="bottom", fontsize=8, color="green")
    ax_traj.set_ylim(0, 90)
    ax_traj.set_yticks([0, 15, 30, 45, 60, 75, 90])
    ax_traj.set_ylabel("Separation from\npointing [deg]", fontsize=10)
    point = scenario.section("pointing")
    ax_traj.set_title(f"Satellite passages: angular separation from pointing "
                      f"(alt {point['altitude_deg']:g} deg, az {point['azimuth_deg']:g} deg); "
                      "thick = in beam", fontsize=11, fontweight="bold")
    ax_traj.grid(True, alpha=0.3)
    ax_traj.tick_params(labelbottom=False)

    # Lower panel: combined waterfall in dB re P0.
    im = ax_wf.imshow(db, origin="lower", aspect="auto", interpolation="none", cmap="plasma",
                      extent=waterfall_extent(scenario), vmin=vmin, vmax=vmax)
    in_band = [c for c in load_signal_catalog(scenario) if carrier_in_band(c.carrier_hz, scenario)]
    for hz in sorted({c.carrier_hz for c in in_band}):
        ids = ", ".join(c.signal_id for c in in_band if c.carrier_hz == hz)
        ax_wf.axhline(hz / 1e6, color="white", ls=":", lw=1.2, alpha=0.9)
        ax_wf.text(0.005, hz / 1e6, f"{ids} {hz / 1e6:.2f} MHz", transform=ax_wf.get_yaxis_transform(),
                   va="bottom", ha="left", fontsize=8, color="white",
                   bbox=dict(facecolor="black", alpha=0.45, edgecolor="none", pad=1.5))
    ax_wf.set_xlim(t[0], t_end)
    ax_wf.set_ylim(scenario.band_min_hz / 1e6, scenario.band_max_hz / 1e6)
    ax_wf.set_xlabel(f"Time [UTC], {scenario.start_utc:%Y-%m-%d}", fontsize=10)
    ax_wf.set_ylabel("Frequency [MHz]", fontsize=10)
    ax_wf.set_title(f"Synthetic waterfall: white noise + satellite emission "
                    f"({total_power.shape[0]} channels x {total_power.shape[1]} samples, "
                    f"{scenario.step_seconds} s)", fontsize=11, fontweight="bold")
    ax_wf.xaxis.set_major_locator(mdates.HourLocator(tz=timezone.utc))
    ax_wf.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=[15, 30, 45], tz=timezone.utc))
    ax_wf.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))
    fig.colorbar(im, ax=ax_wf, orientation="horizontal", label=DB_LABEL, pad=0.09,
                 fraction=0.05, aspect=50)

    # Right: satellite / colour table (GNSSVis table helper).
    ax_table.axis("off")
    legend = [(f"{s.name}", colors[s.norad_id]) for s in geometry.satellites]
    by_norad = {}
    for p in passages:
        by_norad.setdefault(p.norad_id, []).append(
            f"{format_utc(scenario.start_utc + timedelta(seconds=p.entry_offset_seconds))[11:19]}-"
            f"{format_utc(scenario.start_utc + timedelta(seconds=p.exit_offset_seconds))[11:19]}")
    passage_col = ["\n".join(by_norad.get(s.norad_id, [])) or "none" for s in geometry.satellites]
    signal_col = [", ".join(i.replace("GPS_", "") for i in signal_ids_by_norad.get(s.norad_id, []))
                  or "-" for s in geometry.satellites]
    Visualizer.draw_object_table(ax_table, legend,
                                 extra_columns=[("Passage UTC", passage_col),
                                                ("Band", signal_col)],
                                 col_widths=[0.43, 0.12, 0.33, 0.12], fontsize=9,
                                 row_scale=2.6, name_width=24, loc="upper center")
    ax_table.set_title("Satellites (configured NORAD IDs)", fontsize=11, fontweight="bold")
    notes = [
        "Upper panel: separation of each satellite from the",
        "pointing at every sample (thin = outside the beam,",
        "thick = in the beam gate and above the horizon).",
        "Dashed verticals: first included sample and first",
        "excluded sample (exclusive exit) of each passage.",
        "Passage UTC: entry - exclusive exit.",
        "Band: in-band carriers injected during passages.",
        "",
        f"Beam gate: separation <= {half:g} deg; gain",
        f"exp(-4 ln2 (sep/{2 * half:g})^2) inside the gate.",
        f"Colour limits: full range {vmin:.2f} .. {vmax:.2f} dB.",
        "No Doppler, modulation, polarization,",
        "receiver or antenna-temperature model.",
    ]
    ax_table.text(0.0, 0.02, "\n".join(notes), transform=ax_table.transAxes, va="bottom",
                  ha="left", fontsize=8.5, family="monospace", color="0.25")

    seed = scenario.section("randomness").get("seed")
    fig.suptitle(f"{scenario.config.get('scenario_id')}: synthetic, not a measured observation "
                 f"(seed {seed}; relative power, not calibrated)", fontsize=12)
    _save(fig, path)
    return {"vmin_db": vmin, "vmax_db": vmax}
