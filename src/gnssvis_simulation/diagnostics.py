"""Statistical checks and plain-text summaries for intermediate results."""

from __future__ import annotations

from datetime import timedelta
from typing import Dict, List

import numpy as np

from .config import format_utc
from .geometry import GeometryResult, Passage
from .noise import AveragedWhiteNoise

QUANTILES = (0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)


def _batch_se(values_per_batch: np.ndarray) -> float:
    return float(np.std(values_per_batch, ddof=1) / np.sqrt(values_per_batch.size))


def noise_statistics(power: np.ndarray, model: AveragedWhiteNoise,
                     batches: int = 64) -> Dict[str, object]:
    """Compare a (frequency, time) noise array with Gamma(M, P0/M) theory.

    Each scalar check has ``measured``, ``expected``, a standard error ``se``
    and ``z = (measured - expected) / se``. Mean/variance use analytic errors;
    skewness/kurtosis use batch means over time blocks.
    """
    from scipy import stats

    x = np.asarray(power, dtype=np.float64)
    n_f, n_t = x.shape
    n = x.size
    exp = model.expected_moments()
    sigma = exp["std"]

    def check(measured, expected, se):
        return {"measured": float(measured), "expected": float(expected), "se": float(se),
                "z": float((measured - expected) / se) if se > 0 else float("nan")}

    flat = x.ravel()
    blocks = np.array_split(x, batches, axis=1)
    checks = {
        "mean": check(flat.mean(), exp["mean"], sigma / np.sqrt(n)),
        "variance": check(flat.var(ddof=1), exp["variance"],
                          sigma ** 2 * np.sqrt((exp["excess_kurtosis"] + 2.0) / n)),
        "skewness": check(stats.skew(flat), exp["skewness"],
                          _batch_se(np.array([stats.skew(b.ravel()) for b in blocks]))),
        "excess_kurtosis": check(stats.kurtosis(flat), exp["excess_kurtosis"],
                                 _batch_se(np.array([stats.kurtosis(b.ravel()) for b in blocks]))),
    }
    db = 10.0 * np.log10(flat / model.mean_power)
    checks["mean_db_re_p0"] = check(db.mean(), exp["mean_db_re_p0"], exp["std_db"] / np.sqrt(n))

    # Per-channel (row) and per-sample (column) means: flat background?
    row_means, col_means = x.mean(axis=1), x.mean(axis=0)
    row_sd, col_sd = sigma / np.sqrt(n_t), sigma / np.sqrt(n_f)

    def lag1(a, axis):
        a0 = np.moveaxis(a, axis, -1)
        u, v = a0[..., :-1].ravel(), a0[..., 1:].ravel()
        return float(np.corrcoef(u, v)[0, 1]), 1.0 / np.sqrt(u.size)

    rho_t, se_t = lag1(x, 1)
    rho_f, se_f = lag1(x, 0)
    ks = stats.kstest(flat, "gamma", args=(model.shape_parameter, 0.0, model.scale_parameter))
    quant = np.quantile(flat, QUANTILES)
    quant_exp = stats.gamma.ppf(QUANTILES, model.shape_parameter, scale=model.scale_parameter)
    return {
        "n_cells": int(n),
        "checks": checks,
        "channel_means": {"sd_measured": float(row_means.std(ddof=1)), "sd_expected": row_sd,
                          "max_abs_z": float(np.max(np.abs(row_means - exp["mean"])) / row_sd)},
        "sample_means": {"sd_measured": float(col_means.std(ddof=1)), "sd_expected": col_sd,
                         "max_abs_z": float(np.max(np.abs(col_means - exp["mean"])) / col_sd)},
        "lag1_correlation": {"time": check(rho_t, 0.0, se_t), "frequency": check(rho_f, 0.0, se_f)},
        "ks_test": {"statistic": float(ks.statistic), "pvalue": float(ks.pvalue)},
        "quantiles": {f"{q:g}": {"measured": float(m), "expected": float(e)}
                      for q, m, e in zip(QUANTILES, quant, quant_exp)},
        "min": float(flat.min()), "max": float(flat.max()),
        "db_min": float(db.min()), "db_max": float(db.max()),
    }


# -- text rendering -------------------------------------------------------

def format_table(headers: List[str], rows: List[List[object]], numeric: bool = True) -> str:
    """Plain-text table; first column left-aligned, others right-aligned if numeric."""
    cells = [[str(h) for h in headers]] + [[str(c) for c in r] for r in rows]
    widths = [max(len(r[i]) for r in cells) for i in range(len(headers))]
    line = "  ".join("-" * w for w in widths)
    out = ["  ".join(c.ljust(w) for c, w in zip(cells[0], widths)), line]
    out += ["  ".join(c.ljust(w) if i == 0 or not numeric else c.rjust(w)
                      for i, (c, w) in enumerate(zip(r, widths))) for r in cells[1:]]
    return "\n".join(out)


def ascii_timeline(geometry: GeometryResult, width: int = 72) -> str:
    """One row per satellite: '#' any in-beam sample, '-' above horizon, ' ' below."""
    n_t = geometry.selected.shape[1]
    edges = np.linspace(0, n_t, width + 1).astype(int)
    label_w = max(len(f"{s.norad_id} PRN{s.prn:02d}") for s in geometry.satellites)
    lines = []
    hours = (n_t * geometry.step_seconds) / 3600.0
    ruler = [" "] * width
    for h in range(int(hours) + 1):
        col = min(int(round(h / hours * width)), width - 1) if hours else 0
        tag = (geometry.start_utc + timedelta(hours=h)).strftime("%H")
        if col + len(tag) <= width:
            ruler[col:col + len(tag)] = tag
    lines.append(" " * label_w + " |" + "".join(ruler) + "| UTC hour")
    for row, sat in enumerate(geometry.satellites):
        chars = []
        for a, b in zip(edges[:-1], edges[1:]):
            if geometry.selected[row, a:b].any():
                chars.append("#")
            elif geometry.above_horizon[row, a:b].any():
                chars.append("-")
            else:
                chars.append(" ")
        lines.append(f"{sat.norad_id} PRN{sat.prn:02d}".ljust(label_w) + " |" + "".join(chars) + "|")
    minutes = n_t * geometry.step_seconds / width / 60.0
    lines.append(f"legend: '#' in beam (passage), '-' above horizon, ' ' below horizon; "
                 f"1 column = {minutes:g} min")
    return "\n".join(lines)


def passage_rows(geometry: GeometryResult, passages: List[Passage]) -> List[List[object]]:
    rows = []
    index = {n: i for i, n in enumerate(geometry.norad_ids)}
    for p in passages:
        r = index[p.norad_id]
        sl = slice(p.start_index, p.stop_index)
        j = p.start_index + int(np.argmin(geometry.separation_deg[r, sl]))
        sat = geometry.satellites[r]
        rows.append([
            p.norad_id, f"PRN{sat.prn:02d}", sat.block,
            format_utc(geometry.start_utc + timedelta(seconds=p.entry_offset_seconds)),
            format_utc(geometry.start_utc + timedelta(seconds=p.exit_offset_seconds)),
            p.sample_count, p.duration_seconds,
            f"{geometry.separation_deg[r, j]:.3f}",
            format_utc(geometry.start_utc + timedelta(seconds=int(geometry.time_offsets_seconds[j]))),
            f"{geometry.range_m[r, j] / 1e3:.0f}",
            f"{p.entry_clipped}/{p.exit_clipped}".lower(),
        ])
    return rows


# -- emission ---------------------------------------------------------------

def contribution_summary(term, scenario, geometry: GeometryResult) -> Dict[str, object]:
    """Intervals and peak of one satellite/carrier contribution.

    Intervals come from the contribution's own nonzero samples, so they can be
    compared with the geometric passages.
    """
    from .geometry import find_passages

    row = geometry.norad_ids.index(term.norad_id)
    sat = geometry.satellites[row]
    power_t = term.envelope * (term.spectral.max() if term.support.any() else 0.0)
    intervals = find_passages(term.envelope > 0, scenario.step_seconds, term.norad_id)
    summary = {
        "norad_id": term.norad_id, "prn": sat.prn, "block": sat.block,
        "signal_id": term.signal_id, "carrier_hz": term.carrier_hz,
        "active": term.active, "mask_samples": int(term.mask.sum()),
        "support_channels": int(term.support.sum()),
        "intervals": [{
            "entry_utc": format_utc(scenario.start_utc + timedelta(seconds=p.entry_offset_seconds)),
            "exit_utc_exclusive": format_utc(scenario.start_utc + timedelta(seconds=p.exit_offset_seconds)),
            "samples": p.sample_count, "duration_seconds": p.duration_seconds,
            "entry_value": float(power_t[p.start_index]),
            "last_value": float(power_t[p.stop_index - 1]),
        } for p in intervals],
    }
    if term.active:
        j = int(np.argmax(term.envelope))
        i = int(np.argmax(term.spectral))
        summary["peak"] = {
            "value": float(term.spectral[i] * term.envelope[j]),
            "utc": format_utc(scenario.timestamps_utc[j]),
            "time_index": j, "channel_index": i,
            "channel_center_hz": float(scenario.frequencies_hz[i]),
            "range_m": float(geometry.range_m[row, j]),
            "separation_deg": float(geometry.separation_deg[row, j]),
            "range_scaling": float(term.range_scaling[j]),
            "angular_gain": float(term.gain[j]),
            "spectral_weight": float(term.spectral[i]),
        }
    return summary
