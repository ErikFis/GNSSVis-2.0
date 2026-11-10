"""Presentation-quality headless figures using unmodified simulated powers."""

from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .realistic import beam_gain, bpsk_channel_psd

BG = "#091421"
PANEL = "#101f30"
INK = "#e8f0f6"
MUTED = "#a4b6c9"
CYAN = "#51d9d0"
GOLD = "#ffca72"


def style():
    plt.rcParams.update(
        {
            "figure.facecolor": BG,
            "axes.facecolor": PANEL,
            "savefig.facecolor": BG,
            "text.color": INK,
            "axes.labelcolor": MUTED,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.edgecolor": "#34485f",
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titleweight": "bold",
            "axes.titlesize": 15,
            "axes.grid": False,
        }
    )


def plot_realistic_figures(r, output: Path):
    style()
    sc = r.scenario
    g = r.geometry
    hour = sc.time_offsets_seconds / 3600
    f = sc.frequencies_hz / 1e6
    extent = [
        hour[0],
        sc.duration_seconds / 3600,
        sc.band_min_hz / 1e6,
        sc.band_max_hz / 1e6,
    ]
    p0 = float(sc.config["noise"]["mean_power_relative"])
    stamp = sc.start_utc.strftime("%d/%m/%Y %H:%M")
    xlabel = f"Horas desde {stamp} UTC"

    def image(ax, data, title):
        v = 10 * np.log10(np.maximum(data, 1e-10) / p0)
        im = ax.imshow(
            v,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="magma",
            vmin=-2,
            vmax=20,
            interpolation="nearest",
        )
        ax.set(title=title, ylabel="Frequência (MHz)")
        ax.set_yticks([1176.45, 1227.60, 1400, 1575.42])
        ax.set_xticks(np.arange(0, 25, 3))
        return im

    fig = plt.figure(figsize=(15, 9), layout="constrained")
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 3.3, 0.5])
    ax = fig.add_subplot(gs[0])
    below = ax.twinx()
    ax.plot(
        hour,
        r.half_power_mask.sum(axis=0),
        color=CYAN,
        lw=2,
        label="GPS dentro da meia potência",
    )
    ax.fill_between(hour, r.half_power_mask.sum(axis=0), color=CYAN, alpha=0.15)
    below.plot(
        hour, r.solar["separation_deg"], color=GOLD, lw=2, label="Separação do Sol"
    )
    ax.set(
        ylabel="Satélites no feixe",
        xlim=(0, 24),
        ylim=(0, max(4, r.half_power_mask.sum(axis=0).max() + 1)),
    )
    below.set(ylabel="Sol / eixo (°)", ylim=(0, 180))
    below.tick_params(axis="y", labelcolor=GOLD)
    ax.legend(loc="upper left", facecolor=PANEL, edgecolor="none", fontsize=10)
    below.legend(loc="upper right", facecolor=PANEL, edgecolor="none", fontsize=10)
    ax.set_xticks(np.arange(0, 25, 3))
    ax.tick_params(labelbottom=False)
    ax = fig.add_subplot(gs[1])
    im = image(ax, r.total_power, "24 horas · GPS + Sol + ruído")
    ax.set_xlabel(xlabel)
    peak = hour[np.argmax(r.solar_power.sum(axis=0))]
    ax.annotate(
        "Trânsito solar\ncontínuo de banda larga",
        xy=(peak, 1440),
        xytext=(min(peak + 3, 21), 1460),
        color=GOLD,
        fontsize=12,
        arrowprops=dict(arrowstyle="->", color=GOLD),
    )
    for freq, name in [(1176.45, "L5"), (1227.6, "L2"), (1575.42, "L1")]:
        ax.text(0.2, freq + 6, name, color=INK, fontsize=10, weight="bold")
    cax = fig.add_subplot(gs[2])
    cax.axis("off")
    bar = fig.colorbar(
        im, ax=cax, orientation="horizontal", fraction=0.8, aspect=70, pad=0
    )
    bar.set_label(
        "Potência relativa (dB em relação ao fundo médio) · escala fixa −2 a 20 dB"
    )
    fig.suptitle(
        "GNSSVis 2.0 | Geometria e espectro na mesma linha do tempo",
        fontsize=18,
        x=0.06,
        ha="left",
    )
    fig.savefig(output / "01_waterfall_24h.png", dpi=160)
    plt.close(fig)
    fig, axs = plt.subplots(3, 1, figsize=(14, 11), layout="constrained", sharex=True)
    for ax, data, title in zip(
        axs,
        [
            r.noise_power + r.satellite_power,
            r.noise_power + r.solar_power,
            r.total_power,
        ],
        [
            "A · Apenas GPS + fundo",
            "B · Apenas Sol + fundo",
            "C · Superposição em potência linear",
        ],
    ):
        im = image(ax, data, title)
    axs[-1].set_xlabel(xlabel)
    fig.colorbar(
        im,
        ax=axs,
        orientation="horizontal",
        fraction=0.035,
        pad=0.04,
        label="Potência relativa (dB) · mesma escala nos três painéis",
    )
    fig.savefig(output / "02_componentes.png", dpi=150)
    plt.close(fig)
    # L1 detail uses exactly the same data, and never stretches Doppler.
    keep = (f > 1558) & (f < 1590)
    fig, axs = plt.subplots(
        2,
        1,
        figsize=(13, 7),
        layout="constrained",
        gridspec_kw={"height_ratios": [2, 1]},
    )
    ex = [
        0,
        24,
        f[keep][0] - sc.channel_width_hz / 2e6,
        f[keep][-1] + sc.channel_width_hz / 2e6,
    ]
    im = axs[0].imshow(
        10 * np.log10((r.noise_power + r.satellite_power)[keep] / p0),
        origin="lower",
        aspect="auto",
        extent=ex,
        cmap="magma",
        vmin=-2,
        vmax=20,
    )
    axs[0].set(
        title="L1 ampliada · modulação ocupa uma faixa de frequências",
        ylabel="Frequência (MHz)",
    )
    for i in range(len(g.satellites)):
        good = g.above_horizon[i] & (r.half_power_mask[i])
        axs[1].plot(
            hour, np.where(good, r.doppler_hz[i] / 1000, np.nan), lw=1.5, alpha=0.85
        )
    axs[1].axhline(0, color=MUTED, lw=0.6)
    axs[1].set(ylabel="Doppler L1 (kHz)", xlabel=xlabel, xlim=(0, 24))
    axs[1].set_title("Doppler físico durante passagens na meia potência", fontsize=12)
    fig.colorbar(im, ax=axs[0], pad=0.02, label="dB relativos")
    fig.savefig(output / "03_l1_doppler.png", dpi=160)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.4), layout="constrained")
    fx = np.linspace(-25e6, 25e6, 2001)
    for rate, color, label in [
        (1.023e6, CYAN, "BPSK · 1,023 Mchip/s"),
        (10.23e6, GOLD, "BPSK · 10,23 Mchip/s"),
    ]:
        val = bpsk_channel_psd(fx, 50000, 0, rate)
        axs[0].plot(
            fx / 1e6,
            10 * np.log10(np.maximum(val, 1e-6)),
            color=color,
            label=label,
            lw=1.7,
        )
    axs[0].set(
        title="Uma portadora, vários canais",
        xlabel="Desvio da portadora (MHz)",
        ylabel="PSD normalizada (dB)",
        ylim=(-45, 1),
        xlim=(-25, 25),
    )
    axs[0].legend(facecolor=PANEL, edgecolor="none", fontsize=9)
    beam = sc.section("beam")
    angles = np.linspace(-90, 90, 500)
    for freq, color, label in [(1176.45e6, GOLD, "L5"), (1575.42e6, CYAN, "L1")]:
        axs[1].plot(
            angles, beam_gain(angles, 0, freq, beam), lw=2, color=color, label=label
        )
    axs[1].axhline(0.5, color=MUTED, lw=0.8, ls="--")
    axs[1].set(
        title="Feixe suave e cromático",
        xlabel="Desvio horizontal do eixo (°)",
        ylabel="Ganho de potência",
        ylim=(0, 1.05),
        xlim=(-90, 90),
    )
    axs[1].legend(facecolor=PANEL, edgecolor="none")
    fig.savefig(output / "04_modelos.png", dpi=160)
    plt.close(fig)
    # Compare the supplied separable model with the new components on ONE grid.
    from .emission import compute_signal_power

    legacy, _ = compute_signal_power(sc, g)
    fig, axs = plt.subplots(2, 1, figsize=(14, 8), layout="constrained", sharex=True)
    image(
        axs[0],
        r.noise_power + legacy,
        "Modelo anterior · perfis gaussianos fixos + corte geométrico",
    )
    im = image(
        axs[1],
        r.total_power,
        "Novo cenário · espectro espalhado + Doppler + feixe suave + Sol",
    )
    axs[-1].set_xlabel(xlabel)
    fig.colorbar(
        im,
        ax=axs,
        orientation="horizontal",
        fraction=0.035,
        pad=0.04,
        label="Potência relativa (dB) · grades e escalas iguais; amplitudes ilustrativas",
    )
    fig.savefig(output / "05_comparacao.png", dpi=150)
    plt.close(fig)
