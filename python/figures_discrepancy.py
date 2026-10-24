"""Publication figures for the cross-language discrepancy report.

Fig 1: per-quantity relative spread across Octave/Python/Fortran (log scale),
split into citable-reference checks vs observed-only quantities, with orbit
head residuals (absolute metric) at machine level.
Fig 2: aizawa finite-time Lyapunov sensitivity — three panels over d0, x0
perturbation and t_lyap, with the three language values as guide lines.

Run with the workspace venv: .venv/bin/python python/figures_discrepancy.py
Reads results/discrepancy/{per_quantity,sensitivity_aizawa}.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

ROOT = Path(__file__).resolve().parent.parent
DISC = ROOT / "results" / "discrepancy"
FIGS = DISC / "figures"
STYLE = Path(
    "/home/astrid/.local/share/com.ai4s.workbench/runtime/xdg-config/"
    "opencode/skills/publication-figures/openscience.mplstyle"
)

BLUE = "#2a78d6"
AQUA = "#1baf7a"
YELLOW = "#eda100"
INK = "#222222"
GRAY = "#888888"

plt.style.use(str(STYLE) if STYLE.exists() else "default")
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "sans-serif"]


def fig1_spread(rows):
    plot_rows = []
    for r in rows:
        if r["metric"] == "abs_residual":
            plot_rows.append((r, r["spread_abs"], "residual"))
        elif r["metric"] == "rel":
            cls = "observed" if r["observed_only"] else "reference"
            plot_rows.append((r, r["spread"], cls))
    plot_rows.sort(key=lambda t: t[1])
    labels = [f"{r['entry']} · {r['quantity'].replace('_observed', '')}"
              for r, _, _ in plot_rows]
    vals = [max(v, 1e-18) for _, v, _ in plot_rows]
    colors = [{"reference": BLUE, "observed": AQUA, "residual": YELLOW}[c]
              for _, _, c in plot_rows]

    fig, ax = plt.subplots(figsize=(6.4, 8.2))
    y = range(len(vals))
    ax.barh(list(y), vals, color=colors, height=0.72, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlabel("spread across languages (relative; orbit residuals absolute)")
    ax.set_xlim(3e-19, 0.6)
    ax.set_xticks([1e-18, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
    ax.set_xticklabels(["ident.\n(exact)", r"$10^{-16}$", r"$10^{-14}$",
                        r"$10^{-12}$", r"$10^{-10}$", r"$10^{-8}$",
                        r"$10^{-6}$", r"$10^{-4}$", r"$10^{-2}$"], fontsize=7)
    ax.grid(axis="x", alpha=0.25, zorder=0)
    ax.axvline(5e-2, color=GRAY, lw=0.8, ls=":", zorder=2)
    ax.text(5e-2 * 1.15, 0.2, "typical tol", color=GRAY, fontsize=7,
            rotation=90, va="bottom")
    ax.legend(handles=[
        Patch(facecolor=BLUE, label="citable reference"),
        Patch(facecolor=AQUA, label="observed only"),
        Patch(facecolor=YELLOW, label="orbit-head residual (abs)"),
    ], loc="lower right", fontsize=7.5, framealpha=0.9)
    fig.suptitle("Cross-language agreement by quantity (36 checks)", fontsize=9)
    fig.tight_layout()
    out = FIGS / "fig1_spread"
    fig.savefig(f"{out}.png", bbox_inches="tight")
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}.png / .pdf")


def fig2_sensitivity(sens):
    rows = sens["rows"]
    base = sens["baseline_lyap1"]
    lang = sens["lang_values"]
    panels = (
        ("d0", "d0 (initial separation)"),
        ("x0", r"$x_0 + \varepsilon$  ($\varepsilon$ from $10^{-16}$ to $10^{-12}$)"),
        ("t_lyap", r"$t_{\mathrm{lyap}}$  (s)"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
    for ax, (kind, title) in zip(axes, panels):
        sel = [r for r in rows if r["kind"] == kind]
        xs = [r["value"] for r in sel]
        ys = [r["lyap1"] for r in sel]
        ax.semilogx(xs, ys, "o-", color=BLUE, ms=5, lw=1.4, zorder=3)
        ax.axhline(base, color=INK, lw=0.8, zorder=2)
        ax.axhline(lang["fortran"], color=GRAY, lw=0.9, ls="--", zorder=2)
        ax.set_title(title, fontsize=7.5)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"finite-time $\lambda_1$")
    axes[0].set_ylim(0.085, 0.111)
    axes[2].xaxis.set_major_locator(FixedLocator([500, 1000, 2000, 4000, 8000]))
    axes[2].xaxis.set_major_formatter(
        FixedFormatter(["500", "1000", "2000", "4000", "8000"]))
    axes[2].xaxis.set_minor_locator(NullLocator())
    for ax in axes:
        for tick in ax.get_xticklabels():
            tick.set_fontsize(7)
    fig.legend(handles=[
        Line2D([], [], color=INK, lw=0.8, label="Octave ≈ Python (0.10733)"),
        Line2D([], [], color=GRAY, lw=0.9, ls="--", label="Fortran (0.09172)"),
        Line2D([], [], color=BLUE, marker="o", ms=4, lw=1.4,
               label="this experiment (Python)"),
    ], loc="upper right", fontsize=6.5, ncol=3, framealpha=0.9,
        bbox_to_anchor=(0.99, 1.04))
    fig.suptitle("Aizawa: the estimator is invariant to d0 but shifts up to 17% "
                 "for one-ulp x0 changes", fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = FIGS / "fig2_aizawa_sensitivity"
    fig.savefig(f"{out}.png", bbox_inches="tight")
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}.png / .pdf")


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    with open(DISC / "per_quantity.json") as f:
        pq = json.load(f)
    with open(DISC / "sensitivity_aizawa.json") as f:
        sens = json.load(f)
    fig1_spread(pq["rows"])
    fig2_sensitivity(sens)


if __name__ == "__main__":
    main()
