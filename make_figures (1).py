"""Builds the three article figures as vector SVG from lineage_survival_results.json.

Companion code for: Choueiri A. (2026). The Jesus DNA Hypothesis: A Quantitative Reassessment.
Run from the folder containing lineage_survival_results.json (output of lineage_survival_model.py).
Requirements: Python 3.9+, matplotlib; the Inter font is optional.
"""
import json
import math

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

R = json.load(open("lineage_survival_results.json"))

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
ACCENT = "#2a78d6"
GREY = "#b4b2aa"

for f in font_manager.findSystemFonts():
    if "Inter" in f:
        font_manager.fontManager.addfont(f)
plt.rcParams.update({
    "font.family": "Inter",
    "font.size": 8.5,
    "svg.fonttype": "none",
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
W = 6.3  # inches, text-block width


def style(ax):
    ax.grid(axis="both", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, pad=4)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.8)


# ---------- Figure 1: male-line survival by generation ----------
curves = R["survival_curves_gen0_to_80"]
fig, ax = plt.subplots(figsize=(W, 3.15))
style(ax)
gens = list(range(81))
spec = [
    ("growth10_poisson", "10× population growth", GREY, (0, (5, 3)), 1.3),
    ("constant_overdispersed", "High variance in sons", GREY, (0, (1.5, 2)), 1.3),
    ("constant_poisson", "Constant population", ACCENT, "-", 2.0),
]
for key, name, col, ls, lw in spec:
    ys = [p * 100 for p in curves[key]]
    ax.plot(gens, ys, color=col, linestyle=ls, linewidth=lw, solid_capstyle="round")
    ax.annotate(f"{name}, {ys[-1]:.1f}%", (80, ys[-1]), xytext=(6, 0), textcoords="offset points",
                va="center", ha="left", color=INK if col == ACCENT else INK2, fontsize=8)
ax.axvline(65, color=MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
ax.text(64, 82, "2026 at 31-year\ngenerations", ha="right", va="top", color=INK2, fontsize=7.5)
ax.set_yscale("log")
ax.set_ylim(1, 100)
ax.set_xlim(0, 80)
ax.yaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 50, 100]))
ax.yaxis.set_minor_locator(NullLocator())
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
ax.set_xticks(range(0, 81, 10))
ax.set_xlabel("Generations since the founder")
ax.set_ylabel("Male lines still alive (log scale)")
fig.subplots_adjust(left=0.1, right=0.72, top=0.97, bottom=0.14)
fig.savefig("fig1.svg")
plt.close(fig)

# ---------- Figure 2: size of surviving lines vs a 200,000-man line ----------
sc = R["scenarios"]
rows = [
    ("constant_poisson", "Constant population, Poisson"),
    ("constant_overdispersed", "Constant population, high variance"),
    ("growth10_poisson", "10× growth, Poisson"),
    ("growth10_overdispersed", "10× growth, high variance"),
    ("growth30_poisson", "30× growth, Poisson"),
]
fig, ax = plt.subplots(figsize=(W, 2.65))
style(ax)
ax.grid(axis="y", visible=False)
for i, (key, name) in enumerate(rows):
    y = len(rows) - 1 - i
    s = sc[key]["surviving_lines"]
    surv = sc[key]["p_survive_final"] * 100
    ax.plot([s["median"], s["p99"]], [y, y], color=ACCENT, alpha=0.4, linewidth=3.2, solid_capstyle="round")
    ax.plot(s["median"], y, "o", color=ACCENT, markersize=6.5, zorder=4)
    ax.plot([s["max"]] * 2, [y - 0.17, y + 0.17], color=MUTED, linewidth=1.4)
    ax.text(s["max"] * 1.35, y, f"max {s['max']:,}", va="center", color=INK2, fontsize=7.5)
    ax.text(-0.02, y + 0.12, name, transform=ax.get_yaxis_transform(), ha="right", va="center",
            color=INK, fontsize=8)
    ax.text(-0.02, y - 0.2, f"line survives in {surv:.1f}% of runs", transform=ax.get_yaxis_transform(),
            ha="right", va="center", color=INK2, fontsize=7)
ax.axvline(200_000, color=INK2, linewidth=1, linestyle=(0, (3, 2)))
ax.text(200_000 * 0.82, len(rows) - 0.45, "A 200,000-man\nline", ha="right", va="center",
        color=INK, fontsize=7.5)
ax.set_xscale("log")
ax.set_xlim(10, 1_000_000)
ax.set_ylim(-0.6, len(rows) - 0.25)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.xaxis.set_major_locator(FixedLocator([10, 100, 1000, 10_000, 100_000, 1_000_000]))
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: {10:"10",100:"100",1000:"1,000",10_000:"10,000",100_000:"100,000",1_000_000:"1M"}[int(v)]))
ax.set_xlabel("Living male carriers in a surviving line (log scale)")
fig.subplots_adjust(left=0.36, right=0.96, top=0.98, bottom=0.17)
fig.savefig("fig2.svg")
plt.close(fig)

# ---------- Figure 3: known expansions vs the sheet vs neutral range ----------
known = [
    ("Mongol star cluster", "Empire-era expansion · Zerjal et al. 2003", 16_000_000, False),
    ("R-M222 (Niall)", "Irish dynasty · Moore et al. 2006", 3_000_000, False),
    ("Giocangga lineage", "Qing dynasty · Xue et al. 2005", 1_500_000, False),
    ("A 200,000-man line", "The large-line estimate tested here", 200_000, True),
]
med = [sc[k]["surviving_lines"]["median"] for k, _ in rows]
mx = max(sc[k]["surviving_lines"]["max"] for k, _ in rows)


def compact(v):
    return f"{v / 1e6:g}M" if v >= 1e6 else f"{v / 1e3:g}K"


fig, ax = plt.subplots(figsize=(W, 2.55))
style(ax)
ax.grid(axis="y", visible=False)
n = len(known) + 1
for i, (name, sub, v, acc) in enumerate(known):
    y = n - 1 - i
    ax.plot(v, y, "o", color=ACCENT if acc else GREY, markersize=7.5, zorder=4)
    ax.text(v * 1.5, y, compact(v), va="center", color=INK if acc else INK2, fontsize=8,
            fontweight="semibold" if acc else "normal")
    ax.text(-0.02, y + 0.12, name, transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK, fontsize=8)
    ax.text(-0.02, y - 0.2, sub, transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK2, fontsize=7)
lo, hi = min(med), max(med)
ax.plot([lo, hi], [0, 0], color=GREY, linewidth=6, solid_capstyle="round")
ax.plot([mx, mx], [-0.18, 0.18], color=MUTED, linewidth=1.4)
ax.text(mx * 1.4, 0, f"median {lo:g}–{hi:g} (max {mx:,})", va="center", color=INK2, fontsize=7.5)
ax.text(-0.02, 0.12, "Neutral surviving line", transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK, fontsize=8)
ax.text(-0.02, -0.2, "This paper's model, Figure 2", transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK2, fontsize=7)
ax.set_xscale("log")
ax.set_xlim(10, 1e8)
ax.set_ylim(-0.6, n - 0.4)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.xaxis.set_major_locator(FixedLocator([10, 1e3, 1e5, 1e7]))
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: {10: "10", 1e3: "1K", 1e5: "100K", 1e7: "10M"}[v]))
ax.set_xlabel("Estimated living male carriers of a single Y lineage (log scale)")
fig.subplots_adjust(left=0.36, right=0.98, top=0.98, bottom=0.17)
fig.savefig("fig3.svg")
plt.close(fig)
print("ok")


# ---------- Figure 4: illustrative regional split of a surviving line's carriers ----------
regions = [
    ("Lebanon", "Maronite, Orthodox and Druze communities", 30, 40),
    ("Israel and Palestine", "Samaritan and Christian communities", 25, 35),
    ("Syria and Jordan", "Neighbouring populations", 10, 15),
    ("Levantine diaspora", "The Americas, Australia and Europe", 10, 20),
    ("Rest of world", "", 0, 5),
]
fig, ax = plt.subplots(figsize=(W, 2.75))
style(ax)
ax.grid(axis="y", visible=False)
for i, (name, sub, lo_, hi_) in enumerate(regions):
    y = len(regions) - 1 - i
    ax.barh(y, hi_ - lo_, left=lo_, height=0.36, color=GREY, edgecolor="none")
    ax.text(hi_ + 0.8, y, f"{lo_}–{hi_}%" if lo_ else f"<{hi_}%", va="center", color=INK2, fontsize=8)
    ax.text(-0.02, y + (0.12 if sub else 0), name, transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK, fontsize=8)
    if sub:
        ax.text(-0.02, y - 0.2, sub, transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK2, fontsize=7)
    a, b = lo_ * 2000, hi_ * 2000
    c, d = lo_ * 2, hi_ * 2
    ax.text(1.02, y, f"{c}–{d}" if lo_ else f"<{d}", transform=ax.get_yaxis_transform(),
            ha="left", va="center", color=ACCENT, fontsize=8, fontweight="semibold")
    ax.text(1.5, y, f"{a/1000:g}K–{b/1000:g}K" if lo_ else f"<{b/1000:g}K", transform=ax.get_yaxis_transform(),
            ha="left", va="center", color=INK2, fontsize=8)
ax.text(1.02, len(regions) - 0.35, "Men if the line\nhas ~200 (typical)", transform=ax.get_yaxis_transform(), ha="left", va="bottom", color=INK2, fontsize=7)
ax.text(1.5, len(regions) - 0.35, "Men if the line\nhas 200,000", transform=ax.get_yaxis_transform(), ha="left", va="bottom", color=INK2, fontsize=7)
ax.set_xlim(0, 56)
ax.set_ylim(-0.6, len(regions) - 0.3)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks([0, 10, 20, 30, 40, 50])
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
ax.set_xlabel("Assumed share of a surviving line's carriers")
fig.subplots_adjust(left=0.36, right=0.66, top=0.86, bottom=0.17)
fig.savefig("fig4.svg")
plt.close(fig)
print("fig4 ok")

# ---------- Figure 2: intuitive large-line scenarios vs the neutral model ----------
scen = [
    ("Village scenario", "Descendants stay in Galilee villages", 50_000, 100_000, False),
    ("Regional scenario", "Spread across the Levant and diaspora", 200_000, 500_000, True),
    ("Migration scenario", "Strong branching and later migration", 1_000_000, 1_000_000, False),
]
med = [sc[k]["surviving_lines"]["median"] for k, _ in rows]
lo, hi = min(med), max(med)
mx = max(sc[k]["surviving_lines"]["max"] for k, _ in rows)
fig, ax = plt.subplots(figsize=(W, 2.3))
style(ax)
ax.grid(axis="y", visible=False)
n = len(scen) + 1
for i, (name, sub, lo_, hi_, acc) in enumerate(scen):
    y = n - 1 - i
    col = ACCENT if acc else GREY
    if hi_ > lo_:
        ax.plot([lo_, hi_], [y, y], color=col, linewidth=6, solid_capstyle="round")
        lab, x_lab = f"{compact(lo_)}–{compact(hi_)}", hi_ * 1.45
    else:
        ax.plot(lo_, y, "o", color=col, markersize=7.5)
        lab, x_lab = f"~{compact(lo_)}", lo_ * 1.5
    ax.text(x_lab, y, lab, va="center", color=INK if acc else INK2, fontsize=8, fontweight="semibold" if acc else "normal")
    ax.text(-0.02, y + 0.12, name, transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK, fontsize=8)
    ax.text(-0.02, y - 0.2, sub, transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK2, fontsize=7)
ax.plot([lo, hi], [0, 0], color=INK2, linewidth=6, solid_capstyle="round")
ax.plot([mx, mx], [-0.18, 0.18], color=MUTED, linewidth=1.4)
ax.text(mx * 1.4, 0, f"median {lo:g}–{hi:g} (max {mx:,})", va="center", color=INK2, fontsize=7.5)
ax.text(-0.02, 0.12, "Neutral model", transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK, fontsize=8)
ax.text(-0.02, -0.2, "Surviving line, five scenarios", transform=ax.get_yaxis_transform(), ha="right", va="center", color=INK2, fontsize=7)
ax.set_xscale("log")
ax.set_xlim(10, 1e7)
ax.set_ylim(-0.6, n - 0.4)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.xaxis.set_major_locator(FixedLocator([10, 100, 1e3, 1e4, 1e5, 1e6, 1e7]))
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: {10: "10", 100: "100", 1e3: "1K", 1e4: "10K", 1e5: "100K", 1e6: "1M", 1e7: "10M"}[v]))
ax.set_xlabel("Living male carriers of the founder's Y (log scale)")
fig.subplots_adjust(left=0.36, right=0.98, top=0.97, bottom=0.2)
fig.savefig("figS.svg")
plt.close(fig)
print("figS ok")
