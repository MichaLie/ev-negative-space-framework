"""
Render the workflow schematic for the negative-space evidence-gap mapping framework.

The figure is written to the repository's ``figures/`` directory as both PNG and SVG.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# --- Color palette (matches existing figures) ---
NEURO = "#1f77b4"
TUMOR = "#ff7f0e"
CARDIAC = "#2ca02c"

INPUT_BG = "#e8f0fa"     # light blue
ANALYSIS_BG = "#fff4e6"  # light orange
OUTPUT_BG = "#eaf5ea"    # light green
ROBUST_BG = "#f0f0f0"    # light gray
BORDER = "#3a3a3a"
SUBTLE = "#666666"

# --- Figure ---
fig, ax = plt.subplots(figsize=(15, 12))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# Helper for placing rounded panels
def panel(ax, x, y, w, h, fc, ec=BORDER, lw=1.3):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.4,rounding_size=0.8",
        fc=fc, ec=ec, linewidth=lw,
    )
    ax.add_patch(p)
    return p


def card(ax, x, y, w, h, fc="white", ec=BORDER, lw=0.7):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.25,rounding_size=0.5",
        fc=fc, ec=ec, linewidth=lw,
    )
    ax.add_patch(p)
    return p


def arrow(ax, x1, y1, x2, y2, color=BORDER, lw=1.6):
    a = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=18,
        color=color, linewidth=lw, shrinkA=0, shrinkB=0,
    )
    ax.add_patch(a)
    return a


# ============================================================
# PANEL A — Inputs (Query architecture)
# ============================================================
panel(ax, 1.5, 82, 97, 16.5, INPUT_BG)
ax.text(3.5, 96.5, "A  |  Query architecture",
        ha="left", fontsize=12.5, fontweight="bold")

# Three ingredient cards
# Card 1: Core EV clause
card(ax, 3.5, 84.5, 27, 9.5)
ax.text(17, 92.3, "Core EV clause", ha="center",
        fontsize=10, fontweight="bold")
ax.text(17, 88,
        "(extracellular vesicle  OR\nexosome  OR  exosomes)",
        ha="center", va="center", fontsize=8.3,
        style="italic", linespacing=1.5)

# Card 2: Three disease contexts (color-coded)
card(ax, 32, 84.5, 35, 9.5)
ax.text(49.5, 92.3, "3 disease contexts",
        ha="center", fontsize=10, fontweight="bold")
# Color chips with labels — one row each, well separated
chip_y_top = 90.4
chip_h = 1.4
chip_x = 33.5
chip_w = 2.0
row_step = 2.05
ctx_data = [
    ("Neurodegeneration", NEURO,
     "Alzheimer · Parkinson · ALS · FTD · Huntington"),
    ("Tumor metastasis", TUMOR,
     "cancer/tumor  AND  metastasis"),
    ("Cardiac repair", CARDIAC,
     "cardiac/myocardial injury  AND  repair"),
]
for i, (name, col, query) in enumerate(ctx_data):
    y = chip_y_top - i * row_step
    ax.add_patch(Rectangle((chip_x, y - chip_h/2), chip_w, chip_h,
                           fc=col, ec="none"))
    ax.text(chip_x + chip_w + 0.7, y + 0.05, name,
            va="center", fontsize=8.2, fontweight="bold")
    ax.text(chip_x + chip_w + 12.0, y + 0.05, query,
            va="center", fontsize=7.0, style="italic", color=SUBTLE)

# Card 3: Ten pathway axes
card(ax, 69, 84.5, 27.5, 9.5)
ax.text(82.75, 92.3, "10 pathway axes",
        ha="center", fontsize=10, fontweight="bold")
pathways = ["mTOR", "Wnt", "Notch", "NF-κB", "Complement",
            "Sphingolipid/Ceramide", "Autophagy", "Hypoxia/HIF-1",
            "Integrin/Src", "TGF-β"]
left_col = pathways[:5]
right_col = pathways[5:]
col_y_top = 90.6
col_step = 1.25
for i, p in enumerate(left_col):
    ax.text(71.0, col_y_top - i * col_step, "•  " + p,
            fontsize=7.5, va="center")
for i, p in enumerate(right_col):
    ax.text(83.0, col_y_top - i * col_step, "•  " + p,
            fontsize=7.5, va="center")

# Footnote inside panel A (bottom)
ax.text(50, 83.0,
        "PubMed E-utilities  ·  date window Jan 2015 – Feb 2026  ·  30 pathway-context cells (10 × 3)",
        ha="center", fontsize=8, style="italic", color=SUBTLE)

# Arrow A → B
arrow(ax, 50, 81.7, 50, 78.5)


# ============================================================
# PANEL B — Screening & QC
# ============================================================
panel(ax, 1.5, 64, 97, 13, ANALYSIS_BG)
ax.text(3.5, 75.2, "B  |  Screening & quality control",
        ha="left", fontsize=12.5, fontweight="bold")

# Funnel: 929 → 813 → 386
funnel_y = 70.0
def funnel_stage(ax, x, y, n, label):
    card(ax, x - 7, y - 2.6, 14, 5.0, fc="white")
    ax.text(x, y + 0.5, str(n), ha="center", va="center",
            fontsize=15, fontweight="bold", color="#2a2a2a")
    ax.text(x, y - 1.5, label, ha="center", va="center",
            fontsize=7.5, color=SUBTLE, linespacing=1.25)

funnel_stage(ax, 22, funnel_y, "929", "records retrieved")
arrow(ax, 30, funnel_y, 42, funnel_y)
funnel_stage(ax, 50, funnel_y, "813", "unique PMIDs")
arrow(ax, 58, funnel_y, 70, funnel_y)
funnel_stage(ax, 78, funnel_y, "386", "high-relevance\nprimary records")

ax.text(50, 65.2,
        "Filters: title/abstract keyword match  ·  primary research prioritized  ·  mechanistic-language enrichment",
        ha="center", fontsize=8, style="italic", color=SUBTLE)

# Arrows B → C and B → D (split)
arrow(ax, 28, 63.7, 25, 60.5)
arrow(ax, 72, 63.7, 75, 60.5)


# ============================================================
# PANEL C — Full-text adjudication
# ============================================================
panel(ax, 1.5, 36, 47, 23, ANALYSIS_BG)
ax.text(3.5, 56.5, "C  |  Full-text mechanistic adjudication",
        ha="left", fontsize=11.5, fontweight="bold")
ax.text(25, 53.5, "42 papers  ·  top-5 priority axes",
        ha="center", fontsize=9)

# Per-axis horizontal bars
adjud_axes = [
    ("Integrin/Src", 11),
    ("Autophagy", 10),
    ("TGF-β", 8),
    ("Hypoxia/HIF-1", 7),
    ("Wnt", 6),
]
y_start = 50.5
bar_h = 1.5
max_count = 11
bar_max_w = 17  # in figure units
for i, (name, count) in enumerate(adjud_axes):
    y = y_start - i * 1.95
    w = (count / max_count) * bar_max_w
    # Bar
    ax.add_patch(Rectangle((20.5, y - bar_h/2), w, bar_h,
                           fc=NEURO, ec="none", alpha=0.78))
    # Label left of bar
    ax.text(19.5, y, name, ha="right", va="center", fontsize=8.4)
    # Count right of bar
    ax.text(20.5 + w + 0.7, y, str(count), va="center",
            fontsize=8.4, fontweight="bold")

# Footer note
ax.text(25, 38.7,
        "8 structured criteria / paper  ·  directionality proxy (promoting / protective / mixed)\nconfidence tier:  High ≥7/8  ·  Moderate 5–6",
        ha="center", fontsize=7.4, style="italic",
        color=SUBTLE, linespacing=1.4)


# ============================================================
# PANEL D — Cargo-confidence layer
# ============================================================
panel(ax, 51.5, 36, 47, 23, ANALYSIS_BG)
ax.text(53.5, 56.5, "D  |  Cargo-confidence layer",
        ha="left", fontsize=11.5, fontweight="bold")
ax.text(75, 53.5, "31 candidates  (12 miRNAs + 19 proteins)",
        ha="center", fontsize=9)

# Stacked weighted bar of the 5 components
bar_y = 47.0
bar_h = 4.5
bar_x0 = 54.5
bar_total_w = 41
weights = [0.25, 0.25, 0.25, 0.15, 0.10]
labels = ["Identifier\nvalidity",
          "ExoCarta\nsupport",
          "Cross-context\nbreadth",
          "Mechanistic\nenrichment",
          "Biomarker\nreadiness"]
seg_colors = ["#3a7eb5", "#5b9bcd", "#7eb8d8", "#8fc972", "#b8d99c"]

x_cur = bar_x0
for w, lbl, col in zip(weights, labels, seg_colors):
    seg_w = w * bar_total_w
    ax.add_patch(Rectangle((x_cur, bar_y), seg_w, bar_h,
                           fc=col, ec="white", linewidth=1.5))
    ax.text(x_cur + seg_w/2, bar_y + bar_h/2, f"{w:.2f}",
            ha="center", va="center",
            fontsize=9, fontweight="bold", color="white")
    ax.text(x_cur + seg_w/2, bar_y - 1.9, lbl,
            ha="center", va="top", fontsize=6.8,
            linespacing=1.25, color=SUBTLE)
    x_cur += seg_w

# Tiering note
ax.text(75, 41.0,
        "Composite scaled 0–100  ·  tiered High / Moderate / Low",
        ha="center", fontsize=8.2, fontweight="bold")
ax.text(75, 38.0,
        "Weight-sensitivity check: 5 alternative schemes\n(equal · ExoCarta-heavy · breadth-heavy · mechanistic-heavy · biomarker-heavy)",
        ha="center", fontsize=7.2, style="italic",
        color=SUBTLE, linespacing=1.4)

# Arrows C → E and D → E (merge)
arrow(ax, 25, 35.5, 30, 32.0)
arrow(ax, 75, 35.5, 70, 32.0)


# ============================================================
# PANEL E — Outputs
# ============================================================
panel(ax, 1.5, 11.5, 97, 19, OUTPUT_BG)
ax.text(3.5, 28.5, "E  |  Outputs", ha="left",
        fontsize=12.5, fontweight="bold")

# Four output cards
output_cards = [
    {
        "title": "Priority index",
        "body": "(asymmetry ratio)\n×\n(maturity fraction)",
        "footer": "deterministic formula",
    },
    {
        "title": "Ranked roadmap",
        "body": "#1  Integrin/Src   12.19\n#2  Autophagy     5.61\n#3  TGF-β              3.69\n#4  Wnt                  3.21\n#5  Hypoxia/HIF-1  3.05",
        "footer": "+ predefined safety gates",
    },
    {
        "title": "Clinical-trial gap",
        "body": "444 EV trials\n(Feb 2026 snapshot)\n\n0 interventional trials\ntargeting top-priority\npathways",
        "footer": None,
    },
    {
        "title": "Experimental starter kit",
        "body": "donor anchor →\ntarget anchor → EV source\n→ endpoints\n→ safety gates\n→ falsification criteria",
        "footer": "Supplementary Table S18",
    },
]

card_w = 22.5
card_h = 14.5
gap = 1.5
total_w = 4 * card_w + 3 * gap
x_start = (100 - total_w) / 2
for i, oc in enumerate(output_cards):
    x = x_start + i * (card_w + gap)
    card(ax, x, 12.5, card_w, card_h)
    ax.text(x + card_w/2, 25.5, oc["title"],
            ha="center", fontsize=9.5, fontweight="bold")
    ax.text(x + card_w/2, 18.5, oc["body"],
            ha="center", va="center", fontsize=7.3, linespacing=1.55)
    if oc["footer"]:
        ax.text(x + card_w/2, 13.4, oc["footer"],
                ha="center", fontsize=6.9, style="italic", color=SUBTLE)


# ============================================================
# Robustness band (under B–E)
# ============================================================
panel(ax, 1.5, 0.5, 97, 9.5, ROBUST_BG, ec="#9a9a9a", lw=0.8)
ax.text(50, 8.0, "Robustness analyses cross-cutting layers B–E",
        ha="center", fontsize=10, fontweight="bold")
ax.text(50, 5.4,
        "EV nomenclature (4 clauses)  ·  cardiac-only query relaxation  ·  matched-stringency context tiers",
        ha="center", fontsize=8.2)
ax.text(50, 3.4,
        "maturity-threshold sweep (50–150)  ·  12 alternative index formulations  ·  Integrin/Src axis partition",
        ha="center", fontsize=8.2)
ax.text(50, 1.5,
        "All preserve Integrin/Src at rank #1",
        ha="center", fontsize=8, style="italic", color=SUBTLE)

# ============================================================
# Save
# ============================================================
FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)
png_path = FIGURES / "figure1_workflow.png"
svg_path = FIGURES / "figure1_workflow.svg"
plt.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(svg_path, bbox_inches="tight", facecolor="white")
print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")
