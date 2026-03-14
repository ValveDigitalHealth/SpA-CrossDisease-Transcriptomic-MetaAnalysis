#!/usr/bin/env python3
"""
Graphical Abstract: Cross-Disease Transcriptomic Meta-Analysis of SpA
=====================================================================
Pipeline overview diagram showing the analytical workflow and key findings.

Spondyloarthritis AI Computational Biology Institute
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Color palette
TEAL = '#20808D'
DARK_TEAL = '#1B474D'
RUST = '#A84B2F'
GOLD = '#D4A843'
GRAY = '#7A8B8B'
LIGHT_TEAL = '#BCE2E7'
CREAM = '#F7F6F2'
WHITE = '#FFFFFF'
MAROON = '#944454'

fig = plt.figure(figsize=(16, 10), facecolor=WHITE)

# Use full figure for manual layout
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.set_aspect('equal')
ax.axis('off')

# ============================================================
# TITLE
# ============================================================
ax.text(8, 9.6, 'Cross-Disease Transcriptomic Meta-Analysis of Spondyloarthritis',
        fontsize=16, fontweight='bold', ha='center', va='center', color=DARK_TEAL,
        fontfamily='sans-serif')
ax.text(8, 9.25, 'Shared and Disease-Specific Peripheral Blood Signatures Across the SpA Spectrum',
        fontsize=10, ha='center', va='center', color=GRAY, fontfamily='sans-serif')

# ============================================================
# TOP ROW: Data Input
# ============================================================
def draw_box(ax, x, y, w, h, text, subtext, color, text_color='white', fontsize=9, subtextsize=7):
    """Draw a rounded box with centered text."""
    box = FancyBboxPatch((x - w/2, y - h/2), w, h, 
                         boxstyle="round,pad=0.1", 
                         facecolor=color, edgecolor='none', alpha=0.9, zorder=2)
    ax.add_patch(box)
    ax.text(x, y + 0.12, text, fontsize=fontsize, fontweight='bold', ha='center', va='center', 
            color=text_color, zorder=3)
    if subtext:
        ax.text(x, y - 0.15, subtext, fontsize=subtextsize, ha='center', va='center',
                color=text_color, alpha=0.85, zorder=3)

def draw_arrow(ax, x1, y1, x2, y2, color=GRAY):
    """Draw a simple arrow between two points."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.8, 
                               connectionstyle='arc3,rad=0'))

# Row 1: GEO Datasets
row1_y = 8.4
draw_box(ax, 2.0, row1_y, 2.2, 0.7, 'GEO Datasets', '7 datasets, 370 samples', DARK_TEAL)
draw_box(ax, 5.0, row1_y, 2.2, 0.7, 'Disease Groups', 'AS (4), PsA, IBD, jSpA', TEAL)

# Disease labels
diseases = [
    ('AS', 2.5, 7.55, TEAL),
    ('PsA', 4.0, 7.55, RUST),
    ('IBD', 5.5, 7.55, GOLD),
    ('jSpA', 7.0, 7.55, MAROON),
]
for name, x, y, c in diseases:
    pill = FancyBboxPatch((x - 0.45, y - 0.17), 0.9, 0.34, boxstyle="round,pad=0.08",
                          facecolor=c, edgecolor='none', alpha=0.85, zorder=2)
    ax.add_patch(pill)
    ax.text(x, y, name, fontsize=8, fontweight='bold', ha='center', va='center', color=WHITE, zorder=3)

# Arrow from datasets to analysis
draw_arrow(ax, 6.2, row1_y, 8.0, row1_y, GRAY)

# Row 1: Analysis
draw_box(ax, 9.5, row1_y, 2.5, 0.7, 'Differential Expression', 'limma | per-dataset DEG', TEAL)
draw_arrow(ax, 10.8, row1_y, 12.0, row1_y, GRAY)
draw_box(ax, 13.5, row1_y, 2.5, 0.7, 'Meta-Analysis', "Fisher's combined p-value", DARK_TEAL)

# ============================================================
# ROW 2: Core Pipeline (middle)
# ============================================================
row2_y = 6.8

# Central finding box
central_box = FancyBboxPatch((1.5, row2_y - 0.45), 4.5, 0.9,
                             boxstyle="round,pad=0.12",
                             facecolor=RUST, edgecolor=RUST, alpha=0.15, zorder=1)
ax.add_patch(central_box)
central_border = FancyBboxPatch((1.5, row2_y - 0.45), 4.5, 0.9,
                                boxstyle="round,pad=0.12",
                                facecolor='none', edgecolor=RUST, linewidth=2, zorder=2)
ax.add_patch(central_border)

ax.text(3.75, row2_y + 0.15, '40 Meta-Significant AS Genes', 
        fontsize=11, fontweight='bold', ha='center', va='center', color=RUST, zorder=3)
ax.text(3.75, row2_y - 0.18, '35 downregulated  |  5 upregulated  |  RNA processing enrichment',
        fontsize=7.5, ha='center', va='center', color=DARK_TEAL, zorder=3)

# Arrow from meta-analysis down to 40 genes
draw_arrow(ax, 13.5, row1_y - 0.4, 6.0, row2_y + 0.1, GRAY)

# ============================================================
# ROW 3: Downstream Analyses
# ============================================================
row3_y = 5.4
analyses = [
    (2.0, 'Pathway\nEnrichment', 'RNA splicing\nRibosome assembly', TEAL),
    (5.0, 'PPI Network\nAnalysis', 'Hub genes:\nPTBP1, SRSF10', DARK_TEAL),
    (8.0, 'WGCNA', '4 modules per disease\nLow cross-disease\npreservation', GOLD),
    (11.0, 'Immune\nDeconvolution', 'IBD: massive activation\nAS/PsA: subtle shifts', MAROON),
    (14.0, 'Drug\nRepurposing', 'DCK → Cladribine\n(MS repurposing)', RUST),
]

for x, title, sub, color in analyses:
    draw_box(ax, x, row3_y, 2.4, 1.0, title, sub, color, fontsize=9, subtextsize=6.5)

# Arrows from central finding to downstream
for x, _, _, _ in analyses:
    draw_arrow(ax, 3.75, row2_y - 0.5, x, row3_y + 0.55, GRAY)

# ============================================================  
# ROW 4: Key Findings Summary
# ============================================================
row4_y = 3.4

# Big findings banner
banner = FancyBboxPatch((0.5, row4_y - 0.8), 15, 1.6,
                        boxstyle="round,pad=0.15",
                        facecolor=CREAM, edgecolor=DARK_TEAL, linewidth=1.5, zorder=1)
ax.add_patch(banner)

ax.text(8, row4_y + 0.5, 'KEY FINDINGS', fontsize=12, fontweight='bold', 
        ha='center', va='center', color=DARK_TEAL, zorder=2)

findings = [
    (2.5, 'RNA Processing\nCore Signature', 'Splicing machinery\nsystematically\ndownregulated in AS', TEAL),
    (6.0, 'Hub Genes', 'PTBP1, SRSF10,\nEXOSC10, SNIP1, GFM1\nin PPI network', DARK_TEAL),
    (9.5, 'Limited Cross-\nDisease Overlap', 'Distinct modules\nacross AS, IBD, PsA\nsuggesting divergent\npathobiology', GOLD),
    (13.0, 'Cladribine as\nCandidate Drug', 'DCK target identified\nvia network analysis;\napproved for MS', RUST),
]

for x, title, detail, color in findings:
    # Small colored dot
    circle = plt.Circle((x - 1.2, row4_y - 0.05), 0.12, color=color, zorder=3)
    ax.add_patch(circle)
    ax.text(x, row4_y + 0.05, title, fontsize=8, fontweight='bold', ha='center', 
            va='center', color=DARK_TEAL, zorder=3)
    ax.text(x, row4_y - 0.45, detail, fontsize=6.5, ha='center', va='center',
            color=GRAY, zorder=3, linespacing=1.3)

# ============================================================
# ROW 5: Sensitivity & Validation
# ============================================================
row5_y = 1.6

val_box = FancyBboxPatch((0.5, row5_y - 0.6), 15, 1.2,
                         boxstyle="round,pad=0.12",
                         facecolor=LIGHT_TEAL, edgecolor='none', alpha=0.3, zorder=1)
ax.add_patch(val_box)

ax.text(8, row5_y + 0.3, 'VALIDATION & ROBUSTNESS', fontsize=10, fontweight='bold',
        ha='center', va='center', color=DARK_TEAL, zorder=2)

validations = [
    '40-gene signature stable across FDR thresholds (0.01-0.20)',
    'Leave-one-out: genes rank in top 5-10% even with dataset removed',
    'Sensitivity analysis confirms multi-study replication, not single-study bias',
]
for i, v in enumerate(validations):
    ax.text(2.5, row5_y - 0.15 - i * 0.25, f'✓  {v}', fontsize=7.5, ha='left', va='center',
            color=DARK_TEAL, zorder=2)

# Institute attribution
ax.text(8, 0.3, 'Spondyloarthritis AI Computational Biology Institute',
        fontsize=8, ha='center', va='center', color=GRAY, fontstyle='italic')

fig_out = '/home/user/workspace/project1_cross_disease_metaanalysis_spa/figures/graphical_abstract.png'
plt.savefig(fig_out, dpi=300, bbox_inches='tight', facecolor=WHITE)
plt.close()
print(f"Saved: {fig_out}")
