#!/usr/bin/env python3
"""
Inject all figures into manuscript_v3_ARD.docx.
Inserts image XML paragraphs before each figure legend in document.xml.
"""

import os
import shutil
import re

BASE = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
UNPACKED = os.path.join(BASE, "unpacked_ms")
FIGURES_DIR = os.path.join(BASE, "figures")

PAGE_WIDTH_EMU = 5943600

FIGURES = [
    ("graphical_abstract", "graphical_abstract.png", "rId20", "Graphical Abstract"),
    ("fig1", "fig1_volcano_plots.png", "rId21", "Figure 1"),
    ("fig2", "fig2_meta_analysis_summary.png", "rId22", "Figure 2"),
    ("fig3", "fig5_ppi_network.png", "rId23", "Figure 3"),
    ("fig4", "fig7_wgcna_summary.png", "rId24", "Figure 4"),
    ("fig5", "fig8_deconvolution_summary_v2.png", "rId25", "Figure 5"),
    ("fig6", "fig9_sensitivity_analysis.png", "rId26", "Figure 6"),
    ("sfig1", "fig3_deg_comparison.png", "rId27", "Supplementary Figure S1"),
    ("sfig2", "fig4_AS_enrichment.png", "rId28", "Supplementary Figure S2"),
    ("sfig3", "fig6_hub_genes.png", "rId29", "Supplementary Figure S3"),
]


def get_image_dimensions_emu(filepath):
    from PIL import Image
    img = Image.open(filepath)
    w, h = img.width, img.height
    scale = PAGE_WIDTH_EMU / w
    new_w = PAGE_WIDTH_EMU
    new_h = int(h * scale)
    return new_w, new_h


def copy_images():
    media_dir = os.path.join(UNPACKED, "word", "media")
    os.makedirs(media_dir, exist_ok=True)
    for fig_id, filename, rid, label in FIGURES:
        src = os.path.join(FIGURES_DIR, filename)
        dst = os.path.join(media_dir, filename)
        shutil.copy2(src, dst)
        print(f"  Copied {filename} -> word/media/{filename}")


if __name__ == "__main__":
    print("Injecting figures into manuscript...")
    copy_images()
    print("Done!")
