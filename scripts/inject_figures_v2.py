#!/usr/bin/env python3
"""
Inject all figures into manuscript_v3_ARD.docx v2.
Works backwards from end of document to avoid offset shifts.
"""

import os
import shutil

BASE = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
UNPACKED = os.path.join(BASE, "unpacked_ms")
FIGURES_DIR = os.path.join(BASE, "figures")

PAGE_WIDTH_EMU = 5943600

FIGURES_REVERSE = [
    ("Supplementary Figure S3.", "fig6_hub_genes.png", "rId29", "Supplementary Figure S3"),
    ("Supplementary Figure S2.", "fig4_AS_enrichment.png", "rId28", "Supplementary Figure S2"),
    ("Supplementary Figure S1.", "fig3_deg_comparison.png", "rId27", "Supplementary Figure S1"),
    ("Figure 6.", "fig9_sensitivity_analysis.png", "rId26", "Figure 6"),
    ("Figure 5.", "fig8_deconvolution_summary_v2.png", "rId25", "Figure 5"),
    ("Figure 4.", "fig7_wgcna_summary.png", "rId24", "Figure 4"),
    ("Figure 3.", "fig5_ppi_network.png", "rId23", "Figure 3"),
    ("Figure 2.", "fig2_meta_analysis_summary.png", "rId22", "Figure 2"),
    ("Figure 1.", "fig1_volcano_plots.png", "rId21", "Figure 1"),
]


def get_image_dimensions_emu(filepath):
    from PIL import Image
    img = Image.open(filepath)
    w, h = img.width, img.height
    scale = PAGE_WIDTH_EMU / w
    return PAGE_WIDTH_EMU, int(h * scale)


def copy_images():
    media_dir = os.path.join(UNPACKED, "word", "media")
    os.makedirs(media_dir, exist_ok=True)
    for label, filename, rid, alt in FIGURES_REVERSE:
        src = os.path.join(FIGURES_DIR, filename)
        dst = os.path.join(media_dir, filename)
        shutil.copy2(src, dst)
    shutil.copy2(
        os.path.join(FIGURES_DIR, "graphical_abstract.png"),
        os.path.join(media_dir, "graphical_abstract.png")
    )
    print("Copied all images to word/media/")


if __name__ == "__main__":
    print("Step 1: Copying images...")
    copy_images()
    print("Done!")
