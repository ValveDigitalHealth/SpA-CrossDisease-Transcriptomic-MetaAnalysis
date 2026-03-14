#!/usr/bin/env python3
"""
Inject all figures into manuscript_v3_ARD.docx — v2.
Strategy: work backwards from end of document to avoid offset shifts.
Find each figure legend paragraph by its exact bold label text,
then insert the image paragraph directly before it.
"""

import os
import shutil
import re

BASE = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
UNPACKED = os.path.join(BASE, "unpacked_ms")
FIGURES_DIR = os.path.join(BASE, "figures")

PAGE_WIDTH_EMU = 5943600  # 6.5 inches in EMU

# Figure mapping: (search_label, source_file, rId, alt_text)
# IMPORTANT: Order from LAST to FIRST so we work backwards and don't shift earlier positions
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
    new_w = PAGE_WIDTH_EMU
    new_h = int(h * scale)
    return new_w, new_h

def make_image_paragraph(rid, cx, cy, fig_id, description):
    doc_pr_id = int(rid.replace("rId", "")) + 100
    return f'''    <w:p>
      <w:pPr>
        <w:spacing w:before="240" w:after="120"/>
        <w:jc w:val="center"/>
      </w:pPr>
      <w:r>
        <w:drawing>
          <wp:inline distT="0" distB="0" distL="0" distR="0">
            <wp:extent cx="{cx}" cy="{cy}"/>
            <wp:effectExtent l="0" t="0" r="0" b="0"/>
            <wp:docPr id="{doc_pr_id}" name="{fig_id}" descr="{description}"/>
            <a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
              <a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
                <pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">
                  <pic:nvPicPr>
                    <pic:cNvPr id="{doc_pr_id}" name="{fig_id}"/>
                    <pic:cNvPicPr/>
                  </pic:nvPicPr>
                  <pic:blipFill>
                    <a:blip r:embed="{rid}"/>
                    <a:stretch>
                      <a:fillRect/>
                    </a:stretch>
                  </pic:blipFill>
                  <pic:spPr>
                    <a:xfrm>
                      <a:off x="0" y="0"/>
                      <a:ext cx="{cx}" cy="{cy}"/>
                    </a:xfrm>
                    <a:prstGeom prst="rect">
                      <a:avLst/>
                    </a:prstGeom>
                  </pic:spPr>
                </pic:pic>
              </a:graphicData>
            </a:graphic>
          </wp:inline>
        </w:drawing>
      </w:r>
    </w:p>
'''

def copy_images():
    media_dir = os.path.join(UNPACKED, "word", "media")
    os.makedirs(media_dir, exist_ok=True)
    for label, filename, rid, alt in FIGURES_REVERSE:
        src = os.path.join(FIGURES_DIR, filename)
        dst = os.path.join(media_dir, filename)
        shutil.copy2(src, dst)
    # Also copy graphical abstract
    shutil.copy2(
        os.path.join(FIGURES_DIR, "graphical_abstract.png"),
        os.path.join(media_dir, "graphical_abstract.png")
    )
    print("  Copied all images to word/media/")

def add_relationships():
    rels_path = os.path.join(UNPACKED, "word", "_rels", "document.xml.rels")
    with open(rels_path, "r") as f:
        content = f.read()
    
    new_rels = ""
    # Add graphical abstract
    new_rels += '  <Relationship Id="rId20" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/graphical_abstract.png"/>\n'
    for label, filename, rid, alt in FIGURES_REVERSE:
        new_rels += f'  <Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{filename}"/>\n'
    
    content = content.replace("</Relationships>", new_rels + "</Relationships>")
    
    with open(rels_path, "w") as f:
        f.write(content)
    print("  Updated document.xml.rels")

def inject_figures():
    doc_path = os.path.join(UNPACKED, "word", "document.xml")
    with open(doc_path, "r") as f:
        lines = f.readlines()
    
    # First, find the "FIGURE LEGENDS" heading line
    fl_line = None
    for i, line in enumerate(lines):
        if "FIGURE LEGENDS" in line and "SUPPLEMENTARY" not in line:
            fl_line = i
            break
    
    if fl_line is None:
        print("  ERROR: Could not find FIGURE LEGENDS section")
        return
    
    print(f"  Found FIGURE LEGENDS at line {fl_line}")
    
    # Only search within the FIGURE LEGENDS section (from fl_line to end)
    # Find each figure label and its containing <w:p> start
    
    # Process figures in reverse order to maintain line offsets
    for label, filename, rid, alt in FIGURES_REVERSE:
        filepath = os.path.join(FIGURES_DIR, filename)
        cx, cy = get_image_dimensions_emu(filepath)
        img_xml = make_image_paragraph(rid, cx, cy, alt.replace(" ", "_").lower(), alt)
        
        # Find the line containing this label text AFTER fl_line
        target_line = None
        for i in range(fl_line, len(lines)):
            if f">{label}<" in lines[i] or f">{label}</w:t>" in lines[i]:
                target_line = i
                break
        
        if target_line is None:
            print(f"  WARNING: Could not find '{label}' after line {fl_line}")
            continue
        
        # Find the <w:p> that contains this label - walk backwards to find its opening
        p_start = None
        for i in range(target_line, fl_line - 1, -1):
            if "<w:p>" in lines[i]:
                p_start = i
                break
        
        if p_start is None:
            print(f"  WARNING: Could not find <w:p> for '{label}'")
            continue
        
        # Insert image paragraph before this <w:p>
        img_lines = img_xml.split('\n')
        # Filter empty lines and add newline
        img_lines = [l + '\n' for l in img_lines if l.strip()]
        
        for j, img_line in enumerate(img_lines):
            lines.insert(p_start + j, img_line)
        
        print(f"  Inserted {alt} image before line {p_start} ({cx}x{cy} EMU)")
        
        # Update fl_line offset since we inserted lines above other positions
        # Actually, since we process in reverse, all future targets are at higher indices
        # and we insert before them, so no offset needed for earlier figures
    
    # Now add graphical abstract image + legend after "FIGURE LEGENDS" heading
    # Find the updated position of FIGURE LEGENDS
    for i, line in enumerate(lines):
        if "FIGURE LEGENDS" in line and "SUPPLEMENTARY" not in line:
            fl_line = i
            break
    
    # Find the end of the FIGURE LEGENDS heading paragraph (</w:p>)
    fl_p_end = None
    for i in range(fl_line, min(fl_line + 10, len(lines))):
        if "</w:p>" in lines[i]:
            fl_p_end = i
            break
    
    if fl_p_end is not None:
        # Create graphical abstract image + legend
        ga_filepath = os.path.join(FIGURES_DIR, "graphical_abstract.png")
        ga_cx, ga_cy = get_image_dimensions_emu(ga_filepath)
        ga_img_xml = make_image_paragraph("rId20", ga_cx, ga_cy, "graphical_abstract", "Graphical Abstract")
        
        ga_legend_xml = '''    <w:p>
      <w:pPr>
        <w:spacing w:after="120" w:before="240" w:line="480"/>
      </w:pPr>
      <w:r>
        <w:rPr>
          <w:rFonts w:ascii="Calibri" w:cs="Calibri" w:eastAsia="Calibri" w:hAnsi="Calibri"/>
          <w:b/>
          <w:bCs/>
          <w:sz w:val="24"/>
          <w:szCs w:val="24"/>
        </w:rPr>
        <w:t xml:space="preserve">Graphical Abstract.</w:t>
      </w:r>
      <w:r>
        <w:rPr>
          <w:rFonts w:ascii="Calibri" w:cs="Calibri" w:eastAsia="Calibri" w:hAnsi="Calibri"/>
          <w:sz w:val="24"/>
          <w:szCs w:val="24"/>
        </w:rPr>
        <w:t xml:space="preserve"> Overview of the cross-disease transcriptomic meta-analysis pipeline, from dataset curation through differential expression, meta-analysis, network analysis, immune deconvolution, and sensitivity analysis.</w:t>
      </w:r>
    </w:p>
'''
        
        combined = ga_img_xml + ga_legend_xml
        insert_lines = [l + '\n' for l in combined.split('\n') if l.strip()]
        
        insert_pos = fl_p_end + 1
        for j, line in enumerate(insert_lines):
            lines.insert(insert_pos + j, line)
        
        print(f"  Inserted Graphical Abstract after FIGURE LEGENDS heading")
    
    with open(doc_path, "w") as f:
        f.writelines(lines)
    print("  Updated document.xml")

if __name__ == "__main__":
    print("Step 1: Copying images...")
    copy_images()
    
    print("\nStep 2: Adding relationships...")
    add_relationships()
    
    print("\nStep 3: Injecting figures...")
    inject_figures()
    
    print("\nDone!")
