#!/usr/bin/env python3
"""
Inject all figures into manuscript_v3_ARD.docx.
Strategy: 
  - Copy image files to word/media/
  - Add relationships to document.xml.rels
  - Insert image XML paragraphs before each figure legend in document.xml
  - Add graphical abstract at the very beginning (before abstract)
"""

import os
import shutil
import re

BASE = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
UNPACKED = os.path.join(BASE, "unpacked_ms")
FIGURES_DIR = os.path.join(BASE, "figures")

# Figure mapping: (manuscript label, source file, rId)
# Page content width = 12240 - 1440 - 1440 = 9360 DXA = 6.5 inches = 5943600 EMU
PAGE_WIDTH_EMU = 5943600

FIGURES = [
    # Main figures (appear before their legend text)
    ("graphical_abstract", "graphical_abstract.png", "rId20", "Graphical Abstract"),
    ("fig1", "fig1_volcano_plots.png", "rId21", "Figure 1"),
    ("fig2", "fig2_meta_analysis_summary.png", "rId22", "Figure 2"),
    ("fig3", "fig5_ppi_network.png", "rId23", "Figure 3"),  # PPI = manuscript Fig 3
    ("fig4", "fig7_wgcna_summary.png", "rId24", "Figure 4"),  # WGCNA = manuscript Fig 4
    ("fig5", "fig8_deconvolution_summary_v2.png", "rId25", "Figure 5"),  # Deconv = manuscript Fig 5
    ("fig6", "fig9_sensitivity_analysis.png", "rId26", "Figure 6"),  # Sensitivity = manuscript Fig 6
    # Supplementary figures
    ("sfig1", "fig3_deg_comparison.png", "rId27", "Supplementary Figure S1"),
    ("sfig2", "fig4_AS_enrichment.png", "rId28", "Supplementary Figure S2"),
    ("sfig3", "fig6_hub_genes.png", "rId29", "Supplementary Figure S3"),
]

def get_image_dimensions_emu(filepath):
    """Get image dimensions scaled to fit page width, in EMUs."""
    from PIL import Image
    img = Image.open(filepath)
    w, h = img.width, img.height
    # Scale to page width
    scale = PAGE_WIDTH_EMU / w
    new_w = PAGE_WIDTH_EMU
    new_h = int(h * scale)
    return new_w, new_h

def copy_images():
    """Copy all figure images to word/media/."""
    media_dir = os.path.join(UNPACKED, "word", "media")
    os.makedirs(media_dir, exist_ok=True)
    for fig_id, filename, rid, label in FIGURES:
        src = os.path.join(FIGURES_DIR, filename)
        # Use a clean name in media folder
        dst = os.path.join(media_dir, filename)
        shutil.copy2(src, dst)
        print(f"  Copied {filename} -> word/media/{filename}")

def add_relationships():
    """Add image relationships to document.xml.rels."""
    rels_path = os.path.join(UNPACKED, "word", "_rels", "document.xml.rels")
    with open(rels_path, "r") as f:
        content = f.read()
    
    # Add relationships before closing tag
    new_rels = ""
    for fig_id, filename, rid, label in FIGURES:
        rel = f'  <Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{filename}"/>\n'
        new_rels += rel
    
    content = content.replace("</Relationships>", new_rels + "</Relationships>")
    
    with open(rels_path, "w") as f:
        f.write(content)
    print("  Updated document.xml.rels with image relationships")

def make_image_paragraph(rid, cx, cy, fig_id, description):
    """Generate the XML for an inline image paragraph."""
    # Use a unique docPr id based on rId number
    doc_pr_id = int(rid.replace("rId", "")) + 100
    
    return f'''    <w:p>
      <w:pPr>
        <w:jc w:val="center"/>
        <w:spacing w:before="240" w:after="120"/>
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

def inject_figures_into_document():
    """Insert figure images before their respective legends in document.xml."""
    doc_path = os.path.join(UNPACKED, "word", "document.xml")
    with open(doc_path, "r") as f:
        content = f.read()
    
    # We need to add namespace declarations for drawing elements
    # Check if wp namespace already exists
    if 'xmlns:wp=' not in content:
        content = content.replace(
            'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"',
            'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
            'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingml" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
        )
    elif 'xmlns:r=' not in content:
        content = content.replace(
            'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"',
            'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
        )
    
    # Now inject each figure before its legend
    # Strategy: Find the paragraph containing each figure legend text and insert image paragraph before it
    
    # Main figures 1-6: find "Figure N." bold text in the FIGURE LEGENDS section
    figure_insertions = {
        "Figure 1.": ("fig1", "fig1_volcano_plots.png", "rId21"),
        "Figure 2.": ("fig2", "fig2_meta_analysis_summary.png", "rId22"),
        "Figure 3.": ("fig3", "fig5_ppi_network.png", "rId23"),
        "Figure 4.": ("fig4", "fig7_wgcna_summary.png", "rId24"),
        "Figure 5.": ("fig5", "fig8_deconvolution_summary_v2.png", "rId25"),
        "Figure 6.": ("fig6", "fig9_sensitivity_analysis.png", "rId26"),
        "Supplementary Figure S1.": ("sfig1", "fig3_deg_comparison.png", "rId27"),
        "Supplementary Figure S2.": ("sfig2", "fig4_AS_enrichment.png", "rId28"),
        "Supplementary Figure S3.": ("sfig3", "fig6_hub_genes.png", "rId29"),
    }
    
    for label_text, (fig_id, filename, rid) in figure_insertions.items():
        filepath = os.path.join(FIGURES_DIR, filename)
        cx, cy = get_image_dimensions_emu(filepath)
        
        img_xml = make_image_paragraph(rid, cx, cy, fig_id, label_text.rstrip('.'))
        
        # Find the paragraph that contains this label - look for the bold text
        # The pattern is: <w:t ...>Figure N.</w:t> inside a bold run
        # We need to find the <w:p> that contains this and insert before it
        
        # Use a regex to find the paragraph containing this exact figure label
        # The label appears as bold text: <w:t xml:space="preserve">Figure N.</w:t>
        escaped_label = label_text.replace(".", r"\.")
        
        # Find the <w:p> element that contains this figure label
        # Pattern: from <w:p> ... label text ... to the start of that <w:p>
        pattern = rf'(    <w:p>\s*\n\s*<w:pPr>.*?</w:pPr>\s*\n\s*<w:r>\s*\n\s*<w:rPr>.*?<w:b/>.*?</w:rPr>\s*\n\s*<w:t[^>]*>{escaped_label}</w:t>)'
        
        match = re.search(pattern, content, re.DOTALL)
        if match:
            # Insert image paragraph before this paragraph
            insert_point = match.start()
            # Find the <w:p> opening tag
            p_start = content.rfind("    <w:p>", 0, insert_point + 10)
            if p_start >= 0:
                content = content[:p_start] + img_xml + content[p_start:]
                print(f"  Inserted {label_text} image ({cx}x{cy} EMU) before legend")
            else:
                print(f"  WARNING: Could not find <w:p> start for {label_text}")
        else:
            print(f"  WARNING: Could not find legend text for {label_text}")
    
    # Now add graphical abstract - insert it right after the "FIGURE LEGENDS" heading
    ga_filepath = os.path.join(FIGURES_DIR, "graphical_abstract.png")
    ga_cx, ga_cy = get_image_dimensions_emu(ga_filepath)
    ga_xml = make_image_paragraph("rId20", ga_cx, ga_cy, "graphical_abstract", "Graphical Abstract - Study Pipeline Overview")
    
    # Also add a label paragraph for the graphical abstract
    ga_label_xml = '''    <w:p>
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
    
    # Insert GA image + legend right after "FIGURE LEGENDS" heading
    fl_marker = '<w:t xml:space="preserve">FIGURE LEGENDS</w:t>'
    fl_pos = content.find(fl_marker)
    if fl_pos >= 0:
        # Find the end of that paragraph
        end_p = content.find("</w:p>", fl_pos)
        if end_p >= 0:
            insert_after = end_p + len("</w:p>") + 1  # +1 for newline
            content = content[:insert_after] + ga_xml + ga_label_xml + content[insert_after:]
            print(f"  Inserted Graphical Abstract after FIGURE LEGENDS heading")
    
    with open(doc_path, "w") as f:
        f.write(content)
    print("  Updated document.xml with all figure images")

if __name__ == "__main__":
    print("Step 1: Copying images to word/media/...")
    copy_images()
    
    print("\nStep 2: Adding image relationships...")
    add_relationships()
    
    print("\nStep 3: Injecting figures into document.xml...")
    inject_figures_into_document()
    
    print("\nDone! Ready to repack.")
