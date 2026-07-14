import fitz # PyMuPDF
import os

def get_best_fit_rect(quad_rect, panel_w, panel_h):
    """
    Calculates the largest rectangle of dimensions (panel_w, panel_h) 
    that fits entirely inside quad_rect while preserving its aspect ratio.
    Centers it within quad_rect.
    """
    quad_w = quad_rect.x1 - quad_rect.x0
    quad_h = quad_rect.y1 - quad_rect.y0
    
    scale = min(quad_w / panel_w, quad_h / panel_h)
    new_w = panel_w * scale
    new_h = panel_h * scale
    
    # Calculate centering offsets
    dx = (quad_w - new_w) / 2
    dy = (quad_h - new_h) / 2
    
    return fitz.Rect(
        quad_rect.x0 + dx, 
        quad_rect.y0 + dy, 
        quad_rect.x0 + dx + new_w, 
        quad_rect.y0 + dy + new_h
    )

def main():
    print("=== STARTING MAIN FIGURE 4 ASSEMBLY ===")
    print("Assembling high-fidelity 5-panel portrait Main Figure 4 (1440pt x 2160pt)...")
    
    # Create output document
    doc_out = fitz.open()
    # 720pt x 1150pt canvas to add vertical breathing room
    page_out = doc_out.new_page(width=720, height=1150)
    
    # Define Panels
    panels = [
        {
            "id": "A",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\Figure5A_VRS_Histogram.pdf",
            "quad": fitz.Rect(25, 50, 345, 350),
            "label_pos": fitz.Point(25, 35)
        },
        {
            "id": "B",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\VRS_validation_beataml2_scatter.pdf",
            "quad": fitz.Rect(375, 50, 695, 350),
            "label_pos": fitz.Point(375, 35)
        },
        {
            "id": "C",
            "path": r"d:\Proj_AML\04_Figures\06_Drug_Response\drug_sensitivity_heatmap.pdf",
            "quad": fitz.Rect(25, 420, 400, 780),
            "label_pos": fitz.Point(25, 395)
        },
        {
            "id": "E",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\Figure5E_VIALE_Validation.pdf",
            "quad": fitz.Rect(410, 420, 695, 780),
            "label_pos": fitz.Point(410, 395)
        },
        {
            "id": "D",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\Figure5D_Salvage_Boxplots.pdf",
            "quad": fitz.Rect(25, 830, 695, 1130),
            "label_pos": fitz.Point(25, 805)
        }
    ]
    
    # Place each panel on the canvas
    for panel in panels:
        path = panel["path"]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Source panel for {panel['id']} not found: {path}")
            
        print(f"  Placing Panel {panel['id']} ({os.path.basename(path)})...")
        doc_in = fitz.open(path)
        
        # Calculate optimal fit rectangle to preserve aspect ratio
        actual_w = doc_in[0].rect.width
        actual_h = doc_in[0].rect.height
        fit_rect = get_best_fit_rect(panel["quad"], actual_w, actual_h)
        
        # Place vector overlay
        page_out.show_pdf_page(fit_rect, doc_in, 0)
        print(f"    [OK] Placed Panel {panel['id']} at {fit_rect}")
        doc_in.close()
        
    # Draw vector-based panel labels for Panel C and Panel E
    print("  Drawing bold labels for Panel C and Panel E...")
    for panel in panels:
        if panel["id"] == "C":
            # Draw native-looking title for the heatmap
            page_out.insert_text(
                panel["label_pos"], 
                "C. Global Ex Vivo Resistance Profile", 
                fontsize=9.6, 
                fontname="hebo", # Helvetica-Bold
                color=(0, 0, 0.54)
            )
            # Draw subtitle
            subtitle_pos = fitz.Point(panel["label_pos"].x, panel["label_pos"].y + 11)
            page_out.insert_text(
                subtitle_pos, 
                "Targeted drug sensitivity across patient cohort", 
                fontsize=8.4, 
                fontname="helv", 
                color=(0, 0, 0.54)
            )
        elif panel["id"] == "E":
            # Draw native-looking title for the clinical trial validation plot
            page_out.insert_text(
                panel["label_pos"], 
                "E. Clinical Trial Validation (VIALE-A)", 
                fontsize=9.6, 
                fontname="hebo", # Helvetica-Bold
                color=(0, 0, 0.54)
            )
            # Draw subtitle
            subtitle_pos = fitz.Point(panel["label_pos"].x, panel["label_pos"].y + 11)
            page_out.insert_text(
                subtitle_pos, 
                "VRS tracks the 24-gene resistance signature", 
                fontsize=8.4, 
                fontname="helv", 
                color=(0, 0, 0.54)
            )
        
    # Save the consolidated figure
    pdf_out = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures\Figure4_Consolidated.pdf"
    png_out = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures\Figure4_Consolidated.png"
    
    os.makedirs(os.path.dirname(pdf_out), exist_ok=True)
    doc_out.save(pdf_out)
    doc_out.close()
    print(f"  [OK] Successfully compiled primary consolidated Figure 4: {pdf_out}")
    
    # Render high-resolution PNG preview (at 300 DPI)
    print("  Rendering primary consolidated PNG preview...")
    doc_png = fitz.open(pdf_out)
    page_png = doc_png[0]
    pix = page_png.get_pixmap(dpi=300)
    pix.save(png_out)
    doc_png.close()
    print(f"  [OK] Successfully rendered primary consolidated PNG: {png_out}")
    
    print("\n=== MAIN FIGURE 4 ASSEMBLY COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
