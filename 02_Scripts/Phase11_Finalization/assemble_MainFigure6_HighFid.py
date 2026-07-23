import os
import shutil
import fitz

def get_best_fit_rect(quad_rect, panel_w, panel_h):
    quad_w = quad_rect.x1 - quad_rect.x0
    quad_h = quad_rect.y1 - quad_rect.y0
    
    scale = min(quad_w / panel_w, quad_h / panel_h)
    new_w = panel_w * scale
    new_h = panel_h * scale
    
    dx = (quad_w - new_w) / 2
    dy = (quad_h - new_h) / 2
    
    return fitz.Rect(
        quad_rect.x0 + dx, 
        quad_rect.y0 + dy, 
        quad_rect.x0 + dx + new_w, 
        quad_rect.y0 + dy + new_h
    )

def copy_file_safe(src, dst):
    try:
        if not os.path.exists(src):
            print(f"  [WARNING] Source not found: {src}")
            return False
        if os.path.exists(dst) and os.path.samefile(src, dst):
            return True
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        return True
    except Exception as e:
        print(f"  [WARNING] Copy failed from {src} to {dst}: {e}")
        return False

def main():
    print("=== STARTING MAIN FIGURE 6 ASSEMBLY (SALVAGE THERAPEUTICS) ===")
    
    width, height = 800, 1000
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=width, height=height)
    
    panels = [
        {
            "id": "A",
            "path": r"d:\Proj_AML\04_Figures\27_Cluster2_Salvage\Figure_Cluster2_Drug_Classes.pdf",
            "quad": fitz.Rect(10, 10, 390, 320)
        },
        {
            "id": "B",
            "path": r"d:\Proj_AML\04_Figures\27_Cluster2_Salvage\Figure_Cluster_Comparison.pdf",
            "quad": fitz.Rect(410, 10, 790, 320)
        },
        {
            "id": "C",
            "path": r"d:\Proj_AML\04_Figures\27_Cluster2_Salvage\Figure_Cluster2_Top10_Drugs.pdf",
            "quad": fitz.Rect(10, 330, 390, 640)
        },
        {
            "id": "D",
            "path": r"d:\Proj_AML\04_Figures\27_Cluster2_Salvage\Figure_Combination_Therapy.pdf",
            "quad": fitz.Rect(410, 330, 790, 640)
        },
        {
            "id": "E",
            "path": r"d:\Proj_AML\04_Figures\27_Cluster2_Salvage\Figure_Panobinostat_Boxplot.pdf",
            "quad": fitz.Rect(10, 650, 390, 960)
        },
        {
            "id": "F",
            "path": r"d:\Proj_AML\04_Figures\27_Cluster2_Salvage\Figure_Selumetinib_Boxplot.pdf",
            "quad": fitz.Rect(410, 650, 790, 960)
        }
    ]
    
    for panel in panels:
        path = panel["path"]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Component {panel['id']} not found: {path}")
            
        print(f"  Placing Panel {panel['id']}...")
        doc_in = fitz.open(path)
        actual_w = doc_in[0].rect.width
        actual_h = doc_in[0].rect.height
        fit_rect = get_best_fit_rect(panel["quad"], actual_w, actual_h)
        page_out.show_pdf_page(fit_rect, doc_in, 0)
        doc_in.close()
        
    dest_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures\Figure6_Consolidated.pdf"
    dest_png = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures\Figure6_Consolidated.png"
    
    os.makedirs(os.path.dirname(dest_pdf), exist_ok=True)
    doc_out.save(dest_pdf)
    
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png)
    doc_out.close()
    
    latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Figures"
    copy_file_safe(dest_pdf, os.path.join(latex_dir, "Figure6_Consolidated.pdf"))
    copy_file_safe(dest_png, os.path.join(latex_dir, "Figure6_Consolidated.png"))
    
    print("=== MAIN FIGURE 6 ASSEMBLY COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
