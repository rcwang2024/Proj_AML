import subprocess
import os
import shutil
import fitz

def main():
    print("=== STARTING FIGURE S3 GENERATION AND SYNCHRONIZATION ===")
    
    # 1. Run R script to generate component PDFs
    r_script = r"d:\Proj_AML\02_Scripts\Phase11_Finalization\03_generate_FigureS3_HighFid.R"
    print(f"Running component generator: {r_script}")
    subprocess.run(["Rscript", r_script], check=True)
    
    # 2. Canvas assembly settings
    path_in = r"d:\Proj_AML\04_Figures\12_ELN_Comparison\HighFid"
    dest_pdf_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS3.pdf"
    dest_png_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS3.png"
    
    panels = [
        "s3_pA.pdf", # Top-Left
        "s3_pB.pdf", # Top-Right
        "s3_pC.pdf", # Bottom-Left
        "s3_pD.pdf"  # Bottom-Right
    ]
    
    print("\nAssembling high-fidelity Figure S3 canvas (2000pt x 2000pt)...")
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=2000, height=2000)
    
    margin = 50
    w = 925
    h = 925
    
    for i, p_name in enumerate(panels):
        p_path = os.path.join(path_in, p_name)
        if not os.path.exists(p_path):
            raise FileNotFoundError(f"Component not found: {p_path}")
            
        col = i % 2
        row = i // 2
        
        rect = fitz.Rect(
            margin + col*(w + margin),
            margin + row*(h + margin),
            margin + col*(w + margin) + w,
            margin + row*(h + margin) + h
        )
        
        src = fitz.open(p_path)
        page_out.show_pdf_page(rect, src, 0)
        src.close()
        print(f"  [OK] Placed {p_name}")
        
    os.makedirs(os.path.dirname(dest_pdf_1), exist_ok=True)
    doc_out.save(dest_pdf_1)
    
    # Render PNG preview at ultra high res (300 DPI)
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png_1)
    doc_out.close()
    print(f"  [OK] Successfully compiled primary consolidated Figure S3: {dest_pdf_1}")
    
    # 3. Synchronize with target folders
    targets = [
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS3.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS3.png"
        ),
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS3.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS3.png"
        )
    ]
    
    def copy_file_safe(src, dst):
        try:
            if os.path.exists(dst) and os.path.samefile(src, dst):
                print(f"  [INFO] Junction/link detected: {dst} is identical to {src}. No copy needed.")
            else:
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy(src, dst)
                print(f"  [OK] Copied to: {dst}")
        except Exception:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy(src, dst)
            print(f"  [OK] Copied to: {dst} (via backup handler)")

    print("\nSynchronizing outputs across workspace targets...")
    for pdf_tgt, png_tgt in targets:
        copy_file_safe(dest_pdf_1, pdf_tgt)
        copy_file_safe(dest_png_1, png_tgt)
        
    print("\n=== FIGURE S3 COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
