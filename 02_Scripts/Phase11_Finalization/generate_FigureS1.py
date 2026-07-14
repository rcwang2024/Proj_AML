import subprocess
import os
import shutil

def main():
    print("=== STARTING FIGURE S1 GENERATION AND SYNCHRONIZATION ===")
    
    # 1. Run R script to generate raw components
    r_script = r"d:\Proj_AML\02_Scripts\Phase11_Finalization\01_generate_FigureS1.R"
    print(f"Running component generator: {r_script}")
    subprocess.run(["Rscript", r_script], check=True)
    
    # 2. Define primary output files
    src_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS1.pdf"
    src_png = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS1.png"
    
    # 3. Synchronized target locations
    targets = [
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS1.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS1.png"
        ),
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS1.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS1.png"
        )
    ]
    
    # Helper for safe copying (symlinks / junctions support)
    def copy_file_safe(src, dst):
        try:
            if os.path.exists(dst) and os.path.samefile(src, dst):
                print(f"  [INFO] Junction/link detected: {dst} is identical to {src}. No copy needed.")
            else:
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy(src, dst)
                print(f"  [OK] Copied to: {dst}")
        except Exception as e:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy(src, dst)
            print(f"  [OK] Copied to: {dst} (via backup handler)")

    # 4. Perform synchronization copies
    print("\nSynchronizing outputs across workspace targets...")
    for pdf_tgt, png_tgt in targets:
        copy_file_safe(src_pdf, pdf_tgt)
        copy_file_safe(src_png, png_tgt)
        
    print("\n=== FIGURE S1 COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
