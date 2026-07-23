import os
import sys
import subprocess
import time

def run_script(cmd, cwd=None):
    print(f"\n========================================\n[RUN] {' '.join(cmd)}")
    start_time = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    elapsed = time.time() - start_time
    if res.returncode != 0:
        print(f"[FAIL] (took {elapsed:.2f}s)")
        print("----- STDOUT -----")
        print(res.stdout)
        print("----- STDERR -----")
        print(res.stderr)
        return False
    print(f"[SUCCESS] (took {elapsed:.2f}s)")
    if res.stdout.strip():
        print("----- STDOUT -----")
        print(res.stdout.strip())
    return True

def main():
    print("=== STARTING MASTER REGENERATION AND SYNCHRONIZATION PIPELINE ===")
    
    script_dir = r"d:\Proj_AML\02_Scripts\Phase11_Finalization"
    os.chdir(r"d:\Proj_AML")
    
    # 1. Individual R Scripts for Main Figures & Supplementary Figures
    r_scripts = [
        "01_generate_MainFigure1_HighFid.R",
        "03_generate_MainFigure3_HighFid.R",
        "04_generate_MainFigure4_HighFid.R",
        "05_generate_MainFigure5_HighFid.R",
        "generate_MainFigure5_PanelA.R",
        "generate_Figure5B_VRS_Validation_Boxplot.R",
        "generate_Figure5D_Salvage_Boxplots.R",
        "generate_Figure5E_VIALE_Validation.R",
        # Supplementary Figures R scripts
        "01_generate_FigureS1.R",
        "generate_FigureS2_HighFid_Pediatric.R",
        "02_generate_FigureS2_HighFid.R",
        "03_generate_FigureS3_HighFid.R",
        "04_generate_FigureS4_HighFid.R",
        "05_generate_FigureS5_HighFid.R",
        "06_generate_FigureS6_HighFid.R",
        "generate_FigureS10_HighFid.R",
        "generate_FigureS11_HighFid.R",
        "generate_FigureS15_HighFid.R",
        "../Phase14_AMLSG_Validation/05_compile_european_validation_panels.R",
        "../Phase9_ExternalValidation/05_generate_CITEseq_ADT_Validation.R",
        "../Phase9_ExternalValidation/06_generate_CITEseq_longitudinal_dynamics.R",
        "../Phase9_ExternalValidation/07_generate_CITEseq_UMAP_plots.R"
    ]
    
    # 2. Python Assembly Scripts
    py_scripts = [
        "assemble_MainFigure4_SingleCell.py",
        "assemble_MainFigure5_HighFid.py",
        "assemble_MainFigure6_HighFid.py",
        "generate_FigureS2.py",
        "generate_FigureS16.py",
        "assemble_s3_high_fid.py",
        "assemble_s4_high_fid_v2.py",
        "assemble_s5_high_fid.py",
        "assemble_s6_high_fid.py",
        "assemble_new_supplementary_figures.py",
        "assemble_s17_to_s20.py"
    ]
    
    # 3. Synchronization & Packaging Scripts
    sync_and_pack_scripts = [
        "synchronize_main_figures.py",
        "synchronize_supplementary_figures.py",
        "sync_overleaf_repo.py",
        "rebuild_overleaf_package.py"
    ]
    
    # 4. High-Impact Biological Enhancement Scripts
    enhancement_scripts = [
        r"d:\Proj_AML\02_Scripts\Phase10_Enhancements\05_depmap_salvage_validation.R",
        r"d:\Proj_AML\02_Scripts\Phase10_Enhancements\06_proteomic_cohort_validation.R"
    ]
    
    failed = []
    
    print("\n--- PHASE 0: Running High-Impact Biological Enhancement Scripts ---")
    for script_path in enhancement_scripts:
        if not os.path.exists(script_path):
            print(f"[SKIP] Script not found: {script_path}")
            continue
        success = run_script(["Rscript", script_path])
        if not success:
            failed.append(os.path.basename(script_path))
            
    print("\n--- PHASE 1: Running R Panel & Main Figure Scripts ---")
    for script in r_scripts:
        script_path = os.path.join(script_dir, script)
        if not os.path.exists(script_path):
            print(f"[SKIP] Script not found: {script_path}")
            continue
        success = run_script(["Rscript", script_path])
        if not success:
            failed.append(script)
            
    print("\n--- PHASE 2: Running Python Figure Assembly Scripts ---")
    for script in py_scripts:
        script_path = os.path.join(script_dir, script)
        if not os.path.exists(script_path):
            print(f"[SKIP] Script not found: {script_path}")
            continue
        success = run_script(["python", script_path], cwd=script_dir)
        if not success:
            failed.append(script)
            
    print("\n--- PHASE 3: Running Synchronization & Packaging Scripts ---")
    for script in sync_and_pack_scripts:
        script_path = os.path.join(script_dir, script)
        if not os.path.exists(script_path):
            print(f"[SKIP] Script not found: {script_path}")
            continue
        success = run_script(["python", script_path], cwd=script_dir)
        if not success:
            failed.append(script)
            
    print("\n========================================")
    if failed:
        print(f"[WARNING] Pipeline finished with failures in {len(failed)} scripts: {', '.join(failed)}")
        sys.exit(1)
    else:
        print("[SUCCESS] All pipeline scripts executed successfully with exit code 0!")
        sys.exit(0)

if __name__ == "__main__":
    main()
