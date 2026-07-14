import os
import shutil

def main():
    print("=== STARTING MAIN FIGURES MASTER SYNCHRONIZATION ===")
    
    src_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures"
    
    targets = [
        r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Figures",
        r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\02_Main_Figures",
        r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\01_Manuscript_Source\Figures"
    ]
    
    # We synchronize Figure1 to Figure5, both PDF and PNG formats
    fig_basenames = [f"Figure{i}_Consolidated" for i in range(1, 6)]
    extensions = [".pdf", ".png"]
    
    files_to_sync = [base + ext for base in fig_basenames for ext in extensions]
    
    def copy_file_safe(src, dst):
        try:
            # Check if source exists
            if not os.path.exists(src):
                print(f"  [WARNING] Source file not found: {src}")
                return False
                
            # If target exists and is exactly the same, skip to preserve links
            if os.path.exists(dst) and os.path.samefile(src, dst):
                print(f"  [INFO] Junction/link detected: {dst} is identical. Skipping.")
                return True
                
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy(src, dst)
            print(f"  [OK] Copied {os.path.basename(src)} -> {os.path.dirname(dst)}")
            return True
        except Exception as e:
            # Backup handler in case samefile fails on Windows paths
            try:
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy(src, dst)
                print(f"  [OK] Copied {os.path.basename(src)} -> {os.path.dirname(dst)} (via backup)")
                return True
            except Exception as ex:
                print(f"  [ERROR] Failed to copy to {dst}: {ex}")
                return False

    print(f"Master source: {src_dir}")
    print(f"Syncing {len(files_to_sync)} figure assets to {len(targets)} target destinations...\n")
    
    success_count = 0
    total_copies = len(files_to_sync) * len(targets)
    
    for filename in files_to_sync:
        src_path = os.path.join(src_dir, filename)
        for target_dir in targets:
            dst_path = os.path.join(target_dir, filename)
            if copy_file_safe(src_path, dst_path):
                success_count += 1
                
    print(f"\nSynchronization completed: {success_count}/{total_copies} operations successful.")
    print("=== MAIN FIGURES SYNCHRONIZATION COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
