import os
import zipfile

def main():
    print("=== STARTING OVERLEAF SUBMISSION PACKAGE REBUILD ===")
    
    src_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source"
    zip_path = r"d:\Proj_AML\05_Submission\Overleaf_Submission_Package_2026-05-16.zip"
    
    print(f"Source folder: {src_dir}")
    print(f"Creating archive: {zip_path}\n")
    
    # We want to walk through all files and add them, preserving relative paths
    # from inside the src_dir (so that Manuscript.tex is at the root of the zip)
    
    if not os.path.exists(src_dir):
        raise FileNotFoundError(f"Source directory not found: {src_dir}")
        
    # Open zip file
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        file_count = 0
        for root, dirs, files in os.walk(src_dir):
            for file in files:
                # Skip Excel temporary lock files
                if file.startswith("~$"):
                    print(f"  [INFO] Skipping Excel lock file: {file}")
                    continue
                # We skip temporary check/render PNGs that were generated during audit to keep the package clean
                if "check_render" in file or "extracted_from_pdf" in file:
                    print(f"  [INFO] Skipping temporary file: {file}")
                    continue
                    
                full_path = os.path.join(root, file)
                # Compute relative path from src_dir
                rel_path = os.path.relpath(full_path, src_dir)
                
                zipf.write(full_path, rel_path)
                print(f"  [ADD] {rel_path}")
                file_count += 1
                
    print(f"\nSuccessfully archived {file_count} files into: {zip_path}")
    print("=== OVERLEAF SUBMISSION PACKAGE REBUILD COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
