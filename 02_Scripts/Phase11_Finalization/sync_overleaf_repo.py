import os
import shutil

def main():
    print("=== STARTING OVERLEAF REPO SYNCHRONIZATION ===")
    
    src_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source"
    dst_dir = r"d:\Proj_AML\Overleaf_Repo"
    
    print(f"Source: {src_dir}")
    print(f"Destination: {dst_dir}\n")
    
    if not os.path.exists(src_dir):
        raise FileNotFoundError(f"Source folder not found: {src_dir}")
    if not os.path.exists(dst_dir):
        raise FileNotFoundError(f"Destination folder (Overleaf_Repo) not found: {dst_dir}")
        
    copied_count = 0
    skipped_count = 0
    
    for root, dirs, files in os.walk(src_dir):
        # Prevent walking into .git or copying it
        if ".git" in dirs:
            dirs.remove(".git")
            
        for file in files:
            # Skip temp lock files
            if file.startswith("~$"):
                continue
                
            src_file = os.path.join(root, file)
            # Compute relative path
            rel_path = os.path.relpath(src_file, src_dir)
            dst_file = os.path.join(dst_dir, rel_path)
            
            # Skip if destination is identical (same size and mod time or same content)
            # To be safe, check size and exist
            if os.path.exists(dst_file):
                try:
                    if os.path.samefile(src_file, dst_file):
                        # Junction or hard link
                        skipped_count += 1
                        continue
                except Exception:
                    pass
                    
                src_size = os.path.getsize(src_file)
                dst_size = os.path.getsize(dst_file)
                if src_size == dst_size:
                    # Let's check mod times or just assume same if size is equal. 
                    # To be perfectly correct, let's copy if mod time is newer or just overwrite if not identical.
                    src_mtime = os.path.getmtime(src_file)
                    dst_mtime = os.path.getmtime(dst_file)
                    if abs(src_mtime - dst_mtime) < 1.0: # Close enough
                        skipped_count += 1
                        continue
            
            # Make sure destination folder exists
            os.makedirs(os.path.dirname(dst_file), exist_ok=True)
            shutil.copy2(src_file, dst_file)
            print(f"  [SYNC] Copied: {rel_path}")
            copied_count += 1
            
    print(f"\nSynchronization complete: {copied_count} files copied, {skipped_count} files up-to-date.")
    print("=== OVERLEAF REPO SYNCHRONIZATION COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
