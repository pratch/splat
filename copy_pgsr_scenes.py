import os
import subprocess

source_dir = "/ist-nas/users/pratchp/projects/PGSR/output_all/"
dest_dir = "./pgsr_data/"

# Create destination directory if it doesn't exist
os.makedirs(dest_dir, exist_ok=True)

# Use rsync to copy files incrementally
# -a: archive mode (preserves permissions, timestamps, etc.)
# -v: verbose
# --progress: show progress
# --include: only include matching patterns
# --exclude: exclude everything else
cmd = [
    "rsync", "-av", "--progress",
    "--include", "gs_*.splat",
    "--include", "mesh_*.ply",
    "--exclude", "*",
    source_dir,
    dest_dir
]

print(f"Syncing files from {source_dir} to {dest_dir}...")
subprocess.run(cmd)
