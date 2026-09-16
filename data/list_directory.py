import os
from pathlib import Path

# Edit this list if your extracted folder names differ
ROOTS = ["Tab", "MySQL", "SDY1529-DR58_Tab", "SDY1529-DR58_MySQL"]

output_lines = []

for root_name in ROOTS:
    root = Path(root_name)
    if not root.exists():
        continue
    output_lines.append(f"=== {root_name} ===")
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        filenames.sort()
        rel_dir = os.path.relpath(dirpath, root)
        for fname in filenames:
            full_path = Path(dirpath) / fname
            try:
                size = full_path.stat().st_size
            except OSError:
                size = -1
            rel_path = os.path.join(rel_dir, fname) if rel_dir != "." else fname
            output_lines.append(f"{rel_path}\t{size} bytes")
    output_lines.append("")

with open("directory_listing.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(output_lines))

print(f"Wrote {len(output_lines)} lines to directory_listing.txt")
