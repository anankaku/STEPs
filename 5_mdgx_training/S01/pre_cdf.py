import re
import os
from pathlib import Path

# Path to the topology file (used in cpptraj as "parm")
top_path = Path("../../3_antechamber/S01/s01.top")

# Path where make_cdf.in will be written
output_file = Path("make_cdf.in")

# Directory containing Conf*.pdb files
pdb_dir = Path("../../4_mdgx_topology/S01")


def numeric_key(path: Path):
    """
    Extract the first integer from the filename for numeric sorting.
    This prevents incorrect ordering like:
    Conf1, Conf10, Conf2, ...
    """
    m = re.search(r"(\d+)", path.stem)
    return int(m.group(1)) if m else 0


# Collect and numerically sort all Conf*.pdb files
pdb_files = sorted(pdb_dir.glob("Conf*.pdb"), key=numeric_key)

# Define the directory where cpptraj will be executed.
# Here we assume you will run:
#   cd S01
#   cpptraj -i make_cdf.in
run_dir = output_file.parent  # ./S01


def rel_to_run_dir(p: Path) -> str:
    """
    Convert an absolute path to a relative path
    with respect to the cpptraj execution directory.
    This ensures trajin/parm paths are correct.
    """
    return os.path.relpath(p.resolve(), run_dir.resolve()).replace("\\", "/")


# Ensure output directory exists
output_file.parent.mkdir(parents=True, exist_ok=True)

# Write the cpptraj input file
with output_file.open("w") as f:
    # Write topology path (relative to run_dir)
    f.write(f"parm {rel_to_run_dir(top_path)}\n")

    # Write trajin lines with correct relative paths
    for pdb in pdb_files:
        f.write(f"trajin {rel_to_run_dir(pdb)}\n")

    # Output NetCDF trajectory
    # coords.cdf will be created inside run_dir (./S01)
    f.write("trajout coords.cdf netcdf\n")

    # Required cpptraj execution commands
    f.write("run\n")
    f.write("quit\n")

print(f"Generated {output_file} with {len(pdb_files)} structures.")
print(f"Run using:")
print(f"   cd {run_dir}")
print(f"   cpptraj -i {output_file.name}")
