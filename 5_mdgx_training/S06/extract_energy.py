from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional


def extract_last_scf_energy(log_path: Path, keyword: str = "SCF Done") -> Optional[float]:
    """
    Extract the last SCF energy (Hartree) from a Gaussian .log file.

    Example line:
        SCF Done:  E(RB3LYP) =  -228.123456789     A.U. after   10 cycles
    """
    last_energy: Optional[float] = None

    # Stream-read the log file to handle large files efficiently
    with log_path.open("r", errors="replace") as fh:
        for line in fh:
            if keyword not in line:
                continue

            parts = line.split()

            # Preferred robust parsing: read the number right after "="
            try:
                eq_idx = parts.index("=")
                last_energy = float(parts[eq_idx + 1])
            except (ValueError, IndexError):
                # Fallback: attempt fixed-position parsing (less robust)
                try:
                    last_energy = float(parts[4])
                except (ValueError, IndexError):
                    continue

    return last_energy


def numeric_key(path: Path) -> int:
    """
    Natural-sort helper: extract the first integer from filename (stem).
    If none exists, return a large number so it goes to the end.
    """
    m = re.search(r"(\d+)", path.stem)
    return int(m.group(1)) if m else 10**18


def build_energy_dat_from_folder(
    folder: Path,
    pattern: str = "Conf*.log",
    output_file: str = "energy.dat",
    strict: bool = True,
) -> List[float]:
    """
    Scan a folder for Gaussian log files, extract the last SCF energy from each,
    and write energies (Hartree) to output_file.

    Parameters
    ----------
    folder : Path
        Folder containing Gaussian .log files.
    pattern : str
        Glob pattern used to select log files (default: Conf*.log).
    output_file : str
        Output file name (written inside 'folder' unless an absolute path is given).
    strict : bool
        If True, raise an error if any file has no SCF Done line.
        If False, skip problematic files with warnings.

    Returns
    -------
    List[float]
        Extracted energies in the same order as the sorted log list.
    """
    if not folder.exists():
        raise FileNotFoundError(f"[ERROR] Folder not found: {folder}")
    if not folder.is_dir():
        raise NotADirectoryError(f"[ERROR] Not a directory: {folder}")

    log_files = sorted(folder.glob(pattern), key=numeric_key)

    if not log_files:
        raise FileNotFoundError(f"[ERROR] No files matched: {folder}/{pattern}")

    energies: List[float] = []

    for log_path in log_files:
        e = extract_last_scf_energy(log_path)

        if e is None:
            msg = f"[ERROR] 'SCF Done' not found or unparsable in: {log_path.name}"
            if strict:
                raise ValueError(msg)
            print(msg)
            continue

        energies.append(e)
        print(f"[OK] {log_path.name}: {e:.12f} Hartree")

    # Always write output in the same directory as this Python script
    script_dir = Path(__file__).resolve().parent
    out_path = script_dir / output_file

    out_path.write_text("".join(f"{e:.12f}\n" for e in energies))
    print(f"\nWrote {len(energies)} energies -> {out_path.resolve()}")

    return energies


if __name__ == "__main__":
    """
    Command-line interface.

    Basic usage:
        python code.py path/to/log/folder

    Useful options:
        python code.py path/to/log/folder --pattern "*.log"
        python code.py path/to/log/folder --pattern "Conf*.log" --output energy.dat
        python code.py path/to/log/folder --loose
    """
    parser = argparse.ArgumentParser(
        description="Extract last SCF energies from Gaussian log files in a folder."
    )

    # Required positional argument: folder path
    parser.add_argument(
        "folder",
        type=str,
        help="Folder containing Gaussian .log files"
    )

    # Optional: file selection pattern
    parser.add_argument(
        "--pattern",
        type=str,
        default="Conf*.log",
        help="Glob pattern to select log files (default: Conf*.log)"
    )

    # Optional: output file name/path
    parser.add_argument(
        "--output",
        type=str,
        default="energy.dat",
        help="Output file name/path (default: energy.dat)"
    )

    # Optional: non-strict mode
    parser.add_argument(
        "--loose",
        action="store_true",
        help="Skip missing/unparsable SCF energies instead of stopping"
    )

    args = parser.parse_args()

    build_energy_dat_from_folder(
        folder=Path(args.folder),
        pattern=args.pattern,
        output_file=args.output,
        strict=not args.loose,
    )