from __future__ import annotations

from pathlib import Path
from typing import List, Optional


def extract_last_scf_energy(log_path: Path, keyword: str = "SCF Done") -> Optional[float]:
    """
    Return the last SCF energy (Hartree) found in a Gaussian .log file.
    Gaussian line example:
      SCF Done:  E(RB3LYP) =  -228.123456789     A.U. after   10 cycles
    """
    last_energy: Optional[float] = None

    with log_path.open("r", errors="replace") as fh:
        for line in fh:
            if keyword not in line:
                continue

            # robust parse: find token after '='
            # tokens typically: ['SCF','Done:','E(RB3LYP)','=', '-228.123', 'A.U.', ...]
            parts = line.split()
            try:
                eq_idx = parts.index("=")
                last_energy = float(parts[eq_idx + 1])
            except (ValueError, IndexError):
                # fallback: your original approach (5th token) if '=' parse fails
                try:
                    last_energy = float(parts[4])
                except (ValueError, IndexError):
                    continue

    return last_energy


def build_energy_dat(
    num_conform: int,
    log_template: str = "Conf{idx}.log",
    output_file: str = "energy.dat",
    strict: bool = True,
) -> List[float]:
    """
    Extract energies in order Conf1..ConfN and write to output_file.
    If strict=True, raise an error when a file is missing or has no SCF Done line.
    """
    energies: List[float] = []

    for i in range(1, num_conform + 1):
        log_path = Path(log_template.format(idx=i))

        if not log_path.exists():
            msg = f"[ERROR] Missing file: {log_path}"
            if strict:
                raise FileNotFoundError(msg)
            print(msg)
            continue

        e = extract_last_scf_energy(log_path)
        if e is None:
            msg = f"[ERROR] '{'SCF Done'}' not found or unparsable in: {log_path}"
            if strict:
                raise ValueError(msg)
            print(msg)
            continue

        energies.append(e)
        print(f"[OK] {log_path}: {e:.12f} Hartree")

    out_path = Path(output_file)
    out_path.write_text("".join(f"{e:.12f}\n" for e in energies))

    print(f"\nWrote {len(energies)} energies -> {out_path.resolve()}")
    return energies


if __name__ == "__main__":
    build_energy_dat(num_conform=10, output_file="energy.dat", strict=True)
