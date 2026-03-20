#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path
from dataclasses import dataclass

@dataclass
class XYZAtom:
    idx1: int      # 1-based index in the XYZ file
    elem: str
    x: float
    y: float
    z: float
    raw: str       # original line (for writing back)

RANGE_RE = re.compile(r"^\s*(\d+)\s*-\s*(\d+)\s*$")

def parse_indices_1based(expr: str) -> set[int]:
    """
    Parse indices like:
      '1 2 3 10-15, 20'
    into a set of integers (1-based).
    """
    expr = expr.strip()
    if not expr:
        return set()

    parts = re.split(r"[,\s]+", expr)
    out: set[int] = set()
    for p in parts:
        if not p:
            continue
        m = RANGE_RE.match(p)
        if m:
            a = int(m.group(1))
            b = int(m.group(2))
            if a > b:
                a, b = b, a
            out.update(range(a, b + 1))
        else:
            out.add(int(p))
    return out

def read_xyz(path: Path) -> tuple[str, str, list[XYZAtom]]:
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < 3:
        raise RuntimeError("XYZ file too short (needs at least 3 lines).")

    try:
        nat = int(lines[0].strip())
    except ValueError:
        raise RuntimeError("First line of XYZ must be an integer atom count.")

    comment = lines[1] if len(lines) >= 2 else ""
    atom_lines = lines[2:2 + nat]
    if len(atom_lines) != nat:
        raise RuntimeError(
            f"XYZ atom count mismatch: header says {nat}, but found {len(atom_lines)} atom lines."
        )

    atoms: list[XYZAtom] = []
    for i, ln in enumerate(atom_lines, start=1):
        fields = ln.split()
        if len(fields) < 4:
            raise RuntimeError(f"Bad atom line at index {i}: {ln}")
        elem = fields[0]
        x, y, z = float(fields[1]), float(fields[2]), float(fields[3])
        atoms.append(XYZAtom(idx1=i, elem=elem, x=x, y=y, z=z, raw=ln))

    header = str(nat)
    return header, comment, atoms

def write_xyz(path: Path, comment: str, atoms: list[XYZAtom]) -> None:
    out_lines = [str(len(atoms)), comment]
    out_lines += [a.raw for a in atoms]
    path.write_text("\n".join(out_lines) + "\n")

def main() -> None:
    print("\n=== Trim XYZ atoms by index (interactive, 1-based) ===\n")

    xyz_in = Path(input("Input XYZ path: ").strip()).expanduser()
    if not xyz_in.exists():
        raise SystemExit(f"[ERROR] File not found: {xyz_in}")

    header, comment, atoms = read_xyz(xyz_in)
    n = len(atoms)
    print(f"[OK] Read {n} atoms from {xyz_in}")

    print("\nEnter atom indices to REMOVE (1-based).")
    print("You can use ranges, e.g.: 1 2 3 15-20, 33")
    rm_expr = input("Indices to remove: ").strip()
    rm = parse_indices_1based(rm_expr)

    bad = sorted(i for i in rm if i < 1 or i > n)
    if bad:
        raise SystemExit(f"[ERROR] Out of range indices (1..{n}): {bad}")

    removed = [a for a in atoms if a.idx1 in rm]
    kept = [a for a in atoms if a.idx1 not in rm]

    print("\n=== Preview removed atoms ===")
    for a in removed[:30]:
        print(f"{a.idx1:4d} {a.elem:2s}  {a.x: .6f} {a.y: .6f} {a.z: .6f}")
    if len(removed) > 30:
        print(f"... ({len(removed) - 30} more not shown)")

    print("\n=== Summary ===")
    print(f"Removed: {len(removed)}")
    print(f"Kept   : {len(kept)}")

    out_name = input("\nOutput XYZ path (e.g., uncapped.xyz): ").strip()
    xyz_out = Path(out_name).expanduser()

    keep_comment = input("Keep original comment line? [Y/n]: ").strip().lower()
    out_comment = comment if keep_comment in ("", "y", "yes") else ""

    write_xyz(xyz_out, out_comment, kept)
    print(f"\n[DONE] Wrote {len(kept)} atoms to {xyz_out}\n")

if __name__ == "__main__":
    main()
