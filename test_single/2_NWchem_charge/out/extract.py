#!/usr/bin/env python3
from __future__ import annotations

import re
import argparse
from pathlib import Path
from dataclasses import dataclass

DATA_ROW_RE = re.compile(r"^\s*(\d+)\s+([A-Za-z]+)\s+")

@dataclass
class AtomRow:
    idx1: int          # 1-based index (as printed in the RESP table)
    elem: str
    x: float
    y: float
    z: float
    q_esp: float
    q_constr: float

def extract_resp_table(resp_out: Path, debug: bool = False) -> list[AtomRow]:
    """
    Parse the RESP table:

    Atom  Coordinates              Charge
                                   ESP     constr
    1  C   x y z   q_esp   q_constr
    ...

    Returns a list of AtomRow in table order.
    """
    lines = resp_out.read_text(errors="ignore").splitlines()

    # 1) Find header anchor
    anchor_idx = None
    for i, line in enumerate(lines):
        if ("Atom" in line) and ("Coordinates" in line):
            anchor_idx = i
            break
    if anchor_idx is None:
        raise RuntimeError("Could not find table header anchor containing 'Atom' and 'Coordinates'.")

    # 2) Find first data row
    start_idx = None
    for j in range(anchor_idx, min(anchor_idx + 2000, len(lines))):
        if DATA_ROW_RE.match(lines[j]):
            start_idx = j
            break
    if start_idx is None:
        if debug:
            snippet = "\n".join(lines[anchor_idx:anchor_idx+30])
            raise RuntimeError(
                "Found 'Atom/Coordinates' header but no data rows after it.\n"
                "Snippet near header:\n" + snippet
            )
        raise RuntimeError("Found header but no data rows after it.")

    if debug:
        print("\n[DEBUG] Header anchor line:")
        print(lines[anchor_idx])
        print("\n[DEBUG] First data row:")
        print(lines[start_idx])
        print("\n[DEBUG] Next few lines:")
        for k in range(start_idx, min(start_idx + 5, len(lines))):
            print(lines[k])
        print()

    # 3) Parse table rows
    atoms: list[AtomRow] = []
    for line in lines[start_idx:]:
        if not line.strip():
            break
        if not DATA_ROW_RE.match(line):
            break

        fields = line.split()
        # Expect at least: idx, elem, x, y, z, esp, constr
        if len(fields) < 7:
            continue

        try:
            idx1 = int(fields[0])
            elem = fields[1]
            x = float(fields[2]); y = float(fields[3]); z = float(fields[4])
            q_esp = float(fields[-2])
            q_constr = float(fields[-1])
        except ValueError:
            continue

        atoms.append(AtomRow(idx1=idx1, elem=elem, x=x, y=y, z=z, q_esp=q_esp, q_constr=q_constr))

    if not atoms:
        raise RuntimeError(
            "Located table start but parsed zero rows. "
            "Check whether the table format matches: idx elem x y z esp constr."
        )

    # Optional sanity check: idx1 should be 1..N
    # (Not strictly required, but helpful if something is off)
    return atoms


def parse_indices_1based(s: str) -> set[int]:
    parts = re.split(r"[,\s]+", s.strip())
    out: set[int] = set()
    for p in parts:
        if p:
            out.add(int(p))  # keep as 1-based here
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("resp_output", help="Path to NWChem RESP output file (e.g. R20_RESP.output)")
    ap.add_argument("-o", "--out", default="R20charge.dat", help="Output charge file")
    ap.add_argument("--debug", action="store_true", help="Print where the table was detected")
    ap.add_argument("--show-removed", action="store_true", help="Print the atoms you removed (cap atoms)")
    ap.add_argument("--show-kept", action="store_true", help="Print the atoms you kept (uncapped)")
    args = ap.parse_args()

    resp_path = Path(args.resp_output).expanduser()
    if not resp_path.exists():
        raise SystemExit(f"[ERROR] File not found: {resp_path}")

    atoms = extract_resp_table(resp_path, debug=args.debug)
    n = len(atoms)
    print(f"[OK] Parsed {n} atoms from RESP table.")

    cap_input = input("Cap atom indices to REMOVE (1-based, space/comma separated): ").strip()
    cap1 = parse_indices_1based(cap_input)

    bad = sorted(i for i in cap1 if i < 1 or i > n)
    if bad:
        raise SystemExit(f"[ERROR] These indices are out of range (1..{n}): {bad}")

    removed = [a for a in atoms if a.idx1 in cap1]
    kept = [a for a in atoms if a.idx1 not in cap1]

    if args.show_removed:
        print("\n=== Removed (cap) atoms ===")
        for a in removed:
            print(f"{a.idx1:4d} {a.elem:2s}  {a.x: .6f} {a.y: .6f} {a.z: .6f}   ESP {a.q_esp: .6f}   constr {a.q_constr: .6f}")

    if args.show_kept:
        print("\n=== Kept (uncapped) atoms ===")
        for a in kept:
            print(f"{a.idx1:4d} {a.elem:2s}  {a.x: .6f} {a.y: .6f} {a.z: .6f}   ESP {a.q_esp: .6f}   constr {a.q_constr: .6f}")

    keep_charges = [a.q_constr for a in kept]

    print("\n=== Summary ===")
    print(f"Removed atoms : {len(removed)}")
    print(f"Remaining    : {len(kept)}")
    print(f"Sum(all constr)     = {sum(a.q_constr for a in atoms):+.6f}")
    print(f"Sum(remaining constr)= {sum(keep_charges):+.6f}")

    ndp_raw = input("Round decimals? (blank = no rounding; e.g., 4): ").strip()
    if ndp_raw:
        ndp = int(ndp_raw)
        keep_charges = [round(q, ndp) for q in keep_charges]

    out_path = Path(args.out).expanduser()
    out_path.write_text("\n".join(f"{q:.10f}" for q in keep_charges) + "\n")
    print(f"\n[DONE] Wrote {len(keep_charges)} charges to {out_path}\n")


if __name__ == "__main__":
    main()
