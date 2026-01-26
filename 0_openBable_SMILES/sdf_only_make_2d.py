#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List

from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def draw_2d_with_indices_from_sdf(
    in_sdf: Path,
    out_png: Path,
    size: int,
    index_base: int,
) -> None:
    """SDF -> 2D PNG with atom indices (uses SDF atom order). Overwrites out_png."""
    mol = Chem.SDMolSupplier(str(in_sdf), removeHs=False, sanitize=True)[0]
    if mol is None:
        raise ValueError(f"RDKit failed to read SDF: {in_sdf}")

    # Compute 2D coordinates for drawing (does not change atom order)
    AllChem.Compute2DCoords(mol)

    # Annotate each atom with its index (0- or 1-based)
    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", str(atom.GetIdx() + index_base))

    out_png.parent.mkdir(parents=True, exist_ok=True)
    img = Draw.MolToImage(mol, size=(size, size))
    img.save(str(out_png))  # overwrite if exists


def write_atommap_csv_from_sdf(
    in_sdf: Path,
    out_csv: Path,
    index_base: int,
) -> None:
    """SDF -> atommap CSV. Overwrites out_csv."""
    mol = Chem.SDMolSupplier(str(in_sdf), removeHs=False, sanitize=True)[0]
    if mol is None:
        raise ValueError(f"RDKit failed to read SDF: {in_sdf}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:  # overwrite if exists
        w = csv.writer(f)
        w.writerow(["atom_index", "symbol", "atomic_num", "is_h"])
        for atom in mol.GetAtoms():
            idx = atom.GetIdx() + index_base
            sym = atom.GetSymbol()
            w.writerow([idx, sym, atom.GetAtomicNum(), 1 if sym == "H" else 0])


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Batch convert existing SDF files to 2D PNG (with atom indices) and atommap CSV. Overwrites outputs."
    )
    ap.add_argument("--sdfdir", required=True, help="Folder containing .sdf files (e.g., out/sdf)")
    ap.add_argument(
        "--outdir",
        default="out",
        help="Output root folder (should be the same as your pipeline out/ to overwrite png/atommap)",
    )
    ap.add_argument("--png-size", type=int, default=450, help="PNG image size (pixels)")
    ap.add_argument("--index-base", type=int, choices=[0, 1], default=1, help="Atom index base (use 1 to start from 1)")
    args = ap.parse_args()

    sdfdir = Path(args.sdfdir)
    outdir = Path(args.outdir)

    if not sdfdir.exists():
        raise SystemExit(f"[ERROR] sdfdir not found: {sdfdir}")

    sdfs: List[Path] = sorted(sdfdir.glob("*.sdf"))
    if not sdfs:
        raise SystemExit(f"[ERROR] No .sdf files found in: {sdfdir}")

    png_dir = outdir / "png"
    map_dir = outdir / "atommap"
    png_dir.mkdir(parents=True, exist_ok=True)
    map_dir.mkdir(parents=True, exist_ok=True)

    ok, fail = 0, 0
    for sdf in sdfs:
        stem = sdf.stem
        out_png = png_dir / f"{stem}.png"
        out_map = map_dir / f"{stem}_atommap.csv"

        try:
            draw_2d_with_indices_from_sdf(sdf, out_png, size=args.png_size, index_base=args.index_base)
            write_atommap_csv_from_sdf(sdf, out_map, index_base=args.index_base)
            ok += 1
            print(f"[OK] {sdf.name} -> {out_png.name}, {out_map.name} (overwritten if existed)")
        except Exception as e:
            fail += 1
            print(f"[FAIL] {sdf.name}: {e}")

    print(f"[SUMMARY] OK: {ok}  FAIL: {fail}")
    print(f"[OUT] PNG: {png_dir}")
    print(f"[OUT] MAP: {map_dir}")


if __name__ == "__main__":
    main()
