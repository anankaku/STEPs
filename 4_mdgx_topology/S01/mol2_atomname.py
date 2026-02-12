#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def parse_mol2_atom_names(mol2_path: Path) -> Dict[int, str]:
    """
    Parse mol2 atom names from @<TRIPOS>ATOM section.

    Returns
    -------
    Dict[int, str]
        Mapping: 0-based atom index -> atom name (string as in mol2).
    """
    atom_names: Dict[int, str] = {}

    in_atom_block = False
    with mol2_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.upper().startswith("@<TRIPOS>ATOM"):
                in_atom_block = True
                continue

            if s.upper().startswith("@<TRIPOS>") and not s.upper().startswith("@<TRIPOS>ATOM"):
                # Leaving ATOM block
                if in_atom_block:
                    break
                continue

            if in_atom_block:
                # Typical ATOM line:
                # atom_id  atom_name  x  y  z  atom_type  subst_id  subst_name  charge
                # e.g.:
                # 1        C1        0.0 0.0 0.0 c3       1         RES1        -0.123
                parts = s.split()
                if len(parts) < 2:
                    continue
                try:
                    atom_id_1based = int(parts[0])
                except ValueError:
                    continue
                atom_name = parts[1]
                atom_names[atom_id_1based - 1] = atom_name

    if not atom_names:
        raise ValueError(
            f"[ERROR] No atom names parsed from mol2. "
            f"Check if file has @<TRIPOS>ATOM section: {mol2_path}"
        )
    return atom_names


def load_mol_from_mol2(mol2_path: Path) -> Chem.Mol:
    """
    Load RDKit Mol from mol2. Keep it robust with sanitize options.
    """
    # RDKit: MolFromMol2File may fail for some mol2; sanitize=False then sanitize later.
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"[ERROR] RDKit failed to read mol2: {mol2_path}")

    # Try sanitize; if it fails, still proceed (we mainly need connectivity + 2D depiction)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    # Make sure we have explicit Hs if you want them to show; keep as-is by default.
    return mol


def make_2d(mol: Chem.Mol) -> Chem.Mol:
    """
    Generate 2D coordinates for drawing.
    """
    mol2 = Chem.Mol(mol)
    # Compute 2D coords (works even if no conformer)
    AllChem.Compute2DCoords(mol2)
    return mol2


def draw_with_atom_names(
    mol: Chem.Mol,
    atom_names: Dict[int, str],
    out_path: Path,
    size: Tuple[int, int] = (900, 700),
    font_size: float = 0.8,
    show_hs: bool = True,
    legend: str = "",
) -> None:
    """
    Draw molecule to PNG/SVG with atom name labels.
    """
    w, h = size
    ext = out_path.suffix.lower()

    if ext == ".svg":
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    elif ext in [".png", ".jpg", ".jpeg"]:
        drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    else:
        raise ValueError("[ERROR] Output must end with .png or .svg (or .jpg/.jpeg).")

    opts = drawer.drawOptions()
    # Use atom labels we provide (mol2 atom name)
    # RDKit expects labels via opts.atomLabels dict (key: atom idx)
    opts.atomLabels = {i: atom_names.get(i, "") for i in range(mol.GetNumAtoms())}

    # Tuning
    opts.addAtomIndices = False
    opts.addStereoAnnotation = True
    opts.annotationFontScale = float(font_size)

    # Hide H if user wants
    mol_to_draw = mol
    if not show_hs:
        mol_to_draw = Chem.RemoveHs(mol_to_draw)

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol_to_draw, legend=legend)
    drawer.FinishDrawing()

    data = drawer.GetDrawingText()
    if isinstance(data, str):
        out_path.write_text(data, encoding="utf-8")
    else:
        out_path.write_bytes(data)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert MOL2 to 2D depiction with mol2 atom names labeled."
    )
    ap.add_argument(
        "-i", "--mol2",
        required=True,
        help="Path to input .mol2",
    )
    ap.add_argument(
        "-o", "--out",
        required=True,
        help="Output image path (.png or .svg)",
    )
    ap.add_argument(
        "--size",
        default="900,700",
        help="Image size as W,H (default: 900,700)",
    )
    ap.add_argument(
        "--font-scale",
        type=float,
        default=0.8,
        help="Annotation font scale (default: 0.8). Increase if labels are too small.",
    )
    ap.add_argument(
        "--no-h",
        action="store_true",
        help="Do not show explicit hydrogens in the drawing.",
    )
    ap.add_argument(
        "--legend",
        default="",
        help="Optional legend text under molecule.",
    )
    args = ap.parse_args()

    mol2_path = Path(args.mol2).expanduser().resolve()
    out_path = Path(args.out).expanduser().resolve()

    w_str, h_str = args.size.split(",", 1)
    size = (int(w_str.strip()), int(h_str.strip()))

    atom_names = parse_mol2_atom_names(mol2_path)
    mol = load_mol_from_mol2(mol2_path)
    mol = make_2d(mol)

    draw_with_atom_names(
        mol=mol,
        atom_names=atom_names,
        out_path=out_path,
        size=size,
        font_size=args.font_scale,
        show_hs=(not args.no_h),
        legend=args.legend,
    )

    print(f"[OK] Wrote: {out_path}")


if __name__ == "__main__":
    main()
