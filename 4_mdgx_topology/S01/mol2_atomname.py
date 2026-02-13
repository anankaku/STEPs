#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


# -----------------------------
# MOL2 parsing (atom names)
# -----------------------------
def parse_mol2_atom_names(mol2_path: Path) -> Dict[int, str]:
    """
    Parse atom names from the @<TRIPOS>ATOM section of a MOL2 file.

    Returns
    -------
    Dict[int, str]
        Mapping: 0-based atom index -> atom name (exact string from MOL2).
    """
    atom_names: Dict[int, str] = {}
    in_atom_block = False

    with mol2_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            upper = s.upper()
            if upper.startswith("@<TRIPOS>ATOM"):
                # Enter ATOM block
                in_atom_block = True
                continue

            if upper.startswith("@<TRIPOS>") and not upper.startswith("@<TRIPOS>ATOM"):
                # Leave ATOM block when we hit the next section
                if in_atom_block:
                    break
                continue

            if in_atom_block:
                # Typical ATOM line format:
                # atom_id  atom_name  x  y  z  atom_type  subst_id  subst_name  charge
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
            f"[ERROR] Failed to parse atom names. "
            f"Check that the file contains a @<TRIPOS>ATOM section: {mol2_path}"
        )

    return atom_names


# -----------------------------
# RDKit loading + 2D coordinates
# -----------------------------
def load_mol2_rdkit(mol2_path: Path) -> Chem.Mol:
    """
    Load a MOL2 file into an RDKit Mol.

    Notes
    -----
    - We use sanitize=False for robustness (some MOL2 files break RDKit sanitization).
    - We keep removeHs=False to preserve atom indexing as much as possible.
    """
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"[ERROR] RDKit failed to read MOL2: {mol2_path}")

    # Try sanitization; ignore failures because we mainly need connectivity + depiction
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    return mol


def compute_2d_coords(mol: Chem.Mol) -> Chem.Mol:
    """
    Compute 2D coordinates for drawing.
    """
    m = Chem.Mol(mol)
    AllChem.Compute2DCoords(m)
    return m


# -----------------------------
# Drawing
# -----------------------------
def draw_mol_with_names(
    mol: Chem.Mol,
    atom_names: Dict[int, str],
    out_path: Path,
    size: Tuple[int, int] = (900, 700),
    font_scale: float = 0.85,
    hide_h_labels: bool = False,
    add_atom_indices: bool = False,
    legend: str = "",
) -> None:
    """
    Draw a 2D depiction and label atoms with MOL2 atom names.

    Parameters
    ----------
    mol
        RDKit molecule with 2D coordinates.
    atom_names
        Dict mapping RDKit atom index -> MOL2 atom name.
        Assumes RDKit atom ordering matches MOL2 ATOM ordering.
    out_path
        Output path (.png or .svg; .jpg/.jpeg also supported).
    size
        Image size (width, height).
    font_scale
        Scaling factor for annotation text (labels).
    hide_h_labels
        If True, do not show labels on hydrogen atoms (keeps atoms, hides text).
    add_atom_indices
        If True, overlay RDKit atom indices (useful for debugging).
    legend
        Optional legend text under the molecule.
    """
    w, h = size
    ext = out_path.suffix.lower()

    # Choose drawing backend based on output extension
    if ext == ".svg":
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    elif ext in (".png", ".jpg", ".jpeg"):
        drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    else:
        raise ValueError("[ERROR] Output must end with .png or .svg (or .jpg/.jpeg).")

    opts = drawer.drawOptions()
    opts.addAtomIndices = bool(add_atom_indices)
    opts.addStereoAnnotation = True
    opts.annotationFontScale = float(font_scale)

    # Build per-atom label dict (index -> label)
    labels: Dict[int, str] = {}
    n_atoms = mol.GetNumAtoms()

    for i in range(n_atoms):
        label = atom_names.get(i, "")
        if hide_h_labels and mol.GetAtomWithIdx(i).GetAtomicNum() == 1:
            label = ""
        labels[i] = label

    # Override RDKit labels with MOL2 atom names
    opts.atomLabels = labels

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, legend=legend)
    drawer.FinishDrawing()

    # Write output
    data = drawer.GetDrawingText()
    if isinstance(data, str):
        out_path.write_text(data, encoding="utf-8")
    else:
        out_path.write_bytes(data)


# -----------------------------
# CLI
# -----------------------------
def parse_size(s: str) -> Tuple[int, int]:
    """
    Parse W,H string into (width, height).
    """
    try:
        w_str, h_str = s.split(",", 1)
        w = int(w_str.strip())
        h = int(h_str.strip())
        if w <= 0 or h <= 0:
            raise ValueError
        return (w, h)
    except Exception:
        raise ValueError("[ERROR] --size must be formatted as W,H (e.g., 900,700).")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Read MOL2 and draw a 2D depiction with MOL2 atom names labeled."
    )
    ap.add_argument("-i", "--mol2", required=True, help="Path to input .mol2")
    ap.add_argument("-o", "--out", required=True, help="Output image path (.png or .svg)")
    ap.add_argument("--size", default="900,700", help="Image size as W,H (default: 900,700)")
    ap.add_argument("--font-scale", type=float, default=0.85, help="Label font scale (default: 0.85)")
    ap.add_argument("--hide-h-labels", action="store_true", help="Hide labels on hydrogen atoms")
    ap.add_argument("--add-indices", action="store_true", help="Overlay RDKit atom indices (debug)")
    ap.add_argument("--legend", default="", help="Optional legend text under the molecule")
    args = ap.parse_args()

    mol2_path = Path(args.mol2).expanduser().resolve()
    out_path = Path(args.out).expanduser().resolve()
    size = parse_size(args.size)

    # Parse MOL2 atom names and load molecule
    atom_names = parse_mol2_atom_names(mol2_path)
    mol = load_mol2_rdkit(mol2_path)

    # Compute 2D coordinates then draw
    mol2d = compute_2d_coords(mol)
    draw_mol_with_names(
        mol=mol2d,
        atom_names=atom_names,
        out_path=out_path,
        size=size,
        font_scale=args.font_scale,
        hide_h_labels=args.hide_h_labels,
        add_atom_indices=args.add_indices,
        legend=args.legend,
    )

    print(f"[OK] Wrote: {out_path}")


if __name__ == "__main__":
    main()
