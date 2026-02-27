#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, Draw


DIHEDRAL_LINE = re.compile(r"^\s*(GridSample|RandomSample)\s+(.*)$")
ATOMREF = re.compile(r":\d+@([A-Za-z0-9_]+)")


def load_mol2(mol2_path: Path) -> Chem.Mol:
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
    if mol is None:
        raise SystemExit(f"[ERROR] RDKit failed to read mol2: {mol2_path}")
    # Ensure 2D coords exist
    AllChem.Compute2DCoords(mol)
    return mol


def get_atomname_map(mol: Chem.Mol) -> Dict[str, int]:
    """
    Return mapping: Tripos atom name -> RDKit atom index (0-based).
    RDKit stores mol2 atom names in atom prop '_TriposAtomName' (commonly).
    """
    name_to_idx: Dict[str, int] = {}
    missing = 0
    for a in mol.GetAtoms():
        if a.HasProp("_TriposAtomName"):
            nm = a.GetProp("_TriposAtomName").strip()
        elif a.HasProp("atomName"):
            nm = a.GetProp("atomName").strip()
        else:
            nm = ""
        if not nm:
            missing += 1
            continue
        # If duplicates exist, keep first but warn later by counting duplicates
        if nm not in name_to_idx:
            name_to_idx[nm] = a.GetIdx()
    if missing:
        print(f"[WARN] {missing} atoms have no stored mol2 atom name property in RDKit.")
        print("       If labels are missing, consider converting mol2 -> sdf with OpenBabel while preserving atom names.")
    return name_to_idx


def parse_mdgx_dihedrals(mdgx_in_path: Path) -> List[Tuple[str, List[str], str]]:
    """
    Return list of (kind, [a1,a2,a3,a4], original_line)
    where ai are names after '@' (e.g., C2, N1, C10, C11).
    """
    dihs = []
    for line in mdgx_in_path.read_text().splitlines():
        m = DIHEDRAL_LINE.match(line)
        if not m:
            continue
        kind = m.group(1)
        rest = m.group(2)
        atoms = ATOMREF.findall(rest)
        if len(atoms) >= 4:
            dihs.append((kind, atoms[:4], line.rstrip()))
    return dihs


def draw_labeled_2d(mol: Chem.Mol, outpath: Path, name_to_idx: Dict[str, int]) -> None:
    labels = {}
    for a in mol.GetAtoms():
        idx = a.GetIdx()
        nm = a.GetProp("_TriposAtomName").strip() if a.HasProp("_TriposAtomName") else str(idx + 1)
        # 同時顯示 atom name + index（比較不會迷路）
        labels[idx] = f"{nm}({idx+1})"
    img = Draw.MolToImage(mol, size=(900, 700), atomLabels=labels)
    img.save(str(outpath))


def draw_dihedral_highlight(
    mol: Chem.Mol,
    outpath: Path,
    atoms4: List[str],
    name_to_idx: Dict[str, int],
) -> None:
    # highlight atoms; RDKit 允許傳入 highlightAtoms list
    highlight = []
    for nm in atoms4:
        if nm in name_to_idx:
            highlight.append(name_to_idx[nm])
    # label all atoms too
    labels = {}
    for a in mol.GetAtoms():
        idx = a.GetIdx()
        nm = a.GetProp("_TriposAtomName").strip() if a.HasProp("_TriposAtomName") else str(idx + 1)
        labels[idx] = nm
    img = Draw.MolToImage(
        mol,
        size=(900, 700),
        atomLabels=labels,
        highlightAtoms=highlight,
    )
    img.save(str(outpath))


def main():
    if len(sys.argv) != 4:
        print("Usage: python check_dihedral_atomnames.py <mol2> <mdgx_input> <outdir>")
        sys.exit(2)

    mol2_path = Path(sys.argv[1]).resolve()
    mdgx_path = Path(sys.argv[2]).resolve()
    outdir = Path(sys.argv[3]).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    mol = load_mol2(mol2_path)
    name_to_idx = get_atomname_map(mol)
    dihs = parse_mdgx_dihedrals(mdgx_path)

    # 1) Report atom name availability
    print("=== Atom-name check against mdgx dihedrals ===")
    missing_any = False
    for i, (kind, atoms4, line) in enumerate(dihs, 1):
        missing = [nm for nm in atoms4 if nm not in name_to_idx]
        ok = "OK" if not missing else f"MISSING: {missing}"
        if missing:
            missing_any = True
        print(f"[{i:02d}] {kind} {atoms4} -> {ok}")
        print(f"     {line}")

    if missing_any:
        print("\n[WARN] Some atom names referenced by mdgx are missing in the mol2 atom-name list.")
        print("       This usually means your mol2 atom names differ (e.g., C10 vs C1) or names were not preserved.")
        print("       Fix atom names in mol2 or regenerate mol2 with the desired naming scheme.\n")

    # 2) Draw full labeled 2D
    full_img = outdir / "mol2_atomnames_2D.png"
    draw_labeled_2d(mol, full_img, name_to_idx)
    print(f"[OK] Wrote labeled 2D: {full_img}")

    # 3) Draw per-dihedral highlight images
    for i, (_, atoms4, _) in enumerate(dihs, 1):
        img_path = outdir / f"dihedral_{i:02d}_{'-'.join(atoms4)}.png"
        draw_dihedral_highlight(mol, img_path, atoms4, name_to_idx)

    print(f"[OK] Wrote {len(dihs)} dihedral highlight images in: {outdir}")


if __name__ == "__main__":
    main()
