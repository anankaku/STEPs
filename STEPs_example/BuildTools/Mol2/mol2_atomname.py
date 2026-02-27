#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


# -----------------------------
# MOL2 parsing (ATOM/BOND)
# -----------------------------
def parse_mol2_atoms_bonds(mol2_path: Path) -> Tuple[List[dict], List[Tuple[int, int, str]]]:
    atoms: List[dict] = []
    bonds: List[Tuple[int, int, str]] = []

    in_atom = False
    in_bond = False

    with mol2_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            up = s.upper()

            if up.startswith("@<TRIPOS>ATOM"):
                in_atom, in_bond = True, False
                continue
            if up.startswith("@<TRIPOS>BOND"):
                in_atom, in_bond = False, True
                continue
            if up.startswith("@<TRIPOS>"):
                in_atom, in_bond = False, False
                continue

            if in_atom:
                parts = s.split()
                if len(parts) < 6:
                    continue
                atom_id1 = int(parts[0])
                name = parts[1]
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                atype = parts[5]
                atoms.append({"id1": atom_id1, "name": name, "x": x, "y": y, "z": z, "atype": atype})

            elif in_bond:
                parts = s.split()
                if len(parts) < 4:
                    continue
                a1 = int(parts[1])
                a2 = int(parts[2])
                btype = parts[3].lower()
                bonds.append((a1, a2, btype))

    if not atoms:
        raise ValueError(f"[ERROR] No atoms parsed from MOL2: {mol2_path}")

    return atoms, bonds


def guess_element(atom_name: str, atom_type: str) -> str:
    letters = ""
    for ch in atom_name:
        if ch.isalpha():
            letters += ch
        else:
            break
    letters = letters.capitalize()

    two = letters[:2]
    if two in ("Cl", "Br", "Si", "Na", "Li", "Al", "Ca", "Fe", "Zn", "Mg", "Cu", "Mn"):
        return two
    if letters:
        return letters[0]

    t = atom_type.strip()
    if "." in t:
        t = t.split(".", 1)[0]
    t = "".join([c for c in t if c.isalpha()]).capitalize()
    two = t[:2]
    if two in ("Cl", "Br", "Si", "Na", "Li", "Al", "Ca", "Fe", "Zn", "Mg", "Cu", "Mn"):
        return two
    if t:
        return t[0]

    return "C"


def bondtype_to_rdkit(btype: str) -> Chem.BondType:
    if btype in ("1", "single"):
        return Chem.BondType.SINGLE
    if btype in ("2", "double"):
        return Chem.BondType.DOUBLE
    if btype in ("3", "triple"):
        return Chem.BondType.TRIPLE
    if btype == "ar":
        return Chem.BondType.AROMATIC
    return Chem.BondType.SINGLE


def build_rdkit_from_mol2(mol2_path: Path) -> Tuple[Chem.Mol, Dict[int, str]]:
    atoms, bonds = parse_mol2_atoms_bonds(mol2_path)

    rw = Chem.RWMol()
    id1_to_idx: Dict[int, int] = {}
    idx_to_name: Dict[int, str] = {}

    for a in atoms:
        elem = guess_element(a["name"], a["atype"])
        atom = Chem.Atom(elem)
        idx = rw.AddAtom(atom)
        id1_to_idx[a["id1"]] = idx
        idx_to_name[idx] = a["name"]

    for a1_id1, a2_id1, btype in bonds:
        i = id1_to_idx.get(a1_id1)
        j = id1_to_idx.get(a2_id1)
        if i is None or j is None:
            continue
        bt = bondtype_to_rdkit(btype)
        rw.AddBond(i, j, bt)

        if btype == "ar":
            rw.GetAtomWithIdx(i).SetIsAromatic(True)
            rw.GetAtomWithIdx(j).SetIsAromatic(True)
            b = rw.GetBondBetweenAtoms(i, j)
            if b is not None:
                b.SetIsAromatic(True)

    mol = rw.GetMol()

    conf = Chem.Conformer(mol.GetNumAtoms())
    for a in atoms:
        idx = id1_to_idx[a["id1"]]
        conf.SetAtomPosition(idx, Chem.rdGeometry.Point3D(a["x"], a["y"], a["z"]))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    return mol, idx_to_name


def draw_2d_with_atom_names(
    mol: Chem.Mol,
    idx_to_name: Dict[int, str],
    out_path: Path,
    size: Tuple[int, int] = (900, 700),
    font_scale: float = 0.85,
    hide_h_labels: bool = False,
    legend: str = "",
) -> None:
    w, h = size
    ext = out_path.suffix.lower()

    if ext == ".svg":
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    elif ext in (".png", ".jpg", ".jpeg"):
        drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    else:
        raise ValueError("[ERROR] Output must be .png/.svg/.jpg/.jpeg")

    m = Chem.Mol(mol)
    AllChem.Compute2DCoords(m)

    for i in range(m.GetNumAtoms()):
        label = idx_to_name.get(i, "")
        if hide_h_labels and m.GetAtomWithIdx(i).GetAtomicNum() == 1:
            label = ""
        m.GetAtomWithIdx(i).SetProp("atomLabel", label)

    opts = drawer.drawOptions()
    opts.annotationFontScale = float(font_scale)
    opts.addStereoAnnotation = True

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, m, legend=legend)
    drawer.FinishDrawing()

    data = drawer.GetDrawingText()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(data, str):
        out_path.write_text(data, encoding="utf-8")
    else:
        out_path.write_bytes(data)


def parse_size(s: str) -> Tuple[int, int]:
    try:
        w_str, h_str = s.split(",", 1)
        return (int(w_str.strip()), int(h_str.strip()))
    except Exception as exc:
        raise ValueError("[ERROR] --size must be W,H (e.g., 900,700)") from exc


def collect_inputs(single_mol2: str | None, glob_pattern: str | None, input_dir: str | None) -> List[Path]:
    files: List[Path] = []
    if single_mol2:
        files.append(Path(single_mol2).expanduser().resolve())

    if glob_pattern:
        files.extend(sorted(Path.cwd().glob(glob_pattern)))

    if input_dir:
        root = Path(input_dir).expanduser().resolve()
        files.extend(sorted(root.glob("*.mol2")))

    uniq: List[Path] = []
    seen = set()
    for f in files:
        rf = f.resolve()
        if rf not in seen:
            seen.add(rf)
            uniq.append(rf)
    return uniq


def build_output_path(mol2_path: Path, out: str | None, out_dir: str | None, ext: str) -> Path:
    if out:
        return Path(out).expanduser().resolve()

    suffix = ext if ext.startswith(".") else f".{ext}"
    if out_dir:
        target_dir = Path(out_dir).expanduser().resolve()
    else:
        target_dir = mol2_path.parent
    return target_dir / f"{mol2_path.stem}{suffix}"


def main() -> None:
    ap = argparse.ArgumentParser(description="Draw 2D figures from one or multiple MOL2 files.")
    ap.add_argument("-i", "--mol2", help="Single input .mol2 file")
    ap.add_argument("--glob", dest="glob_pattern", help="Glob pattern for batch mode, e.g. '*.mol2'")
    ap.add_argument("--input-dir", help="Directory containing .mol2 files for batch mode")
    ap.add_argument("-o", "--out", help="Output image path (single input only)")
    ap.add_argument("--out-dir", help="Output directory for batch mode")
    ap.add_argument("--ext", default="png", choices=["png", "svg", "jpg", "jpeg"], help="Output extension in batch mode")
    ap.add_argument("--size", default="900,700", help="Image size as W,H (default: 900,700)")
    ap.add_argument("--font-scale", type=float, default=0.85, help="Label font scale")
    ap.add_argument("--hide-h-labels", action="store_true", help="Hide labels on hydrogen atoms")
    ap.add_argument("--legend", default="", help="Optional legend text")
    args = ap.parse_args()

    inputs = collect_inputs(args.mol2, args.glob_pattern, args.input_dir)
    if not inputs:
        raise SystemExit("[ERROR] No input MOL2 files. Use -i, --glob, or --input-dir.")

    if args.out and len(inputs) > 1:
        raise SystemExit("[ERROR] -o/--out only supports single input. Use --out-dir for batch mode.")

    size = parse_size(args.size)

    for mol2_path in inputs:
        out_path = build_output_path(mol2_path, args.out, args.out_dir, args.ext)
        mol, idx_to_name = build_rdkit_from_mol2(mol2_path)
        draw_2d_with_atom_names(
            mol=mol,
            idx_to_name=idx_to_name,
            out_path=out_path,
            size=size,
            font_scale=args.font_scale,
            hide_h_labels=args.hide_h_labels,
            legend=args.legend,
        )
        print(f"[OK] {mol2_path.name} -> {out_path}")


if __name__ == "__main__":
    main()
