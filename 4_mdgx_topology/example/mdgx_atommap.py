#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from io import BytesIO
import math
from PIL import Image, ImageDraw, ImageFont


# =============================================================================
# Parsing mdgx config lines (GridSample / RandomSample)
# =============================================================================

# Match GridSample or RandomSample lines
LINE_PAT = re.compile(r"^\s*(GridSample|RandomSample)\s+(.*)$")

# Extract atom names from patterns like :1@C2
ATOM_PAT = re.compile(r":\d+@([A-Za-z0-9_]+)")

# Extract angular range { min max }
RANGE_PAT = re.compile(r"\{\s*([\-0-9\.]+)\s+([\-0-9\.]+)\s*\}")

# Extract restraint force constant Krst
KRST_PAT = re.compile(r"\bKrst\s+([0-9\.]+)")


def parse_configs(genconf: Path) -> List[dict]:
    """
    Parse GridSample / RandomSample entries from mdgx input file.

    Each config line becomes one dictionary describing:
    - sampling type (GridSample or RandomSample)
    - geometry type (dihedral or angle)
    - involved atom names
    - sampling range
    - restraint force constant
    """
    items = []
    for raw in genconf.read_text(errors="ignore").splitlines():
        m = LINE_PAT.match(raw)
        if not m:
            continue

        kind = m.group(1)
        rest = m.group(2)

        atoms = ATOM_PAT.findall(rest)
        rng = RANGE_PAT.search(rest)
        krst = KRST_PAT.search(rest)

        # GridSample → dihedral (4 atoms)
        # RandomSample → angle (3 atoms)
        if kind == "GridSample":
            atoms = atoms[:4]
            geom = "dihedral"
        else:
            atoms = atoms[:3]
            geom = "angle"

        items.append({
            "kind": kind,
            "geom": geom,
            "atoms": atoms,
            "range": (rng.group(1), rng.group(2)) if rng else None,
            "krst": krst.group(1) if krst else None,
            "raw": raw.strip(),
        })

    return items


# =============================================================================
# Load reference structure (PDB or MOL2)
# =============================================================================

def load_ref_mol(path: Path) -> Chem.Mol:
    """
    Load the reference molecular structure.

    The reference structure is used ONLY for drawing 2D topology
    and atom name mapping. It is not modified.
    """
    suf = path.suffix.lower()
    if suf == ".mol2":
        mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    elif suf in {".pdb", ".ent"}:
        mol = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
    else:
        raise SystemExit(f"[ERROR] Unsupported reference file: {path}")

    if mol is None:
        raise SystemExit(f"[ERROR] RDKit failed to read: {path}")

    return mol


def choose_reference(root: Path) -> Path:
    """
    Automatically select a reference structure from the current directory.

    Priority:
    1) non-Conf*.mol2
    2) non-Conf*.pdb
    """
    mol2s = sorted(
        p for p in root.iterdir()
        if p.is_file() and p.suffix.lower() == ".mol2" and not p.name.startswith("Conf")
    )
    if mol2s:
        return mol2s[0]

    pdbs = sorted(
        p for p in root.iterdir()
        if p.is_file() and p.suffix.lower() in {".pdb", ".ent"} and not p.name.startswith("Conf")
    )
    if pdbs:
        return pdbs[0]

    raise SystemExit("[ERROR] No reference (.mol2 or .pdb) found in current directory.")


# =============================================================================
# Atom name handling
# =============================================================================

def atom_name(a: Chem.Atom) -> str:
    """
    Retrieve atom name with priority:
    1) MOL2 Tripos atom name
    2) PDB atom name
    3) Fallback: element + atom index
    """
    if a.HasProp("_TriposAtomName"):
        nm = a.GetProp("_TriposAtomName").strip()
        if nm:
            return nm

    info = a.GetPDBResidueInfo()
    if info is not None:
        nm = info.GetName().strip()
        if nm:
            return nm

    return f"{a.GetSymbol()}{a.GetIdx()+1}"


def name_to_idx(mol: Chem.Mol) -> Dict[str, int]:
    """
    Build a mapping from atom name → RDKit atom index.

    If duplicate atom names exist, the first occurrence is used.
    """
    mp: Dict[str, int] = {}
    for a in mol.GetAtoms():
        nm = atom_name(a)
        if nm not in mp:
            mp[nm] = a.GetIdx()
    return mp


# =============================================================================
# Drawing: one 2D image per config line
# =============================================================================

def draw_one(
    mol: Chem.Mol,
    mp: Dict[str, int],
    item: dict,
    out_png: Path,
    size: Tuple[int, int] = (1600, 1600),
) -> List[str]:
    """
    Draw a single 2D structure image corresponding to ONE mdgx config line.

    - The full molecule is drawn by RDKit (no atom labels).
    - Only bonds involved in the angle/dihedral are highlighted.
    - Atom names are manually overlaid using PIL to avoid label overlap
      and RDKit version incompatibilities.
    """
    atoms = item["atoms"]
    missing = [n for n in atoms if n not in mp]

    # Build legend text
    rng = item["range"]
    krst = item["krst"]
    legend_parts = [f"{item['kind']} ({item['geom']})", "  ".join(atoms)]
    if rng:
        legend_parts.append(f"range {rng[0]}~{rng[1]} deg")
    if krst:
        legend_parts.append(f"Krst {krst}")
    legend = " | ".join(legend_parts)

    # Compute expanded 2D layout for better readability
    AllChem.Compute2DCoords(mol, bondLength=2.2)

    w, h = size
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    opts = drawer.drawOptions()
    opts.baseFontSize = 24
    opts.legendFontSize = 24

    # Disable RDKit atom labels (we draw labels manually later)
    for a in mol.GetAtoms():
        opts.atomLabels[a.GetIdx()] = ""

    # Highlight only bonds relevant to the angle/dihedral
    highlight_bonds = []
    for a, b in zip(atoms[:-1], atoms[1:]):
        if a in mp and b in mp:
            bond = mol.GetBondBetweenAtoms(mp[a], mp[b])
            if bond:
                highlight_bonds.append(bond.GetIdx())

    if hasattr(opts, "highlightBondWidthMultiplier"):
        opts.highlightBondWidthMultiplier = 18

    # Draw base molecule
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        legend=legend,
        highlightBonds=highlight_bonds,
    )

    # Retrieve pixel coordinates for selected atoms
    idxs = [mp[n] for n in atoms if n in mp]
    pts = [drawer.GetDrawCoords(i) for i in idxs] if idxs else []
    cx = sum(p.x for p in pts) / len(pts) if pts else w / 2
    cy = sum(p.y for p in pts) / len(pts) if pts else h / 2

    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()

    # -------------------------------------------------------------------------
    # Overlay atom labels using PIL (robust across RDKit versions)
    # -------------------------------------------------------------------------
    img = Image.open(BytesIO(png_bytes)).convert("RGBA")
    draw = ImageDraw.Draw(img)

    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 36)
    except Exception:
        font = ImageFont.load_default()

    # Offset parameters to reduce label overlap
    extra_y = [-40, 40, -40, 40]
    off = 25

    for k, name in enumerate(atoms):
        if name not in mp:
            continue

        i = mp[name]
        p = drawer.GetDrawCoords(i)

        # Push label outward from the atom cluster center
        vx, vy = (p.x - cx), (p.y - cy)
        norm = math.hypot(vx, vy) or 1.0
        vx, vy = vx / norm, vy / norm

        x = p.x + vx * off
        y = p.y + vy * off + (extra_y[k] if k < len(extra_y) else 0)

        # Draw outlined text for readability
        for dx, dy in [(-2,0), (2,0), (0,-2), (0,2)]:
            draw.text((x+dx, y+dy), name, font=font, fill=(255,255,255,255))
        draw.text((x, y), name, font=font, fill=(0,0,0,255))

    img.convert("RGB").save(out_png)
    return missing


# =============================================================================
# Main entry point
# =============================================================================

def main():
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage:\n"
            "  python mdgx_angles_2Dmaps.py genconf.in outdir [ref.mol2|ref.pdb]"
        )

    genconf = Path(sys.argv[1]).resolve()
    outdir = Path(sys.argv[2]).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    root = Path.cwd()
    ref = Path(sys.argv[3]).resolve() if len(sys.argv) >= 4 else choose_reference(root)

    items = parse_configs(genconf)
    if not items:
        raise SystemExit("[ERROR] No GridSample/RandomSample lines found in genconf.in")

    mol = load_ref_mol(ref)
    mp = name_to_idx(mol)

    print(f"[INFO] Reference: {ref.name}")
    print(f"[INFO] Found {len(items)} config lines -> generating {len(items)} PNG files")

    any_missing = False
    for i, item in enumerate(items, 1):
        atoms_part = "-".join(item["atoms"]) if item["atoms"] else "NOATOMS"
        out_png = outdir / f"{i:02d}_{item['kind']}_{atoms_part}.png"
        missing = draw_one(mol, mp, item, out_png)
        if missing:
            any_missing = True
            print(f"[WARN] {out_png.name} missing atom names: {missing}")

    print(f"[OK] Output written to: {outdir}")
    if any_missing:
        print("[NOTE] Missing names usually indicate atom-name mismatch "
              "between reference structure and mdgx input.")


if __name__ == "__main__":
    main()
