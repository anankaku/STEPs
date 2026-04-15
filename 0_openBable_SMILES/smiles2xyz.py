#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Any, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def require_exe(name: str) -> str:
    """
    Ensure an external executable exists in PATH.

    Parameters
    ----------
    name : str
        Executable name (e.g. 'obabel').

    Returns
    -------
    str
        Full path to the executable.

    Raises
    ------
    SystemExit
        If the executable is not found in PATH.
    """
    exe = shutil.which(name)
    if not exe:
        raise SystemExit(
            f"[ERROR] Cannot find '{name}' in PATH.\n"
            f"Try: conda install -c conda-forge openbabel\n"
        )
    return exe


def run(cmd: List[str]) -> subprocess.CompletedProcess:
    """
    Run a subprocess command and capture stdout/stderr.

    Parameters
    ----------
    cmd : list of str
        Command to execute.

    Returns
    -------
    subprocess.CompletedProcess
        Completed process object.
    """
    return subprocess.run(cmd, text=True, capture_output=True)


def format_mol_id(raw_id: Any) -> str:
    """
    Format the molecule ID for file naming.

    Rules
    -----
    - If the ID is numeric (e.g. 1, 5, 12), convert it to zero-padded 4 digits:
      0001, 0005, 0012
    - If the ID is non-numeric (e.g. A01), keep it as-is

    Parameters
    ----------
    raw_id : Any
        Raw ID value from the CSV.

    Returns
    -------
    str
        Formatted ID string.

    Raises
    ------
    ValueError
        If the ID is empty or invalid.
    """
    mol_id = str(raw_id).strip()

    if not mol_id or mol_id.lower() == "nan":
        raise ValueError("empty id")

    try:
        return f"{int(float(mol_id)):04d}"
    except ValueError:
        return mol_id


# -----------------------------------------------------------------------------
# OpenBabel interfaces
# -----------------------------------------------------------------------------

def obabel_smiles_to_sdf(
    obabel: str,
    smiles: str,
    out_sdf: Path,
    gen3d: bool,
) -> subprocess.CompletedProcess:
    """
    Convert a SMILES string to an SDF file using OpenBabel.

    The SDF file preserves atom ordering and bond connectivity.
    Optionally generates 3D coordinates.

    Parameters
    ----------
    obabel : str
        Path to OpenBabel executable.
    smiles : str
        Input SMILES string.
    out_sdf : Path
        Output SDF file path.
    gen3d : bool
        Whether to generate 3D coordinates (--gen3d).

    Returns
    -------
    subprocess.CompletedProcess
        OpenBabel execution result.
    """
    cmd = [obabel, f"-:{smiles}", "-osdf", "-O", str(out_sdf)]
    if gen3d:
        cmd.append("--gen3d")
    return run(cmd)


def obabel_sdf_to_xyz(
    obabel: str,
    in_sdf: Path,
    out_xyz: Path,
) -> subprocess.CompletedProcess:
    """
    Convert an SDF file to XYZ format using OpenBabel.

    Atom order in the XYZ file follows the atom order in the SDF.
    Note that XYZ format does not store bond information.

    Parameters
    ----------
    obabel : str
        Path to OpenBabel executable.
    in_sdf : Path
        Input SDF file.
    out_xyz : Path
        Output XYZ file.

    Returns
    -------
    subprocess.CompletedProcess
        OpenBabel execution result.
    """
    cmd = [obabel, str(in_sdf), "-oxyz", "-O", str(out_xyz)]
    return run(cmd)


# -----------------------------------------------------------------------------
# RDKit visualization and mapping
# -----------------------------------------------------------------------------

def draw_2d_with_indices_from_sdf(
    in_sdf: Path,
    out_png: Path,
    size: int,
    index_base: int,
) -> None:
    """
    Generate a 2D depiction from an SDF file and annotate atom indices.

    The atom indices shown in the PNG correspond exactly to the atom
    ordering in the SDF (and therefore the XYZ file generated from it).

    Parameters
    ----------
    in_sdf : Path
        Input SDF file.
    out_png : Path
        Output PNG image path.
    size : int
        Image size in pixels (square).
    index_base : int
        Atom index base (0 or 1).
    """
    mol = Chem.SDMolSupplier(
        str(in_sdf),
        removeHs=False,
        sanitize=True,
    )[0]

    if mol is None:
        raise ValueError("RDKit failed to read SDF (mol is None)")

    AllChem.Compute2DCoords(mol)

    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", str(atom.GetIdx() + index_base))

    img = Draw.MolToImage(mol, size=(size, size))
    img.save(str(out_png))


def write_atommap_csv_from_sdf(
    in_sdf: Path,
    out_csv: Path,
    index_base: int,
) -> None:
    """
    Write a CSV file mapping atom indices to element information.

    This file serves as a stable reference between SDF/XYZ atom order
    and downstream QM or MD workflows.

    Parameters
    ----------
    in_sdf : Path
        Input SDF file.
    out_csv : Path
        Output CSV file.
    index_base : int
        Atom index base (0 or 1).
    """
    mol = Chem.SDMolSupplier(
        str(in_sdf),
        removeHs=False,
        sanitize=True,
    )[0]

    if mol is None:
        raise ValueError("RDKit failed to read SDF (mol is None)")

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["atom_index", "symbol", "atomic_num", "is_h"])

        for atom in mol.GetAtoms():
            idx = atom.GetIdx() + index_base
            symbol = atom.GetSymbol()
            atomic_num = atom.GetAtomicNum()
            is_h = 1 if symbol == "H" else 0
            writer.writerow([idx, symbol, atomic_num, is_h])


# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------

def main() -> None:
    """
    Main entry point.

    Pipeline
    --------
    CSV (SMILES + ID)
        -> OpenBabel: SMILES -> SDF
        -> OpenBabel: SDF -> XYZ
        -> RDKit:     SDF -> 2D PNG (with atom indices)
        -> CSV:       atom index mapping

    File naming
    -----------
    Output files are named using the CSV ID column instead of row number.
    Example:
        id = 1   -> STEPs_0001.sdf
        id = 12  -> STEPs_0012.sdf
        id = A01 -> STEPs_A01.sdf
    """
    ap = argparse.ArgumentParser(
        description=(
            "SMILES -> (OpenBabel) SDF -> XYZ, and "
            "(RDKit) 2D PNG with atom indices using the same SDF atom order. "
            "Output files are named using the CSV ID column."
        )
    )
    ap.add_argument("--csv", required=True, help="Input CSV file")
    ap.add_argument("--smiles-col", default="Capped_SMILES", help="Column name containing SMILES")
    ap.add_argument("--id-col", default="id", help="Column name containing molecule ID")
    ap.add_argument("--outdir", default="out", help="Output directory")
    ap.add_argument("--prefix", default="STEPs", help="Output file name prefix")
    ap.add_argument("--gen3d", action="store_true", help="Use OpenBabel --gen3d for 3D generation")
    ap.add_argument("--png-size", type=int, default=450, help="PNG image size (pixels)")
    ap.add_argument("--index-base", type=int, choices=[0, 1], default=1, help="Atom index base")

    args = ap.parse_args()

    obabel = require_exe("obabel")

    df = pd.read_csv(args.csv)

    if args.smiles_col not in df.columns:
        raise SystemExit(
            f"[ERROR] Missing column '{args.smiles_col}'. "
            f"Found: {list(df.columns)}"
        )

    if args.id_col not in df.columns:
        raise SystemExit(
            f"[ERROR] Missing column '{args.id_col}'. "
            f"Found: {list(df.columns)}"
        )

    outdir = Path(args.outdir)
    sdf_dir = outdir / "sdf"
    xyz_dir = outdir / "xyz"
    png_dir = outdir / "png"
    map_dir = outdir / "atommap"

    for d in (sdf_dir, xyz_dir, png_dir, map_dir):
        d.mkdir(parents=True, exist_ok=True)

    failures: List[Dict[str, Any]] = []
    ok = 0

    # Process each row using the CSV ID for file naming
    for row_i, row in enumerate(df.to_dict(orient="records"), start=1):
        smiles = str(row[args.smiles_col]).strip()

        try:
            mol_id_fmt = format_mol_id(row[args.id_col])
        except ValueError as e:
            failures.append(
                {
                    "row": row_i,
                    "reason": str(e),
                }
            )
            continue

        name = f"{args.prefix}_{mol_id_fmt}"
        # If you want file names to be exactly the ID only, use:
        # name = mol_id_fmt

        out_sdf = sdf_dir / f"{name}.sdf"
        out_xyz = xyz_dir / f"{name}.xyz"
        out_png = png_dir / f"{name}.png"
        out_map = map_dir / f"{name}_atommap.csv"

        if not smiles or smiles.lower() == "nan":
            failures.append(
                {
                    "row": row_i,
                    "id": row[args.id_col],
                    "name": name,
                    "reason": "empty SMILES",
                }
            )
            continue

        try:
            # Step 1: SMILES -> SDF
            p1 = obabel_smiles_to_sdf(obabel, smiles, out_sdf, gen3d=args.gen3d)
            if p1.returncode != 0 or not out_sdf.exists() or out_sdf.stat().st_size == 0:
                raise RuntimeError(
                    "OpenBabel SMILES->SDF failed. "
                    f"stderr: {(p1.stderr or '').strip()[:2000]}"
                )

            # Step 2: SDF -> XYZ
            p2 = obabel_sdf_to_xyz(obabel, out_sdf, out_xyz)
            if p2.returncode != 0 or not out_xyz.exists() or out_xyz.stat().st_size == 0:
                raise RuntimeError(
                    "OpenBabel SDF->XYZ failed. "
                    f"stderr: {(p2.stderr or '').strip()[:2000]}"
                )

            # Step 3: 2D PNG with atom indices
            draw_2d_with_indices_from_sdf(
                out_sdf,
                out_png,
                size=args.png_size,
                index_base=args.index_base,
            )

            # Step 4: Atom index mapping table
            write_atommap_csv_from_sdf(
                out_sdf,
                out_map,
                index_base=args.index_base,
            )

            ok += 1

        except Exception as e:
            failures.append(
                {
                    "row": row_i,
                    "id": row[args.id_col],
                    "name": name,
                    "smiles": smiles,
                    "reason": str(e),
                }
            )

            for p in (out_sdf, out_xyz, out_png, out_map):
                if p.exists():
                    p.unlink()

    print(f"[SUMMARY] OK: {ok}/{len(df)}  FAIL: {len(failures)}")
    print(f"[OUT] SDF: {sdf_dir}")
    print(f"[OUT] XYZ: {xyz_dir}")
    print(f"[OUT] PNG: {png_dir}")
    print(f"[OUT] MAP: {map_dir}")

    if failures:
        fail_path = outdir / "failures.csv"
        with fail_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=failures[0].keys())
            writer.writeheader()
            writer.writerows(failures)
        print(f"[SUMMARY] failure list saved to: {fail_path}")


if __name__ == "__main__":
    main()