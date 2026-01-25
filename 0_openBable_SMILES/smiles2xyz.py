#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import shutil
import subprocess
from pathlib import Path

import pandas as pd


def require_exe(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise SystemExit(
            f"[ERROR] Cannot find '{name}' in PATH.\n"
            f"Try: conda install -c conda-forge openbabel\n"
        )
    return exe


def run_obabel(obabel: str, smiles: str, out_xyz: Path, gen3d: bool) -> subprocess.CompletedProcess:
    cmd = [obabel, f"-:{smiles}", "-O", str(out_xyz)]
    if gen3d:
        cmd.append("--gen3d")
    return subprocess.run(cmd, text=True, capture_output=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--smiles-col", default="Capped_SMILES")
    ap.add_argument("--outdir", default="xyz_out")
    ap.add_argument("--prefix", default="STEPs")
    ap.add_argument("--gen3d", action="store_true")
    args = ap.parse_args()

    obabel = require_exe("obabel")

    df = pd.read_csv(args.csv)
    if args.smiles_col not in df.columns:
        raise SystemExit(f"[ERROR] Missing column '{args.smiles_col}'. Found: {list(df.columns)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    failures = []
    ok = 0

    for i, smiles in enumerate(df[args.smiles_col].astype(str), start=1):
        smiles = smiles.strip()
        name = f"{args.prefix}_{i:04d}"
        out_xyz = outdir / f"{name}.xyz"

        if not smiles or smiles.lower() == "nan":
            failures.append({"row": i, "name": name, "reason": "empty smiles"})
            continue

        proc = run_obabel(obabel, smiles, out_xyz, gen3d=args.gen3d)

        if proc.returncode == 0 and out_xyz.exists() and out_xyz.stat().st_size > 0:
            ok += 1
        else:
            failures.append(
                {
                    "row": i,
                    "name": name,
                    "smiles": smiles,
                    "returncode": proc.returncode,
                    "stderr": (proc.stderr or "").strip()[:2000],
                }
            )
            if out_xyz.exists() and out_xyz.stat().st_size == 0:
                out_xyz.unlink(missing_ok=True)

    print(f"[SUMMARY] OK: {ok}/{len(df)}  FAIL: {len(failures)}")

    if failures:
        fail_path = outdir / "failures.csv"
        with fail_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=failures[0].keys())
            writer.writeheader()
            writer.writerows(failures)
        print(f"[SUMMARY] failure list saved to: {fail_path}")


if __name__ == "__main__":
    main()
