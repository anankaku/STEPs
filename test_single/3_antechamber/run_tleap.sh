#!/usr/bin/env bash
set -eu

# --- helpers ---
abspath () {
  python - <<EOF
import os
print(os.path.abspath("$1"))
EOF
}

# --- interactive inputs (relative path OK) ---
read -rp "prepi path (relative OK): " PREPI_IN
read -rp "mol2  path (relative OK): " MOL2_IN
read -rp "output dir (relative OK): " OUTDIR_IN
read -rp "vdw padding (default 10.0): " PAD
PAD="${PAD:-10.0}"

# --- normalize to absolute paths ---
PREPI="$(abspath "${PREPI_IN}")"
MOL2="$(abspath "${MOL2_IN}")"
OUTDIR="$(abspath "${OUTDIR_IN}")"

mkdir -p "${OUTDIR}"

BASE="$(basename "${MOL2%.*}")"

# --- run tleap ---
tleap -f - <<EOF
source leaprc.gaff2

loadAmberPrep ${PREPI}
x = loadmol2 "${MOL2}"

setBox x vdw ${PAD}
saveAmberParm x ${OUTDIR}/${BASE}.top ${OUTDIR}/${BASE}.crd
quit
EOF

echo "[OK] wrote ${OUTDIR}/${BASE}.{top,crd}"