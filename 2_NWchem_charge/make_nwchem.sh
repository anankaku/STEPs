#!/bin/bash
set -euo pipefail

TEMPLATE="template_base.nw"
OUTROOT="outputs"
mkdir -p "${OUTROOT}"

read -r -p "Enter xyz file path: " XYZ
if [[ ! -f "${XYZ}" ]]; then
  echo "ERROR: file not found: ${XYZ}"
  exit 1
fi

base="$(basename "$(dirname "${XYZ}")")"
job="${base}"
workdir="${OUTROOT}/${job}"
mkdir -p "${workdir}"

# copy template
cp "${TEMPLATE}" "${workdir}/${job}.nw"

# set jobname
sed -i "s/JOBNAME/${job}/g" "${workdir}/${job}.nw"

# insert geometry (skip first 2 lines of xyz)
geom="$(tail -n +3 "${XYZ}")"
sed -i "/GEOMETRY_PLACEHOLDER/r /dev/stdin" "${workdir}/${job}.nw" <<< "${geom}"
sed -i "/GEOMETRY_PLACEHOLDER/d" "${workdir}/${job}.nw"

# keep original xyz for traceability
cp "${XYZ}" "${workdir}/input.xyz"

echo
echo "Created: ${workdir}/${job}.nw"
echo "Next: edit esp constrain block, then run:"
echo "  (cd ${workdir} && nwchem ${job}.nw > ${job}.out)"
