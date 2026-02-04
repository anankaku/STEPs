# Molecular Topology Generation

*(NWChem → AmberTools → AMBER / mdgx compatible)*

This document serves as a technical record for converting **QM-derived electrostatic potentials from NWChem** into **AMBER-compatible force field parameters** using **AmberTools**.
The workflow is designed to support **GAFF2 parameterization and downstream mdgx-based dihedral fitting**, with an emphasis on **deterministic atom identity and reproducibility**.

---

## Overview and Rationale

Although PDB files are commonly used as intermediate representations in molecular workflows, this pipeline intentionally adopts a **MOL2-centered approach**.

**Rationale for using MOL2 instead of PDB**

* PDB files contain atomic coordinates and atom names, but **do not explicitly encode force-field–relevant chemical information**, such as atom types, bond orders, or aromaticity.
* Atom names in PDB (e.g., `C1`, `CA`) are **labels**, not force-field atom types (e.g., `c3`, `ca`, `n` in GAFF2).
* In force field parameterization, AmberTools must infer atom types from **connectivity and chemical context**, which is ambiguous when using geometry-only PDB representations.

By contrast, the **MOL2 format explicitly encodes atomic connectivity**, allowing AmberTools to **deterministically assign GAFF2 atom types** and ensuring consistent RESP charge mapping.
This step is **not redundant**, but a deliberate design choice to define atom identity prior to topology generation.

---

## 1. Raw NWChem ESP Output

The following block contains the electrostatic potential data extracted from the NWChem output file.
The **`ESP constr`** column is used to ensure total charge neutrality under capping constraints.

```text
    Atom              Coordinates                           Charge


                                                  ESP         ESP   
                                                             constr 
 
    1 C     0.279138   -0.142063    0.070878   -0.219403   -0.291620
    2 C     0.221609   -0.013827    0.015578    0.627099    0.693934
    3 O     0.288492    0.086247   -0.000610   -0.607867   -0.611229
    4 N     0.090417   -0.017891   -0.021449   -0.245906   -0.254115
    5 C     0.004338   -0.126924    0.003645    0.168570    0.130571
    6 C    -0.056456   -0.194496   -0.101966    0.322680    0.345772
    7 C    -0.141504   -0.301400   -0.081583   -0.287237   -0.301768
    8 C    -0.167434   -0.342837    0.047862   -0.066070   -0.058558
    9 C    -0.109271   -0.276457    0.154630   -0.145396   -0.145829
   10 C    -0.024508   -0.169419    0.132897   -0.144271   -0.133478
   11 F    -0.031540   -0.154917   -0.229076   -0.264116   -0.265515
   12 C     0.032074    0.102308   -0.077325    0.099868   -0.043612
   13 C    -0.045646    0.178803    0.029938    0.428626    0.611656
   14 O    -0.069751    0.129908    0.138679   -0.556074   -0.571318
   15 N    -0.086195    0.302722   -0.005362   -0.118277   -0.336063
   16 C    -0.053820    0.364495   -0.131786    0.009193    0.164977
   17 C    -0.152205    0.385908    0.092812    0.025476    0.164977
   18 H     0.387003   -0.139299    0.057731    0.072277    0.080966
   19 H     0.237749   -0.229289    0.021000    0.024694    0.051423
   20 H     0.257891   -0.148984    0.177311    0.061465    0.076526
   21 H    -0.185595   -0.350700   -0.166796    0.167425    0.169586
   22 H    -0.233504   -0.426397    0.064995    0.109715    0.108455
   23 H    -0.130848   -0.307642    0.255648    0.107051    0.107697
   24 H     0.017879   -0.114988    0.215669    0.159091    0.152723
   25 H    -0.034103    0.076865   -0.160517    0.037928    0.073867
   26 H     0.115168    0.163290   -0.113780    0.061725    0.073867
   27 H    -0.079935    0.299386   -0.215262    0.000409    0.001018
   28 H     0.052539    0.389166   -0.138438    0.032148    0.001018
   29 H    -0.111123    0.456528   -0.141015    0.050332    0.001018
   30 H    -0.166478    0.326497    0.182838    0.050611    0.001018
   31 H    -0.249262    0.419083    0.055437    0.009252    0.001018
   32 H    -0.090892    0.473129    0.116515    0.028982    0.001018
                                            ------------------------
                                                0.000000    0.000000


```

> **Important:**
> The atom order in this table must exactly match the structure file used in Antechamber.

---

## 2. Charge File Content (`capped_charge.dat`)

The constrained ESP charges are written to a plain-text file for use with `antechamber`.

```text
-0.2916200000
0.6939340000
-0.6112290000
...
0.0010180000
```

**Requirements**

* One charge per atom
* Order must match the structural input
* Total charge must be neutral

---

## 3. Workflow Execution

### Step A: Initial Geometry Conversion (QM → MOL2)

The optimized geometry (`s01.xyz`) from NWChem is converted to **MOL2** format to establish explicit bond connectivity.

```bash
obabel /path/to/your.xyz -O s0X.mol2
```

This step provides:

* Explicit atomic connectivity
* Stable atom ordering
* A chemically interpretable structure for GAFF2 typing

---

### Step B: Charge Assignment and Atom Typing with Antechamber

Two files are generated:

* `s0X.mol2`: fully typed, RESP-charged reference structure
```bash
# Assign RESP charges
antechamber -i geom.mol2 -fi mol2 -c rc -cf capped_charge.dat -o s0X.mol2  -fo mol2  -s 2
```
* `s0X.prepi`: AMBER residue definition used by `tleap`
```bash
# Generate GAFF2 atom types
antechamber -i geom.mol2 -fi mol2 -c rc -cf capped_charge.dat -o s0X.prepi -fo prepi
```

The resulting `s0X.mol2` is used to:

* Verify atom order
* Inspect GAFF2 atom typing
* Serve as a reference structure throughout the pipeline

---

### Step D: Topology Building via TLeap

Create the script `build.tleap` with the following content:

```bash
source leaprc.gaff2

loadAmberPrep s0X.prepi
x = loadmol2 "s0X.mol2"

setBox x vdw 10.0
saveAmberParm x s0X.top s0X.crd
quit
```

Execute:

```bash
tleap -f build.tleap
```

---

## 4. Output Summary

| File         | Description                                   |
| ------------ | --------------------------------------------- |
| `s0X.mol2`   | GAFF2-typed, RESP-charged reference structure |
| `s0X.prepi`  | AMBER residue definition                      |
| `s0X.top`    | AMBER topology file                           |
| `s0X.crd`    | Coordinate file                               |
| Periodic box | 10.0 Å VDW padding                            |

---

## Notes on mdgx / STEPs Compatibility

* `s0X.top` and `s0X.crd` are **direct inputs to mdgx**
* Atom ordering is preserved across:

  * RESP fitting
  * Antechamber
  * TLeap
  * mdgx dihedral scans
* This workflow is fully compatible with:

  * ORCA or NWChem single-point energy evaluation
  * mdgx-based force field fitting
  * STEPs-style parameterization

---

## Summary

This MOL2-centered workflow prioritizes **force-field correctness over geometric convenience**.
By explicitly defining atom connectivity and chemical context prior to topology generation, the pipeline ensures robust, reproducible AMBER force field parameters suitable for systematic mdgx-based fitting.

---
## Running `run_tleap.sh` (Interactive Mode)

### 1. Run the script

From the directory where the script is located:

```bash
bash run_tleap.sh
```

---

### 2. Interactive prompts

After execution, the script will prompt you for the required inputs:

```text
prepi path (relative OK):
```

Enter the path to the `.prepi` file.
Relative paths (relative to your current working directory) are supported.

Example:

```text
../path/to/your.prepi
```

---

```text
mol2  path (relative OK):
```

Enter the path to the `.mol2` file.

Example:

```text
../path/to/your.mol2
```

---

```text
output dir (relative OK):
```

Enter the directory where output files will be written.
The directory will be created automatically if it does not exist.

Example:

```text
../path/you/want
```

---

```text
vdw padding (default 10.0):
```

Enter the van der Waals padding (in Å) used by `setBox`.

* Press **Enter** to accept the default value (`10.0`)
* Or enter a custom value (e.g. `12.0`)

---

### 3. Output files

The script generates the following files in the specified output directory:

```text
<basename>.top
<basename>.crd
```

Example:

```text
../path/to/create/your.top
../path/to/create/your.crd
```

---

### 4. Notes on default behavior

* All input paths (relative or absolute) are internally converted to absolute paths
* `leaprc.gaff2` is sourced automatically
* No `tleap` input files are copied into data directories

This design keeps the workflow clean and reproducible.

---

### 5. First-time setup (only needed once)

If the script is not executable yet, enable execution permission:

```bash
chmod +x run_tleap.sh
```

After this, the script can be run normally with `./run_tleap.sh` or `bash run_tleap.sh`.

---
